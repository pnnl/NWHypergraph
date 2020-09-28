//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include "Log.hpp"
#include <unordered_set>
#include <algorithms/connected_components.hpp>
#include "common.hpp"
#include "mmio_hy.hpp"
#include <edge_list.hpp>
#include <util/AtomicBitVector.hpp>
#include <util/atomic.hpp>
#include <util/intersection_size.hpp>
#include <docopt.h>
#include <execution>

using namespace nw::hypergraph::bench;
using namespace nw::graph;

static constexpr const char USAGE[] =
    R"(hyccrelabel.exe: nw::graph hypergraph connected components benchmark driver.
  Usage:
      hyccrelabel.exe (-h | --help)
      hyccrelabel.exe [-f FILE...] [--version ID...] [-n NUM] [--succession STR] [--relabel] [--clean] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      -f FILE               input file paths (can have multiples)
      -n NUM                number of trials [default: 1]
      --relabel             relabel the graph or not
      -c, --clean           clean the graph or not
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
      --succession STR      successor/predecessor [default: successor]
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
      -v, --verify          verify results
      -V, --verbose         run in verbose mode
)";

template<class T>
inline bool writeMin(T& old, T& next) {
  T    prev;
  bool success;
  do
    prev = old;
  while (prev > next && !(success = nw::graph::cas(old, prev, next)));
  return success;
}


template<class ExecutionPolicy, class Graph, class Transpose>
auto relabelHyperCC(ExecutionPolicy&& exec, Graph& g, Transpose& g_t, const size_t num_realedges, const size_t num_realnodes) {

  auto labeling   = //Afforest(s_adj, s_trans_adj); 
  ccv1(g);//

  std::vector<vertex_id_t> E(num_realedges), N(num_realnodes);
  if (num_realnodes < num_realedges) {
    nw::util::life_timer _("unrelabeling");
    //E.assign(labeling.begin(), labeling.begin() + num_realedges);
    std::for_each(exec, tbb::counting_iterator(0ul), tbb::counting_iterator(num_realedges), [&](auto i) {
      E[i] = labeling[i];
    });
    //N.assign(labeling.begin() + num_realedges, labeling.end());
    std::for_each(exec, tbb::counting_iterator(0ul), tbb::counting_iterator(num_realnodes), [&](auto i) {
      N[i] = labeling[i + num_realedges];
    }); 
  }
  else {
    nw::util::life_timer _("unrelabeling");
    //N.assign(labeling.begin(), labeling.begin() + num_realnodes);
    std::for_each(exec, tbb::counting_iterator(0ul), tbb::counting_iterator(num_realnodes), [&](auto i) {
      N[i] = labeling[i];
    }); 
    //E.assign(labeling.begin() + num_realnodes, labeling.end());
    std::for_each(exec, tbb::counting_iterator(0ul), tbb::counting_iterator(num_realedges), [&](auto i) {
      E[i] = labeling[i + num_realnodes];
    });
  }

  return std::tuple{N, E};
}

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto                     args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;

  std::vector<long> ids     = parse_ids(args["--version"].asStringList());
  std::vector<long> threads = parse_n_threads(args["THREADS"].asStringList());

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times<bool> times;

  // Appease clang.
  //
  // These are captured in lambdas later, and if I use structured bindings
  // they have to be listed as explicit captures (this is according to the 17
  // standard). That's a little bit noisy where it happens, so I just give
  // them real symbols here rather than the local bindings.
  for (auto&& file : files) {
    auto reader = [&](std::string file, bool verbose, size_t& nrealedges, size_t& nrealnodes) {
      //auto aos_a   = load_graph<directed>(file);
      auto aos_a   = read_mm_relabeling<nw::graph::directed>(file, nrealedges, nrealnodes);
      
      // Run relabeling. This operates directly on the incoming edglist.
      if (args["--relabel"].asBool()) {
        auto degrees = aos_a.degrees();
        aos_a.relabel_by_degree<0>(args["--direction"].asString(), degrees);
      }
      // Clean up the edgelist to deal with the normal issues related to
      // undirectedness.
      if (args["--clean"].asBool()) {
        aos_a.swap_to_triangular<0>(args["--succession"].asString());
        aos_a.lexical_sort_by<0>();
        aos_a.uniq();
        aos_a.remove_self_loops();
      }

      adjacency<0> g(aos_a);
      adjacency<1> g_t(aos_a);
      if (verbose) {
        g_t.stream_stats();
        g.stream_stats();
      }
      std::cout << "num_hyperedges = " << nrealedges << " num_hypernodes = " << nrealnodes << std::endl;
      return std::tuple(aos_a, g, g_t);
    };
    size_t num_realedges, num_realnodes;
    auto&& [aos_a, g, g_t]     = reader(file, verbose, num_realedges, num_realnodes);

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        if (verbose) {
          std::cout << "version " << id << std::endl;
        }

        auto verifier = [&](auto&& result) {
          //TODO This returns the subgraph of each component.
          //Does not work for s overlap
          auto&& [N, E] = result;
          if (verbose) {
            //print_top_n(graph, comp);
            std::map<vertex_id_t, edge_list<>> comps;
            std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
              auto&& [edge, node] = elt;
              vertex_id_t key     = E[edge];
              comps[key].push_back(elt);
            });

            for (auto&& j : comps) {
              auto& [k, v] = j;
              v.close_for_push_back();
            }
            std::cout << comps.size() << " subgraphs and" << std::endl;
          }
          std::unordered_set<vertex_id_t> uni_comps(E.begin(), E.end());
          std::cout << uni_comps.size() << " components found" << std::endl;

          if (verify) {
            std::cerr << " v" << id << " failed verification for " << file << " using " << thread << " threads\n";
          }
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, std::forward<decltype(op)>(op), verifier, true); };
        std::vector<vertex_id_t> E, N;

        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { return relabelHyperCC(std::execution::seq, g, g_t, num_realnodes, num_realnodes); });
              break;
            case 1:
              record([&] { return relabelHyperCC(std::execution::par_unseq, g, g_t, num_realnodes, num_realnodes); });
              break;
            default:
              std::cout << "Unknown version v" << id << "\n";
          }
        }
      }
    }
  }

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("cc", file, times, header, "Time(s)");
  }

  return 0;
}
