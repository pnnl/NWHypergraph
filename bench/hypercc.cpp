//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include <unordered_set>
#include <docopt.h>
#include <edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "edge_list_hy.hpp"
#include "s_overlap.hpp"
#include "algorithms/hyper_connected_components.hpp"
#include "algorithms/s_connected_components.hpp"


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(hycc.exe: hypergraph connected components benchmark driver.
  Usage:
      hycc.exe (-h | --help)
      hycc.exe [-f FILE...] [-a FILE...] [--version ID...] [-B NUM] [-n NUM] [--direction DIR] [--relabel NUM] [-cdV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      -f FILE               edge list or matrix market input file paths (can have multiples)
      -a FILE               hypergraph adjacency fils paths (can have multiples)
      -n NUM                number of trials [default: 1]
      -B NUM                number of bins [default: 32]
      --relabel NUM         relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      -c, --clean           clean the graph or not
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";



int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;
  long num_bins   = args["-B"].asLong() ?: 32;

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

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
    auto reader = [&](std::string file, bool verbose) {
      auto aos_a   = load_graph<directed>(file);
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<>(file);
        auto hyperedge_degrees = hyperedges.degrees();
        // Run relabeling. This operates directly on the incoming edglist.
        const long idx = args["--relabel"].asLong();
        if (-1 != idx) {
          //relabel the column with smaller size
          if (0 == idx) {
            auto iperm = hyperedges.sort_by_degree(args["--direction"].asString());
            hypernodes.permute(iperm);
            std::cout << "relabeling hyperedge adjacency by degree..." << std::endl;
          }
          else {
            auto iperm = hypernodes.sort_by_degree(args["--direction"].asString());
            hyperedges.permute(iperm);
            std::cout << "relabeling hypernodes adjacency by degree..." << std::endl;
          }
        }
        //may need to update the degree vector, if we relabel the graph
        hyperedge_degrees = hyperedges.degrees();
        return std::tuple(aos_a, hyperedges, hypernodes, hyperedge_degrees);
      }
      std::vector<index_t> hyperedge_degrees =  aos_a.degrees<0>();
      // Run relabeling. This operates directly on the incoming edglist.
      const long idx = args["--relabel"].asLong();
      if (-1 != idx) {
        //relabel the column with smaller size
        std::vector<index_t> degrees;
        if (0 == idx) {
          degrees = hyperedge_degrees;
          std::cout << "relabeling hyperedges by degree..." << std::endl;
          nw::hypergraph::relabel_by_degree_bipartite<0>(aos_a, args["--direction"].asString(), degrees);
        }
        else {
          degrees = aos_a.degrees<1>();
          std::cout << "relabeling hypernodes by degree..." << std::endl;
          nw::hypergraph::relabel_by_degree_bipartite<1>(aos_a, args["--direction"].asString(), degrees);
        }
      }
      //we may need to get the new degrees 
      //if we relabel the edge list
      hyperedge_degrees = aos_a.degrees<0>();
      adjacency<0> hyperedges(aos_a);
      adjacency<1> hypernodes(aos_a);
      if (verbose) {
        hypernodes.stream_stats();
        hyperedges.stream_stats();
      }
      std::cout << "num_hyperedges = " << aos_a.max()[0] + 1 << " num_hypernodes = " << aos_a.max()[1] + 1 << std::endl;
      return std::tuple(aos_a, hyperedges, hypernodes, hyperedge_degrees);
    };

    auto&&[ aos_a, hyperedges, hypernodes, hyperedgedegrees ] = reader(file, verbose);

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        if (verbose) {
          std::cout << "version " << id << std::endl;
        }

        auto verifier = [&](auto&& result) {
          auto&& [N, E] = result;
          if (verbose) {
            //This returns the subgraph of each component.
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
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, std::forward<decltype(op)>(op), verifier, true); };
        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { return baseline(std::execution::seq, aos_a); });
              break;
            case 1:
              record([&] { return baseline(std::execution::par_unseq, aos_a); });
              break;
            case 2:
              record([&] { return lpCC(hypernodes, hyperedges); });
              break;
            case 3:
              record([&] { return lpCCv2(hypernodes, hyperedges); });
              break;
            case 4:
              record([&] { return lpaNoFrontierCC(std::execution::par_unseq, hypernodes, hyperedges); });
              break;
            case 5:
              record([&] { return lpCC_parallelv2(std::execution::par_unseq, hypernodes, hyperedges); });
              break;
            case 8:
              record([&] { return relabelHyperCC(std::execution::seq, aos_a); });
              break;
            case 9:
              record([&] { return relabelHyperCC(std::execution::par_unseq, aos_a); });
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
