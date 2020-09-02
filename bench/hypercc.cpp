//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

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
      hycc.exe [-f FILE...] [--version ID...] [-n NUM] [--succession STR] [--relabel] [--clean] [--direction DIR] [-dvV] [--log FILE] [--log-header] [--overlap NUM] [THREADS]...

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
      --overlap NUM         s overlap [default: 1]
      -d, --debug           run in debug mode
      -v, --verify          verify results
      -V, --verbose         run in verbose mode
)";



int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;

  size_t s_overlap = args["--overlap"].asLong();

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
      auto hyperedgedegrees = aos_a.degrees<0>();

      // Run relabeling. This operates directly on the incoming edglist.
      if (args["--relabel"].asBool()) {
        //relabel the column with smaller size
        if (aos_a.max()[0] > aos_a.max()[1]) {
          auto hypernodedegrees = aos_a.degrees<1>();
          std::cout << "relabeling hypernodes..." << std::endl;
          //TODO NOT WORKING
          nw::hypergraph::relabel_by_degree_bipartite<1>(aos_a, args["--direction"].asString(), hypernodedegrees);
        }
        else {
          std::cout << "relabeling hyperedges..." << std::endl;
          //TODO NOT WORKING
          nw::hypergraph::relabel_by_degree_bipartite<1>(aos_a, args["--direction"].asString(), hyperedgedegrees);
        }
      }
      // Clean up the edgelist to deal with the normal issues related to
      // undirectedness.
      if (args["--clean"].asBool()) {
        aos_a.swap_to_triangular<0>(args["--succession"].asString());
        aos_a.lexical_sort_by<0>();
        aos_a.uniq();
        aos_a.remove_self_loops();
      }

      adjacency<0> hyperedges(aos_a);
      adjacency<1> hypernodes(aos_a);
      if (verbose) {
        hypernodes.stream_stats();
        hyperedges.stream_stats();
      }
      std::cout << "num_hyperedges = " << aos_a.max()[0] + 1 << " num_hypernodes = " << aos_a.max()[1] + 1 << std::endl;
      return std::tuple(aos_a, hyperedges, hypernodes, hyperedgedegrees);
    };

    // auto&& graphs     = reader(file, verbose);
    // auto&& aos_a      = std::get<0>(graphs);
    // auto&& hyperedges = std::get<1>(graphs);
    // auto&& hypernodes = std::get<2>(graphs);
    // auto&& hyperedgedegrees = std::get<3>(graphs);

    auto&&[ aos_a, hyperedges, hypernodes, hyperedgedegrees ] = reader(file, verbose);
    auto twograph_reader = [&](adjacency<0>& edges, adjacency<1>& nodes, std::vector<nw::graph::index_t>& edgedegrees, size_t s = 1) {
      nw::util::life_timer _("build adj line graph");
      if (ids.end() != std::find(ids.begin(), ids.end(), 6)) {
        //create line graph only when needed by the algorithm
        auto&& linegraph =  to_two_graphv2<undirected>(std::execution::par_unseq, hyperedges, hypernodes, hyperedgedegrees, s_overlap);
        nw::graph::adjacency<0> s_adj(linegraph);
        return s_adj;
      }
      return nw::graph::adjacency<0>(0);
    };
    auto&& s_adj = twograph_reader(hyperedges, hypernodes, hyperedgedegrees, s_overlap);

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
          //TODO
          if (verbose) {
            //print_top_n(graph, comp);
          }
          std::map<vertex_id_t, edge_list<>> comps;
          auto N = std::get<0>(result);
          std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
            auto&& [edge, node] = elt;
            vertex_id_t key     = N[node];
            comps[key].push_back(elt);
          });

          for (auto&& j : comps) {
            auto& [k, v] = j;
            v.close_for_push_back();
          }
          std::cout << comps.size() << " components found" << std::endl;

          if (verify) {
            std::cerr << " v" << id << " failed verification for " << file << " using " << thread << " threads\n";
          }
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
              record([&] { return svCC(std::execution::seq, aos_a, hypernodes, hyperedges); });
              break;
            case 4:
              record([&] { return bfsCC(std::execution::par_unseq, hypernodes, hyperedges); });
              break;
            case 5:
              record([&] { return lpCC(std::execution::par_unseq, hypernodes, hyperedges); });
              break;
            case 6:
              record([&] { return base_two(std::execution::seq, hypernodes, s_adj); });
              break;
            case 7:
              record([&] { return relabelHyperCC(std::execution::seq, aos_a); });
              break;
            case 8:
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
