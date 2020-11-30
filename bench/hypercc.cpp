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
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
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
      const long idx = args["--relabel"].asLong();
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<>(file);
        // Run relabeling. This operates directly on the incoming edglist.
        if (-1 != idx) {
          nw::hypergraph::relabel_by_degree(hyperedges, hypernodes, idx, args["--direction"].asString());
        }
        std::cout << "num_hyperedges = " << hyperedges.size() << " num_hypernodes = " << hypernodes.size() << std::endl;
        return std::tuple(aos_a, hyperedges, hypernodes);
      }
      else {
        // Run relabeling. This operates directly on the incoming edglist.
        if (-1 != idx) {
          std::cout << "relabeling edge_list by degree..." << std::endl;
          if (1 == idx)
            nw::hypergraph::relabel_by_degree<1>(aos_a, args["--direction"].asString());
          else
            nw::hypergraph::relabel_by_degree<0>(aos_a, args["--direction"].asString());
        }
        adjacency<0> hyperedges(aos_a);
        adjacency<1> hypernodes(aos_a);
        std::cout << "num_hyperedges = " << hyperedges.size() << " num_hypernodes = " << hypernodes.size() << std::endl;
        return std::tuple(aos_a, hyperedges, hypernodes);
      }
    };

    auto&&[ aos_a, hyperedges, hypernodes] = reader(file, verbose);
    auto&& hyperedge_degrees = hyperedges.degrees(std::execution::par_unseq);

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
          std::unordered_set<vertex_id_t> uni_comps(E.begin(), E.end());
          std::cout << uni_comps.size() << " components found" ;
          if (verbose) {
            std::cout << ":";
            for(auto it = uni_comps.begin(); it != uni_comps.end(); ++it)
              std::cout << *it << std::endl;
            for(size_t i = 0, e = E.size(); i < e; ++i)
              std::cout << i << " " << E[i] << std::endl;
          }
          std::cout << std::endl;
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
