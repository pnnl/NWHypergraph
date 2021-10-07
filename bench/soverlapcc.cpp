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
#include <containers/edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "s_overlap.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "algorithms/s_connected_components.hpp"


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(scc.exe: s-overlap connected components benchmark driver.
  Usage:
      scc.exe (-h | --help)
      scc.exe [-f FILE...] [--version ID...] [--feature ID...] [--loader-version ID] [-n NUM] [-B NUM] [-s NUM...] [--relabel NUM] [--direction DIR] [--adjoin] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      --loader-version ID   soverlap computation loader kernel version [default: 4] 
                            0)Efficient_Blocked 1)Efficient_Cyclic 2)Naive 3)Map_Blocked 4)Map_Cyclic 
                            5)Ensemble_Blocked 6)Ensemble_Cyclic 7)Map_Frontier_Blocked 8)Map_Frontier_Cyclic 
                            9)HashMap_Frontier_Blocked 10)HashMap_Frontier_Cyclic 11)Efficient_Frontier_Blocked 12)Efficient_Frontier_Cyclic
                            13)HashMap_Blocked 14)HashMap_Cyclic 15)Vector_Blocked 16)Vector_Cyclic
                            17)Static_HashMap_Blocked 18)Static_HashMap_Cyclic
      --feature ID          heuristics in finding soverlap 0)all 1)degree-based pruning 2)skip visited 3)short circuit 4)none [default: 0]
      -f FILE               input file paths (can have multiples and different file format)
      -n NUM                number of trials [default: 1]
      -B NUM                number of bins [default: 32]
      -s NUM                s value of soverlap [default: 1]
      --relabel NUM         relabel the hypergraph - 0(hyperedge)/1(hypernode) [default: -1]
      --direction DIR       hypergraph relabeling direction - ascending/descending [default: ascending]
      --adjoin              adjoin the smaller set among edges and hypernodes with the larger set
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
  long trials  = args["-n"].asLong();
  long num_bins= args["-B"].asLong();
  long loader_version = args["--loader-version"].asLong();
  bool adjoin  = args["--adjoin"].asBool();
  long idx = args["--relabel"].asLong();
  std::string direction = args["--direction"].asString();

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());
  std::bitset features= parse_features(args["--feature"].asStringList());
  std::vector s_values= parse_ids(args["-s"].asStringList());
  if (s_values.empty()) s_values.push_back(1);

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times_WithS<bool> times;

  for (auto&& file : files) {
    size_t nrealedges = 0, nrealnodes = 0;
    auto&& [hyperedges, hypernodes, iperm] = graph_reader<directed>(file, idx, direction, adjoin, nrealedges, nrealnodes);
    auto&& hyperedge_degrees = hyperedges.degrees(std::execution::par_unseq);
    if (debug) {
      hyperedges.stream_indices();
      hypernodes.stream_indices();
      
    
      std::cout << hyperedge_degrees.size() << ": ";
      for (auto d : hyperedge_degrees)
        std::cout << d << " ";
      std::cout << std::endl;
      if (-1 != idx) assert(!iperm.empty());

      if (adjoin) {
        std::cout << "size of the adjoin graph = " << hyperedges.size()
                  << std::endl;
        assert(0 != nrealedges);
        assert(0 != nrealnodes);
      }
    }

    if (verbose) {
      hypernodes.stream_stats();
      hyperedges.stream_stats();
    }

    for (auto&& num_thread : threads) {
      auto _ = set_n_threads(num_thread);
      for (auto&& s : s_values) {
        auto&& s_adj =
          twograph_reader(loader_version, verbose, features, hyperedges,
                          hypernodes, hyperedge_degrees, iperm, nrealedges, nrealnodes, 
                          s, num_thread, num_bins);
        for (auto&& id : ids) {
          auto verifier = [&](auto&& E) {
            // only verify #cc in the result
            std::unordered_map<vertex_id_t, size_t> m;
            for (auto& c : E)
              ++m[c];
            size_t numc = 0;
            if (!debug) {
              for (auto& [k, v] : m) {
                if (1 < v) {
                  ++numc;
                }
              }
            } else {
              std::cout << "Non-singletons: ";
              for (auto& [k, v] : m) {
                if (1 < v) {
                  std::cout << k << " ";
                  ++numc;
                }
              }
              std::cout << std::endl;
            }
            std::unordered_set<vertex_id_t> uni_comps(E.begin(), E.end());
            std::cout << m.size() << " components found" << std::endl;
            std::cout << numc << " non-singleton components found" << std::endl;
          };

          auto record = [&](auto&& op) {
            times.record(file, id, num_thread, s, std::forward<decltype(op)>(op),
                         verifier, true);
          };
          for (int j = 0, e = trials; j < e; ++j) {
            switch (id) {
              case 0:
                record([&] {
                  return linegraph_ccv1(std::execution::par_unseq, hypernodes,
                                        s_adj);
                });
                break;
              case 1:
                record([&] {
                  return linegraph_Afforest(std::execution::par_unseq,
                                            hypernodes, s_adj);
                });
                break;
              default:
                std::cout << "Unknown algorithm version " << id << "\n";
            }
          }
        }
      }  // for s
    }    // for num_thread
  }      // for file

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log_withs("scc", file, times, header, "Time(s)");
  }

  return 0;
}
