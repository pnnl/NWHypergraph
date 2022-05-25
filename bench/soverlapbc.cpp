/**
 * @file soverlapbc.cpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

static constexpr const char USAGE[] =
 R"(sbc.exe : NWhy soverlap betweenness centrality benchmark driver.
  Usage:
      sbc.exe (-h | --help)
      sbc.exe [-f FILE...] [-r NODE | --sources FILE ] [-i NUM] [-n NUM] [-B NUM] [-s NUM] [--relabel NUM] [--direction DIR] [--adjoin] [--loader-version ID] [--seed NUM] [--version ID...] [--log FILE] [--log-header] [-dvV] [THREADS]...

  Options:
      -h, --help            show this screen
      -f FILE               input file paths (can have multiples and different file format)
      -i NUM                number of iteration [default: 1]
      -n NUM                number of trials [default: 1]
      -r NODE               start from node r (default is random)
      -B NUM                number of bins [default: 32]
      -s NUM                s value of soverlap [default: 1]
      --relabel NUM         relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      --direction DIR       hypergraph relabeling direction - ascending/descending [default: ascending]
      --adjoin              adjoin the id spaces of the hyperedges and hypernodes (smaller one comes after the larger one) 
      --sources FILE        sources file
      --seed NUM            random seed [default: 27491095]
      --version ID          algorithm version to run [default: 5]
      --loader-version ID   soverlap computation loader kernel version [default: 14] 
                            0)Efficient_Blocked 1)Efficient_Cyclic 2)Naive 3)Map_Blocked 4)Map_Cyclic 
                            5)Ensemble_Blocked 6)Ensemble_Cyclic 7)Map_Frontier_Blocked 8)Map_Frontier_Cyclic 
                            9)HashMap_Frontier_Blocked 10)HashMap_Frontier_Cyclic 
                            11)Efficient_Frontier_Blocked 12)Efficient_Frontier_Cyclic
                            13)HashMap_Blocked 14)HashMap_Cyclic 15)Vector_Blocked 16)Vector_Cyclic
                            17)Static_HashMap_Blocked 18)Static_HashMap_Cyclic
                            19)Efficient_Blocked_Size 20)Map_Blocked_Size 21)HashMap_Blocked_Size 22)Vector_Blocked_Size
                            23)Static_HashMap_Blocked_Size
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
      -v, --verify          verify results
      -V, --verbose         run in verbose mode
)";

#include <unordered_set>
#include <docopt.h>
#include <nwgraph/edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "s_overlap.hpp"
#include <nwgraph/algorithms/betweenness_centrality.hpp>
#include <nwgraph/experimental/algorithms/betweenness_centrality.hpp>

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;
using namespace nw::util;

using score_t = float;
using accum_t = double;

int main(int argc, char* argv[]) {
  std::vector strings = std::vector<std::string>(argv + 1, argv + argc);
  std::map       args = docopt::docopt(USAGE, strings, true);

  // Read the options
  bool      verify = args["--verify"].asBool();
  bool     verbose = args["--verbose"].asBool();
  bool       debug = args["--debug"].asBool();
  long      trials = args["-n"].asLong() ?: 1;
  long  iterations = args["-i"].asLong() ?: 1;
  long idx         = args["--relabel"].asLong();
  bool adjoin      = args["--adjoin"].asBool();
  std::string direction 
                   = args["--direction"].asString();
  long num_bins    = args["-B"].asLong() ?: 32;
  size_t s         = args["-s"].asLong() ?: 1;
  long loader_version 
                   = args["--loader-version"].asLong() ?: 0;
 
  std::vector     ids = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times times;

  for (auto&& file : files) {
    size_t nrealedges = 0, nrealnodes = 0;
    nw::graph::edge_list<directedness::directed> e;
    using vertex_id_t = vertex_id_t<decltype(e)>;
    if (adjoin) {
      auto&& [h, ht, iperm] = adjoin_graph_reader<vertex_id_t>(
          file, idx, direction, nrealedges, nrealnodes);
      auto&& degrees = h.degrees(std::execution::par_unseq);
      for (auto&& thread : threads) {
        auto _ = set_n_threads(thread);
        auto features = std::bitset<8>();
        auto graph =
            twograph_reader(loader_version, verbose, features, h, ht, degrees,
                            iperm, nrealedges, nrealnodes, s, thread, num_bins);
        // Source should be selected/generated based on line graph
        std::vector<vertex_id_t> sources;
        if (args["--sources"]) {
          sources =
              load_sources_from_file(h, args["--sources"].asString());
          trials = sources.size();
        } else if (args["-r"]) {
          sources.resize(trials);
          std::fill(sources.begin(), sources.end(), args["-r"].asLong());
        } else {
          sources =
              build_random_sources(graph, trials, args["--seed"].asLong());
        }
        for (auto&& id : ids) {
          for (auto&& source : sources) {
            if (verbose) std::cout << "version " << id << std::endl;
            for (int i = 0; i < trials; ++i) {
              std::vector<vertex_id_t> trial_sources(
                  &sources[iterations * i], &sources[iterations * (i + 1)]);
              auto&& [centrality] =
                  times.record(file, id, thread, [&]() -> std::vector<score_t> {
                    switch (id) {
                      case 0:
                        return bc2_v0<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 1:
                        return bc2_v1<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 2:
                        return bc2_v2<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 3:
                        return bc2_v3<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 4:
                        return bc2_v4<score_t, accum_t>(graph, trial_sources,
                                                        thread);
                      case 5:
                        return brandes_bc<score_t, accum_t>(graph, trial_sources,
                                                        thread);
                      default:
                        std::cerr << "Invalid BC version " << id << "\n";
                        return {};
                    }
                  });
              if (verify)
                BCVerifier<score_t, accum_t>(graph, trial_sources, centrality);
            }  // trial
          }    // source
        }      // id
      }        // for num_thread
    } else {
      auto&& [hyperedges, hypernodes, iperm] =
          graph_reader<vertex_id_t>(file, idx, direction);
      auto&& hyperedge_degrees = hyperedges.degrees(std::execution::par_unseq);
      if (debug) {
        hyperedges.stream_indices();
        hypernodes.stream_indices();

        std::cout << hyperedge_degrees.size() << ": ";
        for (auto d : hyperedge_degrees) std::cout << d << " ";
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
      for (auto&& thread : threads) {
        auto _ = set_n_threads(thread);
        auto features = std::bitset<8>();
        auto graph =
            twograph_reader(loader_version, verbose, features, hyperedges,
                            hypernodes, hyperedge_degrees, iperm, nrealedges,
                            nrealnodes, s, thread, num_bins);
        // Source should be selected/generated based on line graph
        std::vector<vertex_id_t> sources;
        if (args["--sources"]) {
          sources =
              load_sources_from_file(hyperedges, args["--sources"].asString());
          trials = sources.size();
        } else if (args["-r"]) {
          sources.resize(trials);
          std::fill(sources.begin(), sources.end(), args["-r"].asLong());
        } else {
          sources =
              build_random_sources(graph, trials, args["--seed"].asLong());
        }
        for (auto&& id : ids) {
          for (auto&& source : sources) {
            if (verbose) std::cout << "version " << id << std::endl;
            for (int i = 0; i < trials; ++i) {
              std::vector<vertex_id_t> trial_sources(
                  &sources[iterations * i], &sources[iterations * (i + 1)]);
              auto&& [centrality] =
                  times.record(file, id, thread, [&]() -> std::vector<score_t> {
                    switch (id) {
                      case 0:
                        return bc2_v0<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 1:
                        return bc2_v1<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 2:
                        return bc2_v2<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 3:
                        return bc2_v3<decltype(graph), score_t, accum_t>(
                            graph, trial_sources);
                      case 4:
                        return bc2_v4<score_t, accum_t>(graph, trial_sources,
                                                        thread);
                      case 5:
                        return brandes_bc<score_t, accum_t>(graph, trial_sources,
                                                        thread);
                      default:
                        std::cerr << "Invalid BC version " << id << "\n";
                        return {};
                    }
                  });
              if (verify)
                BCVerifier<score_t, accum_t>(graph, trial_sources, centrality);
            }  // trial
          }    // source
        }      // id
      }        // for num_thread

    }  // else
  }    // file
  times.print(std::cout);

  if (args["--log"]) {
    auto   file = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("bc", file, times, header, "Time(s)", "Iterations");
  }
  return 0;
}
