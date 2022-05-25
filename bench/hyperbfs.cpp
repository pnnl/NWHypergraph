/**
 * @file hyperbfs.cpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#include "Log.hpp"
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "s_overlap.hpp"
#include "algorithms/hyper_breadth_first_search.hpp"
#include "algorithms/s_breadth_first_search.hpp"
#include <nwgraph/edge_list.hpp>
#include <nwgraph/util/AtomicBitVector.hpp>
#include <nwgraph/util/intersection_size.hpp>
#include <docopt.h>
#include <execution>

using namespace nw::graph;
using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(hyperbfs.exe: hypergraph breadth-first search benchmark driver.
  Usage:
      hyperbfs.exe (-h | --help)
      hyperbfs.exe [-f FILE...] [-r NODE | -s FILE] [-a NUM] [-b NUM] [-B NUM] [-n NUM] [--seed NUM] [--version ID...] [--succession STR] [--relabel] [--clean] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      -f FILE               input file paths (can have multiples)
      -n NUM                number of trials [default: 1]
      -a NUM                alpha parameter [default: 15]
      -b NUM                beta parameter [default: 18]
      -B NUM                number of bins [default: 32]
      -s, --sources FILE    sources file
      --seed NUM            random seed [default: 27491095]
      -r NODE               start from node r (default is random)
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


int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;
  long alpha   = args["-a"].asLong() ?: 15;
  long beta    = args["-b"].asLong() ?: 18;
  long num_bins= args["-B"].asLong() ?: 32;

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  Times<bool> times;

  // Appease clang.
  //
  // These are captured in lambdas later, and if I use structured bindings
  // they have to be listed as explicit captures (this is according to the 17
  // standard). That's a little bit noisy where it happens, so I just give
  // them real symbols here rather than the local bindings.
  for (auto&& file : files) {
    using vertex_id_t = vertex_id_t<nw::graph::bi_edge_list<directedness::directed>>;
    auto reader = [&](std::string file, bool verbose) {
      auto&& aos_a = load_graph(file);
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<vertex_id_t>(file);
        std::cout << "num_hyperedges = " << hyperedges.size() << " num_hypernodes = " << hypernodes.size() << std::endl;
        return std::tuple(hyperedges, hypernodes);
      }
      // Clean up the edgelist to deal with the normal issues related to
      // undirectedness.
      if (args["--clean"].asBool()) {
        nw::graph::swap_to_triangular<0>(aos_a, args["--succession"].asString());
        nw::graph::lexical_sort_by<0>(aos_a);
        nw::graph::uniq(aos_a);
        nw::graph::remove_self_loops(aos_a);
      }

      nw::graph::biadjacency<0> hyperedges(aos_a);
      nw::graph::biadjacency<1> hypernodes(aos_a);
      if (verbose) {
        hypernodes.stream_stats();
        hyperedges.stream_stats();
      }
      std::cout << "num_hyperedges = " << num_vertices(aos_a, 0) << " num_hypernodes = " << num_vertices(aos_a, 1) << std::endl;
      return std::tuple(hyperedges, hypernodes);
    };

    auto&& [ hyperedges, hypernodes ]     = reader(file, verbose);

    //all sources are hyperedges
    std::vector<vertex_id_t> sources;
    if (args["--sources"]) {
      sources = load_sources_from_file(hyperedges, args["--sources"].asString());
      trials  = sources.size();
    } else if (args["-r"]) {
      sources.resize(trials);
      std::fill(sources.begin(), sources.end(), args["-r"].asLong());
    } else {
      sources = build_random_sources(hyperedges, trials, args["--seed"].asLong());
    }

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        for (auto&& source : sources) {
          if (verbose) std::cout << "version " << id << std::endl;

          auto&& [time, parents] = time_op([&] {
            switch (id) {
              case 0:
                return hyperBFS_topdown_parallel_v0(std::execution::par_unseq, source, hypernodes, hyperedges, num_bins);
              case 1:
                return hyperBFS_bottomup_parallel_v0(std::execution::par_unseq, source, hypernodes, hyperedges, num_bins);
              case 2:
                return hyperBFS_topdown_serial_v0(source, hypernodes, hyperedges);
              case 3:
                return hyperBFS_topdown_serial_v1(source, hypernodes, hyperedges);
              case 4:
                return hyperBFS_bottomup_serial_v0(source, hypernodes, hyperedges);
              default:
                std::cerr << "Unknown version " << id << "\n";
                return std::make_tuple(std::vector<vertex_id_t>(), std::vector<vertex_id_t>());
            }
          });

          if (verify) {
            hyperBFSVerifier(hypernodes, hyperedges, source, std::get<0>(parents), std::get<1>(parents));
          }

          times.append(file, id, thread, time, source);
        }
      }
    }
  }

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("hyperbfs", file, times, header, "Time(s)");
  }

  return 0;
}
