/**
 * @file toplexes.cpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#include <docopt.h>
#include <nwgraph/edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "algorithms/toplexes.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(top.exe: hypergraph toplexes benchmark driver.
  Usage:
      top.exe (-h | --help)
      top.exe [-f FILE...] [-a FILE...] [--version ID...] [-B NUM] [-n NUM] [--direction DIR] [--relabel NUM] [-cdV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      -f FILE               edge list or matrix market input file paths (can have multiples)
      -a FILE               hypergraph adjacency fils paths (can have multiples)
      -n NUM                number of trials [default: 1]
      -B NUM                number of bins [default: 32]
      --relabel NUM         relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      -c, --clean           clean the graph or not
      --direction DIR       graph relabeling direction - ascending/descending [default: ascending]
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
    using vertex_id_t = uint32_t;
    auto&&[hyperedges, hypernodes, iperm] = graph_reader<vertex_id_t>(file, args["--relabel"].asLong(), args["--direction"].asString());
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
          std::cout << result.size() << std::endl;
          std::for_each(result.begin(), result.end(), [&](auto& e) {
              
          });
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, std::forward<decltype(op)>(op), verifier, true); };
        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { return toplexes_serial_v0(hyperedges); });
              break;
            case 1:
              record([&] { return toplexes_serial_v1(hyperedges); });
              break;
            case 2:
              record([&] { return toplexes_serial_v2(hyperedges); });
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
    log("top", file, times, header, "Time(s)");
  }

  return 0;
}
