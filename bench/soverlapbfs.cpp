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
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "s_overlap.hpp"
#include "algorithms/hyper_breadth_first_search.hpp"
#include "algorithms/s_breadth_first_search.hpp"
#include <containers/edge_list.hpp>
#include <util/AtomicBitVector.hpp>
#include <util/intersection_size.hpp>
#include <docopt.h>
#include <execution>


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(scc.exe: s-overlap breadth-first search benchmark driver.
  Usage:
      sbfs.exe (-h | --help)
      sbfs.exe [-f FILE...] [--version ID...] [--loader-version ID] [-n NUM] [-a NUM] [-b NUM] [-B NUM] [-s NUM] [--seed NUM] [--relabel NUM] [--direction DIR] [--adjoin] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      --loader-version ID   soverlap computation loader kernel version [default: 14] 
                            0)Efficient_Blocked 1)Efficient_Cyclic 2)Naive 3)Map_Blocked 4)Map_Cyclic 
                            5)Ensemble_Blocked 6)Ensemble_Cyclic 7)Map_Frontier_Blocked 8)Map_Frontier_Cyclic 
                            9)HashMap_Frontier_Blocked 10)HashMap_Frontier_Cyclic 
                            11)Efficient_Frontier_Blocked 12)Efficient_Frontier_Cyclic
                            13)HashMap_Blocked 14)HashMap_Cyclic 15)Vector_Blocked 16)Vector_Cyclic
                            17)Static_HashMap_Blocked 18)Static_HashMap_Cyclic
                            19)Efficient_Blocked_Size 20)Map_Blocked_Size 21)HashMap_Blocked_Size 22)Vector_Blocked_Size
                            23)Static_HashMap_Blocked_Size
      -f FILE               input file paths (can have multiples and different file format)
      -n NUM                number of trials [default: 1]
      -B NUM                number of bins [default: 32]
      -s NUM                s value of s-overlap [default: 1]
      -a NUM                alpha parameter [default: 15]
      -b NUM                beta parameter [default: 18]
      --sources FILE        sources file
      --seed NUM            random seed [default: 27491095]
      -r NODE               start from node r (default is random)
      --relabel NUM         relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      --adjoin              adjoin the id spaces of the hyperedges and hypernodes (smaller one comes after the larger one) 
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
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
  size_t s     = args["-s"].asLong() ?: 1;
  long idx     = args["--relabel"].asLong();
  bool adjoin  = args["--adjoin"].asBool();
  long loader_version 
               = args["--loader-version"].asLong() ?: 0;
  std::string direction 
               = args["--direction"].asString();

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
    size_t nrealedges = 0, nrealnodes = 0;
    auto&& [hyperedges, hypernodes, iperm] = graph_reader<directed>(file, idx, direction, adjoin, nrealedges, nrealnodes);
    auto&& hyperedge_degrees = hyperedges.degrees(std::execution::par_unseq);
    if (verbose) {
      hypernodes.stream_stats();
      hyperedges.stream_stats();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      auto features = std::bitset<8>();
      auto&& graph =
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
        sources = build_random_sources(
            graph, trials,
            args["--seed"].asLong());      
      }
      for (auto&& id : ids) {
        for (auto&& source : sources) {
          if (verbose) std::cout << "version " << id << std::endl;

          auto&& [time, parents] = time_op([&] {
            switch (id) {
              case 0:
                return bfs_v0(graph, source);
              default:
                std::cerr << "Unknown version " << id << "\n";
                return std::vector<vertex_id_t>();
            }
          });

          if (verify) {
            //TODO
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
    log("cc", file, times, header, "Time(s)");
  }

  return 0;
}
