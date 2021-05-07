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
      scc.exe [-f FILE...] [--version ID...] [--feature ID...] [--loader-version ID] [-n NUM] [-B NUM] [-s NUM...] [--relabel NUM] [--direction DIR] [--adjoin NUM] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      --loader-version ID   soverlap computation loader kernal version [default: 4]
      --feature ID          heuristics in finding soverlap 0)all 1)degree-based pruning 2)skip visited 3)short circuit 4)none [default: 0]
      -f FILE               input file paths (can have multiples and different file format)
      -n NUM                number of trials [default: 1]
      -B NUM                number of bins [default: 32]
      -s NUM                s value of soverlap [default: 1]
      --relabel NUM         relabel the hypergraph - 0(hyperedge)/1(hypernode) [default: -1]
      --direction DIR       hypergraph relabeling direction - ascending/descending [default: ascending]
      --adjoin NUM          adjoin hyperedge(0) with hypernode(1) [default: -1]
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
  long adjoin  = args["--adjoin"].asLong();

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

  // Appease clang.
  //
  // These are captured in lambdas later, and if I use structured bindings
  // they have to be listed as explicit captures (this is according to the 17
  // standard). That's a little bit noisy where it happens, so I just give
  // them real symbols here rather than the local bindings.
  for (auto&& file : files) {
    auto reader = [&](std::string file) {
      auto aos_a = load_graph<directed>(file);
      const long idx = args["--relabel"].asLong();
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<>(file);
        // Run relabeling. This operates directly on the incoming edglist.
        if (-1 != idx) {
          nw::hypergraph::relabel_by_degree(hyperedges, hypernodes, idx, args["--direction"].asString());
        }
        std::cout << "num_hyperedges = " << hyperedges.size() << " num_hypernodes = " << hypernodes.size() << std::endl;
        return std::tuple(hyperedges, hypernodes);
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
        return std::tuple(hyperedges, hypernodes);
      }
    };
    auto&&[ hyperedges, hypernodes] = reader(file);
    auto&& hyperedge_degrees = hyperedges.degrees(std::execution::par_unseq);

    if (verbose) {
      hypernodes.stream_stats();
      hyperedges.stream_stats();
    }

    for (auto&& s : s_values) {
    auto&& s_adj = twograph_reader(loader_version, verbose, features, hyperedges, hypernodes, hyperedge_degrees, s, num_bins);

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        auto verifier = [&](auto&& E) {
            //only verify #cc in the result
          std::unordered_map<vertex_id_t, size_t> m;
          for (auto& c : E) {
            ++m[c];
          }
          size_t numc = 0;
          for (auto& [k, v] : m) {
            if (1 < v)
              ++numc;
          }
          std::unordered_set<vertex_id_t> uni_comps(E.begin(), E.end());
          std::cout << m.size() << " components found" << std::endl;
          std::cout << numc << " non-singleton components found" << std::endl;
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, s, std::forward<decltype(op)>(op), verifier, true); };
        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { return linegraph_ccv1(std::execution::par_unseq, hypernodes, s_adj); });
              break;
            case 1:
              record([&] { return linegraph_Afforest(std::execution::par_unseq, hypernodes, s_adj); });
              break;
            default:
              std::cout << "Unknown algorithm version " << id << "\n";
          }
        }
      }
    }//for thread
    }//for s
  }//for file

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log_withs("scc", file, times, header, "Time(s)");
  }

  return 0;
}
