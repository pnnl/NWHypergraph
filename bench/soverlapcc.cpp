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
#include "algorithms/s_connected_components.hpp"


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(scc.exe: s-overlap connected components benchmark driver.
  Usage:
      scc.exe (-h | --help)
      scc.exe [-f FILE...] [--version ID...] [--feature ID...] [--loader-version ID] [-n NUM] [-B NUM] [-s NUM...] [--relabel NUM] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      --loader-version ID   soverlap computation loader kernal version [default: 0]
      --feature ID          heuristics in finding soverlap 0)all 1)degree-based pruning 2)skip visited 3)short circuit 4)none [default: 0]
      -f FILE               input file paths (can have multiples and different file format)
      -n NUM                number of trials [default: 1]
      -B NUM                number of bins [default: 32]
      -s NUM                s value of soverlap [default: 1]
      --relabel NUM         relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
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
  long num_bins= args["-B"].asLong() ?: 32;
  long loader_version = args["--loader-version"].asLong() ?: 0;

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
        return std::tuple(hyperedges, hypernodes, hyperedge_degrees);
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
      std::cout << "num_hyperedges = " << hyperedges.size() << " num_hypernodes = " << hypernodes.size() << std::endl;
      return std::tuple(hyperedges, hypernodes, hyperedge_degrees);
    };
    auto&&[ hyperedges, hypernodes, hyperedgedegrees ] = reader(file, verbose);
    for (auto&& s : s_values) {
    auto twograph_reader = [&](adjacency<0>& edges, adjacency<1>& nodes, std::vector<nw::graph::index_t>& edgedegrees, 
    size_t s = 1, int num_bins = 32) {
      switch (loader_version) {
      case 0:
      {
          nw::graph::edge_list<undirected> &&linegraph = to_two_graph_efficient_parallel_portal<undirected>(verbose, features, std::execution::par_unseq, hyperedges, hypernodes, hyperedgedegrees, s, num_bins);
          //where when an empty edge list is passed in, an adjacency still have two elements
          if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
          nw::graph::adjacency<0> s_adj(linegraph);
          std::cout << "line graph edges = " << linegraph.size() << ", adjacency size = " << s_adj.size() << ", max= " << s_adj.max() << std::endl;
          return s_adj;
      }
      case 1:
      {
          nw::graph::edge_list<undirected> &&linegraph = to_two_graph_naive_parallel_portal<undirected>(verbose, std::execution::par_unseq, hyperedges, hypernodes, s, num_bins);
          //where when an empty edge list is passed in, an adjacency still have two elements
          if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
          nw::graph::adjacency<0> s_adj(linegraph);
          std::cout << "line graph edges = " << linegraph.size() << ", adjacency size = " << s_adj.size() << std::endl;
          return s_adj;
      }
      default:
      {
          std::cerr << "unknown soverlap computation loader" << std::endl;
          return nw::graph::adjacency<0>(0);
      }
      }
    };
    auto&& s_adj = twograph_reader(hyperedges, hypernodes, hyperedgedegrees, s, num_bins);

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        auto verifier = [&](auto&& E) {
            //only verify #cc in the result    
          std::unordered_set<vertex_id_t> uni_comps(E.begin(), E.end());
          std::cout << uni_comps.size() << " components found" << std::endl;
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
