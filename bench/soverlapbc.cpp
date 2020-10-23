//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Jesun Firoz, Xu Tony Liu
//

static constexpr const char USAGE[] =
 R"(sbc.exe : NWhy soverlap betweenness centrality benchmark driver.
  Usage:
      sbc.exe (-h | --help)
      sbc.exe [-f FILE...] [-r NODE | --sources FILE ] [-i NUM] [-n NUM] [-B NUM] [-s NUM] [--relabel NUM] [--loader-version ID] [--seed NUM] [--version ID...] [--log FILE] [--log-header] [-dvV] [THREADS]...

  Options:
      -h, --help              show this screen
      -f FILE                 input file paths (can have multiples and different file format)
      -i NUM                  number of iteration [default: 1]
      -n NUM                  number of trials [default: 1]
      -r NODE                 start from node r (default is random)
      -B NUM                  number of bins [default: 32]
      -s NUM                  s value of soverlap [default: 1]
      --relabel NUM           relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      --sources FILE          sources file
      --seed NUM              random seed [default: 27491095]
      --version ID            algorithm version to run [default: 0]
      --loader-version ID     soverlap computation loader kernal version [default: 0]
      --log FILE              log times to a file
      --log-header            add a header to the log file
      -d, --debug             run in debug mode
      -v, --verify            verify results
      -V, --verbose           run in verbose mode
)";

#include <unordered_set>
#include <docopt.h>
#include <edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "edge_list_hy.hpp"
#include "s_overlap.hpp"
#include "algorithms/betweenness_centrality.hpp"


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
  long num_bins= args["-B"].asLong() ?: 32;
  size_t s     = args["-s"].asLong() ?: 1;
  long loader_version = args["--loader-version"].asLong() ?: 0;
 
  std::vector     ids = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times times;

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
        std::cout << "num_hyperedges = " << hyperedges.size() << " num_hypernodes = " << hypernodes.size() << std::endl;
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

    auto twograph_reader = [&](adjacency<0>& edges, adjacency<1>& nodes, std::vector<nw::graph::index_t>& edgedegrees, 
    size_t s = 1, int num_bins = 32) {
      switch (loader_version) {
      case 0:
      {
          nw::graph::edge_list<undirected> &&linegraph = to_two_graph_efficient_parallel_clean<undirected>(std::execution::par_unseq, hyperedges, hypernodes, hyperedgedegrees, s, num_bins);
          //where when an empty edge list is passed in, an adjacency still have two elements
          if (0 == linegraph.size()) return nw::graph::adjacency<1>(0);
          nw::graph::adjacency<1> s_adj(linegraph);
          std::cout << "line:" << linegraph.size() << " adjacency: " << s_adj.size() << std::endl;
          return s_adj;
      }
      default:
      {
          std::cerr << "unknown soverlap computation loader" << std::endl;
          return nw::graph::adjacency<1>(0);
      }
      }
    };
    auto&& graph = twograph_reader(hyperedges, hypernodes, hyperedgedegrees, s, num_bins);

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        for (auto&& source : sources) {
          if (verbose) std::cout << "version " << id << std::endl;
          for (int i = 0; i < trials; ++i) {
            std::vector<vertex_id_t> trial_sources(&sources[iterations * i], &sources[iterations * (i + 1)]);
            auto&& [centrality] = times.record(file, id, thread, [&]() -> std::vector<score_t> {
            switch (id)
            {
             case 0: return bc2_v0<decltype(graph), score_t, accum_t>(graph, trial_sources);
             case 1: return bc2_v1<decltype(graph), score_t, accum_t>(graph, trial_sources);
             case 2: return bc2_v2<decltype(graph), score_t, accum_t>(graph, trial_sources);
             case 3: return bc2_v3<decltype(graph), score_t, accum_t>(graph, trial_sources);
             case 4: return bc2_v4<score_t, accum_t>(graph, trial_sources, thread);
             case 5: return bc2_v5<score_t, accum_t>(graph, trial_sources, thread);
             default:
              std::cerr << "Invalid BC version " << id << "\n";
              return {};
            }
            });
            if (verify) BCVerifier<score_t, accum_t>(graph, trial_sources, centrality);
          } // trial
        } //source
      } //id
    }//thread
  }//file
  times.print(std::cout);

  if (args["--log"]) {
    auto   file = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("bc", file, times, header, "Time(s)", "Iterations");
  }
  return 0;
}
