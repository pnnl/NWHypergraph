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
    R"(ssssp.exe: NWhy s-overlap single-source shortest paths benchmark driver.
  Usage:
      ssssp.exe (-h | --help)
      ssssp.exe [-f FILE...] [-r NODE | --sources FILE] [-i NUM] [-n NUM] [-B NUM] [-d NUM] [-s NUM] [--relabel NUM] [--loader-version ID] [--seed NUM] [--version ID...] [--log FILE] [--log-header] [-vV] [--debug] [THREADS]...

  Options:
      -h, --help              show this screen
      -f FILE                 input file paths (can have multiples and different file format)
      -i NUM                  number of iteration [default: 1]
      -n NUM                  number of trials [default: 1]
      -B NUM                  number of bins [default: 32]
      -s NUM                  s value of soverlap [default: 1]
      -r NODE                 start from node r
      -d, --delta NUM         value for delta [default: 1]
      --sources FILE          sources file
      --relabel NUM           relabel the graph - 0(hyperedge)/1(hypernode) [default: -1]
      --seed NUM              random seed [default: 27491095]
      --version ID            algorithm version to run [default: 0]
      --loader-version ID     soverlap computation loader kernal version [default: 0]
      --log FILE              log times to a file
      --log-header            add a header to the log file
      --debug                 run in debug mode
      -v, --verify            verify results
      -V, --verbose           run in verbose mode
)";

#include <unordered_set>
#include <docopt.h>
#include <edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "s_overlap.hpp"
#include <algorithms/delta_stepping.hpp>
// #include "algorithms/s_shortest_paths.hpp"


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

using distance_t = std::uint64_t;

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  // Read the options
  bool       verify = args["--verify"].asBool();
  bool      verbose = args["--verbose"].asBool();
  bool        debug = args["--debug"].asBool();
  long       trials = args["-n"].asLong(); // at least one trial
  long     num_bins = args["-B"].asLong();
  long   iterations = args["-i"].asLong(); // at least one iteration
  std::size_t delta = args["--delta"].asLong(); // at least one
  std::size_t s     = args["-s"].asLong();

  long loader_version = args["--loader-version"].asLong();
  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times times;

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

    auto twograph_reader = [&](adjacency<0>& edges, adjacency<1>& nodes, std::vector<nw::graph::index_t>& edgedegrees, 
    size_t s = 1, int num_bins = 32) {
      switch (loader_version) {
      case 0:
      {
          nw::graph::edge_list<undirected, int> &&linegraph = 
          to_two_graph_weighted_efficient_parallel_clean<undirected, int>
          (std::execution::par_unseq, hyperedges, hypernodes, edgedegrees, s, num_bins);
          //where when an empty edge list is passed in, an adjacency still have two elements
          if (0 == linegraph.size()) return nw::graph::adjacency<0, int>(0);
          nw::graph::adjacency<0, int> s_adj(linegraph);
          std::cout << "line:" << linegraph.size() << " adjacency: " << s_adj.size() << std::endl;
          return s_adj;
      }
      default:
      {
          std::cerr << "unknown soverlap computation loader" << std::endl;
          return nw::graph::adjacency<0, int>(0);
      }
      }
    };
    auto&& s_adj = twograph_reader(hyperedges, hypernodes, hyperedge_degrees, s, num_bins);

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    //Source should be selected/generated based on line graph
    std::vector<vertex_id_t> sources;
    if (args["--sources"]) {
      sources = load_sources_from_file(hyperedges, args["--sources"].asString());
      trials  = sources.size();
    } else if (args["-r"]) {
      sources.resize(trials);
      std::fill(sources.begin(), sources.end(), args["-r"].asLong());
    } else {
      sources = build_random_sources(s_adj, trials, args["--seed"].asLong()); // TODO: set source based on line graph vertex max ID
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        for (auto&& source : sources) {
          if (verbose) std::cout << "version " << id << std::endl;

          auto&& [time, distance] = time_op([&] {
	      switch (id) {
	      case 0:
	       	return delta_stepping_v10<distance_t>(s_adj, source, delta); 
		// s_sssp_v0(s_adj, source, delta);
	      default:
		return delta_stepping_v10<distance_t>(s_adj, source, delta);
                // std::cerr << "Unknown version " << id << "\n";
                // return std::make_tuple(std::vector<std::atomic<long unsigned int>, std::allocator<std::atomic<long unsigned int> > >());
		// return true;
		// return std::tuple(0.0);
		// return std::make_tuple(std::vector<vertex_id_t>());
	      }
	    });
	  times.append(file, id, thread, time);
        }
      }
    }
  }
  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("sssp", file, times, header, "Time(s)");
  }
  return 0;
}

