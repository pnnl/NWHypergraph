/**
 * @file soverlapsssp.cpp
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
    R"(ssssp.exe: NWhy s-overlap single-source shortest paths benchmark driver.
  Usage:
      ssssp.exe (-h | --help)
      ssssp.exe [-f FILE...] [-r NODE | --sources FILE] [-i NUM] [-n NUM] [-B NUM] [-d NUM] [-s NUM] [--relabel NUM] [--direction DIR]  [--loader-version ID] [--seed NUM] [--version ID...] [--log FILE] [--log-header] [-vV] [--debug] [THREADS]...

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
      --direction DIR         hypergraph relabeling direction - ascending/descending [default: ascending]
      --seed NUM              random seed [default: 27491095]
      --version ID            algorithm version to run [default: 0]
      --loader-version ID     soverlap computation loader kernal version [default: 0]
                              0)Efficient_Blocked 1)Hashmap_Blocked
      --log FILE              log times to a file
      --log-header            add a header to the log file
      --debug                 run in debug mode
      -v, --verify            verify results
      -V, --verbose           run in verbose mode
)";

#include <unordered_set>
#include <docopt.h>
#include <nwgraph/edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "s_overlap.hpp"
#include <nwgraph/algorithms/delta_stepping.hpp>
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
  long idx          = args["--relabel"].asLong();
  std::string direction 
                    = args["--direction"].asString();

  long loader_version = args["--loader-version"].asLong();
  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times times;

  for (auto&& file : files) {
    using vertex_id_t = vertex_id_t<nw::graph::biadjacency<0>>;
    auto&& [hyperedges, hypernodes, iperm] = weighted_graph_reader<vertex_id_t, int>(file, idx, direction);
    auto&& hyperedge_degrees = hyperedges.degrees(std::execution::par_unseq);

    if (verbose) {
      hypernodes.stream_stats();
      hyperedges.stream_stats();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      auto&& graph =
          weighted_twograph_reader<nw::graph::directedness::undirected, int, nw::graph::biadjacency<0, int>,
          nw::graph::biadjacency<1, int>, 
          vertex_id_t>(loader_version, hyperedges,
                          hypernodes, hyperedge_degrees, s, thread, num_bins);
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

          auto&& [time, distance] = time_op([&] {
	      switch (id) {
	      case 0:
	       	return delta_stepping<distance_t>(graph, source, delta); 
		// s_sssp_v0(s_adj, source, delta);
	      default:
		      return delta_stepping<distance_t>(graph, source, delta);
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

