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
#include "edge_list_hy.hpp"
#include "s_overlap.hpp"
#include "algorithms/hyper_breadth_first_search.hpp"
#include "algorithms/s_breadth_first_search.hpp"
#include <edge_list.hpp>
#include <util/AtomicBitVector.hpp>
#include <util/intersection_size.hpp>
#include <docopt.h>
#include <tbb/iterators.h>
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
  long iterations = args["-i"].asLong() ?: 1;
  long alpha      = args["-a"].asLong() ?: 15;
  long beta       = args["-b"].asLong() ?: 18;
  long num_bins   = args["-B"].asLong() ?: 32;

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
    auto reader = [&](std::string file, bool verbose) {
      nw::graph::edge_list<nw::graph::directed> aos_a = load_graph<nw::graph::directed>(file);
      std::vector<nw::graph::index_t> hyperedgedegrees = aos_a.degrees<0>();

      // Run relabeling. This operates directly on the incoming edglist.
      if (args["--relabel"].asBool()) {
        //relabel the column with smaller size
        if (aos_a.max()[0] > aos_a.max()[1]) {
          std::vector<nw::graph::index_t> hypernodedegrees = aos_a.degrees<1>();
          std::cout << "relabeling hypernodes..." << std::endl;
                    //TODO NOT WORKING
          nw::hypergraph::relabel_by_degree_bipartite<0>(aos_a, args["--direction"].asString(), hypernodedegrees);
        }
        else {
          std::cout << "relabeling hyperedges..." << std::endl;
          //TODO NOT WORKING
          nw::hypergraph::relabel_by_degree_bipartite<0>(aos_a, args["--direction"].asString(), hyperedgedegrees);
        }
      }
      // Clean up the edgelist to deal with the normal issues related to
      // undirectedness.
      if (args["--clean"].asBool()) {
        aos_a.swap_to_triangular<0>(args["--succession"].asString());
        aos_a.lexical_sort_by<0>();
        aos_a.uniq();
        aos_a.remove_self_loops();
      }

      nw::graph::adjacency<0> hyperedges(aos_a);
      nw::graph::adjacency<1> hypernodes(aos_a);
      if (verbose) {
        hypernodes.stream_stats();
        hyperedges.stream_stats();
      }
      std::cout << "num_hyperedges = " << aos_a.max()[0] + 1 << " num_hypernodes = " << aos_a.max()[1] + 1 << std::endl;
      return std::tuple(aos_a, hyperedges, hypernodes, hyperedgedegrees);
    };

    auto&& [ aos_a, hyperedges, hypernodes, hyperedgedegrees ]     = reader(file, verbose);

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

    edge_list<undirected> s_over;
    nw::graph::adjacency<0> s_adj;
    nw::graph::adjacency<1> s_trans_adj;   
    if (std::find(ids.begin(), ids.end(), 4) != ids.end()) {
      s_over = nw::hypergraph::to_two_graphv6<undirected>(std::execution::seq, hyperedges, hypernodes, hyperedgedegrees);    // 1) find 2-graph corresponding to s-overlapped hyper edges
      s_adj  = build_adjacency<0>(s_over);        // 2) convert new edge_list to new_adjacency
      s_trans_adj = build_adjacency<1>(s_over);
    }

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        for (auto&& source : sources) {
          if (verbose) 
            std::cout << "version " << id << std::endl;

          auto&& [time, parents] = time_op([&] {
            switch (id) {
              case 0:
                return hyperBFS_topdown_serial_v0(source, hypernodes, hyperedges);
              case 1:
                return hyperBFS_topdown_serial_v1(source, hypernodes, hyperedges);
              case 2:
                return hyperBFS_bottomup_serial_v0(source, hypernodes, hyperedges);
              case 3:
                return hyperBFS_hybrid_serial_v0(source, hypernodes, hyperedges, aos_a.size());
              case 4:
                return S_BFS_v0(std::execution::seq, source, hypernodes, s_trans_adj, s_adj, num_bins, alpha, beta);
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
