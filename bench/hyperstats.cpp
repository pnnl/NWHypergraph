//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include <docopt.h>
#include <edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(hystats.exe: hypergraph stats driver.
  Usage:
      hystats.exe (-h | --help)
      hystats.exe [-f FILE...] [-a FILE...] [-dV] [--log FILE] [--log-header]

  Options:
      -h, --help            show this screen
      -f FILE               edge list or matrix market input file paths (can have multiples)
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
    auto reader = [&](std::string file, bool verbose) {
      auto aos_a   = load_graph<directed>(file);
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<>(file);
        auto hyperedge_degrees = hyperedges.degrees();
        auto hypernode_degrees = hyperedges.degrees();
        return std::tuple(hyperedges, hypernodes, hyperedge_degrees, hypernode_degrees);
      }
      std::vector<index_t> hyperedge_degrees =  aos_a.degrees<0>();
      std::vector<index_t> hypernode_degrees =  aos_a.degrees<1>();
      // Run relabeling. This operates directly on the incoming edglist.
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
      return std::tuple(hyperedges, hypernodes, hyperedge_degrees, hypernode_degrees);
    };

    auto&&[ hyperedges, hypernodes, hyperedge_degrees, hypernode_degrees ] = reader(file, verbose);

    //hyperedges
    auto M = hyperedges.size();
    size_t hyperE_max_degree = 0;
    size_t hyperE_total_degree = 0;
    std::for_each(hyperedge_degrees.begin(), hyperedge_degrees.end(), [&](auto deg) {
        hyperE_max_degree = std::max(hyperE_max_degree, deg);
        hyperE_total_degree += deg;
    });
    double hyperE_avg_degree = hyperE_total_degree * 1.0 / M;
    //hypernodes
    auto N = hypernodes.size();
    size_t hyperN_max_degree = 0;
    size_t hyperN_total_degree = 0;
    std::for_each(hypernode_degrees.begin(), hypernode_degrees.end(), [&](auto deg) {
        hyperN_max_degree = std::max(hyperN_max_degree, deg);
        hyperN_total_degree += deg;
    });
    double hyperN_avg_degree = hyperN_total_degree * 1.0 / N;
    std::cout << "num_hyperedges = " << M << " num_hypernodes = " << N << std::endl;
    std::cout << "d_edge_max= " << hyperE_max_degree << " d_node_max= " << hyperN_max_degree << std::endl;
    std::cout << "d_edge_total= " << hyperE_total_degree << " d_node_total= " << hyperN_total_degree << std::endl;
    std::cout << "d_edge_avg= " << hyperE_avg_degree << " d_node_avg= " << hyperN_avg_degree << std::endl;
    if (debug)
    {
        hypernodes.stream_indices();
        hyperedges.stream_indices();
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
