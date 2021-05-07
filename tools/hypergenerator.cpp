//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include <fstream>
#include <docopt.h>
#include "common.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "generators/configuration_model.hpp"
#include "io/mmio.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(hygen.exe: random hypergraph generator driver.
  Usage:
      hygen.exe (-h | --help)
      hygen.exe [-f FILE] [--output FILE] [--deg-seqa FILE] [--deg-seqb FILE]

  Options:
      -h, --help            show this screen
      -f FILE               input matrix market file to extract degree sequences
      --output FILE         output random graph as matrix market format
      --degs-edges FILE     degree sequence of edges in FILE
      --degs-nodes FILE     degree sequence of nodes in FILE
)";


int main(int argc, char* argv[]) {
    std::vector<std::string> strings(argv + 1, argv + argc);
    auto args = docopt::docopt(USAGE, strings, true);

    bool debug = args["--debug"].asBool();
    std::string input = args["-f"].asString();
    std::string output = args["--output"].asString();
    std::string deg_seqa = args["--degs-edges"].asString();
    std::string deg_seqa = args["--degs-nodes"].asString();

    auto reader = [&](std::string file) {
      auto aos_a = load_graph<nw::graph::directed>(file);
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<>(file);
        auto hyperedge_degrees = hyperedges.degrees();
        auto hypernode_degrees = hypernodes.degrees();
        return std::tuple(hyperedge_degrees, hypernode_degrees);
      }
      else {
        std::vector<index_t> hyperedge_degrees =  aos_a.degrees<0>();
        std::vector<index_t> hypernode_degrees =  aos_a.degrees<1>();
        return std::tuple(hyperedge_degrees, hypernode_degrees);
      }
    };
    std::vector hyperedge_degrees, hypernode_degrees;
    if ("" != input) {
        auto&&[hyperedge_degrees, hypernode_degrees] = reader(input);
    }

    if ("" != output) {
      std::cout << "Writing mtx file" << std::endl;
      if (0 == idx)
        write_mm_hy(output, hyperedges, M, N);
      else if (1 == idx)
        write_mm_hy(output, hypernodes, N, M);
      else
        std::cerr << "unrecognized flag" << std::endl;
    }

  return 0;
}
