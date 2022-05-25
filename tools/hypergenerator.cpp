/**
 * @file hypergenerator.cpp
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
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "generators/configuration_model.hpp"
#include "io/loader.hpp"

using namespace nw::graph;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(hygen.exe: random unweighted hypergraph generator driver
  Usage:
      hygen.exe (-h | --help)
      hygen.exe [-f FILE] [--output FILE] [--deg-seqa FILE] [--deg-seqb FILE]

  Options:
      -h, --help            show this screen
      -f FILE               input matrix market file to extract degree sequences
      --output FILE         output random hypergraph as matrix market format
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
    std::string deg_seqb = args["--degs-nodes"].asString();
    using vertex_id_t = vertex_id_t<nw::graph::bi_edge_list<nw::graph::directedness::directed>>;
    auto reader = [&](std::string file) {
      auto aos_a = load_graph(file);
       std::vector<vertex_id_t> hyperedge_degrees, hypernode_degrees;
      if (0 == aos_a.size()) {
        auto&& [hyperedges, hypernodes] = load_adjacency<vertex_id_t>(file);
        hyperedge_degrees = hyperedges.degrees();
        hypernode_degrees = hypernodes.degrees();
        return std::tuple(hyperedge_degrees, hypernode_degrees);
      }
      else {
        hyperedge_degrees = degrees<0>(aos_a);
        hypernode_degrees = degrees<1>(aos_a);
        return std::tuple(hyperedge_degrees, hypernode_degrees);
      }
    };

    nw::graph::bi_edge_list<nw::graph::directedness::directed> edges;
    if ("" != input) {
        auto&&[hyperedge_degrees, hypernode_degrees] = reader(input);
        edges = configuration_model(hyperedge_degrees, hypernode_degrees, false);
    }

    if ("" != output) {
      std::cout << "Writing mtx file" << std::endl;
      size_t num_hyperedges = num_vertices(edges, 0), num_hypernodes = num_vertices(edges, 1);
      nw::graph::biadjacency<0> hyperedges(edges);
      write_mm_hy(output, hyperedges, num_hyperedges, num_hypernodes);
    }

  return 0;
}
