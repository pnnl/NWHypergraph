/**
 * @file mmio_to_hyadj.cpp
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
#include "io/loader.hpp"

static constexpr const char USAGE[] =
    R"(mm2adj.exe: convert matrix market file to hypergraph biadjacency file driver.
  Usage:
      mm2adj.exe (-h | --help)
      mm2adj.exe [-i FILE] [-o FILE] [-dV]

  Options:
      -h, --help            show this screen
      -i FILE               matrix market input file path
      -o FILE               hypergraph biadajacency output file path
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";

using namespace nw::hypergraph;

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  std::string input_file = args["-i"].asString();
  std::string output_file = args["-o"].asString();

  auto aos_a = load_graph(input_file);
  if (0 == aos_a.size()) {
    std::cerr << input_file << " is not matrix market file, convert abort" << std::endl;
    exit(1);
  }
  nw::graph::biadjacency<0> E(aos_a);
  nw::graph::biadjacency<1> N(aos_a);
  std::cout << "num_hyperedges:" << E.size() << " num_hypernodes:" << N.size() << std::endl;
  if (false == write_adj_hypergraph(output_file, E, N)){
    std::cerr << "write failed" << std::endl;
  }
  return 0;
}
