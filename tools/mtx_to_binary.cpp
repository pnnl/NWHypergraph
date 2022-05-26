/**
 * @file mtx_to_binary.cpp
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
    R"(mm2bin.exe: convert matrix market file to binary edge list file driver.
  Usage:
      mm2bin.exe (-h | --help)
      mm2bin.exe [-i FILE] [-o FILE] [--unipartite] [-dV]

  Options:
      -h, --help            show this screen
      -i FILE               matrix market input file path
      -o FILE               binary edge list output file path
      --unipartite          unipartite graph or not [default: false]
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";

using namespace nw::hypergraph;

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  bool unipartite = args["--unipartite"].asBool();
  std::string input_file = args["-i"].asString();
  std::string output_file = args["-o"].asString();

  if (unipartite) {
    auto aos_a = load_graph<nw::graph::directedness::undirected>(input_file);
    if (0 == aos_a.size()) {
      std::cerr << input_file << " is not matrix market file, convert abort" << std::endl;
      exit(1);
    }
    aos_a.serialize(output_file);
  }
  else {
    auto aos_a = load_graph(input_file);
    if (0 == aos_a.size()) {
      std::cerr << input_file << " is not matrix market file, convert abort" << std::endl;
      exit(1);
    }
    aos_a.serialize(output_file);
  }
  return 0;
}
