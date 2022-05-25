/**
 * @file hyadj_to_mmio.cpp
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
/*
* The file format is at
* https://github.com/jshun/ligra/tree/ligra-h#input-format-for-ligra-h-applications
**/

static constexpr const char USAGE[] =
    R"(adj2mm.exe: convert hypergraph biadjacency file to matrix market file driver.
  Usage:
      adj2mm.exe (-h | --help)
      adj2mm.exe [-i FILE] [-o FILE] [--transpose] [-dV]

  Options:
      -h, --help            show this screen
      -i FILE               hypergraph biadajacency input file path
      -o FILE               matrix market output file path
      --transpose           swap column 0 with column 1
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
  bool transpose = args["--transpose"].asBool();


  auto&& [hyperedges, hypernodes] = load_adjacency<vertex_id_t<nw::graph::biadjacency<0>>>(input_file);
  if (!transpose)
    write_mm<0>(output_file, hyperedges, "general");
  else
    write_mm<0>(output_file, hypernodes, "general");
  


  return 0;
}
