/**
 * @file mmio_to_mmio.cpp
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
#include <nwgraph/edge_list.hpp>
#include <nwgraph/io/mmio.hpp>

static constexpr const char USAGE[] =
    R"(mm2mm.exe: transpose matrix market file driver.
  Usage:
      mm2mm.exe (-h | --help)
      mm2mm.exe [-i FILE] [-o FILE] [--transpose] [-dV]

  Options:
      -h, --help            show this screen
      -i FILE               matrix market input file path
      -o FILE               matrix market output file path
      --transpose           swap column 0 with column 1
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  std::string input_file = args["-i"].asString();
  std::string output_file = args["-o"].asString();
  bool transpose = args["--transpose"].asBool();

  if (!nw::graph::is_mm(input_file)) {
    std::cerr << "not matrix market file, convert abort" << std::endl;
    exit(1);
  }

  nw::graph::directedness direction = nw::graph::get_mm_symmetry(input_file);

  if (nw::graph::directedness::directed == direction) {
    auto&& aos_a = nw::graph::read_mm<nw::graph::directedness::directed>(input_file);
    std::cout << "general matrix" << std::endl;
    aos_a.stream_stats();

    if (!transpose)
      nw::graph::write_mm<0>(output_file, aos_a, "general");
    else
      nw::graph::write_mm<1>(output_file, aos_a, "general");
  } else {
    auto aos_a = nw::graph::read_mm<nw::graph::directedness::undirected>(input_file);
    std::cout << "symmetric matrix" << std::endl;
    aos_a.stream_stats();

    if (!transpose)
      nw::graph::write_mm<0>(output_file, aos_a, "symmetric");
    else
      nw::graph::write_mm<1>(output_file, aos_a, "symmetric");
  }
  return 0;
}
