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
#include "common.hpp"


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

using namespace nw::hypergraph::tools;
using namespace nw::hypergraph;

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  std::string input_file = args["-i"].asString();
  std::string output_file = args["-o"].asString();
  bool transpose = args["--transpose"].asBool();

  auto aos_a   = load_graph<directed>(input_file);
  if (0 == aos_a.size()) {
    std::cerr << "not matrix market file, convert abort" << std::endl;
    exit(1);
  }

  if (!transpose)
    write_mm<0>(output_file, aos_a);
  else
    write_mm<1>(output_file, aos_a);

  return 0;
}
