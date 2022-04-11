//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2022 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include <map>

#include <docopt.h>
#include "io/loader.hpp"
#include <nwgraph/edge_list.hpp>
#include <nwgraph/io/mmio.hpp>

#include "xtensor/xcsv.hpp"

using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(imdb2mtx.exe: create mtx file from imdb tsv files.
  Usage:
      imdb2mtx.exe (-h | --help)
      imdb2mtx.exe [--title FILE] [--name FILE] [--principal FILE] [-o FILE] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --title FILE          movie title file path
      --name FILE           actor name file path
      --principal FILE      movie to actor file path
      -o FILE               matrix market output file path
      -d, --debug           run in debug mode
      -V, --verbose         run in verbose mode
)";

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  std::string output_file = args["-o"].asString();
  
  nw::graph::bi_edge_list<nw::graph::directedness::directed> edges(0);
  using vertex_id_t = vertex_id_t<decltype(edges)>;

  std::string title_basics_tsv = args["--title"].asString();
  std::string name_basics_tsv = args["--name"].asString();
  std::string title_principals_tsv = args["--principal"].asString();

  std::ifstream                 title_basics_stream(title_basics_tsv);
  auto                          titles     = xt::load_csv<std::string>(title_basics_stream, '\t');
  auto                          titles_shp = titles.shape();
  std::map<std::string, vertex_id_t> titles_map;
  std::map<vertex_id_t, std::string> titles_map_transpose;
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < titles_shp[0]; ++i) {
    if (titles(i, 1) == "movie") {
      titles_map[titles(i, 0)] = i;
      titles_map_transpose[i] = titles(i, 0);
    }
  }

  std::ifstream                 name_basics_stream(name_basics_tsv);
  auto                          names     = xt::load_csv<std::string>(name_basics_stream, '\t');
  auto                          names_shp = names.shape();
  std::map<std::string, vertex_id_t> names_map;
  // this vector store the PrimaryName of the actors
  std::vector<std::string> names_map_transpose(names_shp[0]);
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < names_shp[0]; ++i) {
    names_map[names(i, 0)] = i;
    names_map_transpose[i] = names(i, 1);
  }

  std::ifstream title_principals_stream(title_principals_tsv);
  auto          title_principals = xt::load_csv<std::string>(title_principals_stream, '\t');
  auto          shp              = title_principals.shape();

  //nw::graph::edge_list<nw::graph::directedness::directed> edges(0);
  edges.open_for_push_back();
  //skip the header by starting from i=1
  for (vertex_id_t i = 1; i < shp[0]; ++i) {
    if (title_principals(i, 3) == "actor" || title_principals(i, 3) == "actress") {
      auto title = title_principals(i, 0);
      auto name  = title_principals(i, 2);
      //if the title does not show, then it is not a movie.
      auto it_title = titles_map.find(title);
      if (it_title == titles_map.end()) {
        continue;
      }
      //same if the person is not an actor
      auto it_name = names_map.find(name);
      if (it_name == names_map.end()) {
        continue;
      }
      edges.push_back(it_title->second, it_name->second);
    }
  }
  edges.close_for_push_back();
  edges.stream_stats();


  auto G = nw::graph::biadjacency<0>(edges);
  auto H = nw::graph::biadjacency<1>(edges);

  write_mm<0>(output_file, G, "general");

  return 0;
}