//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Andrew Lumsdaine, Xu Tony Liu
//

#include <unordered_set>
#include <docopt.h>
#include <edge_list.hpp>
#include "Log.hpp"
#include "common.hpp"
#include "s_overlap.hpp"
#include "containers/edge_list_hy.hpp"
#include "containers/compressed_hy.hpp"
#include "algorithms/s_connected_components.hpp"


#include <fstream>
#include <iostream>

#include <deque>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "xtensor/xcsv.hpp"

#include "bfs_edge_range.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;

static constexpr const char USAGE[] =
    R"(imdb.exe: imdb benchmark driver.
  Usage:
      imdb.exe (-h | --help)
      imdb.exe [--title FILE] [--name FILE] [--principal FILE] [--version ID...] [--loader-version ID] [-B NUM] [-s NUM...] [--relabel NUM] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      --loader-version ID   soverlap computation loader kernal version [default: 4]
      --title FILE          movie title file path
      --name FILE           actor name file path
      --principal FILE      movie to actor file path
      -B NUM                number of bins [default: 32]
      -s NUM                s value of soverlap [default: 1]
      --relabel NUM         relabel the hypergraph - 0(hyperedge)/1(hypernode) [default: -1]
      --direction DIR       hypergraph relabeling direction - ascending/descending [default: ascending]
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
  long num_bins= args["-B"].asLong();
  long loader_version = args["--loader-version"].asLong();
  long idx = args["--relabel"].asLong();
  std::string direction = args["--direction"].asString();

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());
  auto _ = set_n_threads(threads[0]);
  std::vector s_values= parse_ids(args["-s"].asStringList());
  if (s_values.empty()) s_values.push_back(1);

  std::string title_basics_tsv = args["--title"].asString();
  std::string name_basics_tsv = args["--name"].asString();
  std::string title_principals_tsv = args["--principal"].asString();

  nw::util::timer               t0("load titles");
  std::ifstream                 title_basics_stream(title_basics_tsv);
  auto                          titles     = xt::load_csv<std::string>(title_basics_stream, '\t');
  auto                          titles_shp = titles.shape();
  std::map<std::string, vertex_id_t> titles_map;
  for (vertex_id_t i = 0; i < titles_shp[0]; ++i) {
    if (titles(i, 1) == "movie") {
      titles_map[titles(i, 0)] = i;
    }
  }
  t0.stop();
  std::cout << t0 << std::endl;

  nw::util::timer               t1("load names");
  std::ifstream                 name_basics_stream(name_basics_tsv);
  auto                          names     = xt::load_csv<std::string>(name_basics_stream, '\t');
  auto                          names_shp = names.shape();
  std::map<std::string, vertex_id_t> names_map;
  for (vertex_id_t i = 0; i < names_shp[0]; ++i) {
    names_map[names(i, 0)] = i;
  }
  t1.stop();
  std::cout << t1 << std::endl;

  nw::util::timer t2("load hypergraph");
  std::ifstream title_principals_stream(title_principals_tsv);
  auto          title_principals = xt::load_csv<std::string>(title_principals_stream, '\t');
  auto          shp              = title_principals.shape();

  t2.stop();
  std::cout << t2 << std::endl;

  nw::util::timer t3("build hypergraph");

  nw::graph::edge_list<nw::graph::directedness::directed> edges;
  edges.open_for_push_back();

  size_t title_counter = 0;
  size_t name_counter  = 0;
  for (vertex_id_t i = 1; i < shp[0]; ++i) {
    if (title_principals(i, 3) == "actor" || title_principals(i, 3) == "actress") {

      auto title = title_principals(i, 0);
      auto name  = title_principals(i, 2);
#if 0
      if (name == "nm0837064") {
        auto idx = titles_map[title];
        auto aa  = titles(idx, 2);
        std::cout << title << " " << aa << std::endl;
      }
      if (titles_map.find(title) == titles_map.end()) {
	      titles_map[title] = title_counter++;
      }
      if (names_map.find(name) == names_map.end()) {
	      names_map[name] = name_counter++;
      }
#endif
      edges.push_back(titles_map[title], names_map[name]);
    }
  }
  edges.close_for_push_back(false);

  t3.stop();
  std::cout << t3 << std::endl;
  edges.stream_stats();

  nw::util::timer t4("build biadjacencies");

  auto G = nw::graph::adjacency<0>(edges);
  auto H = nw::graph::adjacency<1>(edges);

  t4.stop();
  std::cout << t4 << std::endl;

  nw::util::timer t5("build s_overlap");

  index_t s = *std::min_element(s_values.begin(), s_values.end());
  auto&& degrees = H.degrees();
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  tbb::parallel_for(
      nw::graph::cyclic<decltype(H), vertex_id_t>(std::forward<decltype(H)>(H), num_bins),
      [&](auto& i) {
        int worker_index = tbb::task_arena::current_thread_index();
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          if (degrees[hyperE] < s) continue;
          std::map<size_t, vertex_id_t> K;
          for (auto&& [hyperN] : hyperE_ngh) {
            for (auto&& [anotherhyperE] : G[hyperN]) {
              if (degrees[anotherhyperE] < s) continue;
              if (hyperE < anotherhyperE) ++K[anotherhyperE];
            }
          }
          for (auto&& [anotherhyperE, val] : K) {
            if (val >= s)
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      },
      tbb::auto_partitioner());
  nw::graph::edge_list<nw::graph::directedness::undirected> s_overlap(0);
  s_overlap.open_for_push_back();
  // do this in serial
  std::for_each(tbb::counting_iterator<int>(0),
                tbb::counting_iterator<int>(num_bins), [&](auto i) {
                  std::for_each(two_graphs[i].begin(), two_graphs[i].end(),
                                [&](auto&& e) { s_overlap.push_back(e); });
                });
  s_overlap.close_for_push_back(false);
#if 0
  nw::graph::edge_list<nw::graph::directedness::undirected, size_t> s_overlap;
  s_overlap.open_for_push_back();

  for (size_t i = 0; i < H.size(); ++i) {

    if ((i % 8192) == 0) {
      std::cout << i << std::endl;
    }

    for (auto&& [k] : H[i]) {
      for (auto&& [j] : G[k]) {
        if (j > i) {
          s_overlap.push_back(i, j, k);
        }
      }
    }
  }
  s_overlap.close_for_push_back();
#endif
  t5.stop();
  std::cout << t5 << std::endl;

  nw::util::timer t6("build s_overlap adjacency");

  auto L = nw::graph::adjacency<0>(s_overlap);

  t6.stop();
  std::cout << t6 << std::endl;

  // Kevin Bacon is nm0000102
  // David Suchet is nm0837064
  // Kyra Sedgwick is nm0001718

  size_t kevin_bacon   = names_map["nm0000102"];
  size_t david_suchet  = names_map["nm0837064"];
  size_t kyra_sedgwick = names_map["nm0001718"];

  std::vector<size_t> distance(L.size());
  std::vector<size_t> parents(L.size());
  std::vector<size_t> together_in(L.size());

  for (auto&& [u, v] : nw::graph::bfs_edge_range(L, kevin_bacon)) {
    distance[v]    = distance[u] + 1;
    parents[v]     = u;
  }

  std::cout << "Kevin Bacon has a Bacon number of " << distance[kevin_bacon] << std::endl;

  std::cout << "Kyra Sedgwick has a bacon number of " << distance[kyra_sedgwick] << std::endl;
  size_t d = distance[kyra_sedgwick];
  while (kyra_sedgwick != kevin_bacon) {
    std::cout << names(kyra_sedgwick, 1) << " starred with " << names(parents[kyra_sedgwick], 1) << std::endl;
    kyra_sedgwick = parents[kyra_sedgwick];
    if (d-- == 0) {
      break;
    }
  }

  std::cout << "David Suchet has a Bacon number of " << distance[david_suchet] << std::endl;
  d = distance[david_suchet];
  while (david_suchet != kevin_bacon) {
    std::cout << names(david_suchet, 1) << " starred with " << names(parents[david_suchet], 1) << std::endl;
    david_suchet = parents[david_suchet];
    if (d-- == 0) {
      break;
    }
  }

  return 0;
}





int main1(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong();
  long num_bins= args["-B"].asLong();
  long loader_version = args["--loader-version"].asLong();
  bool adjoin  = args["--adjoin"].asBool();
  long idx = args["--relabel"].asLong();
  std::string direction = args["--direction"].asString();

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());
  std::bitset features= parse_features(args["--feature"].asStringList());
  std::vector s_values= parse_ids(args["-s"].asStringList());
  if (s_values.empty()) s_values.push_back(1);

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times_WithS<bool> times;

  for (auto&& file : files) {
    size_t nrealedges = 0, nrealnodes = 0;
    //std::vector<index_t> hyperedge_degrees;
    auto&& [hyperedges, hypernodes, iperm] = graph_reader<directed>(file, idx, direction, adjoin, nrealedges, nrealnodes);
    auto&& hyperedge_degrees = hyperedges.degrees(std::execution::par_unseq);
    if (verbose) {
      hyperedges.stream_indices();
      hypernodes.stream_indices();
      
    
      std::cout << hyperedge_degrees.size() << ": ";
      for (auto d : hyperedge_degrees)
        std::cout << d << " ";
      std::cout << std::endl;
    }
    if (-1 != idx)
      assert(!iperm.empty());
    
    if (adjoin) {
      std::cout << "size of the merged adjacency = " << hyperedges.size() << std::endl;
      assert(0 != nrealedges);
      assert(0 != nrealnodes);
    }

   

    if (verbose) {
      //hypernodes.stream_stats();
      //hyperedges.stream_stats();
    }

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& s : s_values) {
      auto&& s_adj =
          twograph_reader(loader_version, verbose, features, hyperedges,
                          hypernodes, hyperedge_degrees, iperm, nrealedges, nrealnodes, s, num_bins);

      for (auto&& thread : threads) {
        auto _ = set_n_threads(thread);
        for (auto&& id : ids) {
          auto verifier = [&](auto&& E) {
            // only verify #cc in the result
            std::unordered_map<vertex_id_t, size_t> m;
            for (auto& c : E) {
              ++m[c];
            }
            size_t numc = 0;
            for (auto& [k, v] : m) {
              if (1 < v) ++numc;
            }
            std::unordered_set<vertex_id_t> uni_comps(E.begin(), E.end());
            std::cout << m.size() << " components found" << std::endl;
            std::cout << numc << " non-singleton components found" << std::endl;
          };

          auto record = [&](auto&& op) {
            times.record(file, id, thread, s, std::forward<decltype(op)>(op),
                         verifier, true);
          };
          for (int j = 0, e = trials; j < e; ++j) {
            switch (id) {
              case 0:
                record([&] {
                  return linegraph_ccv1(std::execution::par_unseq, hypernodes,
                                        s_adj);
                });
                break;
              case 1:
                record([&] {
                  return linegraph_Afforest(std::execution::par_unseq,
                                            hypernodes, s_adj);
                });
                break;
              default:
                std::cout << "Unknown algorithm version " << id << "\n";
            }
          }
        }
      }  // for thread
    }    // for s
  }      // for file

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log_withs("imdb", file, times, header, "Time(s)");
  }

  return 0;
}
