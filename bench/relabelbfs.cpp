//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include "Log.hpp"
#include "io/mmio.hpp"
#include "common.hpp"
#include <edge_list.hpp>
#include <algorithms/bfs.hpp>
#include <util/AtomicBitVector.hpp>
#include <docopt.h>
#include "algorithms/relabel_x.hpp"

using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;
using namespace nw::graph;

static constexpr const char USAGE[] =
    R"(hyperbfsrelabel.exe: nw::graph hypergraph breadth-first search benchmark driver.
  Usage:
      hyperbfsrelabel.exe (-h | --help)
      hyperbfsrelabel.exe -f FILE... [-r NODE | -s FILE] [-a NUM] [-b NUM] [-B NUM] [-n NUM] [--seed NUM] [--version ID...] [--succession STR] [--relabel] [--clean] [--direction DIR]  [--log FILE] [--log-header] [-dvV] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      -f FILE               input file paths (can have multiples)
      -n NUM                number of trials [default: 1]
      -a NUM                alpha parameter [default: 15]
      -b NUM                beta parameter [default: 18]
      -B NUM                number of bins [default: 32]
      -r NODE               start from node r (default is random)
      -s, --sources FILE    sources file
      --seed NUM            random seed [default: 27491095]
      --relabel             relabel the graph or not
      -c, --clean           clean the graph or not
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
      --succession STR      successor/predecessor [default: successor]
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
      -v, --verify          verify results
      -V, --verbose         run in verbose mode
)";



template<typename Graph>
size_t BU_step(Graph& g, std::vector<vertex_id_t>& parent, 
nw::graph::AtomicBitVector<>& front, nw::graph::AtomicBitVector<>& next) {
  nw::util::life_timer _(__func__);
  size_t num = g.max() + 1;    // number of hypernodes/hyperedges
  size_t awake_count = 0;
  next.clear();
    //or use a parallel for loop
  std::for_each(counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num), [&](auto u) {
    if (std::numeric_limits<vertex_id_t>::max() == parent[u]) {
      std::for_each(g.begin()[u].begin(), g.begin()[u].end(), [&](auto&& x) {
        auto v = std::get<0>(x);
        if (front.atomic_get(v)) {
          //v has not been claimed
          parent[u] = v;      
          next.atomic_set(u);
          ++awake_count;
          return;
        }
      });
    }
  });
  return awake_count;
}

template<typename Graph>
size_t TD_step(Graph& g, std::vector<vertex_id_t>& parent,
std::vector<vertex_id_t>& front, std::vector<vertex_id_t>& next) {
  nw::util::life_timer _(__func__);
  size_t scout_count = 0;
  std::for_each(front.begin(), front.end(), [&](auto& u) {
    std::for_each(g.begin()[u].begin(), g.begin()[u].end(), [&](auto&& x) {
        auto v = std::get<0>(x);
        if (std::numeric_limits<vertex_id_t>::max() == parent[v]) {
          //v has not been claimed
          parent[v] = u;
          next.push_back(v);
          scout_count += g[v].size();
        }
    });
  });
  return scout_count;
}
inline void queue_to_bitmap(std::vector<vertex_id_t>& queue, nw::graph::AtomicBitVector<>& bitmap) {
  std::for_each(std::execution::par_unseq, queue.begin(), queue.end(), [&](auto&& u) { 
    bitmap.atomic_set(u); 
  });
}
inline void bitmap_to_queue(nw::graph::AtomicBitVector<>& bitmap, std::vector<vertex_id_t>& queue) {
  tbb::parallel_for(bitmap.non_zeros(nw::graph::pow2(15)), [&](auto&& range) {
    for (auto &&i = range.begin(), e = range.end(); i != e; ++i) {
      queue.push_back(*i);
    }
  });
}

template<class ExecutionPolicy, class Graph, class Transpose>
auto relabelhyperBFS_hybrid(ExecutionPolicy&& exec, vertex_id_t source_hyperedge, 
Graph& g, Transpose& g_t, const size_t num_realedges, const size_t num_realnodes,
size_t numpairs, int alpha = 15, int beta = 18) {
  size_t n = g.max() + 1;
  std::vector<vertex_id_t> parent(n);
  parent[source_hyperedge] = source_hyperedge;
  std::vector<vertex_id_t> frontier, nextfrontier;
  frontier.reserve(n);
  nextfrontier.reserve(n);
  frontier.push_back(source_hyperedge);
  nw::graph::AtomicBitVector front(n), cur(n);
  size_t edges_to_check = numpairs;
  size_t scout_count = g[source_hyperedge].size();
  while (!frontier.empty()) {
    if (scout_count > edges_to_check / alpha) {
      size_t awake_count, old_awake_count;
      awake_count = frontier.size();
      do {
        old_awake_count = awake_count;
        awake_count = BU_step(g, parent, front, cur);
        std::swap(front, cur);
        cur.clear();
      } while ((awake_count >= old_awake_count) || (awake_count > n / beta));
      bitmap_to_queue(front, frontier);
      scout_count = 1;
    }
    else {
      edges_to_check -= scout_count;
      scout_count = TD_step(g, parent, frontier, nextfrontier);
      std::swap(frontier, nextfrontier);
      nextfrontier.clear();
    }
  }//while

  std::vector<vertex_id_t> E, N;
  if (num_realnodes < num_realedges) {
    E.assign(parent.begin(), parent.begin() + num_realedges);
    N.assign(parent.begin() + num_realedges, parent.end());
    //for all the parent of hyperN, substract the offset
    for (size_t i = 0, e = N.size(); i != e; ++i) {
      N[i] -= num_realedges;
    }
  }
  else {
    N.assign(parent.begin(), parent.begin() + num_realnodes);
    E.assign(parent.begin() + num_realnodes, parent.end());
    //for all the parent of hyperE, substract the offset
    for (size_t i = 0, e = E.size(); i != e; ++i) {
      E[i] -= num_realnodes;
    }
  }

  return std::tuple{N, E};
}

template<typename GraphN, typename GraphE>
bool hyperBFSVerifier(GraphN& hypernodes, GraphE& hyperedges, vertex_id_t source_hyperedge, 
std::vector<vertex_id_t>& parentE, std::vector<vertex_id_t>& parentN) {
  size_t num_hyperedges = hyperedges.max() + 1;
  size_t num_hypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> depthE(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> depthN(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  depthE[source_hyperedge] = 0;
  std::vector<vertex_id_t> to_visitN, to_visitE;
  to_visitE.reserve(num_hyperedges);
  to_visitN.reserve(num_hypernodes);
  to_visitE.push_back(source_hyperedge);
  auto edges = hyperedges.begin();
  auto nodes = hypernodes.begin();
  //mark depth information for the hypergraph
  while (false == (to_visitE.empty() && to_visitN.empty())) {
    for (auto it = to_visitE.begin(); it != to_visitE.end(); it++) {
      vertex_id_t hyperE = *it;
      for (auto u : edges[hyperE]) {
        vertex_id_t hyperN = std::get<0>(u);
        if (depthN[hyperN] == std::numeric_limits<vertex_id_t>::max()) {
          depthN[hyperN] = depthE[hyperE] + 1;
          to_visitN.push_back(hyperN);
        }
      }
    }
    to_visitE.clear();
    for (auto it = to_visitN.begin(); it != to_visitN.end(); it++) {
      vertex_id_t hyperN = *it;
      for (auto u : nodes[hyperN]) {
        vertex_id_t hyperE = std::get<0>(u);
        if (depthE[hyperE] == std::numeric_limits<vertex_id_t>::max()) {
          depthE[hyperE] = depthN[hyperN] + 1;
          to_visitE.push_back(hyperE);
        }
      }
    }
    to_visitN.clear();
  }//while
  //verify parent information
  for (vertex_id_t hyperE = 0; hyperE < num_hyperedges; ++hyperE) {
    if ((depthE[hyperE] != std::numeric_limits<vertex_id_t>::max()) && std::numeric_limits<vertex_id_t>::max() != (parentE[hyperE])) {
      if (hyperE == source_hyperedge) {
        //verify source
        if (!((parentE[hyperE] == hyperE) && (depthE[hyperE] == 0))) {
          std::cout << "Source wrong " << hyperE << " " << parentE[hyperE] << " " << depthE[hyperE] << std::endl;
          return false;
        }
        continue;
      }
      bool parent_found = false;
      for (auto u : edges[hyperE]) {
        vertex_id_t hyperN = std::get<0>(u);
        if (hyperN == parentE[hyperE]) {
          //if(it != edges[v].end()) {
          if (depthN[hyperN] != depthE[hyperE] - 1) {
            std::cout << "Wrong depths for " << hyperE << " & " << hyperN << std::endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        std::cout << "Couldn't find edge from " << parentE[hyperE] << " to " << hyperE << std::endl;
        return false;
      }
    } else if (depthE[hyperE] != parentE[hyperE]) {
      std::cout << "Reachability mismatch " << hyperE << " " << depthE[hyperE] << " " << parentE[hyperE] << std::endl;
      return false;
    }
  }
  for (vertex_id_t hyperN = 0; hyperN < num_hypernodes; ++hyperN) {
    if ((depthN[hyperN] != std::numeric_limits<vertex_id_t>::max()) && std::numeric_limits<vertex_id_t>::max() != (parentN[hyperN])) {
      bool parent_found = false;
      for (auto u : nodes[hyperN]) {
        vertex_id_t hyperE = std::get<0>(u);
        if (hyperE == parentN[hyperN]) {
          //if(it != edges[v].end()) {
          if (depthE[hyperE] != depthN[hyperN] - 1) {
            std::cout << "Wrong depths for " << hyperN << " & " << hyperE << std::endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        std::cout << "Couldn't find edge from " << parentN[hyperN] << " to " << hyperN << std::endl;
        return false;
      }
    } else if (depthN[hyperN] != parentN[hyperN]) {
      std::cout << "Reachability mismatch " << hyperN << " " << depthN[hyperN] << " " << parentN[hyperN] << std::endl;
      return false;
    }
  }
  return true;
}


int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong();
  long alpha      = args["-a"].asLong() ?: 15;
  long beta       = args["-b"].asLong() ?: 18;
  long num_bins   = args["-B"].asLong() ?: 32;

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  std::vector ids     = parse_ids(args["--version"].asStringList());
  std::vector threads = parse_n_threads(args["THREADS"].asStringList());

  Times<bool> times;

  // Appease clang.
  //
  // These are captured in lambdas later, and if I use structured bindings
  // they have to be listed as explicit captures (this is according to the 17
  // standard). That's a little bit noisy where it happens, so I just give
  // them real symbols here rather than the local bindings.
  for (auto&& file : files) {
    auto reader = [&](std::string file, bool verbose, size_t& nrealedges, size_t& nrealnodes) {
      //auto aos_a   = load_graph<directed>(file);
      auto aos_a   = read_mm_relabeling<nw::graph::directed>(file, nrealedges, nrealnodes);
      if (0 == aos_a.size()) {
        return read_and_relabel_adj_hypergraph(file, nrealedges, nrealnodes);
      }

      nw::graph::adjacency<0> g(aos_a);
      if (verbose) {
        g.stream_stats();
      }
      return g;
    };
    size_t num_realedges, num_realnodes;
    auto&& g    = reader(file, verbose, num_realedges, num_realnodes);
    nw::graph::adjacency<1> g_t(0, 0);
    std::cout << "size of the merged adjacency = " << g.size() << std::endl;
    std::cout << "num_hyperedges = " << num_realedges << " num_hypernodes = " << num_realnodes << std::endl;
    //all sources are hyperedges
    std::vector<vertex_id_t> sources;
    if (args["-r"]) {
      sources.resize(trials);
      std::fill(sources.begin(), sources.end(), args["-r"].asLong());
    } else {
      sources = {0};//build_random_sources(g, trials, args["--seed"].asLong());
    }
    if (debug) {
      g.stream_indices();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        for (auto&& source : sources) {
          if (verbose) 
            std::cout << "version " << id << std::endl;
          using Graph = nw::graph::adjacency<0>;
          using Transpose = nw::graph::adjacency<1>;
          using ExecutionPolicy = decltype(std::execution::par_unseq);
          auto&& [time, parents] = time_op([&] {
            switch (id) {
              case 0: {
                auto v11 = nw::graph::bfs_v11<Graph, Transpose>;
                //bfs_v11(graph, gx, source, num_bins, alpha, beta);
                using BFSV11 = decltype(v11);
                auto&& [N, E] = nw::hypergraph::relabel_x_parallel<ExecutionPolicy, BFSV11, vertex_id_t>(std::execution::par_unseq, num_realedges, num_realnodes, v11, g, g_t, source, num_bins, alpha, beta);
                if (num_realnodes < num_realedges) {
                  //for all the parent of hyperN, substract the offset
                  std::for_each(std::execution::par_unseq, tbb::counting_iterator(0ul), tbb::counting_iterator(N.size()), [&](auto i) {
                    N[i] -= num_realedges;
                  }); 
                }
                else {
                  //for all the parent of hyperE, substract the offset
                   std::for_each(std::execution::par_unseq, tbb::counting_iterator(0ul), tbb::counting_iterator(E.size()), [&](auto i) {
                    E[i] -= num_realnodes;
                  }); 
                }
                return std::tuple(N, E);
              }
              default:
                std::cerr << "Unknown version " << id << "\n";
                return std::make_tuple(std::vector<vertex_id_t>(), std::vector<vertex_id_t>());
            }
          });

          if (verify) {
            hyperBFSVerifier(g, g_t, source, std::get<0>(parents), std::get<1>(parents));
          }

          times.append(file, id, thread, time, source);
        }
      }
    }
  }

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("hyperbfs", file, times, header, "Time(s)");
  }

  return 0;
}
