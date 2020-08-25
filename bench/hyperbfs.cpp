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
#include "common.hpp"
#include "edge_list.hpp"
#include "util/AtomicBitVector.hpp"
#include "util/atomic.hpp"
#include "util/intersection_size.hpp"
#include <docopt.h>

using namespace nw::graph;
using namespace nw::util;
using namespace nw::hypergraph::bench;

static constexpr const char USAGE[] =
    R"(hyperbfs.exe: nw::graph hypergraph breadth-first search benchmark driver.
  Usage:
      hyperbfs.exe (-h | --help)
      hyperbfs.exe -f FILE... [-r NODE | -s FILE] [-a NUM] [-b NUM] [-B NUM] [-n NUM] [--seed NUM] [--version ID...] [--succession STR] [--relabel] [--clean] [--direction DIR]  [--log FILE] [--log-header] [-dvV] [THREADS]...

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


template<class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph(ExecutionPolicy&& exec, HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1) {
  life_timer _(__func__);
  //auto n_nbs = adjacency<0>(H);
  //auto e_nbs = adjacency<1>(H);
  edge_list<undirected> two_graph(0);
  two_graph.open_for_push_back();
  auto counter = 0;
  for (size_t i = 0; i < e_nbs.size(); ++i) {
    for (size_t j = i + 1; j < e_nbs.size(); ++j) {
      ++counter;
      size_t count = nw::graph::intersection_size(e_nbs[i], e_nbs[j]);    //       if intersection_size(n_nbs(i), n_nbs(j)) >= s
      if (count >= s) {
        two_graph.push_back(i, j);    //         add (i,j) to new edge_list
      }
    }
  }
  std::cout << counter << " intersections performed" << std::endl;
  two_graph.close_for_push_back();
  return two_graph;
}

template<class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv2(ExecutionPolicy&& exec, HyperEdge& e_nbs, HyperNode& n_nbs, std::vector<index_t>& hyperedgedegrees, size_t s = 1) {
  life_timer _(__func__);
  //auto n_nbs = adjacency<0>(H);
  //auto e_nbs = adjacency<1>(H);
  edge_list<undirected> two_graph(0);
  two_graph.open_for_push_back();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();
  auto counter = 0;
  //O(n*average degree of hyperedges*average degree of hypernodes*average degree of hyperedges)
  for (size_t i = 0; i < e_nbs.size(); ++i) { //O(n)
    auto hyperE = i;
    if (hyperedgedegrees[hyperE] < s) continue;
    std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) { //O(average degree of hyperedges)
      auto hyperN = std::get<0>(x);
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& y) { //O(average degree of hypernodes)
        //so we check compid of each hyperedge
        auto anotherhyperE = std::get<0>(y);
        if (hyperE > anotherhyperE) return;
        ++counter;
        //O(average degree of hyperedges)
        size_t count = nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]);    //       if intersection_size(n_nbs(i), n_nbs(j)) >= s
        if (count >= s) {
          two_graph.push_back(hyperE, anotherhyperE);    //         add (i,j) to new edge_list
        }
      });
    });
  }
  std::cout << counter << " intersections performed" << std::endl;

  two_graph.close_for_push_back();
  return two_graph;
}

template<class HyperGraph>
auto base_two(HyperGraph& H, size_t s = 1) {
  return base_two(std::execution::seq, H, s);
}

template<class ExecutionPolicy, class HyperNode, class SGraph, class STranspose>
auto base_two(ExecutionPolicy&& exec, HyperNode& hypernodes, SGraph& s_adj, STranspose& s_tran_adj) {
  auto E = ccv1(exec, s_adj);//Afforest(s_adj, s_trans_adj);                 // 3) run whatever on new_adjacency
  auto nhypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> N(nhypernodes);
  //for each hypernode, find N[i]
  for (size_t hyperN = 0; hyperN < nhypernodes; ++hyperN) {
    //get the id of the first neighbor (hyperedge) of hypernode hyperN
    auto hyperE = std::get<0>(*hypernodes[hyperN].begin());
    N[hyperN] = E[hyperE];
  }
  return std::tuple(N, E);
}

template<class ExecutionPolicy, class HyperGraph>
auto to_relabel_graph(ExecutionPolicy&& exec, HyperGraph& aos_a) {

  auto n_nbs = aos_a.max()[1] + 1;
  auto e_nbs = aos_a.max()[0] + 1;
  edge_list<undirected> relabel_graph(0);
  relabel_graph.open_for_push_back();
  //we relabel the smaller set (either hypernode or hyperedge)
  if (n_nbs < e_nbs) {
    std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
      auto&& [edge, node] = elt;
      node = node + e_nbs;
      relabel_graph.push_back(edge, node);    //         add (i,j) to new edge_list
    });
  }
  else {
    std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
      auto&& [edge, node] = elt;
      edge = edge + n_nbs;
      relabel_graph.push_back(edge, node);    //         add (i,j) to new edge_list
    });
  }

  relabel_graph.close_for_push_back();
  return relabel_graph;
}

template<class ExecutionPolicy, class HyperGraph>
auto relabelHyperCC(ExecutionPolicy&& exec, HyperGraph& aos_a) {
  auto relabel_g = to_relabel_graph(exec, aos_a);    // 1) find 2-graph corresponding to s-overlapped hyper edges
  //relabel_g. template lexical_sort_by<0>();
  //relabel_g.uniq();
  //relabel_g.remove_self_loops();

  //relabel_g.stream_edges();

  auto s_adj  = adjacency<0>(relabel_g);        // 2) convert new edge_list to new_adjacency
  //auto s_trans_adj = adjacency<1>(relabel_g);
  //s_adj.stream_indices();

  auto labeling   = //Afforest(s_adj, s_trans_adj); 
  ccv1(exec, s_adj);//

  auto n_nbs = aos_a.max()[1] + 1;
  auto e_nbs = aos_a.max()[0] + 1;
  std::vector<vertex_id_t> N, E;
  if (n_nbs < e_nbs) {
    E.assign(labeling.begin(), labeling.begin() + e_nbs);
    N.assign(labeling.begin() + e_nbs, labeling.end());
  }
  else {
    N.assign(labeling.begin(), labeling.begin() + n_nbs);
    E.assign(labeling.begin() + n_nbs, labeling.end());
  }
  return std::tuple(N, E);
}

/*
* Topdown bfs in serial.
*/
template<typename GraphN, typename GraphE>
auto hyperBFS_topdown_serial_v0(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> parentE(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> parentN(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE;
  //initial node frontier includes every node
  frontierN.reserve(num_hypernodes);
  frontierE.reserve(num_hyperedges);
  frontierE.push_back(source_hyperedge);
  parentE[source_hyperedge] = source_hyperedge; //parent of root is itself

  while (false == (frontierE.empty() && frontierN.empty())) {
    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      std::for_each(hyperedges.begin()[hyperE].begin(), hyperedges.begin()[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) {
          frontierN.push_back(hyperN);
          parentN[hyperN] = hyperE;
        }
      });
    });
    frontierE.clear();
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      std::for_each(hypernodes.begin()[hyperN].begin(), hypernodes.begin()[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
          frontierE.push_back(hyperE);
          parentE[hyperE] = hyperN;
        }
      });
    });
    frontierN.clear();
  }    //while

  return std::tuple{parentN, parentE};
}

/*
* Topdown bfs in serial without using bitmap
*/
template<typename GraphN, typename GraphE>
auto hyperBFS_topdown_serial_v1(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> parentE(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> parentN(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> frontierN, frontierE;
  //initial node frontier includes every node
  frontierN.reserve(num_hypernodes);
  frontierE.reserve(num_hyperedges);
  frontierE.push_back(source_hyperedge);
  parentE[source_hyperedge] = source_hyperedge; //parent of root is itself

  while (false == (frontierE.empty() && frontierN.empty())) {
    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      std::for_each(hyperedges.begin()[hyperE].begin(), hyperedges.begin()[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        if (std::numeric_limits<vertex_id_t>::max() == parentN[hyperN]) {
          frontierN.push_back(hyperN);
          parentN[hyperN] = hyperE;
        }
      });
    });
    frontierE.clear();
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      std::for_each(hypernodes.begin()[hyperN].begin(), hypernodes.begin()[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        if (std::numeric_limits<vertex_id_t>::max() == parentE[hyperE]) {
          frontierE.push_back(hyperE);
          parentE[hyperE] = hyperN;
        }
      });
    });
    frontierN.clear();
  }    //while

  return std::tuple{parentN, parentE};
}

/*
* Bottomup bfs in serial.
*/
template<typename GraphN, typename GraphE>
auto hyperBFS_bottomup_serial_v0(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> parentE(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> parentN(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  parentE[source_hyperedge] = source_hyperedge;

  for (vertex_id_t hyperE = 0; hyperE < num_hyperedges; ++hyperE) {
    if (std::numeric_limits<vertex_id_t>::max() == parentE[hyperE]) {
      //all neighbors of hyperedges are hypernode
      std::for_each(hyperedges.begin()[hyperE].begin(), hyperedges.begin()[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        parentE[hyperN] = hyperN;
        return;    //return once found the first valid parent
      });
    }
  }
  for (vertex_id_t hyperN = 0; hyperN < num_hypernodes; ++hyperN) {
    if (std::numeric_limits<vertex_id_t>::max() == parentN[hyperN]) {
    //all neighbors of hypernodes are hyperedges
    std::for_each(hypernodes.begin()[hyperN].begin(), hypernodes.begin()[hyperN].end(), [&](auto&& x) {
      //so we check compid of each hyperedge
      auto hyperE = std::get<0>(x);
      parentE[hyperE] = hyperN;
      return;
    });
    }
  }

  return std::tuple{parentN, parentE};
}

template<typename Graph>
size_t BU_step(Graph& g, std::vector<vertex_id_t>& parent, 
nw::graph::AtomicBitVector<>& front, nw::graph::AtomicBitVector<>& next) {
  life_timer _(__func__);
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
  life_timer _(__func__);
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
/*
* 
* To save time computing the number of edges exiting the frontier, this
* implementation precomputes the degrees in bulk at the beginning by storing
* them in parent array as negative numbers. Thus the encoding of parent is:
* parent[x] < 0 implies x is unvisited and parent[x] = -out_degree(x)
* parent[x] >= 0 implies x been visited
*/
template<typename Graph>
std::vector<vertex_id_t> init_parent(Graph& g) {
  auto num = g.max() + 1;
  std::vector<vertex_id_t> parent(num);
  std::for_each(counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num), [&](auto u) {
    parent[u] = 0 != g[u].size() ? -g[u].size() : -1;
  });
  return parent;
}
/*
* Direction-optimizing bfs in serial.
*/
template<typename GraphN, typename GraphE>
auto hyperBFS_hybrid_serial_v0(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges,
size_t numpairs, int alpha = 15, int beta = 18) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> parentE(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> parentN(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE;
  frontierN.reserve(num_hypernodes);
  frontierE.reserve(num_hyperedges);
  frontierE.push_back(source_hyperedge);
  parentE[source_hyperedge] = source_hyperedge;

  size_t edges_to_check = numpairs;
  size_t scout_count = hyperedges[source_hyperedge].size();
  while (false == (frontierE.empty() && frontierN.empty())) {
    if (scout_count > edges_to_check / alpha) {
      //BOTTOMUP
      size_t awake_count, old_awake_count;
      queue_to_bitmap(frontierE, visitedE);
      awake_count = frontierE.size();
      while (true) {
        old_awake_count = awake_count;
        awake_count = BU_step(hyperedges, parentE, visitedN, visitedE);
        if ((awake_count < old_awake_count) && (awake_count <= numpairs / beta)) {
          bitmap_to_queue(visitedN, frontierN);
          break;
        }
        old_awake_count = awake_count;
        awake_count = BU_step(hypernodes, parentN, visitedE, visitedN);
        if ((awake_count < old_awake_count) && (awake_count <= numpairs / beta)) {
          bitmap_to_queue(visitedE, frontierE);
          break;
        } 
      } //while
      scout_count = 1;  
    }
    else {
      //TOPDOWN
      edges_to_check -= scout_count;
      scout_count = TD_step(hyperedges, parentE, frontierE, frontierN);
      scout_count = TD_step(hypernodes, parentN, frontierN, frontierE);
    }
  }    //while
  return std::tuple{parentN, parentE};
}

/*
* Direction-optimizing bfs in serial.
*/
template<typename GraphN, typename GraphE>
auto hyperBFS_hybrid_serial_v1(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges,
size_t numpairs, int alpha = 15, int beta = 18) {
  life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> parentE(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> parentN(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE;
  //initial node frontier includes every node
  frontierN.reserve(num_hypernodes);
  frontierE.reserve(num_hyperedges);
  frontierE.push_back(source_hyperedge);
  parentE[source_hyperedge] = 0;

  size_t scout_count = hyperedges[source_hyperedge].size();
  while (false == (frontierE.empty() && frontierN.empty())) {
    if (scout_count > numpairs / alpha) {
      //TOPDOWN for frontierE
      scout_count = 0;
      std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
        //all neighbors of hyperedges are hypernode
        auto levelE = parentE[hyperE];
        std::for_each(hyperedges.begin()[hyperE].begin(), hyperedges.begin()[hyperE].end(), [&](auto&& x) {
          auto hyperN = std::get<0>(x);
          if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) {
            frontierN.push_back(hyperN);
            parentN[hyperN] = hyperE;
          }
        });
        scout_count += hyperedges[hyperE].size();
      });
      frontierE.clear();
    }
    else {
      //BOTTOMUP for frontierE
      scout_count = 0;
      std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
        //all neighbors of hyperedges are hypernode
        auto levelE = parentE[hyperE];
        std::for_each(hyperedges.begin()[hyperE].begin(), hyperedges.begin()[hyperE].end(), [&](auto&& x) {
          auto hyperN = std::get<0>(x);
          if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) {
            frontierN.push_back(hyperN);
            parentN[hyperN] = hyperE;
            return;
          }
        });
        ++scout_count;
      });
      frontierE.clear();
    }
    //split direction optimization processes for frontierE and frontierN
    if (scout_count > numpairs / alpha) {
      //TOPDOWN for frontierN
      scout_count = 0;
      std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
        //all neighbors of hypernodes are hyperedges
        auto levelN = parentN[hyperN];
        std::for_each(hypernodes.begin()[hyperN].begin(), hypernodes.begin()[hyperN].end(), [&](auto&& x) {
          //so we check compid of each hyperedge
          auto hyperE = std::get<0>(x);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
            parentE[hyperE] = hyperN;
          }
        });
        scout_count += hypernodes[hyperN].size();
      });
      frontierN.clear();
    }
    else {
      //BOTTOMUP for frontierN
      std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
        //all neighbors of hypernodes are hyperedges
        auto levelN = parentN[hyperN];
        std::for_each(hypernodes.begin()[hyperN].begin(), hypernodes.begin()[hyperN].end(), [&](auto&& x) {
          //so we check compid of each hyperedge
          auto hyperE = std::get<0>(x);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
            parentE[hyperE] = hyperN;
            return;
          }
        });
        ++scout_count;
      });
      frontierN.clear();
    } //if
  } //while

  return std::tuple{parentN, parentE};
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

int main() {
  return 0;
}
int main2(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;
  long iterations = args["-i"].asLong() ?: 1;
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
    auto reader = [&](std::string file, bool verbose) {
      nw::graph::edge_list<nw::graph::directed> aos_a = load_graph<nw::graph::directed>(file);
      std::vector<nw::graph::index_t> hyperedgedegrees = aos_a.degrees<0>();

      // Run relabeling. This operates directly on the incoming edglist.
      if (args["--relabel"].asBool()) {
        //relabel the column with smaller size
        if (aos_a.max()[0] > aos_a.max()[1]) {
          std::vector<nw::graph::index_t> hypernodedegrees = aos_a.degrees<1>();
          std::cout << "relabeling hypernodes..." << std::endl;
          aos_a.relabel_by_degree_bipartite<1>(args["--direction"].asString(), hypernodedegrees);
        }
        else {
          std::cout << "relabeling hyperedges..." << std::endl;
          aos_a.relabel_by_degree_bipartite<0>(args["--direction"].asString(), hyperedgedegrees);
        }
      }
      // Clean up the edgelist to deal with the normal issues related to
      // undirectedness.
      if (args["--clean"].asBool()) {
        aos_a.swap_to_triangular<0>(args["--succession"].asString());
        aos_a.lexical_sort_by<0>();
        aos_a.uniq();
        aos_a.remove_self_loops();
      }

      adjacency<0> hyperedges(aos_a);
      adjacency<1> hypernodes(aos_a);
      if (verbose) {
        hypernodes.stream_stats();
        hyperedges.stream_stats();
      }
      std::cout << "num_hyperedges = " << aos_a.max()[0] + 1 << " num_hypernodes = " << aos_a.max()[1] + 1 << std::endl;
      return std::tuple(aos_a, hyperedges, hypernodes, hyperedgedegrees);
    };

    auto&& graphs     = reader(file, verbose);
    auto&& aos_a      = std::get<0>(graphs);
    auto&& hyperedges = std::get<1>(graphs);
    auto&& hypernodes = std::get<2>(graphs);
    auto&& hyperedgedegrees = std::get<3>(graphs);

    //all sources are hyperedges
    std::vector<vertex_id_t> sources;
    if (args["--sources"]) {
      sources = load_sources_from_file(hyperedges, args["--sources"].asString());
      trials  = sources.size();
    } else if (args["-r"]) {
      sources.resize(trials);
      std::fill(sources.begin(), sources.end(), args["-r"].asLong());
    } else {
      sources = build_random_sources(hyperedges, trials, args["--seed"].asLong());
    }

    edge_list<undirected> s_over;
    edge_list<undirected> s_over1 = to_two_graph(std::execution::seq, hyperedges, hypernodes);
    if (std::find(ids.begin(), ids.end(), 6) != ids.end())
      s_over = to_two_graphv2(std::execution::seq, hyperedges, hypernodes, hyperedgedegrees);    // 1) find 2-graph corresponding to s-overlapped hyper edges


    auto s_adj  = adjacency<0>(s_over);        // 2) convert new edge_list to new_adjacency
    auto s_trans_adj = adjacency<1>(s_over);

    if (debug) {
      hypernodes.stream_indices();
      hyperedges.stream_indices();
    }

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        for (auto&& source : sources) {
          if (verbose) 
            std::cout << "version " << id << std::endl;

          auto&& [time, parents] = time_op([&] {
            switch (id) {
              case 0:
                return hyperBFS_topdown_serial_v0(source, hypernodes, hyperedges);
              case 1:
                return hyperBFS_topdown_serial_v1(source, hypernodes, hyperedges);
              case 2:
                return hyperBFS_bottomup_serial_v0(source, hypernodes, hyperedges);
              case 3:
                return hyperBFS_hybrid_serial_v0(source, hypernodes, hyperedges, aos_a.size());
              default:
                std::cerr << "Unknown version " << id << "\n";
                return std::make_tuple(std::vector<vertex_id_t>(), std::vector<vertex_id_t>());
            }
          });

          if (verify) {
            hyperBFSVerifier(hypernodes, hyperedges, source, std::get<0>(parents), std::get<1>(parents));
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
