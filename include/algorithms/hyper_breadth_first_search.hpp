//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2020
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//


#pragma once
#include <nwgraph/util/timer.hpp>
#include <nwgraph/util/AtomicBitVector.hpp>
#include <nwgraph/util/atomic.hpp>
#include <nwgraph/adaptors/vertex_range.hpp>

namespace nw {
namespace hypergraph {

template<class T>
inline bool writeMin(T& old, T& next) {
  T    prev;
  bool success = false;
  do
    prev = old;
  while (prev > next && !(success = nw::graph::cas(old, prev, next)));
  return success;
}
/*
* Topdown bfs in serial.
*/
template<typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
auto hyperBFS_topdown_serial_v0(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
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
        if (std::numeric_limits<vertex_id_t>::max() == parentN[hyperN]) {
          if (visitedN.atomic_get(hyperN) == 0 &&
              visitedN.atomic_set(hyperN) == 0) {
            frontierN.push_back(hyperN);
            parentN[hyperN] = hyperE;
          }
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
          if (visitedE.atomic_get(hyperE) == 0 &&
              visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
            parentE[hyperE] = hyperN;
          }
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
template<typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
auto hyperBFS_topdown_serial_v1(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
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

template<class ExecutionPolicy, typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
auto hyperBFS_topdown_parallel_v0(ExecutionPolicy&& ep, const vertex_id_t source_hyperedge, GraphN& nodes, GraphE& edges, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t num_hyperedges = num_vertices(edges);
  size_t num_hypernodes = num_vertices(nodes);

  std::vector<vertex_id_t> parentE(num_hyperedges);
  std::vector<vertex_id_t> parentN(num_hypernodes);

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE;
  frontierN.reserve(num_hypernodes);
  frontierE.reserve(num_hyperedges);
  //initial edge frontier includes every node
  //or use a parallel for loop
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hyperedges), [&](auto i) {
    parentE[i] = std::numeric_limits<vertex_id_t>::max();
  });
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hypernodes), [&](auto i) {
    parentN[i] = std::numeric_limits<vertex_id_t>::max();
  });
  frontierE.push_back(source_hyperedge);
  parentE[source_hyperedge] = source_hyperedge; //parent of root is itself


  std::vector<vertex_id_t> frontier[num_bins];
  auto traverse = [&](auto& g, auto& cur, auto& bitmap, auto& parents) {
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0ul, cur.size()), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      for (auto i = r.begin(), e = r.end(); i < e; ++i) {
        vertex_id_t x = cur[i];
        //all neighbors of hyperedges are hypernode, vice versa
        std::for_each(g[x].begin(), g[x].end(), [&](auto&& j) {
          auto y = std::get<0>(j);
          if (std::numeric_limits<vertex_id_t>::max() == parents[y]) {
            if (writeMin(parents[y], x)) {
              if (0 == bitmap.atomic_get(y) && 0 == bitmap.atomic_set(y)) {
                frontier[worker_index].push_back(y);
              }
            }
          }
        });
      } //for
    }, tbb::auto_partitioner());
  };

  std::vector<size_t> size_array(num_bins);
  auto curtonext = [&](auto& next) {
    size_t size = 0;
    for (int i = 0; i < num_bins; ++i) {
      //calculate the size of each thread-local frontier
      size_array[i] = size;
      //accumulate the total size of all thread-local frontiers
      size += frontier[i].size();
    }
    //resize next frontier
    next.resize(size); 
    std::for_each(ep, nw::graph::counting_iterator(0), nw::graph::counting_iterator(num_bins), [&](auto i) {
      //copy each thread-local frontier to next frontier based on their size offset
      auto begin = std::next(next.begin(), size_array[i]);
      std::copy(ep, frontier[i].begin(), frontier[i].end(), begin);
      frontier[i].clear();
    });
  };

  while (false == (frontierE.empty() && frontierN.empty())) {
    traverse(edges, frontierE, visitedN, parentN);
    curtonext(frontierN);
    visitedN.clear();
    frontierE.clear();
    traverse(nodes, frontierN, visitedE, parentE);
    curtonext(frontierE);
    visitedE.clear();
    frontierN.clear();
  } 
  return std::tuple(parentN, parentE);
}

/*
* Bottomup bfs in serial.
*/
template<typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
auto hyperBFS_bottomup_serial_v0(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
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
            writeMin(parentE[hyperE], hyperN);//parentE[hyperE] = hyperN;
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
             writeMin(parentN[hyperN], hyperE);//parentN[hyperN] = hyperE;
            return;
        });
    }
  }

  return std::tuple{parentN, parentE};
}

template<class ExecutionPolicy, typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
auto hyperBFS_bottomup_parallel_v0(ExecutionPolicy&& ep, const vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> parentE(num_hyperedges);
  std::vector<vertex_id_t> parentN(num_hypernodes);
  //initial edge frontier includes every node
  //or use a parallel for loop
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hyperedges), [&](auto i) {
    parentE[i] = std::numeric_limits<vertex_id_t>::max();
  });
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hypernodes), [&](auto i) {
    parentN[i] = std::numeric_limits<vertex_id_t>::max();
  });
  parentE[source_hyperedge] = source_hyperedge; //parent of root is itself

  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0ul, num_hyperedges), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::this_task_arena::current_thread_index();
    for (vertex_id_t hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
      if (std::numeric_limits<vertex_id_t>::max() == parentE[hyperE]) {
      //all neighbors of hyperedges are hypernode
        std::for_each(hyperedges.begin()[hyperE].begin(), hyperedges.begin()[hyperE].end(), [&](auto&& x) {
            auto hyperN = std::get<0>(x);
            writeMin(parentE[hyperE], hyperN);
            return;    //return once found the first valid parent
        });
      }
    }
  }, tbb::auto_partitioner());
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0ul, num_hypernodes), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::this_task_arena::current_thread_index();
    for (vertex_id_t hyperN = r.begin(), e = r.end(); hyperN < e; ++hyperN) {
    if (std::numeric_limits<vertex_id_t>::max() == parentN[hyperN]) {
    //all neighbors of hypernodes are hyperedges
        std::for_each(hypernodes.begin()[hyperN].begin(), hypernodes.begin()[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
            auto hyperE = std::get<0>(x);
            writeMin(parentN[hyperN], hyperE);
            return;    //return once found the first valid parent
        });
    }
    }
  }, tbb::auto_partitioner());

  return std::tuple(parentN, parentE);
}

template<typename Graph, typename vertex_id_t = vertex_id_t<Graph>>
size_t BU_step(Graph& g, std::vector<vertex_id_t>& parent, 
nw::graph::AtomicBitVector<>& front, nw::graph::AtomicBitVector<>& next) {
  nw::util::life_timer _(__func__);
  size_t num = g.max() + 1;    // number of hypernodes/hyperedges
  size_t awake_count = 0;
  next.clear();
    //or use a parallel for loop
  std::for_each(nw::graph::counting_iterator<vertex_id_t>(0), nw::graph::counting_iterator<vertex_id_t>(num), [&](auto u) {
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

template<typename Graph, typename vertex_id_t = vertex_id_t<Graph>>
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
template<class ExecutionPolicy, class Vector>
inline void queue_to_bitmap(ExecutionPolicy&& ep, Vector& queue, nw::graph::AtomicBitVector<>& bitmap) {
  std::for_each(ep, queue.begin(), queue.end(), [&](auto&& u) { 
    bitmap.atomic_set(u); 
  });
}
template<class Vector>
inline void bitmap_to_queue(nw::graph::AtomicBitVector<>& bitmap, Vector& queue) {
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
template<typename Graph, typename vertex_id_t = vertex_id_t<Graph>>
std::vector<vertex_id_t> init_parent(Graph& g) {
  auto num = g.max() + 1;
  std::vector<vertex_id_t> parent(num);
  std::for_each(nw::graph::counting_iterator<vertex_id_t>(0), nw::graph::counting_iterator<vertex_id_t>(num), [&](auto u) {
    parent[u] = 0 != g[u].size() ? -g[u].size() : -1;
  });
  return parent;
}
/*
* Direction-optimizing bfs in serial. TODO Cannot converge
*/
template<typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
auto hyperBFS_hybrid_serial_v1(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges,
size_t numpairs, int alpha = 15, int beta = 18) {
  nw::util::life_timer _(__func__);
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
      queue_to_bitmap(std::execution::par_unseq, frontierE, visitedE);
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
* Direction-optimizing bfs in serial, verified working version
*/
template<typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
auto hyperBFS_hybrid_serial_v0(vertex_id_t source_hyperedge, GraphN& hypernodes, GraphE& hyperedges,
size_t numpairs, int alpha = 15, int beta = 18) {
  nw::util::life_timer _(__func__);
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


template<typename GraphN, typename GraphE, typename vertex_id_t = vertex_id_t<GraphN>>
bool hyperBFSVerifier(GraphN& nodes, GraphE& edges, vertex_id_t source_hyperedge, 
std::vector<vertex_id_t>& parentN, std::vector<vertex_id_t>& parentE) {
  size_t num_hyperedges = num_vertices(edges);
  size_t num_hypernodes = num_vertices(nodes);
  std::vector<vertex_id_t> depthE(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> depthN(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  depthE[source_hyperedge] = 0;
  std::vector<vertex_id_t> to_visitN, to_visitE;
  to_visitE.reserve(num_hyperedges);
  to_visitN.reserve(num_hypernodes);
  to_visitE.push_back(source_hyperedge);

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
  for (vertex_id_t hyperE = 0; hyperE < num_hyperedges; ++hyperE) {
    if ((depthE[hyperE] != std::numeric_limits<vertex_id_t>::max()) &&
        (std::numeric_limits<vertex_id_t>::max() != parentE[hyperE])) {
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
          if (!((depthE[hyperE] == depthN[hyperN] - 1) || (depthE[hyperE] - 1 == depthN[hyperN]))) {
            std::cout << "Wrong depths for hyperedge " << hyperE << " & " << hyperN << std::endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        std::cout << "Couldn't find edge from node " << parentE[hyperE] << " to hyperedge " << hyperE << std::endl;
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
          if (!((depthE[hyperE] == depthN[hyperN] - 1) || (depthE[hyperE] - 1 == depthN[hyperN]))) {
            std::cout << "Wrong depths for node " << hyperN << " & " << hyperE << std::endl;
            return false;
          }
          parent_found = true;
          break;
        }
      }
      if (!parent_found) {
        std::cout << "Couldn't find edge from hyperedge " << parentN[hyperN] << " to node " << hyperN << std::endl;
        return false;
      }
    } else if (depthN[hyperN] != parentN[hyperN]) {
      std::cout << "Reachability mismatch " << hyperN << " " << depthN[hyperN] << " " << parentN[hyperN] << std::endl;
      return false;
    }
  }
  return true;
}


}//namespace hypergraph
}//namespace nw