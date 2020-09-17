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
#include <util/timer.hpp>
#include <util/AtomicBitVector.hpp>
#include <util/atomic.hpp>
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

template <typename T>
inline bool compare_and_swap(T& x, T old_val, T new_val) {
  return __sync_bool_compare_and_swap(&x, *(&old_val), *(&new_val));
}

/*
* baseline has been verified with the Python version
*/
template<class eputionPolicy, typename Graph>
auto baseline(eputionPolicy&& ep, Graph& aos_a) {
  nw::util::life_timer _(__func__);

  size_t num_hyperedges = aos_a.max()[0] + 1;    // number of hyperedges
  size_t num_hypernodes = aos_a.max()[1] + 1;    // number of hypernodes

  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> E(num_hyperedges, std::numeric_limits<vertex_id_t>::max());

  std::for_each(ep, aos_a.begin(), aos_a.end(), [&](auto&& elt) {
    auto&& [edge, node] = elt;
    if (E[edge] == std::numeric_limits<vertex_id_t>::max()) E[edge] = edge;
    if (E[edge] == N[node]) return;
    if (N[node] == std::numeric_limits<vertex_id_t>::max()) {
      N[node] = E[edge];
    } else if (N[node] > E[edge]) {
      auto temp = N[node];
      std::replace(ep, N.begin(), N.end(), temp, E[edge]);
      std::replace(ep, E.begin(), E.end(), temp, E[edge]);
    } else if (N[node] < E[edge]) {
      auto temp = E[edge];
      std::replace(ep, N.begin(), N.end(), temp, N[node]);
      std::replace(ep, E.begin(), E.end(), temp, N[node]);
    }
  });
  // for each u in adjacency
  //   if (vertex is unlabeled) label[u] = u
  //   for each neighbor v of u {
  //     if label[u] < label[v] {
  //       label[v] = label[u]
  //     } else if (label[v] < label[u]) {
  //       label[u] = label[v]
  //      }
  //   }
  // }

  // for each hyperedge he in he_adjacency
  //   if (he is unlabeled) label[he] = he
  //   for each hypernode hn in he
  //     for each hyperedge hf in hn
  //       if label[hf] < label[he] update
  //       if label[he] < label[hf] update

  // for every hyper_edge in original list of hyperedges -- or use edge_range(adjacency)
  // get vertex (or edge) number
  // look corresponding component
  // push_back edge to comps[component].push_back(hyper_edge)

  // OR
  // sort aos_a according to component number
  // std::sort( ... , [])
  
  return std::tuple(N, E);
}

template<class ExecutionPolicy, typename Graph, typename GraphN, typename GraphE>
auto svCC(ExecutionPolicy&& ep, Graph& aos_a, GraphN& cn, GraphE& ce) {
  nw::util::life_timer _(__func__);

  size_t num_hyperedges = ce.max() + 1;    // number of hyperedges
  size_t num_hypernodes = cn.max() + 1;    // number of hypernodes
  std::vector<vertex_id_t> N(num_hypernodes), E(num_hyperedges);
  std::iota(E.begin(), E.end(), 0);
  std::iota(N.begin(), N.end(), 0);
  /*
  for (size_t k = 0; k < E.size(); ++k) {
    E[k] = k;
  }
  */

  bool changed = true;
  for (size_t num_iter = 0; num_iter < aos_a.size(); ++num_iter) {
    if (false == changed) break;
    changed = false;

    for (vertex_id_t edge = 0; edge < num_hyperedges; ++edge) {
      for (auto&& elt : ce[edge]) {
        auto&& [node] = elt;

        if (N[node] > E[edge]) {
          N[node] = E[edge];
          while (N[node] != N[N[node]]) {
            N[node] = N[N[node]];
          }
          changed = true;
        } else if (N[node] < E[edge]) {
          E[edge] = N[node];
          while (E[node] != E[E[node]]) {
            E[node] = E[E[node]];
          }
          changed = true;
        }
      }
    }

    for (vertex_id_t node = 0; node < num_hypernodes; ++node) {
      for (auto&& elt : cn[node]) {
        auto&& [edge] = elt;

        if (N[node] > E[edge]) {
          N[node] = E[edge];
          while (N[node] != N[N[node]]) {
            N[node] = N[N[node]];
          }
          changed = true;
        } else if (N[node] < E[edge]) {
          E[edge] = N[node];
          while (E[node] != E[E[node]]) {
            E[node] = E[E[node]];
          }
          changed = true;
        }
      }
    }
    /*
   for (vertex_id_t n = 0; n < num_hyperedges; n++) {
    while (E[n] != E[E[n]]) {
      E[n] = E[E[n]];
    }
  }

  for (vertex_id_t n = 0; n < num_hypernodes; n++) {
    while (N[n] != N[N[n]]) {
      N[n] = N[N[n]];
    }
  }
  */
  }    //for

  return std::tuple{N, E};
}


template<typename GraphN, typename GraphE>
auto lpCC(GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  frontierN.reserve(num_hypernodes);
  //initial node frontier includes every node
  std::iota(frontierE.begin(), frontierE.end(), 0);
  std::iota(E.begin(), E.end(), 0);

  auto nodes = hypernodes.begin();
  auto edges = hyperedges.begin();
  while (false == (frontierE.empty() && frontierN.empty())) {
    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      auto labelE = E[hyperE];
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        auto labelN = N[hyperN];
        if (labelE < labelN) {
          writeMin(N[hyperN], labelE);
          if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) frontierN.push_back(hyperN);
        }
      });
    });
    //reset bitmap for N
    visitedN.clear();
    frontierE.clear();
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      auto labelN = N[hyperN];
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        auto labelE = E[hyperE];
        if (labelN < labelE) {
          writeMin(E[hyperE], labelN);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
          }
        }
      });
    });
    //reset bitmap for E
    visitedE.clear();
    frontierN.clear();
  }    //while

  return std::tuple{N, E};
}

template<typename GraphN, typename GraphE>
auto lpCCv2(GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> N(num_hypernodes);

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN(num_hypernodes), frontierE;
  frontierE.reserve(num_hyperedges);
  //initial node frontier includes every node
  std::iota(frontierN.begin(), frontierN.end(), 0);
  std::iota(N.begin(), N.end(), 0);

  auto nodes = hypernodes.begin();
  auto edges = hyperedges.begin();
  while (false == (frontierN.empty() && frontierE.empty())) {
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      auto labelN = N[hyperN];
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        auto labelE = E[hyperE];
        if (labelN < labelE) {
          writeMin(E[hyperE], labelN);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
          }
        }
      });
    });
    //reset bitmap for E
    visitedE.clear();
    frontierN.clear();
    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      auto labelE = E[hyperE];
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        auto labelN = N[hyperN];
        if (labelE < labelN) {
          writeMin(N[hyperN], labelE);
          if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) 
            frontierN.push_back(hyperN);
        }
      });
    });
    //reset bitmap for N
    visitedN.clear();
    frontierE.clear();
  }    //while

  return std::tuple{N, E};
}

inline bool updateAtomic(std::vector<vertex_id_t>& dest, std::vector<vertex_id_t>& source, std::vector<vertex_id_t>& prevDest,
                         vertex_id_t d, vertex_id_t s) {    //atomic Update
  auto origID = dest[d];
  return (writeMin(dest[d], source[s]) && origID == prevDest[d]);
}


template <typename T>
bool hook(T u, T v, std::vector<T>& compu, std::vector<T>& compv) {
  T p1 = compu[u];
  T p2 = compv[v];
  T high, low;
  bool success = false;
  while (p1 != p2) {
    if (p1 > p2) {
      high = p1;
      low = p2;
      volatile T p_high = compu[high];
      if (p_high == low || (success = compare_and_swap(compu[high], high, low))) break;
      p1 = compu[high];
      p2 = compv[low];
    }
    else {
      high = p2;
      low = p1;
      volatile T p_high = compv[high];
      if (p_high == low || (success = compare_and_swap(compv[high], high, low))) break;
      p1 = compu[low];
      p2 = compv[high];
    }
  }    // while
  return success;
}

template<class ExecutionPolicy, typename GraphN, typename GraphE>
auto lpCC_parallelv1(ExecutionPolicy&& ep, GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes);

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  tbb::concurrent_vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  frontierN.reserve(num_hypernodes);
  //initial edge frontier includes every node
  //or use a parallel for loop
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hyperedges), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  });
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hypernodes), [&](auto i) {
    N[i] = std::numeric_limits<vertex_id_t>::max();
  });
  
  auto nodes = hypernodes.begin();
  auto edges = hyperedges.begin();
  while (false == (frontierE.empty() && frontierN.empty())) {
    std::for_each(ep, frontierE.begin(), frontierE.end(), [&](auto hyperE) {
      //all neighbors of hyperedges are hypernode
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        if (E[hyperE] < N[hyperN]){ 
          if (writeMin(N[hyperN], E[hyperE])) {
            if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) 
              frontierN.push_back(hyperN);
          }
        }
      });
    });
    //reset bitmap for N
    visitedN.clear();
    frontierE.clear();
    std::for_each(ep, frontierN.begin(), frontierN.end(), [&](auto hyperN) {
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        if (N[hyperN] < E[hyperE]) {
          if (writeMin(E[hyperE], N[hyperN])) {
            if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) 
              frontierE.push_back(hyperE);
          }
        }
      });
    });
    //reset bitmap for E
    visitedE.clear();
    frontierN.clear();
  }    //while
  return std::tuple{N, E};
}


template<class ExecutionPolicy, typename GraphN, typename GraphE>
auto lpCC_parallelv2(ExecutionPolicy&& ep, GraphN& hypernodes, GraphE& hyperedges, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes);

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  frontierN.reserve(num_hypernodes);
  //initial edge frontier includes every node
  //or use a parallel for loop
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hyperedges), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  });
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hypernodes), [&](auto i) {
    N[i] = std::numeric_limits<vertex_id_t>::max();
  });
  
  auto nodes = hypernodes.begin();
  auto edges = hyperedges.begin();
  std::vector<vertex_id_t> frontier[num_bins];
  while (false == (frontierE.empty() && frontierN.empty())) {

    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0ul, frontierE.size()), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      for (auto i = r.begin(), e = r.end(); i < e; ++i) {
        vertex_id_t hyperE = frontierE[i];
        vertex_id_t labelE = E[hyperE];
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) {
          auto hyperN = std::get<0>(x);
          if (labelE < N[hyperN]) {
            if (writeMin(N[hyperN], labelE)) {
              if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0)
                frontier[worker_index].push_back(hyperN);
            }
          }
        });
      } //for
    }, tbb::auto_partitioner());
    for (int i = 0; i < num_bins; ++i) {
      frontierN.insert(frontierN.end(), frontier[i].begin(), frontier[i].end());
      frontier[i].clear();
    }
    //std::sort(frontierN.begin(), frontierN.end());
    //std::unique(frontierN.begin(), frontierN.end());
    //reset bitmap for N
    visitedN.clear();
    frontierE.clear();
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0ul, frontierN.size()), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      for (auto i = r.begin(), e = r.end(); i < e; ++i) {
        vertex_id_t hyperN = frontierN[i];
        vertex_id_t labelN = N[hyperN];
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&x) {
          //so we check compid of each hyperedge
          auto hyperE = std::get<0>(x);
          if (labelN < E[hyperE]) {
            if (writeMin(E[hyperE], labelN)) {
              if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0)
                frontier[worker_index].push_back(hyperE);
            }
          }
        });
      } //for
    }, tbb::auto_partitioner());
    for (int i = 0; i < num_bins; ++i) {
      frontierE.insert(frontierE.end(), frontier[i].begin(), frontier[i].end());
      frontier[i].clear();
    }
    //std::sort(frontierE.begin(), frontierE.end());
    //std::unique(frontierE.begin(), frontierE.end());
    //reset bitmap for E
    visitedE.clear();
    frontierN.clear();
  }    //while
  return std::tuple{N, E};
}


template<class ExecutionPolicy, typename GraphN, typename GraphE>
auto lpaNoFrontierCC(ExecutionPolicy&& ep, GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes);


  nw::graph::AtomicBitVector activeN(num_hypernodes);
  nw::graph::AtomicBitVector activeE(num_hyperedges);
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hyperedges), [&](auto i) {
    E[i] = i;
    activeE.atomic_set(i);
  });
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hypernodes), [&](auto i) {
    N[i] = std::numeric_limits<vertex_id_t>::max();
    activeN.atomic_set(i);
  });

  auto nodes  = hypernodes.begin();
  auto edges  = hyperedges.begin();
  auto change = true;
  while (change) {
    change = false;
    tbb::parallel_for(activeE.non_zeros(num_hyperedges), [&](auto&& range) {
      for (auto&& i = range.begin(), e = range.end(); i != e; ++i) {
        auto hyperE = *i;
        auto labelE = E[hyperE];
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
          auto hyperN = std::get<0>(x);
          auto labelN = N[hyperN];
          if (labelE < labelN) {
            writeMin(N[hyperN], labelE);
            change = true;
          }
          else
            activeN.atomic_set(hyperN);
        });
      }
    });
    tbb::parallel_for(activeN.non_zeros(num_hypernodes), [&](auto&& range) {
      for (auto&& i = range.begin(), e = range.end(); i != e; ++i) {
        auto hyperN = *i;
        auto labelN = N[hyperN];
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
          //so we check compid of each hyperedge
          auto hyperE = std::get<0>(x);
          auto labelE = E[hyperE];
          if (labelN < labelE) {
            writeMin(E[hyperE], labelN);
            change = true;
          } else {
              activeE.atomic_set(hyperE);
          }
        });
      }
    });
  }    //while

  return std::tuple{N, E};
}


}//namespace hypergraph
}//namespace nw