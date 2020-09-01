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
  bool success;
  do
    prev = old;
  while (prev > next && !(success = nw::graph::cas(old, prev, next)));
  return success;
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
inline bool updateAtomic(std::vector<vertex_id_t>& dest, std::vector<vertex_id_t>& source, std::vector<vertex_id_t>& prevDest,
                         vertex_id_t d, vertex_id_t s) {    //atomic Update
  auto origID = dest[d];
  return (writeMin(dest[d], source[s]) && origID == prevDest[d]);
}
// Verifies CC result by performing a BFS from a vertex in each component
// - Asserts search does not reach a vertex with a different component label
// - If the graph is directed, it performs the search as if it was undirected
// - Asserts every vertex is visited (degree-0 vertex should have own label)
template<class ExecutionPolicy, typename GraphN, typename GraphE>
auto bfsCC(ExecutionPolicy&& ep, GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());
  std::vector<vertex_id_t> prevE(num_hyperedges);
  std::vector<vertex_id_t> prevN(num_hypernodes);

  nw::graph::AtomicBitVector   visitedN(num_hypernodes);
  nw::graph::AtomicBitVector   visitedE(num_hyperedges);
  std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  frontierN.reserve(num_hypernodes);
  //frontierE.resize(num_hyperedges);
  //initial node frontier includes every node
  std::iota(frontierE.begin(), frontierE.end(), 0);
  std::iota(E.begin(), E.end(), 0);

  /*
  //or use a parallel for loop
  for (vertex_id_t i = 0; i < num_hyperedges; ++i) {
  //std::for_each(ep, std::begin(0), std::end(num_hypernodes), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  }
  //);
  */

  std::vector<size_t> feworkload, fnworkload;

  while (false == (frontierE.empty() && frontierN.empty())) {
    feworkload.push_back(frontierE.size());

    std::for_each(frontierE.begin(), frontierE.end(), [&](auto& hyperE) {
      //all neighbors of hyperedges are hypernode
      auto labelE = E[hyperE];
      std::for_each(hyperedges[hyperE].begin(), hyperedges[hyperE].end(), [&](auto&& x) {
        auto hyperN = std::get<0>(x);
        auto labelN = N[hyperN];
        if (labelE < labelN) {
          //updateAtomic(N, E, prevN, hyperN, hyperE);
          N[hyperN] = labelE;
          auto old  = N[hyperN];
          while (old == N[hyperN]) {
            auto s = nw::graph::cas(N[hyperN], old, labelE);
            if (s) {
              if (visitedN.atomic_get(hyperN) == 0) {
                // visitedN.atomic_set(hyperN);
                frontierN.push_back(hyperN);
              }
              break;
            }
            old = N[hyperN];
            if (N[hyperN] < labelE) {
              break;
            }
          }
        }
        /*
          if (N[hyperN] == E[hyperE]) return;
          else if (E[hyperE] < N[hyperN]) {
            writeMin(N[hyperN], E[hyperE]);
            if (visitedN.atomic_get(hyperN) == 0 && visitedN.atomic_set(hyperN) == 0) {
              frontierN.push_back(hyperN);
            }
          }
*/
      });
    });
    //reset bitmap for N
    visitedN.clear();
    std::for_each(ep, frontierE.begin(), frontierE.end(), [&](auto& i) {
      //all neighbors of hyperedges are hypernode
      prevE[i] = E[i];
    });
    frontierE.clear();
    fnworkload.push_back(frontierN.size());
    std::for_each(frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
      //all neighbors of hypernodes are hyperedges
      auto labelN = N[hyperN];
      std::for_each(hypernodes[hyperN].begin(), hypernodes[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        auto labelE = E[hyperE];
        if (labelN < labelE) {
          updateAtomic(E, N, prevE, hyperE, hyperN);
          //while (writeMin(E[hyperE], labelN) && prevE[hyperN] == labelE);
          if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
            frontierE.push_back(hyperE);
          }
        }
        /*
          if (E[hyperE] == N[hyperN]) return;
          else if (N[hyperN] < E[hyperE]) {
            writeMin(E[hyperE], N[hyperN]);
            if (visitedE.atomic_get(hyperE) == 0 && visitedE.atomic_set(hyperE) == 0) {
              frontierE.push_back(hyperE);
            }
          }
          */
      });
    });
    //reset bitmap for E

    std::for_each(ep, frontierN.begin(), frontierN.end(), [&](auto& i) {
      //all neighbors of hyperedges are hypernode
      prevN[i] = N[i];
    });
    visitedE.clear();
    frontierN.clear();
  }    //while

  std::cout << "feworkload:";
  for (auto i : feworkload)
    std::cout << i << " ";
  std::cout << std::endl;
  std::cout << "fnworkload:";
  for (auto i : fnworkload)
    std::cout << i << " ";
  std::cout << std::endl;
  return std::tuple{N, E};
}


template<class ExecutionPolicy, typename GraphN, typename GraphE>
auto lpCC(ExecutionPolicy&& ep, GraphN& hypernodes, GraphE& hyperedges) {
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
  //std::iota(frontierE.begin(), frontierE.end(), 0);
  //std::iota(E.begin(), E.end(), 0);
  
  //or use a parallel for loop
  std::for_each(ep, counting_iterator<vertex_id_t>(0), counting_iterator<vertex_id_t>(num_hyperedges), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  });
  
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

template<class ExecutionPolicy, typename GraphN, typename GraphE>
auto lpaNoFrontierCC(ExecutionPolicy&& ep, GraphN& hypernodes, GraphE& hyperedges) {
  nw::util::life_timer _(__func__);
  size_t     num_hypernodes = hypernodes.max() + 1;    // number of hypernodes
  size_t     num_hyperedges = hyperedges.max() + 1;    // number of hyperedges

  std::vector<vertex_id_t> E(num_hyperedges);
  std::vector<vertex_id_t> N(num_hypernodes, std::numeric_limits<vertex_id_t>::max());

  nw::graph::AtomicBitVector visitedN(num_hypernodes);
  nw::graph::AtomicBitVector visitedE(num_hyperedges);
  //std::vector<vertex_id_t> frontierN, frontierE(num_hyperedges);
  //frontierN.reserve(num_hypernodes);
  //frontierE.resize(num_hyperedges);
  //initial node frontier includes every node
  //std::iota(frontierE.begin(), frontierE.end(), 0);
  std::iota(E.begin(), E.end(), 0);
  std::iota(N.begin(), N.end(), 0);
  /*
  //or use a parallel for loop
  for (vertex_id_t i = 0; i < num_hyperedges; ++i) {
  //std::for_each(ep, std::begin(0), std::end(num_hypernodes), [&](auto i) {
    frontierE[i] = i;
    E[i] = i;
  }
  //);
  */
  auto nodes  = hypernodes.begin();
  auto edges  = hyperedges.begin();
  auto change = true;
  //while (false == (frontierE.empty() && frontierN.empty())) {
  while (false == change) {
    //std::for_each(ep, frontierE.begin(), frontierE.end(), [&](auto& hyperE) {

    for (vertex_id_t hyperE = 0; hyperE < num_hyperedges; ++hyperE) {
      //find each active hyperedge
      if (visitedE.atomic_get(hyperE)) {
        //push every hypernode's label to all neighbors of hyperedges
        auto labelE = E[hyperE];
        visitedE.atomic_set(hyperE);
        std::for_each(ep, edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) {
          auto hyperN = std::get<0>(x);
          auto labelN = N[hyperN];
          if (labelE == labelN) {
            //set to unactive
            visitedN.atomic_set(hyperN);
            return;
          } else if (labelE < labelN) {
            //updateAtomic(N, E, prevN, hyperN, hyperE);
            writeMin(N[hyperN], labelE);
            //while (writeMin(N[hyperN], labelE) && prevN[hyperN] == labelN);
            if (!visitedN.atomic_get(hyperN)) {
              visitedN.atomic_set(hyperN);
              change = true;
              //frontierN.push_back(hyperN);
            }
          } else {
            //once we found all the previous update are all false, we append all the false-modified
            //hypernodes to frontierN
            //safe to write without atomic action
            E[hyperE] = labelN;
            //writeMin(E[hyperE], labelN);
            labelE = labelN;
            visitedE.atomic_set(hyperE);
            /*
            std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& y) {
              auto hyperfalseN = std::get<0>(y);
              if (!visitedN.atomic_get(hyperfalseN)){
                visitedN.atomic_set(hyperfalseN);
                frontierN.push_back(hyperfalseN);
              }
            });
            */
          }
        });
      }
    }
    //reset bitmap for N
    //visitedN.clear();
    //frontierE.clear();
    for (vertex_id_t hyperN = 0; hyperN < num_hypernodes; ++hyperN) {
      if (visitedN.atomic_get(hyperN)) {
        //std::for_each(ep, frontierN.begin(), frontierN.end(), [&](auto& hyperN) {
        //all neighbors of hypernodes are hyperedges
        auto labelN = N[hyperN];

        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& x) {
          //so we check compid of each hyperedge
          auto hyperE = std::get<0>(x);
          auto labelE = E[hyperE];
          if (labelN == E[hyperE]) {
            visitedE.atomic_set(hyperE);
            return;
          } else if (labelN < labelE) {
            writeMin(E[hyperE], labelN);
            //updateAtomic(E, N, prevE, hyperE, hyperN);
            //while (writeMin(N[hyperN], labelE) && prevN[hyperN] == labelN);
            if (!visitedE.atomic_get(hyperE)) {
              visitedE.atomic_set(hyperE);
              change = true;
              //frontierE.push_back(hyperE);
            }
          } else {
            N[hyperE] = labelE;
            //writeMin(N[hyperN], labelE);
            labelN = labelE;
            /*
            std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& y) {
              auto hyperfalseE = std::get<0>(y);
              if (!visitedE.atomic_get(hyperfalseE)){
                visitedE.atomic_set(hyperfalseE);
                frontierE.push_back(hyperfalseE);
              }
            });
            */
          }
        });
      }
    }
    //reset bitmap for E
    //visitedE.clear();
    //frontierN.clear();
  }    //while
  return std::tuple{N, E};
}


}//namespace hypergraph
}//namespace nw