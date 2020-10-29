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
#include <util/intersection_size.hpp>
#include <util/parallel_for.hpp>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>
#include <util/AtomicBitVector.hpp>
#include <mutex>
#include <bitset>
#include <map>
#include "tbb/task_scheduler_init.h"

using namespace nw::graph;

namespace nw {
namespace hypergraph {

enum heuristics {
ALL = 0,
DEG_PRUNING,
SKIP_VISITED,
SHORT_CIRCUIT,
NONE,
UPPER_TRIANGULAR
};

void print_heuristics_verbose(heuristics h) {
  switch (h) {
  case ALL:
    std::cout << "All are turned on" << std::endl; return;
  case DEG_PRUNING:
    std::cout << "Degree pruning is turned on" << std::endl; return;
  case SKIP_VISITED:
    std::cout << "Skip visited is turned on" << std::endl; return;
  case UPPER_TRIANGULAR:
    std::cout << "Upper triangular is turned on" << std::endl; return;
  case SHORT_CIRCUIT:
    std::cout << "Short circuit is turned on" << std::endl; return;
  case NONE:
    std::cout << "None heuristic is turned on" << std::endl; return;
  default:
    std::cout << "Unknow heuristic is passed" << std::endl; return;
  }
}
/// Basic helper used for all of the inner set intersections.
///
/// This wraps `std::set_intersection` to produce the size of the set rather
/// than the set itself, and also handles the fact that our iterator value types
/// are tuples where we only care about the first element for ordering.
///
/// @tparam           A The type of the first iterator.
/// @tparam           B The type of the second iterator.
/// @tparam           C The type of the third iterator.
/// @tparam           D The type of the fourth iterator.
/// @tparam ExecutionPolicy The type of the parallel execution policy.
///
/// @param            i The beginning of the first range.
/// @param           ie The end of the first range.
/// @param            j The beginning of the second range.
/// @param           je The end of the second range.
/// @param           ep The parallel execution policy.
///
/// @returns            The size of the intersected set.

/*
*
we can pass in s  and check s against n  there when we increment n and as soon we see n==s we can immediately return without doing the full intersection.
and return true/false .
We don't have to use std::set_intersection, instead write our own as this one above and make the return type to be bool
*/
template <class A, class B, class C, class D>
bool is_intersection_size_s(A i, B&& ie, C j, D&& je, size_t s = 1) {
  // Custom comparator because we know our iterator operator* produces tuples
  // and we only care about the first value.
  static constexpr auto lt = [](auto&& x, auto&& y) { return std::get<0>(x) < std::get<0>(y); };

  // Use our own trivial loop for the intersection size when the execution
  // policy is sequential, otherwise rely on std::set_intersection.
  //
  // @todo We really don't need set intersection. You'd hope that it would be
  //       efficient with the output counter, but it just isn't. Parallelizing
  //       the intersection size seems non-trivial though.
    std::size_t n = 0;
    while (i != ie && j != je) {
      if (lt(*i, *j)) {
        ++i;
      } else if (lt(*j, *i)) {
        ++j;
      } else {
        ++n;
        if (n >= s) return true;
        ++i;
        ++j;
      }
    }
    return false;
}

/*
* Squeeze the edge lists such that the ids are consecutive in the new edge list
*/
template<directedness edge_directedness = nw::graph::undirected>
auto squeeze_edgelist(std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> &two_graphs) {
    nw::util::life_timer _(__func__);
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    vertex_id_t index = 0;
    std::unordered_map<vertex_id_t, vertex_id_t> relabel_map;
    int num_bins = two_graphs.size();
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto &&e) {
        auto &&[x, y] = e;
        if (relabel_map.end() == relabel_map.find(x)) {
          //if x has not been relabeled
          relabel_map[x] = index;
          ++index;
        }
        auto newx = relabel_map[x];
        if (relabel_map.end() == relabel_map.find(y)) {
          //if y has not been relabeled
          relabel_map[y] = index;
          ++index;
        }
        auto newy = relabel_map[y];
        //std::cout << x << " " << y << " into " << newx << " " << newy << std::endl;
        result.push_back(newx, newy);
      });
    });
    result.close_for_push_back();

    return result;
}
/*
* Squeeze the edge lists such that the ids are consecutive in the new edge list.
* Here the edge list is weighted.
*/
template<directedness edge_directedness = nw::graph::undirected, class T>
auto squeeze_weighted_edgelist(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>> &two_graphs) {
    nw::util::life_timer _(__func__);
    nw::graph::edge_list<edge_directedness, T> result(0);
    result.open_for_push_back();
    //do this in serial
    vertex_id_t index = 0;
    std::unordered_map<vertex_id_t, vertex_id_t> relabel_map;
    int num_bins = two_graphs.size();
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto &&e) {
        auto &&[x, y, w] = e;
        if (relabel_map.end() == relabel_map.find(x)) {
          //if x has not been relabeled
          relabel_map[x] = index;
          ++index;
        }
        auto newx = relabel_map[x];
        if (relabel_map.end() == relabel_map.find(y)) {
          //if y has not been relabeled
          relabel_map[y] = index;
          ++index;
        }
        auto newy = relabel_map[y];
        //std::cout << x << " " << y << " into " << newx << " " << newy << std::endl;
        result.push_back(newx, newy, w);
      });
    });
    result.close_for_push_back();

    return result;
}

/*
* Parallel naive version, clean version, for perf testing
*/
template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_naive_parallel(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1, int num_bins = 32) {
  size_t M = e_nbs.size();
  std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  {
  nw::util::life_timer _(__func__);
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::task_arena::current_thread_index();
    for (auto i = r.begin(), e = r.end(); i != e; ++i) {
      for (size_t j = i + 1; j < M; ++j) {
        if (nw::graph::intersection_size(e_nbs[i], e_nbs[j]) >= s) {
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(i), std::forward<vertex_id_t>(j)));
        }
      }
    }
  }, tbb::auto_partitioner());
  }
  return squeeze_edgelist(two_graphs);
}

/*
* Parallel naive version with statistic counters, for benchmarking
*/
template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_naive_parallel_with_counter(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1, int num_bins = 32) {

  size_t M = e_nbs.size();
  std::atomic<size_t> counter = 0, nedges = 0;
  std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  {
  nw::util::life_timer _(__func__);
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::task_arena::current_thread_index();
    for (auto i = r.begin(), e = r.end(); i != e; ++i) {
      for (size_t j = i + 1; j < M; ++j) {
        ++counter;
        if (nw::graph::intersection_size(e_nbs[i], e_nbs[j]) >= s) {
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(i), std::forward<vertex_id_t>(j)));
          ++nedges;
        }
      }
    }
  }, tbb::auto_partitioner());
  std::cout << counter << " intersections performed, " 
  << nedges << " edges added" << std::endl;
  }
  return squeeze_edgelist(two_graphs);;
}

/*
* TODO Have not squeeze.
*/
template<directedness edge_directedness = nw::graph::undirected, class HyperEdge, class HyperNode>
auto to_two_graph_naive_serial(HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();
  size_t counter = 0, nedges = 0;
  for (size_t i = 0; i < e_nbs.size(); ++i) {
    for (size_t j = i + 1; j < e_nbs.size(); ++j) {
      ++counter;
      if (nw::graph::intersection_size(e_nbs[i], e_nbs[j]) >= s) {
        two_graph.push_back(i, j);
        ++nedges;
      }
    }
  }
  std::cout << counter << " intersections performed, " 
  << nedges << " edges added" << std::endl;
  two_graph.close_for_push_back();
  return two_graph;
}

/*
* clean without counter. All features on. Fastest version.
*/
template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_clean(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE >= anotherhyperE) continue;
            if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;  
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      } //for
    }, tbb::auto_partitioner());
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            //so we check compid of each hyperedge        
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) continue;
            //filter edges deg(e) < s
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            //avoid duplicate intersections
            if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;         
            //O(average degree of hyperedges)
            if (is_intersection_size_s(edges[hyperE].begin(), edges[hyperE].end(),
            edges[anotherhyperE].begin(), edges[anotherhyperE].end(), s)) {
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          }//each neighbor of hyperN
        }//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    }
    return squeeze_edgelist(two_graphs);
  }//else
}

/*
* clean without counter. All features on. Fastest version.
* Store s value as edge weight.
*/
template<directedness edge_directedness = nw::graph::undirected, class T, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_weighted_efficient_parallel_clean(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) {
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) {
            auto anotherhyperE = std::get<0>(y);
            if (hyperE >= anotherhyperE) return;
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    nw::graph::edge_list<edge_directedness, T> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        auto [i, j] = e;
        result.push_back(i, j, 1);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>> two_graphs(num_bins);
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) {
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) {
            auto anotherhyperE = std::get<0>(y);
            //so we check compid of each hyperedge        
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;
            //filter edges deg(e) < s
            if (hyperedgedegrees[anotherhyperE] < s) return;
            //avoid duplicate intersections
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;         
            //O(average degree of hyperedges)
            size_t s_value = nw::graph::intersection_size(edges[hyperE], edges[anotherhyperE]);
            if (s <= s_value) {
              auto e = std::make_tuple<vertex_id_t, vertex_id_t, T>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE), std::forward<T>(s_value));
              two_graphs[worker_index].push_back(e);
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    }
    return squeeze_weighted_edgelist<nw::graph::undirected, T>(two_graphs);
  }//else
}


/*
* Clean without counter, optional feature.
*/
template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_optional_features_clean(std::bitset<8>& features, ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            if (hyperE >= anotherhyperE) return;
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;

            //filter edges deg(e) < s
            if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[anotherhyperE] < s) return;

            //avoid duplicate intersections
            if (features[heuristics::SKIP_VISITED] && visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;
            
            //O(average degree of hyperedges)
            //short circuit
            auto r = features[heuristics::SHORT_CIRCUIT] ? (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) :
            (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]));
            if (r) {
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    }
    return squeeze_edgelist(two_graphs);
  }//else
}

/*
* With counter, all features on.
*/
template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_with_counter(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nedges = 0;
    //avoid intersection when s=1
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            if (hyperE >= anotherhyperE) return;
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            ++nedges;
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    std::cout << nedges << " edges added" << std::endl;
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nvisited = 0, nbelowS = 0, nduplicates = 0, nintersections = 0, nedges = 0;
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            ++nvisited;
            
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;
            ++nbelowS;

            //filter edges deg(e) < s
            if (hyperedgedegrees[anotherhyperE] < s) return;

            ++nduplicates;
            //avoid duplicate intersections
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;
            
            ++nintersections;
            //O(average degree of hyperedges)
            //short circuit
            if (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) {
              ++nedges;
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    std::cout << nvisited << " visits, "
    << nbelowS << " hyperedges has less S nodes, "
    << nduplicates << " duplicate intersections avoid, "
    << nintersections << " intersections performed, " 
    << nedges << " edges added" << std::endl;
    }
    return squeeze_edgelist(two_graphs);
  }//else
}


/*
* With counter, optional feature.
*/
template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_optional_features_with_counter(std::bitset<8>& features, ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nedges = 0;
    //avoid intersection when s=1
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            if (hyperE >= anotherhyperE) return;
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            ++nedges;
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    std::cout << nedges << " edges added" << std::endl;
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nvisited = 0, nbelowS = 0, nduplicates = 0, nintersections = 0, nedges = 0;
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            ++nvisited;
            
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;
            ++nbelowS;

            //filter edges deg(e) < s
            if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[anotherhyperE] < s) return;

            ++nduplicates;
            //avoid duplicate intersections
            if (features[heuristics::SKIP_VISITED] && visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;
            
            ++nintersections;
            //O(average degree of hyperedges)
            //short circuit
            auto r = features[heuristics::SHORT_CIRCUIT] ? (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) :
            (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]));
            if (r) {
              ++nedges;
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    std::cout << nvisited << " visits, "
    << nbelowS << " hyperedges has less S nodes, "
    << nduplicates << " duplicate intersections avoid, "
    << nintersections << " intersections performed, " 
    << nedges << " edges added" << std::endl;
  }
    return squeeze_edgelist(two_graphs);
  }//else
}

template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_with_map_parallel(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) 
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) 
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();
    return result;
  }
  return squeeze_edgelist(two_graphs);
}

template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_naive_parallel_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1, int num_bins = 32) {
  if(!verbose)
    return to_two_graph_naive_parallel(ep, e_nbs, n_nbs, s, num_bins);
  else
    return to_two_graph_naive_parallel_with_counter(ep, e_nbs, n_nbs, s, num_bins);
}

template<directedness edge_directedness = nw::graph::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_portal(bool verbose, std::bitset<8>& features, ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  //if none feature selected, then default to turn all on
  if (features.none()) features[heuristics::ALL] = true;
  if (features[heuristics::NONE]) {
    //if none is set, manully override each heuristic to false
    for (int i = heuristics::DEG_PRUNING; i != heuristics::NONE; ++i){
      features[i] = false;
    }
  }

  if (features[heuristics::ALL]) {
    if (!verbose) 
      return to_two_graph_efficient_parallel_clean(ep, e_nbs, n_nbs, hyperedgedegrees, s, num_bins);
    else {
      print_heuristics_verbose(heuristics::ALL);
      return to_two_graph_efficient_parallel_with_counter(ep, e_nbs, n_nbs, hyperedgedegrees, s, num_bins);
    }
  }
  else {
    if (!verbose) {
      //clean version without counter
      for (int i = heuristics::DEG_PRUNING; i <= heuristics::NONE; ++i) {
        if (features[i]) print_heuristics_verbose((heuristics)i);
      }
      return to_two_graph_efficient_parallel_optional_features_clean(features, ep, e_nbs, n_nbs, hyperedgedegrees, s, num_bins);
    }
    else {
      //with counter
      for (int i = heuristics::DEG_PRUNING; i <= heuristics::NONE; ++i) {
        if (features[i]) print_heuristics_verbose((heuristics)i);
      }
      return to_two_graph_efficient_parallel_optional_features_with_counter(features, ep, e_nbs, n_nbs, hyperedgedegrees, s, num_bins);
    }
  }
}

}//namespace hypergraph
}//namespace nw
