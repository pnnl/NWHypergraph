//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2021
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//
#pragma once
#include <bitset>
#include <util/intersection_size.hpp>
#include <adaptors/cyclic_neighbor_range.hpp>
#include <adaptors/vertex_range.hpp>

#include "util/slinegraph_helper.hpp"
#include <tbb/task_arena.h>
#include <tbb/blocked_range2d.h>

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

/**
* This function prints the heuristics for efficient algorithm.
*
* @param[in] h heuristic enum ID
*
*/
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

/// Short circuit for the set intersection.
///
/// This wraps `std::set_intersection` to produce the size of the set rather
/// than the set itself, and also handles the fact that our iterator value types
/// are tuples where we only care about the first element for ordering.
/// we can pass in s and check s against n there when we increment n and as soon we see n==s we can immediately return without doing the full intersection.
/// and return true/false .
/// We don't have to use std::set_intersection, instead write our own as this one above and make the return type to be bool
///
/// @tparam           A The type of the first iterator.
/// @tparam           B The type of the second iterator.
/// @tparam           C The type of the third iterator.
/// @tparam           D The type of the fourth iterator.
///
/// @param            i The beginning of the first range.
/// @param           ie The end of the first range.
/// @param            j The beginning of the second range.
/// @param           je The end of the second range.
/// @param            s The value of s to short circuit.
///
/// @returns            true/false whether the first range and the second range 
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
/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. This version is for s = 1.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
* @param[in] range_begin the begin of the range to split
* @param[in] range_end the end of the range to split
*
*/
template <class HyperEdge, class HyperNode>
void to_two_graph_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, size_t range_begin,
    size_t range_end) {
  nw::util::life_timer _(__func__);
  auto M = edges.size();
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(range_begin, range_end, (range_end - range_begin) / num_bins),
      [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
          std::fill(visitedE.begin(), visitedE.end(), false);
          // all neighbors of hyperedges are hypernode
          for (auto&& [hyperN] : edges[hyperE]) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperE >= anotherhyperE) continue;
              //if hyperE-anotherhyperE has been inserted, abort
              //this is to avoid duplicate edges inserted into slinegraph
              if (visitedE[anotherhyperE])
                continue;
              else
                visitedE[anotherhyperE] = true;
              two_graphs[worker_index].push_back(
                  std::make_pair<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
            }
          }
        }  // for
      },
      tbb::auto_partitioner());
}

/**
* Efficient computation of a s-line graph of a hypergraph for s > 1. 
* It uses blocked range as workload distribution strategy.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
* @param[in] range_begin the begin of the range to split
* @param[in] range_end the end of the range to split
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
*/
template <class HyperEdge, class HyperNode>
void to_two_graph_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, size_t range_begin,
    size_t range_end, std::vector<index_t>& hyperedgedegrees, size_t s) {
  nw::util::life_timer _(__func__);
  auto M = edges.size();
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(range_begin, range_end, (range_end - range_begin) / num_bins),
      [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
          std::fill(visitedE.begin(), visitedE.end(), false);
          if (hyperedgedegrees[hyperE] < s) continue;
          // all neighbors of hyperedges are hypernode
          for (auto&& [hyperN] : edges[hyperE]) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              // so we check compid of each hyperedge
              // travese upper triangluar with lhs > rhs
              // avoid self edge with lhs == rhs
              if (hyperE >= anotherhyperE) continue;
              // filter edges deg(e) < s
              if (hyperedgedegrees[anotherhyperE] < s) continue;
              // avoid duplicate intersections
              if (visitedE[anotherhyperE])
                continue;
              else
                visitedE[anotherhyperE] = true;
              // O(average degree of hyperedges)
              if (is_intersection_size_s(edges[hyperE].begin(),
                                         edges[hyperE].end(),
                                         edges[anotherhyperE].begin(),
                                         edges[anotherhyperE].end(), s)) {
                two_graphs[worker_index].push_back(
                    std::make_pair<vertex_id_t, vertex_id_t>(
                        std::forward<vertex_id_t>(hyperE),
                        std::forward<vertex_id_t>(anotherhyperE)));
              }
            }  // each neighbor of hyperN
          }    // each neighbor of hyperE
        }      // for each hyperE
      },
      tbb::auto_partitioner());
}
/**
* Efficient computation of a weighted s-line graph of a hypergraph for s = 1. 
* It uses blocked range as workload distribution strategy.
*
* @tparam T the type of weight in the 1-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
* @param[in] range_begin the begin of the range to split
* @param[in] range_end the end of the range to split
*
*/
template <class T, class HyperEdge, class HyperNode>
void to_weighted_two_graph_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, size_t range_begin,
    size_t range_end) {
  nw::util::life_timer _(__func__);
  auto M = edges.size();
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(range_begin, range_end, (range_end - range_begin) / num_bins),
      [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
          std::fill(visitedE.begin(), visitedE.end(), false);
          // all neighbors of hyperedges are hypernode
          for (auto&& [hyperN] : edges[hyperE]) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperE >= anotherhyperE) continue;
              //if hyperE-anotherhyperE has been inserted, abort
              //this is to avoid duplicate edges inserted into slinegraph
              if (visitedE[anotherhyperE])
                continue;
              else
                visitedE[anotherhyperE] = true;
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t, T>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE), 
                      1));
            }
          }
        }  // for
      },
      tbb::auto_partitioner());
}
/**
* Efficient computation of a weighted s-line graph of a hypergraph for s > 1. 
* The weight is the number of overlapping vertices.
* Can not use short circuit (i.e. is_intersection_size_s) any more.
* It uses blocked range as workload distribution strategy.
*
* @tparam T the type of weight in the 1-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
* @param[in] range_begin the begin of the range to split
* @param[in] range_end the end of the range to split
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
*/
template <class T, class HyperEdge, class HyperNode>
void to_weighted_two_graph_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, size_t range_begin,
    size_t range_end, std::vector<index_t>& hyperedgedegrees, size_t s) {
  nw::util::life_timer _(__func__);
  auto M = edges.size();
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(range_begin, range_end, (range_end - range_begin) / num_bins),
      [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
          std::fill(visitedE.begin(), visitedE.end(), false);
          if (hyperedgedegrees[hyperE] < s) continue;
          // all neighbors of hyperedges are hypernode
          for (auto&& [hyperN] : edges[hyperE]) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              // so we check compid of each hyperedge
              // travese upper triangluar with lhs > rhs
              // avoid self edge with lhs == rhs
              if (hyperE >= anotherhyperE) continue;
              // filter edges deg(e) < s
              if (hyperedgedegrees[anotherhyperE] < s) continue;
              // avoid duplicate intersections
              if (visitedE[anotherhyperE])
                continue;
              else
                visitedE[anotherhyperE] = true;
              // O(average degree of hyperedges)
              size_t s_value = nw::graph::intersection_size(edges[hyperE], edges[anotherhyperE]);
              if (s_value >= s) {
                two_graphs[worker_index].push_back(
                    std::make_tuple<vertex_id_t, vertex_id_t>(
                        std::forward<vertex_id_t>(hyperE),
                        std::forward<vertex_id_t>(anotherhyperE), 
                        s_value));
              }
            }  // each neighbor of hyperN
          }    // each neighbor of hyperE
        }      // for each hyperE
      },
      tbb::auto_partitioner());
}
/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses experimental 2 dimensional blocked range as workload distribution strategy. This version is for s = 1.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
* @param[in] rowrange_begin the begin of the 1D range to split
* @param[in] rowrange_end the end of the 1D range to split
* @param[in] row_grainsize the grain size of the 1D range
* @param[in] colrange_begin the begin of the 2D range to split
* @param[in] colrange_end the end of the 2D range to split
* @param[in] col_grainsize he grain size of the 2D range
*/
template<class HyperEdge, class HyperNode>
void to_two_graph_block_range2d(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs, 
HyperEdge& edges, HyperNode& nodes, int num_bins, 
size_t rowrange_begin, size_t rowrange_end, size_t row_grainsize,
size_t colrange_begin, size_t colrange_end, size_t col_grainsize) {
    //std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    auto M = edges.size();
    tbb::parallel_for(tbb::blocked_range2d<vertex_id_t>(rowrange_begin, rowrange_end, row_grainsize, colrange_begin, colrange_end, col_grainsize), 
    [&](tbb::blocked_range2d<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      //for each thread, create a bitmap for visited edges
      for (auto i = r.cols().begin(), ie = r.cols().end(); i != ie; ++i) {
        for (auto j = r.rows().begin(), e = r.rows().end(); j != e; ++j) {
            auto hyperE = num_bins * j + i;
            std::fill(visitedE.begin(), visitedE.end(), false);
                //all neighbors of hyperedges are hypernode
            for (auto &&[hyperN] : edges[hyperE]) {
                for (auto &&[anotherhyperE] : nodes[hyperN]) {
                    //if (hyperE >= anotherhyperE) continue;
                    if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;  
                    two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
                }
            }
            
        }//for cols
      }//for rows
      
    }, tbb::auto_partitioner());
    //return two_graphs;
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses experimental 2 dimensional blocked range as workload distribution strategy. This version is for s > 1.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
* @param[in] rowrange_begin the begin of the 1D range to split
* @param[in] rowrange_end the end of the 1D range to split
* @param[in] row_grainsize the grain size of the 1D range
* @param[in] colrange_begin the begin of the 2D range to split
* @param[in] colrange_end the end of the 2D range to split
* @param[in] col_grainsize he grain size of the 2D range
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
*/
template<class HyperEdge, class HyperNode>
void to_two_graph_block_range2d(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs, 
HyperEdge& edges, HyperNode& nodes, int num_bins, 
size_t rowrange_begin, size_t rowrange_end, size_t row_grainsize,
size_t colrange_begin, size_t colrange_end, size_t col_grainsize,
std::vector<index_t>& hyperedgedegrees, size_t s) {
    //std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    auto M = edges.size();
    tbb::parallel_for(tbb::blocked_range2d<vertex_id_t>(rowrange_begin, rowrange_end, row_grainsize, colrange_begin, colrange_end, col_grainsize), 
    [&](tbb::blocked_range2d<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto i = r.cols().begin(), ie = r.cols().end(); i != ie; ++i) {
        for (auto j = r.rows().begin(), e = r.rows().end(); j != e; ++j) {
            auto hyperE = num_bins * j + i;
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
                        edges[anotherhyperE].begin(), edges[anotherhyperE].end(), s)) 
                        two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
                }//each neighbor of hyperN
            }//each neighbor of hyperE
        } //for cols
      }//for rows
     
    }, tbb::auto_partitioner());
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses cyclic neighbor range as workload distribution strategy. 
* This version is for s = 1.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <class HyperEdge, class HyperNode>
void to_two_graph_cyclic(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  tbb::parallel_for(
      nw::graph::cyclic_neighbor_range(edges, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          std::fill(visitedE.begin(), visitedE.end(), false);
          // all neighbors of hyperedges are hypernode
          for (auto&& [hyperN] : hyperE_ngh) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperE >= anotherhyperE) continue;
              if (visitedE[anotherhyperE])
                continue;
              else
                visitedE[anotherhyperE] = true;
              two_graphs[worker_index].push_back(
                  std::make_pair<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
            }
          }
        }
      },
      tbb::auto_partitioner());
}
/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses cyclic neighbor range as workload distribution strategy. 
* This version is for s > 1.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] num_bins the number of bins to divide the workload
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
*
*/
template <class HyperEdge, class HyperNode>
void to_two_graph_cyclic(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, std::vector<index_t>& hyperedgedegrees, size_t s) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  tbb::parallel_for(
      nw::graph::cyclic_neighbor_range(edges, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          if (hyperedgedegrees[hyperE] < s) continue;
          std::fill(visitedE.begin(), visitedE.end(), false);
          // all neighbors of hyperedges are hypernode
          for (auto&& [hyperN] : hyperE_ngh) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              // so we check compid of each hyperedge
              // travese upper triangluar with lhs > rhs
              // avoid self edge with lhs == rhs
              if (hyperE >= anotherhyperE) continue;
              // filter edges deg(e) < s
              if (hyperedgedegrees[anotherhyperE] < s) continue;
              // avoid duplicate intersections
              if (visitedE[anotherhyperE])
                continue;
              else
                visitedE[anotherhyperE] = true;
              // O(average degree of hyperedges)
              if (is_intersection_size_s(edges[hyperE].begin(),
                                         edges[hyperE].end(),
                                         edges[anotherhyperE].begin(),
                                         edges[anotherhyperE].end(), s))
                two_graphs[worker_index].push_back(
                    std::make_pair<vertex_id_t, vertex_id_t>(
                        std::forward<vertex_id_t>(hyperE),
                        std::forward<vertex_id_t>(anotherhyperE)));
            }  // each neighbor of hyperN
          }    // each neighbor of hyperE
        }
      },
      tbb::auto_partitioner());
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* clean without counter. All features on. 
* For various bin sizes experiments only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] max_bins the maximum number of bins to experiment with
* @returns edge list of 1-line graph or an empty edge list if s > 1
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_vary_size (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int max_bins = 1) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  if (1 == s) {
    linegraph_t two_graphs(num_threads);
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, max_bins, 0, M);
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    auto M = edges.size();
    //block size: 1, 4, 16, 64, 256, 1024, ... until reach 4^num_bins
    for (size_t i = 1; max_bins > 0; i *= 4, --max_bins) {
      linegraph_t two_graphs(num_threads);
      size_t num_bins = M / i; // num_bins: M, M/4, M/16, M/256, ...
      to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, 0, M, hyperedgedegrees, s); 
    }
    //abandon results
    return nw::graph::edge_list<edge_directedness>(0);
  }//else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Clean without counter. All features on. Fastest version.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze(two_graphs);
  } else {
    // when s > 1
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M, hyperedgedegrees, s);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  }  // else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses 2D blocked range as workload distribution strategy. 
* Clean without counter. All features on. 
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_2d (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();

  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
        to_two_graph_block_range2d(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins,
        0, M / num_bins, num_bins, 0, num_bins, num_bins);
        to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, M / num_bins * num_bins, M);
    }
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_threads), [&](auto i) {
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
    {
        nw::util::life_timer _(__func__);
        to_two_graph_block_range2d(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins,
        0, M / num_bins, num_bins, 0, num_bins, num_bins,
        hyperedgedegrees, s);
        to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, 
        M / num_bins * num_bins, M, hyperedgedegrees, s);
    }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Clean without counter. All features on. Fastest version.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_cyclic (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins);
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, hyperedgedegrees, s);
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}



/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Counter version. All features on. For benchmark purpose.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_cyclic_with_counter (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();

  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
      nw::util::life_timer _(__func__);
      std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
      tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE >= anotherhyperE) continue;
            if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;  
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
        }
      }, tbb::auto_partitioner());
      std::cout << "#visits for each thread:" << std::endl;
      for (auto &v : num_visits)
        std::cout << v << " ";
      std::cout << std::endl;
      std::cout << "#edges for each thread:" << std::endl;
      for (auto &v : num_edges)
        std::cout << v << " ";
      std::cout << std::endl;
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    {
      nw::util::life_timer _(__func__);
      std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
      tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          if (hyperedgedegrees[hyperE] < s) continue;
          std::fill(visitedE.begin(), visitedE.end(), false);
            //all neighbors of hyperedges are hypernode
          for (auto &&[hyperN] : hyperE_ngh) {
            for (auto &&[anotherhyperE] : nodes[hyperN]) {
              //so we check compid of each hyperedge        
              //travese upper triangluar with lhs > rhs
              //avoid self edge with lhs == rhs
              ++num_visits[worker_index];
              if (hyperE >= anotherhyperE) continue;
              //filter edges deg(e) < s
              if (hyperedgedegrees[anotherhyperE] < s) continue;
                    //avoid duplicate intersections
              if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;         
                    //O(average degree of hyperedges)
              if (is_intersection_size_s(edges[hyperE].begin(), edges[hyperE].end(),
                edges[anotherhyperE].begin(), edges[anotherhyperE].end(), s)) {
                ++num_edges[worker_index];
                two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
              }
            }//each neighbor of hyperN
          }//each neighbor of hyperE
        }
      }, tbb::auto_partitioner());
      std::cout << "#visits for each thread:" << std::endl;
      for (auto &v : num_visits)
        std::cout << v << " ";
      std::cout << std::endl;
      std::cout << "#edges for each thread:" << std::endl;
      for (auto &v : num_edges)
        std::cout << v << " ";
      std::cout << std::endl;
    }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

/**
* Efficient computation of a weighted s-line graph of a hypergraph. 
* Store s value as edge weight.
* It uses blocked range as workload distribution strategy. 
* Clean version. All heuristics on.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam T the type of the weight in the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of each bin after dividing the workload
* @returns the edge list of the weighted s-line graph
*/
template<directedness edge_directedness = undirected, class T, class HyperEdge, class HyperNode>
auto to_weighted_two_graph_efficient_blocked (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();
  //create an array of line graphs for each thread
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_weighted_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze<undirected, T>(two_graphs);
  }
  else {
    //when s > 1
    to_weighted_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M, hyperedgedegrees, s);
    return create_edgelist_with_squeeze<undirected, T>(two_graphs);
  }//else
}


/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Clean without counter, optional feature.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] features a fixed-size sequence of 8 bits to control which heuristics are enabled
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of each bin after dividing the workload
* @returns the edge list of the weighted s-line graph
*/
template <directedness edge_directedness = undirected, class HyperEdge,
          class HyperNode>
auto to_two_graph_efficient_blocked_optional_features_clean(
    std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads,
    int bin_size = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
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
            if (features[heuristics::SKIP_VISITED] && visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
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
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Counter version. All features on. For benchmark purpose.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_with_counter (HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
    //avoid intersection when s=1
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            ++num_visits[worker_index];
            if (hyperE >= anotherhyperE) return;
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_threads), [&](auto i) {
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
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
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
            ++num_visits[worker_index];
            
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;
            //filter edges deg(e) < s
            if (hyperedgedegrees[anotherhyperE] < s) return;
            //avoid duplicate intersections
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;
            
            //O(average degree of hyperedges)
            //short circuit
            if (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) {
              ++num_edges[worker_index];
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}


/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Counter version, optional feature. For benchmark only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] features a fixed-size sequence of 8 bits to control which heuristics are enabled
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of each bin after dividing the workload
* @returns the edge list of the weighted s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_optional_features_with_counter(std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nedges = 0;
    //avoid intersection when s=1
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
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
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_threads), [&](auto i) {
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
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nvisited = 0, nbelowS = 0, nduplicates = 0, nintersections = 0, nedges = 0;
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
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
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Clean without counter. All features on. Fastest version. Do not squeeze.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_without_sequeeze(
  HyperEdge& edges, 
  HyperNode& nodes, 
  std::vector<index_t>& hyperedgedegrees, 
  size_t s,
  int num_threads, 
  int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes,
                         M / bin_size, 0, M);
    return two_graphs;
  } else {
    // when s > 1
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes,
                         M / bin_size, 0, M, hyperedgedegrees, s);
    return two_graphs;
  }  // else
}

/**
* Efficient computation of a weighted s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Clean without counter. All features on. Without squeeze.
* Store the number of overlapping vertices as edge weight.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam T the type of the weight in the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the weighted s-line graph
*/
template<directedness edge_directedness = nw::graph::undirected, class T, class HyperEdge, class HyperNode>
auto to_weighted_two_graph_efficient_blocked_without_squeeze(
  HyperEdge& edges, 
  HyperNode& nodes, 
  std::vector<index_t>& hyperedgedegrees, 
  size_t s,
  int num_threads, 
  int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_weighted_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M);
    return two_graphs;
  }
  else {
    //when s > 1
    to_weighted_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M, hyperedgedegrees, s);
    return two_graphs;
  }//else
}

/**
* Portal of efficient computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] e_nbs adjacency for hyperedges
* @param[in] n_nbs adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template <directedness edge_directedness = undirected, class HyperEdge,
          class HyperNode>
auto to_two_graph_efficient_cyclic_portal(
    bool verbose, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<index_t>& hyperedgedegrees, size_t s, int numb_threads,
    int num_bins = 32) {
  if (!verbose)
    return to_two_graph_efficient_parallel_cyclic(
        e_nbs, n_nbs, hyperedgedegrees, s, numb_threads, num_bins);
  else
    return to_two_graph_efficient_parallel_cyclic_with_counter(
        e_nbs, n_nbs, hyperedgedegrees, s, numb_threads, num_bins);
}

/**
* Portal of efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] features a fixed-size sequence of 8 bits to control which heuristics are enabled
* @param[in] e_nbs adjacency for hyperedges
* @param[in] n_nbs adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the s-line graph
*/
template <directedness edge_directedness = undirected, class HyperEdge,
          class HyperNode>
auto to_two_graph_efficient_blocked_portal(
    bool verbose, std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads,
    int bin_size = 32) {
  // if none feature selected, then default to turn all on
  if (features.none()) features[heuristics::ALL] = true;
  if (features[heuristics::NONE]) {
    // if none is set, manully override each heuristic to false
    for (int i = heuristics::DEG_PRUNING; i != heuristics::NONE; ++i) {
      features[i] = false;
    }
  }

  if (features[heuristics::ALL]) {
    if (!verbose)
      return to_two_graph_efficient_blocked(e_nbs, n_nbs, hyperedgedegrees, s,
                                            num_threads, bin_size);
    else {
      print_heuristics_verbose(heuristics::ALL);
      return to_two_graph_efficient_blocked_with_counter(
          e_nbs, n_nbs, hyperedgedegrees, s, num_threads, bin_size);
    }
  } else {
    if (!verbose) {
      // clean version without counter
      for (int i = heuristics::DEG_PRUNING; i <= heuristics::NONE; ++i) {
        if (features[i]) print_heuristics_verbose((heuristics)i);
      }
      return to_two_graph_efficient_blocked_optional_features_clean(
          features, e_nbs, n_nbs, hyperedgedegrees, s, num_threads, bin_size);
    } else {
      // with counter
      for (int i = heuristics::DEG_PRUNING; i <= heuristics::NONE; ++i) {
        if (features[i]) print_heuristics_verbose((heuristics)i);
      }
      return to_two_graph_efficient_blocked_optional_features_with_counter(
          features, e_nbs, n_nbs, hyperedgedegrees, s, num_threads, bin_size);
    }
  }
}

}//namespace hypergraph
}//namespace nw