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

#include <nwgraph/util/intersection_size.hpp>
#include <nwgraph/adaptors/cyclic_neighbor_range.hpp>
#include <nwgraph/adaptors/vertex_range.hpp>

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

namespace efficient {
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
template <class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
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
template <class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
void to_two_graph_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, size_t range_begin,
    size_t range_end, std::vector<vertex_id_t>& hyperedgedegrees, size_t s) {
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
template <class T, class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
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
          for (auto&& e : edges[hyperE]) {
            auto hyperN = target(edges, e);
            for (auto&& n : nodes[hyperN]) {
              auto anotherhyperE = target(nodes, n);
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
template <class T, class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
void to_weighted_two_graph_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, size_t range_begin,
    size_t range_end, std::vector<vertex_id_t>& hyperedgedegrees, size_t s) {
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
          for (auto&& e : edges[hyperE]) {
            auto hyperN = target(edges, e);
            for (auto&& n : nodes[hyperN]) {
              auto anotherhyperE = target(nodes, n);
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
template<class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
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
template<class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
void to_two_graph_block_range2d(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs, 
HyperEdge& edges, HyperNode& nodes, int num_bins, 
size_t rowrange_begin, size_t rowrange_end, size_t row_grainsize,
size_t colrange_begin, size_t colrange_end, size_t col_grainsize,
std::vector<vertex_id_t>& hyperedgedegrees, size_t s) {
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
template <class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
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
template <class HyperEdge, class HyperNode, class vertex_id_t = nw::graph::vertex_id_t<HyperEdge>>
void to_two_graph_cyclic(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, int num_bins, std::vector<vertex_id_t>& hyperedgedegrees, size_t s) {
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

} //namespace efficient

}//namespace hypergraph
}//namespace nw