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
#include "nwgraph/adaptors/cyclic_neighbor_range.hpp"

#include "util/slinegraph_helper.hpp"
#include "to_two_graph_efficient.hpp"


namespace nw {
namespace hypergraph {

namespace frontier {

/**
* Computation of a s-line graph of a hypergraph for s > 1. 
* It uses blocked range as workload distribution strategy.
* This is based on frontier, which stores all the eligible hyperedge pairs to set intersect.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bins after dividing the workload
*
*/
template <class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
void to_two_graph_frontier_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, std::vector<vertex_id_t>& hyperedgedegrees,
    size_t s, int num_threads, int bin_size = 32) {
  nw::util::life_timer _(__func__);
  auto M = edges.size();
  using frontier_t = std::vector<std::tuple<vertex_id_t, vertex_id_t>>;
  // 1. create thread local frontier and consolidate
  std::vector<frontier_t> thread_local_frontiers(num_threads);
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(0, M, bin_size),
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
              thread_local_frontiers[worker_index].push_back({hyperE, anotherhyperE});       
            }  // each neighbor of hyperN
          }    // each neighbor of hyperE
        }      // for each hyperE
      },
      tbb::auto_partitioner());
  // 2. consolidate thread_local_frontier to a global frontier
  frontier_t frontier;
  std::vector<size_t> size_array(num_threads);
  size_t size = 0;
  for (int i = 0; i < num_threads; ++i) {
    // calculate the size of each thread-local frontier
    size_array[i] = size;
    // accumulate the total size of all thread-local frontiers
    size += thread_local_frontiers[i].size();
  }
  // resize next frontier
  frontier.resize(size);
  std::for_each(std::execution::par_unseq, nw::graph::counting_iterator(0),
                nw::graph::counting_iterator(num_threads), [&](auto i) {
                  // copy each thread-local frontier to next frontier based
                  // on their size offset
                  auto begin = std::next(frontier.begin(), size_array[i]);
                  std::copy(std::execution::par_unseq,
                            thread_local_frontiers[i].begin(),
                            thread_local_frontiers[i].end(), begin);
                  // thread_local_frontiers[i].clear();
                });
  // 3. process the edge pairs in frontier
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, size, bin_size),
      [&](tbb::blocked_range<size_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto i = r.begin(), e = r.end(); i != e; ++i) {
          auto&& [hyperE, anotherhyperE] = frontier[i];
          // O(average degree of hyperedges)
          if (efficient::is_intersection_size_s(
                  edges[hyperE].begin(), edges[hyperE].end(),
                  edges[anotherhyperE].begin(), edges[anotherhyperE].end(),
                  s)) {
            two_graphs[worker_index].push_back(
                std::make_pair<vertex_id_t, vertex_id_t>(
                    std::forward<vertex_id_t>(hyperE),
                    std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      },
      tbb::auto_partitioner());
}

/**
* Computation of a s-line graph of a hypergraph for s > 1. 
* It uses cyclic range as workload distribution strategy.
* This is based on frontier, which stores all the eligible hyperedge pairs to set intersect.
*
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
void to_two_graph_frontier_cyclic(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, std::vector<vertex_id_t>& hyperedgedegrees,
    size_t s, int num_threads, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  auto M = edges.size();
  using frontier_t = std::vector<std::tuple<vertex_id_t, vertex_id_t>>;
  // 1. create thread local frontier and consolidate
  std::vector<frontier_t> thread_local_frontiers(num_threads);
  tbb::parallel_for(
    nw::graph::cyclic_neighbor_range(edges, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          if (hyperedgedegrees[hyperE] < s) continue;
          std::vector<bool> visitedE(M, false);
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
              thread_local_frontiers[worker_index].push_back({hyperE, anotherhyperE});       
            }  // each neighbor of hyperN
          }    // each neighbor of hyperE
        }
      },
      tbb::auto_partitioner());
  // 2. consolidate thread_local_frontier to a global frontier
  frontier_t frontier;
  std::vector<size_t> size_array(num_threads);
  size_t size = 0;
  for (int i = 0; i < num_threads; ++i) {
    // calculate the size of each thread-local frontier
    size_array[i] = size;
    // accumulate the total size of all thread-local frontiers
    size += thread_local_frontiers[i].size();
  }
  // resize next frontier
  frontier.resize(size);
  std::for_each(std::execution::par_unseq, nw::graph::counting_iterator(0),
                nw::graph::counting_iterator(num_threads), [&](auto i) {
                  // copy each thread-local frontier to next frontier based
                  // on their size offset
                  auto begin = std::next(frontier.begin(), size_array[i]);
                  std::copy(std::execution::par_unseq,
                            thread_local_frontiers[i].begin(),
                            thread_local_frontiers[i].end(), begin);
                  // thread_local_frontiers[i].clear();
                });
  // 3. process the edge pairs in frontier
  tbb::parallel_for(
      nw::graph::cyclic(frontier, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, anotherhyperE] = *j;
          // O(average degree of hyperedges)
          if (efficient::is_intersection_size_s(
                  edges[hyperE].begin(), edges[hyperE].end(),
                  edges[anotherhyperE].begin(), edges[anotherhyperE].end(),
                  s)) {
            two_graphs[worker_index].push_back(
                std::make_pair<vertex_id_t, vertex_id_t>(
                    std::forward<vertex_id_t>(hyperE),
                    std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      },
      tbb::auto_partitioner());
}

}//namespace frontier

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* This is based on frontier, which stores all the eligible hyperedge pairs to set intersect.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bins after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_frontier_blocked(HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  } else {
    // when s > 1
    frontier::to_two_graph_frontier_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, hyperedgedegrees, s, num_threads, bin_size);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  }  // else
}

/**
* Computation of a s-line graph of a hypergraph for s > 1. 
* It uses cyclic range as workload distribution strategy.
* This is based on frontier, which stores all the eligible hyperedge pairs to set intersect.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_frontier_cyclic(HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, num_bins);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  } else {
    // when s > 1
    frontier::to_two_graph_frontier_cyclic(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, hyperedgedegrees, s, num_threads, num_bins);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  }  // else
}

}//namespace hypergraph
}//namespace nw