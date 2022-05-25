/**
 * @file slinegraph_adjoin.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#pragma once
#include <nwgraph/adaptors/cyclic_range_adaptor.hpp>
#include <tbb/task_arena.h>
#include "algorithms/to_two_graph_efficient.hpp"
#include "util/slinegraph_helper.hpp"

namespace nw {
namespace hypergraph {

/**
* Frontier-based computation of a s-line graph of a hypergraph for s > 1. 
* Frontier contains either the original ids of the hypergraph or
* the relabeled ids of the adjoin graph.
* It uses blocked range as workload distribution strategy.
* The container can either be a map or a hashmap.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of hyperedge incidence
* @tparam HypergraphT the type of hypernode incidence
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns edge list of the s-line graph
*
*/
template <directedness edge_directedness = nw::graph::directedness::undirected,
          class Hypergraph, class HypergraphT, class vertex_id_t = vertex_id_t<Hypergraph>>
auto to_two_graph_efficient_frontier_blocked(Hypergraph& h,
                               HypergraphT& ht, std::vector<vertex_id_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s,
                               int num_threads, int num_bins = 32) {
  size_t M = frontier.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    //avoid intersection when s=1
    {
        nw::util::life_timer _(__func__);
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0ul, M, M / num_bins),
            [&](const tbb::blocked_range<size_t>& r) {
              int worker_index = tbb::this_task_arena::current_thread_index();
              std::vector<bool> visitedE(h.size(), false);
              for (size_t i = r.begin(), e = r.end(); i != e; ++i) {
                auto hyperE = frontier[i];
                std::fill(visitedE.begin(), visitedE.end(), false);
                // all neighbors of hyperedges are hypernode
                for (auto&& [hyperN] : h[hyperE]) {
                  for (auto&& [anotherhyperE] : ht[hyperN]) {
                    if (hyperE >= anotherhyperE) continue;
                    if (visitedE[anotherhyperE])
                      continue;
                    else
                      visitedE[anotherhyperE] = true;
                    two_graphs[worker_index].push_back(
                        std::make_tuple<vertex_id_t, vertex_id_t>(
                            std::forward<vertex_id_t>(hyperE),
                            std::forward<vertex_id_t>(anotherhyperE)));
                  }
                }
              }  // for
            },
            tbb::auto_partitioner());
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          tbb::blocked_range<size_t>(0ul, M, M / num_bins),
          [&](const tbb::blocked_range<size_t>& r) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            std::vector<bool> visitedE(h.size(), false);
            for (size_t i = r.begin(), e = r.end(); i != e; ++i) {
              auto hyperE = frontier[i];
              std::fill(visitedE.begin(), visitedE.end(), false);
              if (degrees[hyperE] < s) continue;
              // all neighbors of hyperedges are hypernode
              for (auto&& [hyperN] : h[hyperE]) {
                for (auto&& [anotherhyperE] : ht[hyperN]) {
                  // so we check compid of each hyperedge
                  // travese upper triangluar with lhs > rhs
                  // avoid self edge with lhs == rhs
                  if (hyperE >= anotherhyperE) continue;
                  // filter edges deg(e) < s
                  if (degrees[anotherhyperE] < s) continue;
                  // avoid duplicate intersections
                  if (visitedE[anotherhyperE])
                    continue;
                  else
                    visitedE[anotherhyperE] = true;
                  // O(average degree of hyperedges)
                  if (efficient::is_intersection_size_s(h[hyperE].begin(),
                                             h[hyperE].end(),
                                             h[anotherhyperE].begin(),
                                             h[anotherhyperE].end(), s)) {
                    two_graphs[worker_index].push_back(
                        std::make_tuple<vertex_id_t, vertex_id_t>(
                            std::forward<vertex_id_t>(hyperE),
                            std::forward<vertex_id_t>(anotherhyperE)));
                  }
                }  // each neighbor of hyperN
              }    // each neighbor of hyperE
            }      // for each hyperE
          },
          tbb::auto_partitioner());
    }
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  }//else
}

/*
* clean without counter. All heuristics on. 
* Operates on Adjoin hypergraph.
*/
template <directedness edge_directedness = nw::graph::directedness::undirected,
          class Hypergraph, class HypergraphT, class vertex_id_t = vertex_id_t<Hypergraph>>
auto to_two_graph_efficient_frontier_cyclic(Hypergraph& h,
                               HypergraphT& ht, std::vector<vertex_id_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s,
                               int num_threads, int num_bins = 32) {
  size_t M = frontier.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          nw::graph::cyclic(frontier, num_bins),
          [&](auto& i) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            std::vector<bool> visitedE(h.size(), false);
            for (auto&& j = i.begin(); j != i.end(); ++j) {
              auto hyperE = frontier[*j];
              std::fill(visitedE.begin(), visitedE.end(), false);
              // all neighbors of hyperedges are hypernode
              for (auto&& [hyperN] : h[hyperE]) {
                for (auto&& [anotherhyperE] : ht[hyperN]) {
                  if (hyperE >= anotherhyperE) continue;
                  if (visitedE[anotherhyperE])
                    continue;
                  else
                    visitedE[anotherhyperE] = true;
                  two_graphs[worker_index].push_back(
                      std::make_tuple<vertex_id_t, vertex_id_t>(
                          std::forward<vertex_id_t>(hyperE),
                          std::forward<vertex_id_t>(anotherhyperE)));
                }
              }
            }
          },
          tbb::auto_partitioner());
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    // when s > 1
    // create an array of line graphs for each thread
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          nw::graph::cyclic(frontier, num_bins),
          [&](auto& i) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            std::vector<bool> visitedE(h.size(), false);
            for (auto&& j = i.begin(); j != i.end(); ++j) {
              auto hyperE = frontier[*j];
              if (degrees[hyperE] < s) continue;
              std::fill(visitedE.begin(), visitedE.end(), false);
              // all neighbors of hyperedges are hypernode
              for (auto&& [hyperN] : h[hyperE]) {
                for (auto&& [anotherhyperE] : ht[hyperN]) {
                  // so we check compid of each hyperedge
                  // travese upper triangluar with lhs > rhs
                  // avoid self edge with lhs == rhs
                  if (hyperE >= anotherhyperE) continue;
                  // filter edges deg(e) < s
                  if (degrees[anotherhyperE] < s) continue;
                  // avoid duplicate intersections
                  if (visitedE[anotherhyperE])
                    continue;
                  else
                    visitedE[anotherhyperE] = true;
                  // O(average degree of hyperedges)
                  if (efficient::is_intersection_size_s(h[hyperE].begin(),
                                             h[hyperE].end(),
                                             h[anotherhyperE].begin(),
                                             h[anotherhyperE].end(), s))
                    two_graphs[worker_index].push_back(
                        std::make_tuple<vertex_id_t, vertex_id_t>(
                            std::forward<vertex_id_t>(hyperE),
                            std::forward<vertex_id_t>(anotherhyperE)));
                }  // each neighbor of hyperN
              }    // each neighbor of hyperE
            }
          },
          tbb::auto_partitioner());
    }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}



}//namespace hypergraph
}//namespace nw