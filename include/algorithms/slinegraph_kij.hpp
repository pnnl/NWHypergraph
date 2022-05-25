/**
 * @file slinegraph_kij.hpp
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
#include <nwgraph/adaptors/vertex_range.hpp>

#include "util/slinegraph_helper.hpp"
#include <tbb/task_arena.h> //for tbb::this_task_arena::current_thread_index()

namespace nw {
namespace hypergraph {

/**
* Spgemm-kij style soverlap computation algorithm.
* Here the outer loop iterates over vertices.
* For each vertex, check its eligible incident hyperedges.
* Accumulate the overlap of eligible hyperedge pairs.
* At the end, generate edge list of s-line graph based on overlap information.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution policy
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] ep execution policy
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns edge list of the s-line graph
*/
template <directedness edge_directedness = nw::graph::directedness::undirected, class ExecutionPolicy,
          class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_spgemm_kij_cyclic(ExecutionPolicy&& ep, HyperEdge& edges,
                                    HyperNode& nodes,
                                    std::vector<vertex_id_t>& hyperedgedegrees,
                                    size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();
  std::vector<std::map<size_t, size_t>> soverlap(M, std::map<size_t, size_t>());
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
      nw::util::life_timer _(__func__);
      efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs),
                                     edges, nodes, num_bins);
    }
    return create_edgelist_without_squeeze(two_graphs);
  } else {
    // when s > 1
    // create an array of line graphs for each thread
    {
      nw::util::life_timer _(__func__);
      std::for_each(
          ep, nw::graph::counting_iterator<vertex_id_t>(0ul),
          nw::graph::counting_iterator<vertex_id_t>(N), [&](auto hyperN) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            // add all the eligible hyperedges to frontier
            std::vector<vertex_id_t> frontier;
            for (auto&& [hyperE] : nodes[hyperN]) {
              if (hyperedgedegrees[hyperE] >= s) frontier.push_back(hyperE);
            }  // each neighbor of hyperE
            size_t n = frontier.size();
            if (1 >= n) return;
            //here i, j implies that i < j 
            //since the neighbor list of hyperN are sorted in ascending order
            for (size_t i = 0; i < n - 1; ++i) {
              auto hyperE = frontier[i];
              for (size_t j = i + 1; j < n; ++j) {
                auto anotherhyperE = frontier[j];
                if (soverlap[hyperE][anotherhyperE] < s)
                  ++soverlap[hyperE][anotherhyperE];
              }  // each neighbor of hyperN
            }
          });
      //filter overlap information based on s value
      std::for_each(
          ep, nw::graph::counting_iterator<vertex_id_t>(0ul),
          nw::graph::counting_iterator<vertex_id_t>(M), [&](auto hyperE) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (auto&& [anotherhyperE, val] : soverlap[hyperE]) {
              if (s <= val) {
                two_graphs[worker_index].push_back(
                    std::make_pair<vertex_id_t, vertex_id_t>(
                        std::forward<vertex_id_t>(hyperE),
                        std::forward<vertex_id_t>(anotherhyperE)));
              }
            }
          });
    }
    return create_edgelist_with_squeeze(two_graphs);
  }  // else
}

}//namespace hypergraph
}//namespace nw