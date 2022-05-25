/**
 * @file slinegraph_frontier.hpp
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
#include "nwgraph/adaptors/cyclic_neighbor_range.hpp"

#include "util/slinegraph_helper.hpp"
#include "to_two_graph_efficient.hpp"
#include "to_two_graph_frontier.hpp"

namespace nw {
namespace hypergraph {

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
template <directedness edge_directedness = nw::graph::directedness::undirected,
          class HyperEdge, class HyperNode,
          class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_frontier_blocked(HyperEdge& edges, HyperNode& nodes,
                                   std::vector<vertex_id_t>& hyperedgedegrees,
                                   size_t s, int num_threads,
                                   int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
      nw::util::life_timer _(__func__);
      efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M);
    }
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  } else {
    {
      nw::util::life_timer _(__func__);
      // when s > 1
      frontier::to_two_graph_frontier_blocked(
          std::forward<linegraph_t>(two_graphs), edges, nodes, hyperedgedegrees,
          s, num_threads, bin_size);
    }
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
    {
      nw::util::life_timer _(__func__);
      efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs),
                                     edges, nodes, num_bins);
    }
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  } else {
    {
      nw::util::life_timer _(__func__);
      // when s > 1
      frontier::to_two_graph_frontier_cyclic(
          std::forward<linegraph_t>(two_graphs), edges, nodes, hyperedgedegrees,
          s, num_threads, num_bins);
    }
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  }  // else
}

}//namespace hypergraph
}//namespace nw