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
#include "to_two_graph_efficient.hpp"
#include "experimental/slinegraph_efficient.hpp"
#include "util/slinegraph_helper.hpp"

namespace nw {
namespace hypergraph {

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
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_blocked (HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze(two_graphs);
  } else {
    // when s > 1
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
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
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_parallel_2d (HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();

  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
        efficient::to_two_graph_block_range2d(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins,
        0, M / num_bins, num_bins, 0, num_bins, num_bins);
        efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, M / num_bins * num_bins, M);
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    {
        nw::util::life_timer _(__func__);
        efficient::to_two_graph_block_range2d(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins,
        0, M / num_bins, num_bins, 0, num_bins, num_bins,
        hyperedgedegrees, s);
        efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, 
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
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_cyclic (HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins);
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, hyperedgedegrees, s);
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
template <directedness edge_directedness = nw::graph::directedness::undirected, class T, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_weighted_two_graph_efficient_blocked(
    HyperEdge& edges, HyperNode& nodes, std::vector<vertex_id_t>& hyperedgedegrees,
    size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();
  // create an array of line graphs for each thread
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    efficient::to_weighted_two_graph_blocked(
        std::forward<linegraph_t>(two_graphs), edges, nodes, M / bin_size, 0,
        M);
    return create_edgelist_without_squeeze<nw::graph::directedness::undirected, vertex_id_t, T>(two_graphs);
  } else {
    // when s > 1
    efficient::to_weighted_two_graph_blocked(
        std::forward<linegraph_t>(two_graphs), edges, nodes, M / bin_size, 0, M,
        hyperedgedegrees, s);
    return create_edgelist_with_squeeze<nw::graph::directedness::undirected, vertex_id_t, T>(two_graphs);
  }  // else
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
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_blocked_without_sequeeze(
  HyperEdge& edges, 
  HyperNode& nodes, 
  std::vector<vertex_id_t>& hyperedgedegrees, 
  size_t s,
  int num_threads, 
  int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes,
                         M / bin_size, 0, M);
    return two_graphs;
  } else {
    // when s > 1
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes,
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
template<directedness edge_directedness = nw::graph::directedness::undirected, class T, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_weighted_two_graph_efficient_blocked_without_squeeze(
  HyperEdge& edges, 
  HyperNode& nodes, 
  std::vector<vertex_id_t>& hyperedgedegrees, 
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
template <directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_cyclic_portal(
    bool verbose, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int numb_threads,
    int num_bins = 32) {
  if (!verbose)
    return to_two_graph_efficient_cyclic(
        e_nbs, n_nbs, hyperedgedegrees, s, numb_threads, num_bins);
  else
    return to_two_graph_efficient_cyclic_with_counter(
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
template <directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_blocked_portal(
    bool verbose, std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads,
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