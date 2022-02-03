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
#include <numeric>

#include "experimental/slinegraph_map.hpp"
#include "util/slinegraph_helper.hpp"
#include "to_two_graph_efficient.hpp"
#include "to_two_graph_map.hpp"
#include <tbb/task_arena.h>
#include <tbb/blocked_range2d.h>

namespace nw {
namespace hypergraph {

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a map to store overlapping counts.
*
* @tparam edge_directedness the type of edge directedness in the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] bin_size the size of bins after dividing the workload
*
*/
template <directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_map_blocked(HyperEdge& edges, HyperNode& nodes,
                              std::vector<vertex_id_t>& hyperedgedegrees, size_t s,
                              int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  using container_t = std::map<size_t, size_t>;
  if (1 < s) {
    map::to_two_graph_map_blocked<container_t>(std::forward<linegraph_t>(two_graphs), edges,
                                  nodes, hyperedgedegrees, s, bin_size);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs),
                                    edges, nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy.
* Use a map to store overlapping counts.
*
* @tparam edge_directedness the type of edge directedness in the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_bins the number of bins to divide the workload
*
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_map_cyclic(HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  using container_t = std::map<size_t, size_t>;
  if (1 < s) {
    map::to_two_graph_map_cyclic<container_t>(
        std::forward<linegraph_t>(two_graphs), edges, nodes, hyperedgedegrees,
        s, num_bins);
    return create_edgelist_with_squeeze(two_graphs);
  } else {
    efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes,
                                   num_bins);
    return create_edgelist_without_squeeze(two_graphs);
  }
}


/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] bin_size the size of bins after dividing the workload
*
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_hashmap_blocked(HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  using container_t = std::unordered_map<size_t, size_t>;
  if (1 < s) {
    map::to_two_graph_map_blocked<container_t>(
        std::forward<linegraph_t>(two_graphs), edges, nodes, hyperedgedegrees,
        s, bin_size);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  }
  else {
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] bin_size the size of bins after dividing the workload
*
*/
template <directedness edge_directedness = nw::graph::directedness::undirected, class T, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_weighted_two_graph_hashmap_blocked(
    HyperEdge& edges, HyperNode& nodes, std::vector<vertex_id_t>& hyperedgedegrees,
    size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>;
  linegraph_t two_graphs(num_threads);
  using container_t = std::unordered_map<size_t, size_t>;
  if (1 < s) {
    map::to_weighted_two_graph_map_blocked<container_t>(
        std::forward<linegraph_t>(two_graphs), edges, nodes, hyperedgedegrees,
        s, bin_size);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    efficient::to_weighted_two_graph_blocked(
        std::forward<linegraph_t>(two_graphs), edges, nodes, M / bin_size, 0,
        M);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy.
* Using a hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] num_bins the number of bins to divide the workload
*
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_hashmap_cyclic(HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  using container_t = std::unordered_map<size_t, size_t>;
  if (1 < s) {
    map::to_two_graph_map_cyclic<container_t>(
        std::forward<linegraph_t>(two_graphs), edges, nodes, hyperedgedegrees,
        s, num_bins);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  }
  else {
    efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}


/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a vector to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] bin_size the size of bins after dividing the workload
*
*/
template <directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_vector_blocked(HyperEdge& edges, HyperNode& nodes,
                                 std::vector<vertex_id_t>& hyperedgedegrees,
                                 size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 < s) {
    map::to_two_graph_vector_blocked(std::forward<linegraph_t>(two_graphs),
                                     edges, nodes, hyperedgedegrees, s,
                                     num_threads, bin_size);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs),
                                    edges, nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}
/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy.
* Using a vector to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] num_bins the number of bins to divide the workload
*
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_vector_cyclic(HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
    size_t M = edges.size();
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 < s) {
    map::to_two_graph_vector_cyclic(std::forward<linegraph_t>(two_graphs),
                                     edges, nodes, hyperedgedegrees, s,
                                     num_threads, num_bins);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    efficient::to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs),
                                    edges, nodes, num_bins);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses 2D blocked range as workload distribution strategy. 
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
auto to_two_graph_with_map_parallel2d(HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range2d<vertex_id_t>(0, M / N, M / num_bins, 0, N, M / num_bins), [&](tbb::blocked_range2d<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      for (auto i = r.rows().begin(), ie = r.rows().end(); i != ie; ++i) {
        for (auto hyperE = r.cols().begin(), e = r.cols().end(); hyperE != e; ++hyperE) {
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
        }//for cols
      }//for rows
    }, tbb::auto_partitioner());
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range2d<vertex_id_t>(0, M / N, num_bins, 0, N, num_bins), [&](tbb::blocked_range2d<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      for (auto i = r.rows().begin(), ie = r.rows().end(); i != ie; ++i) {
        for (auto hyperE = r.cols().begin(), e = r.cols().end(); hyperE != e; ++hyperE) {
          std::map<size_t, size_t> K;
          for (auto &&[hyperN] : edges[hyperE]) {
            for (auto &&[anotherhyperE] : nodes[hyperN]) {
              if (hyperE < anotherhyperE) ++K[anotherhyperE];
            }
          }
          for (auto &&[anotherhyperE, val] : K) {
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }//for cols
      }//for rows
    }, tbb::auto_partitioner());
    return create_edgelist_without_squeeze(two_graphs);
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a pre-allocated hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution polify for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] ep the execution policy
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] bin_size the size of bins after dividing the workload
* @returns the edge list of the s-line graph
*
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_static_hashmap_blocked(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    std::vector<std::unordered_map<size_t, size_t>> soverlap(num_threads);
    std::for_each(ep, soverlap.begin(), soverlap.end(), [&](auto& K) {
      K.reserve(M);
    });
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        auto& K = soverlap[worker_index];
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          //if soverlap is not used, continue;
          if (0 == val) continue;
          if (val >= s) 
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          //reset soverlap information for next iteration use
          K[anotherhyperE] = 0; //equal to val = 0;
        }
      }
    }, tbb::auto_partitioner());
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false); 
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE >= anotherhyperE) continue;
            //if hyperE-anotherhyperE has been inserted, abort
            //this is to avoid duplicate edges inserted into slinegraph
            if (visitedE[anotherhyperE])
              continue;
            else
              visitedE[anotherhyperE] = true;
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      }
    }, tbb::auto_partitioner());
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy.
* Using a pre-allocated to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution polify for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] ep the execution policy
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of bins to divide the workload
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_static_hashmap_cyclic(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    std::vector<std::unordered_map<size_t, size_t>> soverlap(num_threads);
    std::for_each(ep, soverlap.begin(), soverlap.end(), [&](auto& K) {
      K.reserve(M);
    });
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        if (hyperedgedegrees[hyperE] < s) continue;
        auto& K = soverlap[worker_index];
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (0 == val) continue;
          if (val >= s) 
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          K[anotherhyperE] = 0;
        }
      }
    }, tbb::auto_partitioner());
  }
  else {
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      std::vector<bool> visitedE(M, false); 
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::fill(visitedE.begin(), visitedE.end(), false);
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE >= anotherhyperE) continue;
            //if hyperE-anotherhyperE has been inserted, abort
            //this is to avoid duplicate edges inserted into slinegraph
            if (visitedE[anotherhyperE])
              continue;
            else
              visitedE[anotherhyperE] = true;
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      }
    }, tbb::auto_partitioner());  
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Use a map to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_map_blocked_portal(bool verbose, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_map_blocked(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_map_blocked_with_counter(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Use a map to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_map_cyclic_portal(bool verbose, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_map_cyclic(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_map_cyclic_with_counter(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Use a hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_hashmap_blocked_portal(bool verbose, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_hashmap_blocked(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_hashmap_blocked_with_counter(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Use a hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_hashmap_cyclic_portal(bool verbose, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_hashmap_cyclic(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_hashmap_cyclic_with_counter(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Use a vector to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_vector_blocked_portal(bool verbose, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_vector_blocked(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_vector_blocked_with_counter(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Use a vector to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_vector_cyclic_portal(bool verbose, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_vector_cyclic(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_vector_cyclic_with_counter(edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Use a pre-allocated hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution polify for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bin after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_static_hashmap_blocked_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_static_hashmap_blocked(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_static_hashmap_blocked_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}
/**
* Portal of computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Use a pre-allocated hashmap to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution polify for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] verbose flag to control stats collection
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class ExecutionPolicy, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_static_hashmap_cyclic_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_static_hashmap_cyclic(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_static_hashmap_cyclic_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/**
* Count the overlapping vertices between each hyperedge pair in the hypergraph. 
* It uses cyclic range as workload distribution strategy.
*
* @tparam edge_directedness the type of the edge directedness of the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of vertex incidence
* @param[in] edges adjacency of hyperedges
* @param[in] nodes adjacency of vertices
* @param[in] num_bins the number of bins to divide the workload
* @returns the overlapping information of the hypergraph
*/
template<class HyperEdge, class HyperNode>
auto to_two_graph_count_neighbors_cyclic(HyperEdge& edges, HyperNode& nodes, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  std::vector<std::unordered_map<size_t, size_t>> two_graphs(M, std::unordered_map<size_t, size_t>());
  tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
    int worker_index = tbb::this_task_arena::current_thread_index();    
    for (auto&& j = i.begin(); j != i.end(); ++j) {
      auto&& [hyperE, hyperE_ngh] = *j;
      for (auto&& x : hyperE_ngh) {
        auto hyperN = std::get<0>(x);
        for (auto&& y : nodes[hyperN]) {
          auto anotherhyperE = std::get<0>(y);
          if (hyperE < anotherhyperE) ++two_graphs[hyperE][anotherhyperE];
        }
      }
    }
  }, tbb::auto_partitioner());
  return two_graphs;
}

/**
* Count the overlapping vertices between each hyperedge pair in the hypergraph. 
* It uses cyclic range as workload distribution strategy.
*
* @tparam edge_directedness the type of the edge directedness of the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of vertex incidence
* @param[in] edges adjacency of hyperedges
* @param[in] nodes adjacency of vertices
* @param[in] min_s the smallest s value to filter the overlapping vertices
* @param[in] num_bins the number of bins to divide the workload
* @returns the overlapping information of the hypergraph
*/
template<class HyperEdge, class HyperNode>
auto to_two_graph_count_neighbors_cyclic(HyperEdge& edges, HyperNode& nodes, size_t min_s, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  std::vector<std::unordered_map<size_t, size_t>> two_graphs(M, std::unordered_map<size_t, size_t>());
  tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
    int worker_index = tbb::this_task_arena::current_thread_index();    
    for (auto&& j = i.begin(); j != i.end(); ++j) {
      auto&& [hyperE, hyperE_ngh] = *j;
      if (min_s > hyperE_ngh.size()) continue;
      for (auto&& x : hyperE_ngh) {
        auto hyperN = std::get<0>(x);
        for (auto&& y : nodes[hyperN]) {
          auto anotherhyperE = std::get<0>(y);
          if (hyperE < anotherhyperE) ++two_graphs[hyperE][anotherhyperE];
        }
      }
    }
  }, tbb::auto_partitioner());
  return two_graphs;
}

/**
* Count the overlapping vertices between each hyperedge pair in the hypergraph. 
* It uses blocked range as workload distribution strategy.
*
* @tparam edge_directedness the type of the edge directedness of the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of vertex incidence
* @param[in] edges adjacency of hyperedges
* @param[in] nodes adjacency of vertices
* @param[in] bin_size the size of bins after dividing the workload
* @returns the overlapping information of the hypergraph
*/
template<class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_count_neighbors_blocked(HyperEdge& edges, HyperNode& nodes, int bin_size = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  std::vector<std::unordered_map<size_t, size_t>> two_graphs(M, std::unordered_map<size_t, size_t>());
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::this_task_arena::current_thread_index();  
    for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
      for (auto&& x : edges[hyperE]) {
        auto hyperN = std::get<0>(x);
        for (auto&& y : nodes[hyperN]) {
          auto anotherhyperE = std::get<0>(y);
          if (hyperE < anotherhyperE) ++two_graphs[hyperE][anotherhyperE];
        }
      }
    }
  }, tbb::auto_partitioner());
  return two_graphs;
}

/**
* Count the overlapping vertices between each hyperedge pair in the hypergraph. 
* It uses blocked range as workload distribution strategy.
*
* @tparam edge_directedness the type of the edge directedness of the s-line graph
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of vertex incidence
* @param[in] edges adjacency of hyperedges
* @param[in] nodes adjacency of vertices
* @param[in] min_s the smallest s value to filter the overlapping vertices
* @param[in] bin_size the size of bins after dividing the workload
* @returns the overlapping information of the hypergraph
*/
template<class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_count_neighbors_blocked(HyperEdge& edges, HyperNode& nodes, size_t min_s, int bin_size = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  std::vector<std::unordered_map<size_t, size_t>> two_graphs(M, std::unordered_map<size_t, size_t>());
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::this_task_arena::current_thread_index();  
    for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
      if (min_s > edges[hyperE].size()) continue;
      for (auto&& x : edges[hyperE]) {
        auto hyperN = std::get<0>(x);
        for (auto&& y : nodes[hyperN]) {
          auto anotherhyperE = std::get<0>(y);
          if (hyperE < anotherhyperE) ++two_graphs[hyperE][anotherhyperE];
        }
      }
    }
  }, tbb::auto_partitioner());
  return two_graphs;
}

/**
* Populate a unweighted s-line graph based on the overlapping counts of each hyperedge. 
*
* @tparam edge_directedness the type of the edge directedness of the s-line graph
* @param[in] neighbor_map the overlapping counts of each hyperedge
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected>
auto populate_linegraph_from_neighbor_map(std::vector<std::unordered_map<size_t, size_t>>& neighbor_map, std::size_t s) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> linegraph(0);
  using vertex_id_t = vertex_id_t<decltype(linegraph)>;
  vertex_id_t index = 0;
  std::unordered_map<vertex_id_t, vertex_id_t> relabel_map;
  linegraph.open_for_push_back();
  if (1 < s) {
    //if s > 1, we squeeze the IDs
  for (std::size_t hyperE = 0, e = neighbor_map.size(); hyperE < e; ++hyperE) {
    for (auto &&[anotherhyperE, val] : neighbor_map[hyperE]) {
      if (val >= s) {
        auto x = hyperE;
        auto y = anotherhyperE;
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
        linegraph.push_back(newx, newy);
      }
    }
  }
  }
  else {
    //if 1 = s, no need to squeeze
    for (std::size_t hyperE = 0, e = neighbor_map.size(); hyperE < e; ++hyperE) {
      for (auto &&[anotherhyperE, val] : neighbor_map[hyperE]) {
        if (val >= s) {
          //std::cout << hyperE << "-" << anotherhyperE << std::endl;
          linegraph.push_back(std::make_tuple<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }
  }
  linegraph.close_for_push_back();
  return linegraph;
}

/**
* Populate a weighted/unweighted s-line graph based on the overlapping counts of each hyperedge. 
*
* @tparam edge_directedness the type of the edge directedness of the s-line graph
* @tparam Hypergraph the type of hyperedge incidence
* @tparam Attributes the type of s-line graph weights
* @param[in] g adjacency of hyperedges
* @param[in] neighbor_map the overlapping counts of each hyperedge
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] weighted the flag a weighted/unweighted s-line graph
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class Hypergraph, typename... Attributes>
auto populate_linegraph_from_neighbor_map(Hypergraph& g, std::vector<std::unordered_map<size_t, size_t>>& neighbor_map, 
std::size_t s, bool weighted = false) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness, Attributes...> linegraph;
  using vertex_id_t = vertex_id_t<decltype(linegraph)>;
  size_t M = g.size();
  //n is the number of hyperedges, m is the number of hypernodes
  //time complexity of counting neighbors is same as the efficient: O(n*deg(edges)*deg(nodes)*deg(edges))
  //time complexity of extract slinegraph from the neighbor counts: O(n*deg(edges)) -> worst is O(n^2)
  //space complexity: O(n*total_deg(H)) -> worst is O(n^2)
  //total_deg(H)=sum of deg(each hyperedge)
  //total_deg(H) >> n?
  if (!weighted) {
    for (size_t hyperE = 0; hyperE < M; ++hyperE) {
      for (auto &&[anotherhyperE, val] : neighbor_map[hyperE]) {
        if (val >= s) 
        //std::cout << hyperE << "-" << anotherhyperE << std::endl;
          linegraph.push_back(std::make_tuple<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE), 1));
      }
    }
  }
  else {
    for (size_t hyperE = 0; hyperE < M; ++hyperE) {
      for (auto &&[anotherhyperE, val] : neighbor_map[hyperE]) {
        if (val >= s) 
        //std::cout << hyperE << "-" << anotherhyperE << std::endl;
          linegraph.push_back(std::make_tuple<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE), val));
      }
    }
  }
  linegraph.close_for_push_back(false);
  return linegraph;
}

}//namespace hypergraph
}//namespace nw