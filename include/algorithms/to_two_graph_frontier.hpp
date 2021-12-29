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
#include <unordered_set>
#include <unordered_map>
#include <adaptors/cyclic_neighbor_range.hpp>
#include <adaptors/cyclic_range_adapter.hpp>
#include "util/slinegraph_helper.hpp"

namespace nw {
namespace hypergraph {
namespace frontier {
/**
* Frontier-based computation of a s-line graph of a hypergraph for s > 1. 
* Frontier contains either the original ids of the hypergraph or
* the relabeled ids of the adjoin graph.
* It uses blocked range as workload distribution strategy.
* The container can either be a map or a hashmap.
*
* @tparam Container the type of container to store overlaping counts
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <class Container = std::unordered_map<size_t, size_t>, class Hypergraph, class HypergraphT>
void to_two_graph_frontier_hashmap_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    Hypergraph& h, HypergraphT& ht, std::vector<index_t>& degrees,
    std::vector<vertex_id_t>& frontier, size_t s, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  auto M = frontier.size();
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0ul, M, M / num_bins),
      [&](const tbb::blocked_range<size_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (size_t i = r.begin(), e = r.end(); i != e; ++i) {
          auto hyperE = frontier[i];
          if (degrees[hyperE] < s) continue;
          Container overlaps;
          for (auto&& [hyperN] : h[hyperE]) {
            for (auto&& [anotherhyperE] : ht[hyperN]) {
              if (degrees[anotherhyperE] < s) continue;
              if (hyperE < anotherhyperE) ++overlaps[anotherhyperE];
            }
          }
          for (auto& [anotherhyperE, val] : overlaps)
            if (val >= s)
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
        }
      },
      tbb::auto_partitioner());
}

/**
* Frontier-based computation of a s-line graph of a hypergraph for s > 1. 
* Frontier contains either the original ids of the hypergraph or
* the relabeled ids of the adjoin graph.
* It uses cyclic range as workload distribution strategy.
* The container can either be a map or a hashmap.
*
* @tparam Container the type of container to store overlaping counts
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <class Container = std::unordered_map<size_t, size_t>,
          class Hypergraph, class HypergraphT>
void to_two_graph_frontier_hashmap_cyclic(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    Hypergraph& h, HypergraphT& ht, std::vector<index_t>& degrees,
    std::vector<vertex_id_t>& frontier, size_t s, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  auto M = frontier.size();
  tbb::parallel_for(
      nw::graph::cyclic(frontier, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto hyperE = frontier[*j];
          if (degrees[hyperE] < s) continue;
          Container K;
          for (auto&& [hyperN] : h[hyperE]) {
            for (auto&& [anotherhyperE] : ht[hyperN]) {
              if (degrees[anotherhyperE] < s) continue;
              if (hyperE < anotherhyperE) ++K[anotherhyperE];
            }
          }
          for (auto&& [anotherhyperE, val] : K) {
            if (val >= s)
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      },
      tbb::auto_partitioner());
}

/**
* Frontier-based computation of a s-line graph of a hypergraph for s = 1. 
* Frontier contains either the original ids of the hypergraph or
* the relabeled ids of the adjoin graph.
* It uses blocked range as workload distribution strategy.
* The container can either be a set or a unordered_set.
*
* @tparam Container the type of container to store overlaping counts
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <class Container = std::set<size_t>, class Hypergraph, class HypergraphT>
void to_two_graph_frontier_set_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    Hypergraph& h, HypergraphT& ht, std::vector<vertex_id_t>& frontier,
    int num_bins = 32) {
  nw::util::life_timer _(__func__);
  auto M = frontier.size();
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0ul, M, M / num_bins),
      [&](const tbb::blocked_range<size_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (size_t i = r.begin(), e = r.end(); i != e; ++i) {
          auto hyperE = frontier[i];
          Container overlaps;
          for (auto&& [hyperN] : h[hyperE]) {
            for (auto&& [anotherhyperE] : ht[hyperN]) {
              if (hyperE < anotherhyperE) overlaps.insert(anotherhyperE);
            }
          }
          for (auto& anotherhyperE : overlaps)
            two_graphs[worker_index].push_back(
                std::make_tuple<vertex_id_t, vertex_id_t>(
                    std::forward<vertex_id_t>(hyperE),
                    std::forward<vertex_id_t>(anotherhyperE)));
        }
      },
      tbb::auto_partitioner());
}

/**
* Frontier-based computation of a s-line graph of a hypergraph for s = 1. 
* Frontier contains either the original ids of the hypergraph or
* the relabeled ids of the adjoin graph.
* It uses cyclic range as workload distribution strategy.
* The container can either be a set or a unordered_set.
*
* @tparam Container the type of container to store overlaping counts
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <class Container = std::set<size_t>, class Hypergraph,
          class HypergraphT>
void to_two_graph_frontier_set_cyclic(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    Hypergraph& h, HypergraphT& ht, std::vector<vertex_id_t>& frontier,
    int num_bins = 32) {
  nw::util::life_timer _(__func__);
  auto M = frontier.size();
  tbb::parallel_for(
      nw::graph::cyclic(frontier, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto hyperE = frontier[*j];
          Container overlaps;
          for (auto&& [hyperN] : h[hyperE]) {
            for (auto&& [anotherhyperE] : ht[hyperN]) {
              if (hyperE < anotherhyperE) overlaps.insert(anotherhyperE);
            }
          }
          for (auto& anotherhyperE : overlaps)
            two_graphs[worker_index].push_back(
                std::make_tuple<vertex_id_t, vertex_id_t>(
                    std::forward<vertex_id_t>(hyperE),
                    std::forward<vertex_id_t>(anotherhyperE)));
        }
      },
      tbb::auto_partitioner());
}

}//namespace frontier
}//namespace hypergraph
}//namespace nw