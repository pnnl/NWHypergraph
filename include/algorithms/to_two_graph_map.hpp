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
#include <map>
#include <unordered_map>
#include <adaptors/cyclic_neighbor_range.hpp>
#include <adaptors/vertex_range.hpp>
#include "util/slinegraph_helper.hpp"

namespace nw {
namespace hypergraph {
namespace map {
/**
* Computation of a s-line graph of a hypergraph for s > 1. 
* It uses blocked range as workload distribution strategy.
* The container can either be a map or a hashmap.
*
* @tparam Container the type of container to store overlaping counts
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] bin_size the size of bins after dividing the workload
*
*/
template <class Container = std::unordered_map<size_t, size_t>, class HyperEdge,
          class HyperNode>
void to_two_graph_map_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, std::vector<index_t>& hyperedgedegrees,
    size_t s, int bin_size = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(0, M, bin_size),
      [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
          if (hyperedgedegrees[hyperE] < s) continue;
          Container K;
          for (auto&& [hyperN] : edges[hyperE]) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperedgedegrees[anotherhyperE] < s) continue;
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
* Computation of a s-line graph of a hypergraph for s > 1. 
* It uses cyclic range as workload distribution strategy.
* The container can either be a map or a hashmap.
*
* @tparam Container the type of container to store overlaping counts
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
template <class Container = std::unordered_map<size_t, size_t>, class HyperEdge,
          class HyperNode>
void to_two_graph_map_cyclic(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, std::vector<index_t>& hyperedgedegrees,
    size_t s, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  tbb::parallel_for(
      nw::graph::cyclic_neighbor_range(edges, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          if (hyperedgedegrees[hyperE] < s) continue;
          Container K;
          for (auto&& [hyperN] : hyperE_ngh) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperedgedegrees[anotherhyperE] < s) continue;
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
* Computation of a weighted s-line graph of a hypergraph for s > 1. 
* It uses blocked range as workload distribution strategy.
* The container can either be a map or a hashmap.
*
* @tparam Container the type of container to store overlaping counts
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[out] two_graphs thread local edge list of s-line graph
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] bin_size the size of bins after dividing the workload
*
*/
template <class Container = std::unordered_map<size_t, size_t>, class T, class HyperEdge,
          class HyperNode>
void to_weighted_two_graph_map_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, std::vector<index_t>& hyperedgedegrees,
    size_t s, int bin_size = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(0, M, bin_size),
      [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
          if (hyperedgedegrees[hyperE] < s) continue;
          Container K;
          for (auto&& [hyperN] : edges[hyperE]) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperedgedegrees[anotherhyperE] < s) continue;
              if (hyperE < anotherhyperE) ++K[anotherhyperE];
            }
          }
          for (auto&& [anotherhyperE, val] : K) {
            if (val >= s)
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE), val));
          }
        }
      },
      tbb::auto_partitioner());
}

/**
* Computation of a s-line graph of a hypergraph for s > 1. 
* It uses blocked range as workload distribution strategy.
* The container to store overlapping counts is a pre-allocated vector.
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
template <class HyperEdge, class HyperNode>
void to_two_graph_vector_blocked(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, std::vector<index_t>& hyperedgedegrees,
    size_t s, int num_threads, int bin_size = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  // create a thread-local vector to store the soverlap information
  std::vector<std::vector<size_t>> soverlap(num_threads,
                                            std::vector<size_t>(M, 0));
  tbb::parallel_for(
      tbb::blocked_range<vertex_id_t>(0, M, bin_size),
      [&](tbb::blocked_range<vertex_id_t>& r) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        auto& K = soverlap[worker_index];
        for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
          if (hyperedgedegrees[hyperE] < s) continue;
          // std::fill(K.begin(), K.end(), 0);
          for (auto&& [hyperN] : edges[hyperE]) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperedgedegrees[anotherhyperE] < s) continue;
              if (hyperE < anotherhyperE) ++K[anotherhyperE];
            }
          }
          for (vertex_id_t anotherhyperE = 0; anotherhyperE < M;
               ++anotherhyperE) {
            // if this slot is unused, no need to reset its value to 0
            if (0 == K[anotherhyperE]) continue;
            if (K[anotherhyperE] >= s) {
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
            }
            K[anotherhyperE] = 0;  // set to 0 for next use
          }
        }
      },
      tbb::auto_partitioner());
}

/**
* Computation of a s-line graph of a hypergraph for s > 1. 
* It uses cyclic range as workload distribution strategy.
* The container to store overlapping counts is a pre-allocated vector.
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
template <class HyperEdge, class HyperNode>
void to_two_graph_vector_cyclic(
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>&& two_graphs,
    HyperEdge& edges, HyperNode& nodes, std::vector<index_t>& hyperedgedegrees,
    size_t s, int num_threads, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  // create a thread-local vector to store the soverlap information
  std::vector<std::vector<size_t>> soverlap(num_threads,
                                            std::vector<size_t>(M, 0));
  tbb::parallel_for(
      nw::graph::cyclic_neighbor_range(edges, num_bins),
      [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        auto& K = soverlap[worker_index];
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          if (hyperedgedegrees[hyperE] < s) continue;
          for (auto&& [hyperN] : hyperE_ngh) {
            for (auto&& [anotherhyperE] : nodes[hyperN]) {
              if (hyperedgedegrees[anotherhyperE] < s) continue;
              if (hyperE < anotherhyperE) ++K[anotherhyperE];
            }
          }
          for (vertex_id_t anotherhyperE = 0; anotherhyperE < M;
               ++anotherhyperE) {
            // if this slot is unused, no need to reset its value to 0
            if (0 == K[anotherhyperE]) continue;
            if (K[anotherhyperE] >= s) {
              two_graphs[worker_index].push_back(
                  std::make_tuple<vertex_id_t, vertex_id_t>(
                      std::forward<vertex_id_t>(hyperE),
                      std::forward<vertex_id_t>(anotherhyperE)));
            }
            K[anotherhyperE] = 0;  // set to 0 for next use
          }
        }
      },
      tbb::auto_partitioner());
}

}//namespace map
}//namespace hypergraph
}//namespace nw