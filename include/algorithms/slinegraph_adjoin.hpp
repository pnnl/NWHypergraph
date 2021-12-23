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
#include <adaptors/cyclic_neighbor_range.hpp>
#include <adaptors/cyclic_range_adapter.hpp>
#include <adaptors/vertex_range.hpp>

#include "to_two_graph_efficient.hpp"
#include "util/slinegraph_helper.hpp"
#include <tbb/task_arena.h> //for tbb::this_task_arena::current_thread_index()
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>

#include "to_two_graph_frontier.hpp"

using namespace nw::hypergraph;

namespace nw {
namespace hypergraph {

/**
* Computation of an s-line graph of a hypergraph. 
* It uses a frontier to contain the workload.
* It uses blocked range as workload distribution strategy.
* Using a map to store overlapping counts for s > 1,
* and a set to avoid duplicate edges for s = 1.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <directedness edge_directedness = undirected,
          class Hypergraph, class HypergraphT>
auto to_two_graph_map_frontier_blocked(Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s,
                               int num_threads, int num_bins = 32) {
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  auto M = frontier.size();
  if (1 < s) {
    using container_t = std::map<size_t, size_t>;
    frontier::to_two_graph_frontier_hashmap_blocked<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, degrees, frontier, s,
        num_bins);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    using container_t = std::set<size_t>;
    frontier::to_two_graph_frontier_set_blocked<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, frontier, num_bins);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/**
* Computation of an s-line graph of a hypergraph. 
* It uses a frontier to contain the workload.
* It uses cyclic range as workload distribution strategy.
* Using a map to store overlapping counts for s > 1,
* and a set to avoid duplicate edges for s = 1.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <directedness edge_directedness = undirected, class Hypergraph,
          class HypergraphT>
auto to_two_graph_map_frontier_cyclic(Hypergraph& h, HypergraphT& ht,
                                      std::vector<index_t>& degrees,
                                      std::vector<vertex_id_t>& frontier,
                                      size_t s, int num_threads,
                                      int num_bins = 32) {
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  auto M = frontier.size();
  if (1 < s) {
    using container_t = std::map<size_t, size_t>;
    frontier::to_two_graph_frontier_hashmap_cyclic<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, degrees, frontier, s,
        num_bins);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    using container_t = std::set<size_t>;
    frontier::to_two_graph_frontier_set_cyclic<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, frontier, num_bins);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/*
* clean without counter. All heuristics on. Fastest version. Using blocked_range.
*/
template <directedness edge_directedness = undirected,
          class Hypergraph, class HypergraphT>
auto to_two_graph_efficient_frontier_blocked(Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
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
template <directedness edge_directedness = undirected,
          class Hypergraph, class HypergraphT>
auto to_two_graph_efficient_frontier_cyclic(Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
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

/**
* Computation of an s-line graph of a hypergraph. 
* It uses a frontier to contain the workload.
* It uses blocked range as workload distribution strategy.
* Using a hashmap to store overlapping counts for s > 1,
* and a hashset to avoid duplicate edges for s = 1.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <directedness edge_directedness = undirected, class Hypergraph,
          class HypergraphT>
auto to_two_graph_hashmap_frontier_blocked(Hypergraph& h, HypergraphT& ht,
                                           std::vector<index_t>& degrees,
                                           std::vector<vertex_id_t>& frontier,
                                           size_t s, int num_threads,
                                           int num_bins = 32) {
  using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  auto M = frontier.size();
  if (1 < s) {
    using container_t = std::unordered_map<size_t, size_t>;
    frontier::to_two_graph_frontier_hashmap_blocked<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, degrees, frontier, s,
        num_bins);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    using container_t = std::unordered_set<size_t>;
    frontier::to_two_graph_frontier_set_blocked<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, frontier, num_bins);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/**
* Computation of an s-line graph of a hypergraph. 
* It uses a frontier to contain the workload.
* It uses cyclic range as workload distribution strategy.
* Using a hashmap to store overlapping counts for s > 1,
* and a hashset to avoid duplicate edges for s = 1.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] frontier the worklist of all the hyperedges to check
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*
*/
template <directedness edge_directedness = undirected,
          class Hypergraph, class HypergraphT>
auto to_two_graph_hashmap_frontier_cyclic(Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s,
                               int num_threads, int num_bins = 32) {
    using linegraph_t =
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  auto M = frontier.size();
  if (1 < s) {
    using container_t = std::unordered_map<size_t, size_t>;
    frontier::to_two_graph_frontier_hashmap_cyclic<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, degrees, frontier, s,
        num_bins);
    return create_edgelist_with_squeeze<edge_directedness>(two_graphs);
  } else {
    using container_t = std::unordered_set<size_t>;
    frontier::to_two_graph_frontier_set_cyclic<container_t>(
        std::forward<linegraph_t>(two_graphs), h, ht, frontier, num_bins);
    return create_edgelist_without_squeeze<edge_directedness>(two_graphs);
  }
}

/**
* Build a frontier based on the number of edges if this is no id permutation
* (meaning frontier contains the original ids of the hypergraph)
* or based on the number of realedges and realnodes
* (meaning frontier contains the relabeled ids of the adjoin graph)
*
* @tparam ExecutionPolicy the type of the execution policy
* @param[in] ep execution policy
* @param[in] iperm permutation array of the original hypergraph
* @param[in] nedges the number of hyperedges in the original hypergraph
* @param[in] nrealedges the number of real (hyper)edges in the adjoin graph
* @param[in] nrealnodes the number of real nodes in the adjoin graph
* @returns a frontier of the hyperedge ids
*/
template <class ExecutionPolicy>
auto build_frontier(ExecutionPolicy& ep, std::vector<vertex_id_t>& iperm,
                    size_t nedges, size_t nrealedges, size_t nrealnodes) {
  //frontier will contain the number of edges no matter what the ids of them become (with adjoin or relabel)
  std::vector<vertex_id_t> frontier; 
  if (iperm.empty()) {
    //without relabel by degree
    size_t start, end;
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, without adjoin or relabel by degree
      start = 0ul;
      end = nedges;
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, without relabel by degree
      if (nrealedges > nrealnodes) {
        //if nrealedges is greater than nrealnodes,
        //then in the adjacency, hyperedges are in the front from index 0 to nrealedges
        start = 0ul;
        end = nrealedges;
      } else {
        //if nrealedges is smaller than nrealnodes,
        // then in the adjacency, hyperedges are in the back from index nrealnodes to nrealnode + nrealedges
        start = nrealnodes;
        end = nrealnodes + nrealedges;
      }
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
  }
  else {
    //with relabel by degree
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, with relabel by degree
      auto start = 0ul;
      auto end = nedges;
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, with relabel by degree
      size_t start, end;
      if (nrealedges < nrealnodes) {
        start = nrealnodes;
        end = nedges;
      } else {
        start = 0ul;
        end = nrealedges;
      }
      frontier.resize(end - start);
      //std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
      std::for_each(ep, nw::graph::counting_iterator<vertex_id_t>(0ul),
                    nw::graph::counting_iterator<vertex_id_t>(end - start),
                    [&](auto i) { frontier[i] = iperm[i + start]; });
    }
  }
  return frontier;
}

/**
* Portal for frontier-based map soverlap computation algorithms using cyclic range.
* It can handle:
* 1) original hypergraph
* 2) original hypergraph with relabeling by degree on either hyperedges/hypernodes
* 3) adjoin hypergraph 
* 4) adjoin hypergraph with relabeling by degree on both hyperedges/hypernodes
* Portal will build frontier using the ids in the Hypergraph 'h'.
* These ids could be different from the origial ids of the hyperedges.
* For adjoin hypergraph, the new ids in h is incremented by the nrealedges/nrealnodes (whichever is larger).
* For orignal hypergraph, the ids are the ids of the hyperedges.
*
* @tparam ExecutionPolicy the type of the execution policy
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] ep execution policy
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] iperm permutation array of the original hypergraph
* @param[in] nrealedges the number of real (hyper)edges in the adjoin graph
* @param[in] nrealnodes the number of real nodes in the adjoin graph
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_map_frontier_cyclic_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s, int num_threads, int num_bins = 32) {
  auto frontier = build_frontier(ep, iperm, h.size(), nrealedges, nrealnodes);
  return to_two_graph_map_frontier_cyclic(h, ht, degrees, frontier, s, num_threads, num_bins);
}

/**
* Portal for frontier-based map soverlap computation algorithms using blocked range.
* It can handle:
* 1) original hypergraph
* 2) original hypergraph with relabeling by degree on either hyperedges/hypernodes
* 3) adjoin hypergraph 
* 4) adjoin hypergraph with relabeling by degree on both hyperedges/hypernodes
* Portal will build frontier using the ids in the Hypergraph 'h'.
* These ids could be different from the origial ids of the hyperedges.
* For adjoin hypergraph, the new ids in h is incremented by the nrealedges/nrealnodes (whichever is larger).
* For orignal hypergraph, the ids are the ids of the hyperedges.
*
* @tparam ExecutionPolicy the type of the execution policy
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] ep execution policy
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] iperm permutation array of the original hypergraph
* @param[in] nrealedges the number of real (hyper)edges in the adjoin graph
* @param[in] nrealnodes the number of real nodes in the adjoin graph
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_map_frontier_blocked_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s, int num_threads, int num_bins = 32) {
  auto frontier = build_frontier(ep, iperm, h.size(), nrealedges, nrealnodes);
  return to_two_graph_map_frontier_blocked(h, ht, degrees, frontier, s, num_threads, num_bins);
}

/*
* Portal for frontier based efficient soverlap computation algorithms.
* It can handle:
* 1) original hypergraph
* 2) original hypergraph with relabeling by degree on either hyperedges/hypernodes
* 3) adjoin hypergraph 
* 4) adjoin hypergraph with relabeling by degree on both hyperedges/hypernodes
* Portal will build frontier using the ids in the Hypergraph 'h'.
* These ids could be different from the origial ids of the hyperedges.
* For adjoin hypergraph, the new ids in h is incremented by the nrealedges/nrealnodes (whichever is larger).
* For orignal
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_efficient_frontier_cyclic_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s, int num_threads, int num_bins = 32) {
  auto frontier = build_frontier(ep, iperm, h.size(), nrealedges, nrealnodes);
  return to_two_graph_efficient_frontier_cyclic(h, ht, degrees, frontier, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_efficient_frontier_blocked_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s, int num_threads, int num_bins = 32) {
  auto frontier = build_frontier(ep, iperm, h.size(), nrealedges, nrealnodes);
  return to_two_graph_efficient_frontier_blocked(h, ht, degrees, frontier, s, num_threads, num_bins);
}

/**
* Portal for frontier-based hashmap soverlap computation algorithms using blocked range.
* It can handle:
* 1) original hypergraph
* 2) original hypergraph with relabeling by degree on either hyperedges/hypernodes
* 3) adjoin hypergraph 
* 4) adjoin hypergraph with relabeling by degree on both hyperedges/hypernodes
* Portal will build frontier using the ids in the Hypergraph 'h'.
* These ids could be different from the origial ids of the hyperedges.
* For adjoin hypergraph, the new ids in h is incremented by the nrealedges/nrealnodes (whichever is larger).
* For orignal hypergraph, the ids are the ids of the hyperedges.
*
* @tparam ExecutionPolicy the type of the execution policy
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] ep execution policy
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] iperm permutation array of the original hypergraph
* @param[in] nrealedges the number of real (hyper)edges in the adjoin graph
* @param[in] nrealnodes the number of real nodes in the adjoin graph
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_hashmap_frontier_blocked_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s, int num_threads, int num_bins = 32) {
  auto frontier = build_frontier(ep, iperm, h.size(), nrealedges, nrealnodes);
  return to_two_graph_hashmap_frontier_blocked(h, ht, degrees, frontier, s, num_threads, num_bins);
}

/**
* Portal for frontier-based hashmap soverlap computation algorithms using cyclic range.
* It can handle:
* 1) original hypergraph
* 2) original hypergraph with relabeling by degree on either hyperedges/hypernodes
* 3) adjoin hypergraph 
* 4) adjoin hypergraph with relabeling by degree on both hyperedges/hypernodes
* Portal will build frontier using the ids in the Hypergraph 'h'.
* These ids could be different from the origial ids of the hyperedges.
* For adjoin hypergraph, the new ids in h is incremented by the nrealedges/nrealnodes (whichever is larger).
* For orignal hypergraph, the ids are the ids of the hyperedges.
*
* @tparam ExecutionPolicy the type of the execution policy
* @tparam edge_directedness the type of the edge directedness
* @tparam Hypergraph the type of the hyperedge incidence, or the adjoin graph
* @tparam HypergraphT the type of the hypernode incidence, or the transpose of adjoin graph
* @param[in] ep execution policy
* @param[in] h adjacency for hyperedges, or the adjacency of the adjoin graph
* @param[in] ht adjacency for hypernodes, or the dual of the adjoin graph
* @param[in] degrees the degrees of hyperedges
* @param[in] iperm permutation array of the original hypergraph
* @param[in] nrealedges the number of real (hyper)edges in the adjoin graph
* @param[in] nrealnodes the number of real nodes in the adjoin graph
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] num_bins the number of bins to divide the workload
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_hashmap_frontier_cyclic_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s, int num_threads, int num_bins = 32) {
  auto frontier = build_frontier(ep, iperm, h.size(), nrealedges, nrealnodes);
  return to_two_graph_hashmap_frontier_cyclic(h, ht, degrees, frontier, s, num_threads, num_bins);
}


}//namespace hypergraph
}//namespace nw