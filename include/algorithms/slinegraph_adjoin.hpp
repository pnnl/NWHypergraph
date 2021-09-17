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

#include "util/slinegraph_helper.hpp"
#include <tbb/task_arena.h> //for tbb::this_task_arena::current_thread_index()
#include <unordered_set>
#include <set>
#include <unordered_map>
#include <map>

namespace nw {
namespace hypergraph {

/*
 * This to_two_graph operates on hypergraph based frontier.
 * The hypergraph can be relabeled by degree.
 * */
template <directedness edge_directedness = undirected, class ExecutionPolicy,
          class Hypergraph, class HypergraphT>
auto to_two_graph_map_frontier_blocked(ExecutionPolicy&& ep, Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s = 1,
                               int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(
      num_bins);
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0ul, frontier.size()),
        [&](const tbb::blocked_range<size_t>& r) {
          int worker_index = tbb::this_task_arena::current_thread_index();
          for (size_t i = r.begin(), e = r.end(); i != e; ++i) {
            auto hyperE = frontier[i];
            if (degrees[hyperE] < s) continue;
            std::map<size_t, size_t> overlaps;
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
  } else {
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          tbb::blocked_range<size_t>(0ul, frontier.size()),
          [&](const tbb::blocked_range<size_t>& r) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (size_t i = r.begin(), e = r.end(); i != e; ++i) {
              auto hyperE = frontier[i];
              std::set<size_t> overlaps;
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    // do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0),
                  nw::graph::counting_iterator<int>(num_bins), [&](auto i) {
                    std::for_each(two_graphs[i].begin(), two_graphs[i].end(),
                                  [&](auto&& e) { result.push_back(e); });
                  });
    result.close_for_push_back();
    return result;
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/*
 * This to_two_graph operates on hypergraph based frontier.
 * The hypergraph can be relabeled by degree.
 * */
template <directedness edge_directedness = undirected, class ExecutionPolicy,
          class Hypergraph, class HypergraphT>
auto to_two_graph_map_frontier_cyclic(ExecutionPolicy&& ep, Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s = 1,
                               int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(
      num_bins);
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(
        nw::graph::cyclic(frontier, num_bins),
        [&](auto& i) {
          int worker_index = tbb::this_task_arena::current_thread_index();
          for (auto&& j = i.begin(); j != i.end(); ++j) {
            auto hyperE = frontier[*j];
            if (degrees[hyperE] < s) continue;
            std::map<size_t, size_t> K;
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
  } else {
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          nw::graph::cyclic(frontier, num_bins),
          [&](auto& i) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (auto&& j = i.begin(); j != i.end(); ++j) {
              auto hyperE = frontier[*j];
              std::set<size_t> overlaps;
              /*
              std::vector<size_t> overlaps(N, 0);
              nsoverlaps
              set, insertion is log(nsoverlap)
              total insertion is nsoverlap*log(nsoverlap),
              traverse over set is nsoverlap*log(nsoverlap)
              total is 2*nsoverlap*log(nsoverlap),

              vector, insertion is 1
              total insertion is nsoverlap
              traverse over vector is N
              nsoverlap*2log(nsoverlap) vs nsoverlap + N
              
              as s increase, nsoverlap drops.
              when s is large, nsoverlap  << N
              when s = 1, nsoverlap is the average degree < N

              hyperE nsoverlap << nsoverlap + N
              */
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    // do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0),
                  nw::graph::counting_iterator<int>(num_bins), [&](auto i) {
                    std::for_each(two_graphs[i].begin(), two_graphs[i].end(),
                                  [&](auto&& e) { result.push_back(e); });
                  });
    result.close_for_push_back();
    return result;
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/*
* clean without counter. All heuristics on. Fastest version. Using blocked_range.
*/
template <directedness edge_directedness = undirected, class ExecutionPolicy,
          class Hypergraph, class HypergraphT>
auto to_two_graph_efficient_frontier_blocked(ExecutionPolicy&& ep, Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s = 1,
                               int num_bins = 32) {
  size_t M = frontier.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_bins);
  if (1 == s) {
    //avoid intersection when s=1
    {
        nw::util::life_timer _(__func__);
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0ul, M),
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          tbb::blocked_range<size_t>(0ul, M),
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
                  if (is_intersection_size_s(h[hyperE].begin(),
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
template <directedness edge_directedness = undirected, class ExecutionPolicy,
          class Hypergraph, class HypergraphT>
auto to_two_graph_efficient_frontier_cyclic(ExecutionPolicy&& ep, Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s = 1,
                               int num_bins = 32) {
  size_t M = frontier.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_bins);
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
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
                  if (is_intersection_size_s(h[hyperE].begin(),
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

/*
 * This to_two_graph operates on hypergraph based frontier. Use an unordered_map instead of map.
 * The hypergraph can be relabeled by degree.
 * */
template <directedness edge_directedness = undirected, class ExecutionPolicy,
          class Hypergraph, class HypergraphT>
auto to_two_graph_hashmap_frontier_blocked(ExecutionPolicy&& ep, Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s = 1,
                               int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(
      num_bins);
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(
        tbb::blocked_range(0ul, frontier.size()),
        [&](auto& i) {
          int worker_index = tbb::this_task_arena::current_thread_index();
          for (auto j = i.begin(); j != i.end(); ++j) {
            auto hyperE = frontier[j];
            if (degrees[hyperE] < s) continue;
            std::unordered_map<size_t, size_t> K;
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
  } else {
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          tbb::blocked_range(0ul, frontier.size()),
          [&](auto& i) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (auto j = i.begin(); j != i.end(); ++j) {
              auto hyperE = frontier[j];
              std::unordered_set<size_t> overlaps;
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    // do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0),
                  nw::graph::counting_iterator<int>(num_bins), [&](auto i) {
                    std::for_each(two_graphs[i].begin(), two_graphs[i].end(),
                                  [&](auto&& e) { result.push_back(e); });
                  });
    result.close_for_push_back();
    return result;
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/*
 * This to_two_graph operates on hypergraph based frontier. Use an unordered_map instead of map.
 * The hypergraph can be relabeled by degree.
 * */
template <directedness edge_directedness = undirected, class ExecutionPolicy,
          class Hypergraph, class HypergraphT>
auto to_two_graph_hashmap_frontier_cyclic(ExecutionPolicy&& ep, Hypergraph& h,
                               HypergraphT& ht, std::vector<index_t>& degrees,
                               std::vector<vertex_id_t>& frontier, size_t s = 1,
                               int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(
      num_bins);
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(
        nw::graph::cyclic(frontier, num_bins),
        [&](auto& i) {
          int worker_index = tbb::this_task_arena::current_thread_index();
          for (auto&& j = i.begin(); j != i.end(); ++j) {
              auto hyperE = frontier[*j];
            if (degrees[hyperE] < s) continue;
            std::unordered_map<size_t, size_t> K;
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
  } else {
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          nw::graph::cyclic(frontier, num_bins),
          [&](auto& i) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (auto&& j = i.begin(); j != i.end(); ++j) {
              auto hyperE = frontier[*j];
              std::unordered_set<size_t> overlaps;
              /*
              std::vector<size_t> overlaps(N, 0);
              nsoverlaps
              set, insertion is log(nsoverlap)
              total insertion is nsoverlap*log(nsoverlap),
              traverse over set is nsoverlap*log(nsoverlap)
              total is 2*nsoverlap*log(nsoverlap),

              vector, insertion is 1
              total insertion is nsoverlap
              traverse over vector is N
              nsoverlap*2log(nsoverlap) vs nsoverlap + N
              
              as s increase, nsoverlap drops.
              when s is large, nsoverlap  << N
              when s = 1, nsoverlap is the average degree < N

              hyperE nsoverlap << nsoverlap + N
              */
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    // do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0),
                  nw::graph::counting_iterator<int>(num_bins), [&](auto i) {
                    std::for_each(two_graphs[i].begin(), two_graphs[i].end(),
                                  [&](auto&& e) { result.push_back(e); });
                  });
    result.close_for_push_back();
    return result;
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/*
* Portal for frontier based map soverlap computation algorithms.
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
auto to_two_graph_map_frontier_cyclic_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s = 1, int num_bins = 32) {
  //frontier will contain the number of edges no matter what the ids of them become (with adjoin or relabel)
  std::vector<vertex_id_t> frontier; 
  if (iperm.empty()) {
    //without relabel by degree
    size_t start, end;
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, without adjoin or relabel by degree
      start = 0ul;
      end = h.size();
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
      auto end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, with relabel by degree
      size_t start, end;
      if (nrealedges < nrealnodes) {
        start = nrealnodes;
        end = h.size();
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
  return to_two_graph_map_frontier_cyclic(ep, h, ht, degrees, frontier, s, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_map_frontier_blocked_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s = 1, int num_bins = 32) {
  //frontier will contain the number of edges no matter what the ids of them become (with adjoin or relabel)
  std::vector<vertex_id_t> frontier; 
  if (iperm.empty()) {
    //without relabel by degree
    size_t start, end;
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, without adjoin or relabel by degree
      start = 0ul;
      end = h.size();
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
      auto end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, with relabel by degree
      size_t start, end;
      if (nrealedges < nrealnodes) {
        start = nrealnodes;
        end = h.size();
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
  return to_two_graph_map_frontier_blocked(ep, h, ht, degrees, frontier, s, num_bins);
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
size_t s = 1, int num_bins = 32) {
  //frontier will contain the number of edges no matter what the ids of them become (with adjoin or relabel)
  std::vector<vertex_id_t> frontier; 
  if (iperm.empty()) {
    //without relabel by degree
    size_t start, end;
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, without adjoin or relabel by degree
      start = 0ul;
      end = h.size();
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
      auto end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, with relabel by degree
      size_t start, end;
      if (nrealedges < nrealnodes) {
        start = nrealnodes;
        end = h.size();
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
  return to_two_graph_efficient_frontier_cyclic(ep, h, ht, degrees, frontier, s, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class Hypergraph, class HypergraphT>
auto to_two_graph_efficient_frontier_blocked_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s = 1, int num_bins = 32) {
  //frontier will contain the number of edges no matter what the ids of them become (with adjoin or relabel)
  std::vector<vertex_id_t> frontier; 
  if (iperm.empty()) {
    //without relabel by degree
    size_t start, end;
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, without adjoin or relabel by degree
      start = 0ul;
      end = h.size();
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
      auto end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, with relabel by degree
      size_t start, end;
      if (nrealedges < nrealnodes) {
        start = nrealnodes;
        end = h.size();
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
  return to_two_graph_efficient_frontier_blocked(ep, h, ht, degrees, frontier, s, num_bins);
}

/*
* Portal for frontier based hashmap soverlap computation algorithms.
* It replace a map with an unordered_map.
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
auto to_two_graph_hashmap_frontier_blocked_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s = 1, int num_bins = 32) {
  std::vector<vertex_id_t> frontier; //frontier will contain the number of edges no matter the ids of it
  if (iperm.empty()) {
    //without relabel by degree
    size_t start, end;
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, without adjoin or relabel by degree
      start = 0ul;
      end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, without relabel by degree
      if (nrealedges > nrealnodes) {
        start = 0ul;
        end = nrealedges;
      } else {
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
      auto end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, with relabel by degree
      size_t start, end;
      if (nrealedges < nrealnodes) {
        start = nrealnodes;
        end = h.size();
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
  return to_two_graph_hashmap_frontier_blocked(ep, h, ht, degrees, frontier, s, num_bins);
}

/*
* Portal for frontier based hashmap soverlap computation algorithms.
* It replace a map with an unordered_map.
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
auto to_two_graph_hashmap_frontier_cyclic_portal(ExecutionPolicy&& ep, Hypergraph& h, HypergraphT& ht,
std::vector<index_t>& degrees, std::vector<vertex_id_t>& iperm,
size_t nrealedges, size_t nrealnodes, 
size_t s = 1, int num_bins = 32) {
  std::vector<vertex_id_t> frontier; //frontier will contain the number of edges no matter the ids of it
  if (iperm.empty()) {
    //without relabel by degree
    size_t start, end;
    if (0 == nrealnodes && 0 == nrealedges) {
      // for original hypergraph, without adjoin or relabel by degree
      start = 0ul;
      end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, without relabel by degree
      if (nrealedges > nrealnodes) {
        start = 0ul;
        end = nrealedges;
      } else {
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
      auto end = h.size();
      frontier.resize(end - start);
      std::copy(ep, counting_iterator<vertex_id_t>(start), counting_iterator<vertex_id_t>(end), frontier.begin());
    }
    else {
      // for adjoin hypergraph, with relabel by degree
      size_t start, end;
      if (nrealedges < nrealnodes) {
        start = nrealnodes;
        end = h.size();
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
  return to_two_graph_hashmap_frontier_cyclic(ep, h, ht, degrees, frontier, s, num_bins);
}


}//namespace hypergraph
}//namespace nw