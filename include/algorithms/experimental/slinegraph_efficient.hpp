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
#include "algorithms/to_two_graph_efficient.hpp"
#include "util/slinegraph_helper.hpp"

namespace nw {
namespace hypergraph {

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* clean without counter. All features on. 
* For various bin sizes experiments only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] max_bins the maximum number of bins to experiment with
* @returns edge list of 1-line graph or an empty edge list if s > 1
*/
template<directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge, class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_blocked_vary_size (HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int max_bins = 1) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  if (1 == s) {
    linegraph_t two_graphs(num_threads);
    {
      nw::util::life_timer _(__func__);
      efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs),
                                      edges, nodes, max_bins, 0, M);
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    auto M = edges.size();
    //block size: 1, 4, 16, 64, 256, 1024, ... until reach 4^num_bins
    for (size_t i = 1; max_bins > 0; i *= 4, --max_bins) {
      linegraph_t two_graphs(num_threads);
      size_t num_bins = M / i;  // num_bins: M, M/4, M/16, M/256, ...
      {
        nw::util::life_timer _(__func__);
        efficient::to_two_graph_blocked(std::forward<linegraph_t>(two_graphs),
                                        edges, nodes, num_bins, 0, M,
                                        hyperedgedegrees, s);
      }
    }
    //abandon results
    return nw::graph::edge_list<edge_directedness>(0);
  }//else
}


/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Counter version. All features on. For benchmark purpose.
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
auto to_two_graph_efficient_cyclic_with_counter (HyperEdge& edges, HyperNode& nodes, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();

  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
      nw::util::life_timer _(__func__);
      std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
      tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE >= anotherhyperE) continue;
            if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;  
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
        }
      }, tbb::auto_partitioner());
      std::cout << "#visits for each thread:" << std::endl;
      for (auto &v : num_visits)
        std::cout << v << " ";
      std::cout << std::endl;
      std::cout << "#edges for each thread:" << std::endl;
      for (auto &v : num_edges)
        std::cout << v << " ";
      std::cout << std::endl;
    }
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    {
      nw::util::life_timer _(__func__);
      std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
      tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
          auto&& [hyperE, hyperE_ngh] = *j;
          if (hyperedgedegrees[hyperE] < s) continue;
          std::fill(visitedE.begin(), visitedE.end(), false);
            //all neighbors of hyperedges are hypernode
          for (auto &&[hyperN] : hyperE_ngh) {
            for (auto &&[anotherhyperE] : nodes[hyperN]) {
              //so we check compid of each hyperedge        
              //travese upper triangluar with lhs > rhs
              //avoid self edge with lhs == rhs
              ++num_visits[worker_index];
              if (hyperE >= anotherhyperE) continue;
              //filter edges deg(e) < s
              if (hyperedgedegrees[anotherhyperE] < s) continue;
                    //avoid duplicate intersections
              if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;         
                    //O(average degree of hyperedges)
              if (efficient::is_intersection_size_s(edges[hyperE].begin(), edges[hyperE].end(),
                edges[anotherhyperE].begin(), edges[anotherhyperE].end(), s)) {
                ++num_edges[worker_index];
                two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
              }
            }//each neighbor of hyperN
          }//each neighbor of hyperE
        }
      }, tbb::auto_partitioner());
      std::cout << "#visits for each thread:" << std::endl;
      for (auto &v : num_visits)
        std::cout << v << " ";
      std::cout << std::endl;
      std::cout << "#edges for each thread:" << std::endl;
      for (auto &v : num_edges)
        std::cout << v << " ";
      std::cout << std::endl;
    }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Clean without counter, optional feature.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] features a fixed-size sequence of 8 bits to control which heuristics are enabled
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of each bin after dividing the workload
* @returns the edge list of the weighted s-line graph
*/
template <directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_blocked_optional_features_clean(
    const std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads,
    int bin_size = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            if (hyperE >= anotherhyperE) return;
            if (features[heuristics::SKIP_VISITED] && visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;

            //filter edges deg(e) < s
            if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[anotherhyperE] < s) return;

            //avoid duplicate intersections
            if (features[heuristics::SKIP_VISITED] && visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;
            
            //O(average degree of hyperedges)
            //short circuit
            auto r = features[heuristics::SHORT_CIRCUIT] ? (efficient::is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) :
            (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]));
            if (r) {
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Counter version. All features on. For benchmark purpose.
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
auto to_two_graph_efficient_blocked_with_counter (HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
    //avoid intersection when s=1
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            ++num_visits[worker_index];
            if (hyperE >= anotherhyperE) return;
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            ++num_visits[worker_index];
            
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;
            //filter edges deg(e) < s
            if (hyperedgedegrees[anotherhyperE] < s) return;
            //avoid duplicate intersections
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;
            
            //O(average degree of hyperedges)
            //short circuit
            if (efficient::is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) {
              ++num_edges[worker_index];
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

/**
* Efficient computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Counter version, optional feature. For benchmark only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] features a fixed-size sequence of 8 bits to control which heuristics are enabled
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of each bin after dividing the workload
* @returns the edge list of the weighted s-line graph
*/
template <directedness edge_directedness = nw::graph::directedness::undirected, class HyperEdge,
          class HyperNode, class vertex_id_t = vertex_id_t<HyperEdge>>
auto to_two_graph_efficient_blocked_optional_features_with_counter(
    const std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<vertex_id_t>& hyperedgedegrees, size_t s, int num_threads,
    int bin_size = 32) {
  size_t M = e_nbs.size();
  size_t N = n_nbs.size();
  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nedges = 0;
    //avoid intersection when s=1
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            if (hyperE >= anotherhyperE) return;
            if (visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;  
            ++nedges;
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          });
        });
      } //for
    }, tbb::auto_partitioner());
    std::cout << nedges << " edges added" << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    std::atomic<size_t> nvisited = 0, nbelowS = 0, nduplicates = 0, nintersections = 0, nedges = 0;
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      std::vector<bool> visitedE(M, false);
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            ++nvisited;
            
            //travese upper triangluar with lhs > rhs
            //avoid self edge with lhs == rhs
            if (hyperE >= anotherhyperE) return;
            ++nbelowS;

            //filter edges deg(e) < s
            if (features[heuristics::DEG_PRUNING] && hyperedgedegrees[anotherhyperE] < s) return;

            ++nduplicates;
            //avoid duplicate intersections
            if (features[heuristics::SKIP_VISITED] && visitedE[anotherhyperE]) return; else visitedE[anotherhyperE] = true;
            
            ++nintersections;
            //O(average degree of hyperedges)
            //short circuit
            auto r = features[heuristics::SHORT_CIRCUIT] ? (efficient::is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) :
            (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]));
            if (r) {
              ++nedges;
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
            }
          });//each neighbor of hyperN
        });//each neighbor of hyperE
      } //for each hyperE
     
    }, tbb::auto_partitioner());
    std::cout << nvisited << " visits, "
    << nbelowS << " hyperedges has less S nodes, "
    << nduplicates << " duplicate intersections avoid, "
    << nintersections << " intersections performed, " 
    << nedges << " edges added" << std::endl;
  }
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}

}//namespace hypergraph
}//namespace nw