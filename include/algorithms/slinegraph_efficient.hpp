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

using namespace nw::hypergraph::efficient;

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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_vary_size (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int max_bins = 1) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  if (1 == s) {
    linegraph_t two_graphs(num_threads);
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, max_bins, 0, M);
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    auto M = edges.size();
    //block size: 1, 4, 16, 64, 256, 1024, ... until reach 4^num_bins
    for (size_t i = 1; max_bins > 0; i *= 4, --max_bins) {
      linegraph_t two_graphs(num_threads);
      size_t num_bins = M / i; // num_bins: M, M/4, M/16, M/256, ...
      to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, 0, M, hyperedgedegrees, s); 
    }
    //abandon results
    return nw::graph::edge_list<edge_directedness>(0);
  }//else
}

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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze(two_graphs);
  } else {
    // when s > 1
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
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
* @param[in] bin_size the number of bins to divide the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_parallel_2d (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();

  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    {
    nw::util::life_timer _(__func__);
    //avoid intersection when s=1
        to_two_graph_block_range2d(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins,
        0, M / num_bins, num_bins, 0, num_bins, num_bins);
        to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, M / num_bins * num_bins, M);
    }
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_threads), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    //create an array of line graphs for each thread
    {
        nw::util::life_timer _(__func__);
        to_two_graph_block_range2d(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins,
        0, M / num_bins, num_bins, 0, num_bins, num_bins,
        hyperedgedegrees, s);
        to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, 
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_cyclic (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins);
    return create_edgelist_without_squeeze(two_graphs);
  }
  else {
    //when s > 1
    to_two_graph_cyclic(std::forward<linegraph_t>(two_graphs), edges, nodes, num_bins, hyperedgedegrees, s);
    return create_edgelist_with_squeeze(two_graphs);
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_cyclic_with_counter (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
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
              if (is_intersection_size_s(edges[hyperE].begin(), edges[hyperE].end(),
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
template<directedness edge_directedness = undirected, class T, class HyperEdge, class HyperNode>
auto to_weighted_two_graph_efficient_blocked (HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();
  //create an array of line graphs for each thread
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_weighted_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes, M / bin_size, 0, M);
    return create_edgelist_without_squeeze<undirected, T>(two_graphs);
  }
  else {
    //when s > 1
    to_weighted_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges,
                               nodes, M / bin_size, 0, M, hyperedgedegrees, s);
    return create_edgelist_with_squeeze<undirected, T>(two_graphs);
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
template <directedness edge_directedness = undirected, class HyperEdge,
          class HyperNode>
auto to_two_graph_efficient_blocked_optional_features_clean(
    std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads,
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
            auto r = features[heuristics::SHORT_CIRCUIT] ? (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_with_counter (HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_threads), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
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
            if (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_optional_features_with_counter(std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
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
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(nw::graph::counting_iterator<int>(0), nw::graph::counting_iterator<int>(num_threads), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
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
            auto r = features[heuristics::SHORT_CIRCUIT] ? (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
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
template<directedness edge_directedness = nw::graph::undirected, class HyperEdge, class HyperNode>
auto to_two_graph_efficient_blocked_without_sequeeze(
  HyperEdge& edges, 
  HyperNode& nodes, 
  std::vector<index_t>& hyperedgedegrees, 
  size_t s,
  int num_threads, 
  int bin_size = 32) {
  size_t M = edges.size();
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_threads);
  if (1 == s) {
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes,
                         M / bin_size, 0, M);
    return two_graphs;
  } else {
    // when s > 1
    to_two_graph_blocked(std::forward<linegraph_t>(two_graphs), edges, nodes,
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
template<directedness edge_directedness = nw::graph::undirected, class T, class HyperEdge, class HyperNode>
auto to_weighted_two_graph_efficient_blocked_without_squeeze(
  HyperEdge& edges, 
  HyperNode& nodes, 
  std::vector<index_t>& hyperedgedegrees, 
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
template <directedness edge_directedness = undirected, class HyperEdge,
          class HyperNode>
auto to_two_graph_efficient_cyclic_portal(
    bool verbose, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<index_t>& hyperedgedegrees, size_t s, int numb_threads,
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
template <directedness edge_directedness = undirected, class HyperEdge,
          class HyperNode>
auto to_two_graph_efficient_blocked_portal(
    bool verbose, std::bitset<8>& features, HyperEdge& e_nbs, HyperNode& n_nbs,
    std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads,
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