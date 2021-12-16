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

#include "util/slinegraph_helper.hpp"
#include "../to_two_graph_efficient.hpp"
#include "../to_two_graph_map.hpp"
#include <tbb/task_arena.h>

using namespace nw::graph;

namespace nw {
namespace hypergraph {

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Counter version. For benchmark purpose.
* Use a map to store overlapping counts.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] bin_size the size of bins after dividing the workload
* @returns the edge list of the s-line graph
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_map_blocked_with_counter(HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) {
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
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          ++num_edges[worker_index];
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy. 
* Counter version. For benchmark purpose.
* Use a map to store overlapping counts.
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
auto to_two_graph_map_cyclic_with_counter(HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        if (hyperedgedegrees[hyperE] < s) continue;
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) {
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
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          ++num_edges[worker_index];
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
}


/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a hashmap to store overlapping counts. Benchmark only.
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_blocked_with_counter(HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        std::unordered_map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) {
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
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::unordered_map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          ++num_edges[worker_index];
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy.
* Using a hashmap to store overlapping counts. Benchmark only.
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_cyclic_with_counter(HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        if (hyperedgedegrees[hyperE] < s) continue;
        std::unordered_map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) {
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
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::unordered_map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          ++num_edges[worker_index];
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  } 
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a pre-allocated hashmap to store overlapping counts. Benchmark only.
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
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_static_hashmap_blocked_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
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
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (0 == val) continue;
          if (val >= s) {
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
          val = 0;
        }
      }
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::unordered_map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          ++num_edges[worker_index];
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
}
/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy.
* Using a pre-allocated hashmap to store overlapping counts. Benchmark only.
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
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_static_hashmap_cyclic_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
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
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (0 ==  val) continue;
          if (val >= s) {
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
          val = 0;
        }
      }
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::unordered_map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          ++num_edges[worker_index];
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy.
* Using a vector to store overlapping counts. Benchmark only.
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_vector_blocked_with_counter(HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
  if (1 < s) {
    nw::util::life_timer _(__func__);
    //create a thread-local vector to store the soverlap information
    std::vector<std::vector<size_t>> soverlap(num_threads, std::vector<size_t>(M, 0));
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      auto& K = soverlap[worker_index];
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        //extract edges out of soverlap information
        for (vertex_id_t anotherhyperE = 0; anotherhyperE < M; ++anotherhyperE) {
          //if this slot is unused, no need to reset its value to 0
          if (0 == K[anotherhyperE]) continue;
          if (K[anotherhyperE] >= s) {
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
          //after query of soverlap, reset K[anotherhyperE] to 0 for next use
          K[anotherhyperE] = 0;
        }
      } //for
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    //create a thread-local vector to store the soverlap information
    std::vector<std::vector<bool>> soverlap(num_threads, std::vector<bool>(M, false));    
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      auto& K = soverlap[worker_index];
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::fill(K.begin(), K.end(), false);
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE >= anotherhyperE) continue;
            //check hyperE-anotherhyperE has been inserted or not
            //to avoid duplicate edges in the slinegraph
            if (K[anotherhyperE]) {
              ++num_counts[worker_index];
              continue; 
            }
            else 
              K[anotherhyperE] = true;
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
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses cyclic range as workload distribution strategy.
* Using a fix-size vector to store overlapping counts. Benchmark only.
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
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_vector_cyclic_with_counter(HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  std::vector<size_t> num_visits(num_threads, 0), num_counts(num_threads, 0), num_edges(num_threads, 0);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      std::vector<size_t> K(M, 0);
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        if (hyperedgedegrees[hyperE] < s) continue;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) {
              ++num_counts[worker_index];
              ++K[anotherhyperE];
            }
          }
        }
        //extract edges out of soverlap information
        for (vertex_id_t anotherhyperE = 0; anotherhyperE < M; ++anotherhyperE) {
          //if this slot is unused, no need to reset its value to 0
          if (0 == K[anotherhyperE]) continue;
          if (K[anotherhyperE] >= s) {
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
          //after query of soverlap, reset K[anotherhyperE] to 0 for next use
          K[anotherhyperE] = 0;
        }
      }
    }, tbb::auto_partitioner());
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_with_squeeze(two_graphs);
  }
  else {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      std::vector<bool> K(M, false);
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::fill(K.begin(), K.end(), false);
        std::unordered_map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            ++num_visits[worker_index];
            if (hyperE >= anotherhyperE) continue;
            if (K[anotherhyperE]) {
              ++num_counts[worker_index];
              continue;
            }
            else
              K[anotherhyperE] = true;
            ++num_edges[worker_index];
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      } //for
    }, tbb::auto_partitioner());  
    std::cout << "#visits for each thread:" << std::endl;
    for (auto &v : num_visits)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#counts for each thread:" << std::endl;
    for (auto &v : num_counts)
      std::cout << v << " ";
    std::cout << std::endl;
    std::cout << "#edges for each thread:" << std::endl;
    for (auto &v : num_edges)
      std::cout << v << " ";
    std::cout << std::endl;
    return create_edgelist_without_squeeze(two_graphs);
  }
}

/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Use a map for overlapping information. For various bin sizes experiments only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution policy for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] ep the execution policy
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] max_bins the maximum number of bins to experiment with
* @returns edge list of 1-line graph or an empty edge list if s > 1
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_map_blocked_vary_size(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int max_bins = 32) {

  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    //block size: 1, 4, 16, 64, 256, 1024, ... until reach 4^num_bins
    for (size_t i = 1; max_bins > 0; i *= 4, --max_bins) {
      size_t num_bins = M / i; // num_bins: M, M/4, M/16, M/256, ...
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          tbb::blocked_range<vertex_id_t>(0, M, M / num_bins),
          [&](tbb::blocked_range<vertex_id_t>& r) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
              if (hyperedgedegrees[hyperE] < s) continue;
              std::map<size_t, size_t> K;
              for (auto&& [hyperN] : edges[hyperE]) {
                for (auto&& [anotherhyperE] : nodes[hyperN]) {
                  if (hyperedgedegrees[anotherhyperE] < s) continue;
                  if (hyperE < anotherhyperE) ++K[anotherhyperE];
                }
              }
              for (auto&& [anotherhyperE, val] : K) {
                if (val >= s)
                  two_graphs[worker_index].push_back(
                      std::make_pair<vertex_id_t, vertex_id_t>(
                          std::forward<vertex_id_t>(hyperE),
                          std::forward<vertex_id_t>(anotherhyperE)));
              }
            }
          },
          tbb::auto_partitioner());
    }
    return nw::graph::edge_list<edge_directedness>(0);
  }
  else {
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, max_bins), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        //std::map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE < anotherhyperE) 
              two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
      }
    }, tbb::auto_partitioner());
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
}
/**
* Computation of a s-line graph of a hypergraph. 
* It uses blocked range as workload distribution strategy. 
* Use a hashmap for overlapping information. For various bin sizes experiments only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution policy for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] ep the execution policy
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] max_bins the maximum number of bins to experiment with
* @returns edge list of 1-line graph or an empty edge list if s > 1
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_blocked_vary_size(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int max_bins = 32) {

  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    // block size: 1, 4, 16, 64, 256, 1024, ... until reach 4^num_bins
    for (size_t i = 1; max_bins > 0; i *= 4, --max_bins) {
      size_t num_bins = M / i;  // num_bins: M, M/4, M/16, M/256, ...
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(
          num_threads);
      nw::util::life_timer _(__func__);
      tbb::parallel_for(
          tbb::blocked_range<vertex_id_t>(0, M, M / num_bins),
          [&](tbb::blocked_range<vertex_id_t>& r) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
              if (hyperedgedegrees[hyperE] < s) continue;
              std::unordered_map<size_t, size_t> K;
              for (auto&& [hyperN] : edges[hyperE]) {
                for (auto&& [anotherhyperE] : nodes[hyperN]) {
                  if (hyperedgedegrees[anotherhyperE] < s) continue;
                  if (hyperE < anotherhyperE) ++K[anotherhyperE];
                }
              }
              for (auto&& [anotherhyperE, val] : K) {
                if (val >= s)
                  two_graphs[worker_index].push_back(
                      std::make_pair<vertex_id_t, vertex_id_t>(
                          std::forward<vertex_id_t>(hyperE),
                          std::forward<vertex_id_t>(anotherhyperE)));
              }
            }
          },
          tbb::auto_partitioner());
    }
    return nw::graph::edge_list<edge_directedness>(0);
  }
  else {
    std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, M / max_bins), [&](tbb::blocked_range<vertex_id_t>& r) {
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
* It uses blocked range as workload distribution strategy. 
* Use a vector for overlapping information. For various bin sizes experiments only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution policy for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] ep the execution policy
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] max_bins the maximum number of bins to experiment with
* @returns edge list of 1-line graph or an empty edge list if s > 1
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_vector_blocked_vary_size(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int max_bins = 32) {

  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    // block size: 1, 4, 16, 64, 256, 1024, ... until reach 4^num_bins
    for (size_t i = 1; max_bins > 0; i *= 4, --max_bins) {
      size_t num_bins = M / i;  // num_bins: M, M/4, M/16, M/256, ...
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(
          num_threads);
      nw::util::life_timer _(__func__);
      // create a thread-local vector to store the soverlap information
      std::vector<std::vector<size_t>> soverlap(num_threads,
                                                std::vector<size_t>(M, 0));
      tbb::parallel_for(
          tbb::blocked_range<vertex_id_t>(0, M, M / num_bins),
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
                      std::make_pair<vertex_id_t, vertex_id_t>(
                          std::forward<vertex_id_t>(hyperE),
                          std::forward<vertex_id_t>(anotherhyperE)));
                }
                K[anotherhyperE] = 0;  // set to 0 for next use
              }
            }
          },
          tbb::auto_partitioner());
    }
    return nw::graph::edge_list<edge_directedness>(0);
  }
  else {
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
    {
    nw::util::life_timer _(__func__);
    //create a thread-local vector to store the soverlap information
    std::vector<std::vector<bool>> soverlap(num_threads, std::vector<bool>(M, false));
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, M / max_bins), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      auto& visitedE = soverlap[worker_index];
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
* It uses blocked range as workload distribution strategy. 
* Use a pre-allocated hashmap for overlapping information. 
* For various bin sizes experiments only.
*
* @tparam edge_directedness the type of the edge directedness
* @tparam ExecutionPolicy the type of the execution policy for std::for_each
* @tparam HyperEdge the type of hyperedge incidence
* @tparam HyperNode the type of hypernode incidence
* @param[in] ep the execution policy
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] hyperedgedegrees the degrees of hyperedges
* @param[in] s the number of overlapping vertices between each hyperedge pair
* @param[in] num_threads the number of threads
* @param[in] max_bins the maximum number of bins to experiment with
* @returns edge list of 1-line graph or an empty edge list if s > 1
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_static_hashmap_blocked_vary_size(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int max_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    // block size: 1, 4, 16, 64, 256, 1024, ... until reach 4^num_bins
    for (size_t i = 1; max_bins > 0; i *= 4, --max_bins) {
      size_t num_bins = M / i;  // num_bins: M, M/4, M/16, M/256, ...
      std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(
          num_threads);
      nw::util::life_timer _(__func__);
      std::vector<std::unordered_map<size_t, size_t>> soverlap(num_threads);
      std::for_each(ep, soverlap.begin(), soverlap.end(),
                    [&](auto& K) { K.reserve(M); });
      tbb::parallel_for(
          tbb::blocked_range<vertex_id_t>(0, M, M / num_bins),
          [&](tbb::blocked_range<vertex_id_t>& r) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
              if (hyperedgedegrees[hyperE] < s) continue;
              auto& K = soverlap[worker_index];
              for (auto&& [hyperN] : edges[hyperE]) {
                for (auto&& [anotherhyperE] : nodes[hyperN]) {
                  if (hyperedgedegrees[anotherhyperE] < s) continue;
                  if (hyperE < anotherhyperE) ++K[anotherhyperE];
                }
              }
              for (auto&& [anotherhyperE, val] : K) {
                // if soverlap is not used, continue;
                if (0 == val) continue;
                if (val >= s)
                  two_graphs[worker_index].push_back(
                      std::make_pair<vertex_id_t, vertex_id_t>(
                          std::forward<vertex_id_t>(hyperE),
                          std::forward<vertex_id_t>(anotherhyperE)));
                // reset soverlap information for next iteration use
                K[anotherhyperE] = 0;  // equal to val = 0;
              }
            }
          },
          tbb::auto_partitioner());
    }
    return nw::graph::edge_list<edge_directedness>(0);
  }
  else {
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, max_bins), [&](tbb::blocked_range<vertex_id_t>& r) {
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


}//namespace hypergraph
}//namespace nw