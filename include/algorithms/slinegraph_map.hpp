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
#include <adaptors/cyclic_neighbor_range.hpp>
#include <adaptors/vertex_range.hpp>

#include "util/slinegraph_helper.hpp"
#include <tbb/task_arena.h>
#include <tbb/blocked_range2d.h>

namespace nw {
namespace hypergraph {


template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_map_blocked(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
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
      }
    }, tbb::auto_partitioner());
  }
  else {
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_map_cyclic(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
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
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) 
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());
  }
  else {
    {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        //std::map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_map_blocked_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_map_cyclic_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_blocked(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        std::unordered_map<size_t, size_t> K;
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
      }
    }, tbb::auto_partitioner());
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_cyclic(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
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
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
          if (val >= s) 
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_blocked_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_cyclic_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_vector_blocked(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    //create a thread-local vector to store the soverlap information
    std::vector<std::vector<size_t>> soverlap(num_threads, std::vector<size_t>(M, 0));
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      auto& K = soverlap[worker_index];
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        //std::fill(K.begin(), K.end(), 0);
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (vertex_id_t anotherhyperE = 0; anotherhyperE < M; ++anotherhyperE) {
          //if this slot is unused, no need to reset its value to 0
          if (0 == K[anotherhyperE]) continue;
          if (K[anotherhyperE] >= s) {
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
          K[anotherhyperE] = 0; //set to 0 for next use
        }
      }
    }, tbb::auto_partitioner());
  }
  else {
    {
    nw::util::life_timer _(__func__);
    //create a thread-local vector to store the soverlap information
    std::vector<std::vector<bool>> soverlap(num_threads, std::vector<bool>(M, false));
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, bin_size), [&](tbb::blocked_range<vertex_id_t>& r) {
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_vector_cyclic(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    //create a thread-local vector to store the soverlap information
    std::vector<std::vector<size_t>> soverlap(num_threads, std::vector<size_t>(M, 0));
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();   
      auto& K = soverlap[worker_index]; 
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        if (hyperedgedegrees[hyperE] < s) continue;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperedgedegrees[anotherhyperE] < s) continue;
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (vertex_id_t anotherhyperE = 0; anotherhyperE < M; ++anotherhyperE) {
          //if this slot is unused, no need to reset its value to 0
          if (0 == K[anotherhyperE]) continue;
          if (K[anotherhyperE] >= s) {
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
          K[anotherhyperE] = 0; //set to 0 for next use
        }
      }
    }, tbb::auto_partitioner());
  }
  else {
    {
    nw::util::life_timer _(__func__);
    //create a thread-local vector to store the soverlap information
    std::vector<std::vector<bool>> soverlap(num_threads, std::vector<bool>(M, false));
    tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::this_task_arena::current_thread_index();
      auto& visitedE = soverlap[worker_index];
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_vector_blocked_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_vector_cyclic_with_counter(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_with_map_parallel2d(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_static_hashmap_blocked(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int bin_size = 32) {
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
  return create_edgelist_with_squeeze(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_static_hashmap_cyclic(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
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
  return create_edgelist_with_squeeze(two_graphs);
}

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
  return create_edgelist_with_squeeze(two_graphs);
}

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
  return create_edgelist_with_squeeze(two_graphs);
}

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

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_map_blocked_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_map_blocked(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_map_blocked_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_map_cyclic_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_map_cyclic(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_map_cyclic_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_blocked_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_hashmap_blocked(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_hashmap_blocked_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_hashmap_cyclic_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_hashmap_cyclic(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_hashmap_cyclic_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_vector_blocked_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_vector_blocked(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_vector_blocked_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_vector_cyclic_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_vector_cyclic(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_vector_cyclic_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_static_hashmap_blocked_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_static_hashmap_blocked(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_static_hashmap_blocked_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_static_hashmap_cyclic_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s, int num_threads, int num_bins = 32) {
  if (!verbose) 
    return to_two_graph_static_hashmap_cyclic(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
  else
    return to_two_graph_static_hashmap_cyclic_with_counter(ep, edges, nodes, hyperedgedegrees, s, num_threads, num_bins);
}

/*
* counts the neighbor hyperedges of the hyperedges
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

/*
* counts the neighbor hyperedges of the hyperedges
*/
template<class HyperEdge, class HyperNode>
auto to_two_graph_count_neighbors_cyclic(HyperEdge& edges, HyperNode& nodes, size_t min_s, int num_bins = 32) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  std::vector<std::map<size_t, size_t>> two_graphs(M, std::map<size_t, size_t>());
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

/*
* counts the neighbor hyperedges of the hyperedges
*/
template<class HyperEdge, class HyperNode>
auto to_two_graph_count_neighbors_blocked(HyperEdge& edges, HyperNode& nodes) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  std::vector<std::unordered_map<size_t, size_t>> two_graphs(M, std::unordered_map<size_t, size_t>());
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
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

/*
* counts the neighbor hyperedges of the hyperedges
*/
template<class HyperEdge, class HyperNode>
auto to_two_graph_count_neighbors_blocked(HyperEdge& edges, HyperNode& nodes, size_t min_s) {
  nw::util::life_timer _(__func__);
  size_t M = edges.size();
  std::vector<std::unordered_map<size_t, size_t>> two_graphs(M, std::unordered_map<size_t, size_t>());
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
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

template<directedness edge_directedness = undirected>
auto populate_linegraph_from_neighbor_map(std::vector<std::unordered_map<size_t, size_t>>& neighbor_map, std::size_t s) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> linegraph(0);
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
  linegraph.close_for_push_back(false);
  return linegraph;
}

template<directedness edge_directedness = undirected, class Hypergraph, typename... Attributes>
auto populate_linegraph_from_neighbor_map(Hypergraph& g, std::vector<std::unordered_map<size_t, size_t>>& neighbor_map, 
std::size_t s, bool weighted = false) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness, Attributes...> linegraph;
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