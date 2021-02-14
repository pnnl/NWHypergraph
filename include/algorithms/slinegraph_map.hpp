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
#include <cyclic_range_adapter.hpp>

#include "util/slinegraph_helper.hpp"
#include "tbb/task_scheduler_init.h"

namespace nw {
namespace hypergraph {


template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_with_map_parallel(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();    
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
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();    
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : edges[hyperE]) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();
    return result;
  }
  return squeeze_edgelist(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_with_map_cyclic(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::task_arena::current_thread_index();    
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
    nw::util::life_timer _(__func__);
    tbb::parallel_for(nw::graph::cyclic(edges, num_bins), [&](auto& i) {
      int worker_index = tbb::task_arena::current_thread_index();    
      for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::map<size_t, size_t> K;
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE < anotherhyperE) ++K[anotherhyperE];
          }
        }
        for (auto &&[anotherhyperE, val] : K) {
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
        }
      }
    }, tbb::auto_partitioner());  
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();
    return result;
  }
  return squeeze_edgelist(two_graphs);
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_with_map_parallel2d(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
  size_t M = edges.size();
  size_t N = nodes.size();
  if (1 < s) {
    nw::util::life_timer _(__func__);
    tbb::parallel_for(tbb::blocked_range2d<vertex_id_t>(0, M / N, num_bins, 0, N, num_bins), [&](tbb::blocked_range2d<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
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
      int worker_index = tbb::task_arena::current_thread_index();
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
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();
    return result;
  }
  return squeeze_edgelist(two_graphs);
}

/*
* counts the neighbor hyperedges of the hyperedges
*/
template<class HyperEdge, class HyperNode>
auto to_two_graph_count_neighbors_cyclic(HyperEdge& edges, HyperNode& nodes, int num_bins = 32) {
  size_t M = edges.size();
  std::vector<std::map<size_t, size_t>> two_graphs(M, std::map<size_t, size_t>());
  nw::util::life_timer _(__func__);
  tbb::parallel_for(nw::graph::cyclic(edges, num_bins), [&](auto& i) {
    int worker_index = tbb::task_arena::current_thread_index();    
    for (auto&& j = i.begin(); j != i.end(); ++j) {
      auto&& [hyperE, hyperE_ngh] = *j;
      for (auto &&[hyperN, w] : edges[hyperE]) {
        for (auto &&[anotherhyperE, anotherw] : nodes[hyperN]) {
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
  size_t M = edges.size();
  std::vector<std::map<size_t, size_t>> two_graphs(M, std::map<size_t, size_t>());
  nw::util::life_timer _(__func__);
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::task_arena::current_thread_index();  
    for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
      for (auto &&[hyperN, w] : edges[hyperE]) {
        for (auto &&[anotherhyperE, anotherw] : nodes[hyperN]) {
          if (hyperE < anotherhyperE) ++two_graphs[hyperE][anotherhyperE];
        }
      }
    }
  }, tbb::auto_partitioner());
  return two_graphs;
}

template<directedness edge_directedness = undirected, class Hypergraph, typename... Attributes>
auto populate_linegraph_from_neighbor_map(Hypergraph& g, std::vector<std::map<size_t, size_t>>& neighbor_map, 
std::size_t s, bool weighted = false) {
  nw::graph::edge_list<edge_directedness, Attributes...> linegraph;
  size_t M = g.size();
  //n is the number of hyperedges, m is the number of hypernodes
  //time complexity of counting neighbors is same as the efficient: O(n*deg(edges)*deg(nodes)*deg(edges))
  //time complexity of extract slinegraph from the neighbor counts: O(n*deg(edges)) -> worst is O(n^2)
  //space complexity: O(n*total_deg(H)) -> worst is O(n^2)
  //total_deg(H)=sum of deg(each hyperedge)
  //total_deg(H) >> n?
  for (size_t hyperE = 0; hyperE < M; ++hyperE) {
    for (auto &&[anotherhyperE, val] : neighbor_map[hyperE]) {
      if (val >= s) {
        //std::cout << hyperE << "-" << anotherhyperE << std::endl;
        if (!weighted)
          val = 1;
        linegraph.push_back(std::make_tuple<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE), val));
      }
    }
  }
  linegraph.close_for_push_back(false);
  return linegraph;
}

}//namespace hypergraph
}//namespace nw