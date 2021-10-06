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
#include <util/intersection_size.hpp>
#include "util/slinegraph_helper.hpp"
#include "tbb/task_arena.h"

namespace nw {
namespace hypergraph {


/*
* Parallel naive version, clean version, for perf testing
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_naive_parallel(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
size_t s, int num_threads, int num_bins = 32) {
  size_t M = e_nbs.size();
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  {
  nw::util::life_timer _(__func__);
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, M / num_bins), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::this_task_arena::current_thread_index();
    for (auto i = r.begin(), e = r.end(); i != e; ++i) {
      for (size_t j = i + 1; j < M; ++j) {
        if (nw::graph::intersection_size(e_nbs[i], e_nbs[j]) >= s) {
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(i), std::forward<vertex_id_t>(j)));
        }
      }
    }
  }, tbb::auto_partitioner());
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/*
* Parallel naive version with statistic counters, for benchmarking
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_naive_parallel_with_counter(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
size_t s, int num_threads, int num_bins = 32) {

  size_t M = e_nbs.size();
  std::vector<size_t> num_visits(num_threads, 0), num_edges(num_threads, 0);
  std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>> two_graphs(num_threads);
  {
  nw::util::life_timer _(__func__);
  tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M, M / num_bins), [&](tbb::blocked_range<vertex_id_t>& r) {
    int worker_index = tbb::this_task_arena::current_thread_index();
    for (auto i = r.begin(), e = r.end(); i != e; ++i) {
      for (size_t j = i + 1; j < M; ++j) {
        ++num_visits[worker_index];
        if (nw::graph::intersection_size(e_nbs[i], e_nbs[j]) >= s) {
          two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(i), std::forward<vertex_id_t>(j)));
          ++num_edges[worker_index];
        }
      }
    }
  }, tbb::auto_partitioner());
  //std::cout << "#visits for each thread:" << std::endl;
  std::cout << "Total #visits: ";
  size_t total_visits = 0;
  for (auto &v : num_visits) {
    total_visits += v;
    //std::cout << v << " ";
  }
  std::cout << total_visits << std::endl;
  size_t total_edges = 0;
  std::cout << "Total #edges: ";
  //std::cout << "#edges for each thread:" << std::endl;
  for (auto &v : num_edges) {
    total_edges += v;
    //std::cout << v << " ";
  }
  std::cout << total_edges << std::endl;
  }
  return create_edgelist_with_squeeze(two_graphs);
}

/*
* TODO Have not squeeze.
*/
template<directedness edge_directedness = undirected, class HyperEdge, class HyperNode>
auto to_two_graph_naive_serial(HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();
  size_t counter = 0, nedges = 0;
  for (size_t i = 0; i < e_nbs.size(); ++i) {
    for (size_t j = i + 1; j < e_nbs.size(); ++j) {
      ++counter;
      if (nw::graph::intersection_size(e_nbs[i], e_nbs[j]) >= s) {
        two_graph.push_back(i, j);
        ++nedges;
      }
    }
  }
  std::cout << counter << " intersections performed, " 
  << nedges << " edges added" << std::endl;
  two_graph.close_for_push_back();
  return two_graph;
}


template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_naive_parallel_portal(bool verbose, ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
size_t s, int num_threads, int num_bins = 32) {
  if(!verbose)
    return to_two_graph_naive_parallel(ep, e_nbs, n_nbs, s, num_threads, num_bins);
  else
    return to_two_graph_naive_parallel_with_counter(ep, e_nbs, n_nbs, s, num_threads, num_bins);
}

}//namespace hypergraph
}//namespace nw