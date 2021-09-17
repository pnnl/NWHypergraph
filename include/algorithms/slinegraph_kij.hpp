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
#include <adaptors/vertex_range.hpp>

#include "util/slinegraph_helper.hpp"
#include <tbb/task_arena.h> //for tbb::this_task_arena::current_thread_index()
#include <util/atomic.hpp>

namespace nw {
namespace hypergraph {

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_spgemm_kij_cyclic(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();
  //std::vector<std::map<size_t, size_t>> soverlap(M, std::map<size_t, size_t>());
  //std::atomic<int>* soverlap = new std::atomic<int>[M * M];
  //std::atomic<int> * a1 = new std::atomic<int>[ n ];
  //std::atomic<size_t> soverlap[M * M];
  // Create
  //std::vector<int> soverlap(M * M, 0);
  //std::vector<std::vector<int>> soverlap(M, std::vector<int>(M, 0));
  size_t** soverlap = new size_t*[M];
  for(size_t i = 0; i < M; ++i) {
    soverlap[i] = new size_t[M];
    std::fill(soverlap[i], soverlap[i] + M, 0);
  }

  std::cout << "!!!!!" << std::endl;
    //std::vector<std::atomic<int>> soverlap(M * M);
//std::cout << "!!!" << soverlap.size() << std::endl;
    // Initialize.
 //for (size_t i = 0; i < M*M; ++i)
    //soverlap[i].store(0, std::memory_order_relaxed);
  //for(auto& e : soverlap)
   //     e.store(0, std::memory_order_relaxed);
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_bins);
  if (1 == s) {
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE >= anotherhyperE) continue;
            if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;  
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
        }
      }, tbb::auto_partitioner());
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
    //create an array of line graphs for each thread
    {
      nw::util::life_timer _(__func__);
      std::for_each(
          ep, nw::graph::counting_iterator<vertex_id_t>(0ul),
            nw::graph::counting_iterator<vertex_id_t>(N), [&](auto hyperN) {
            // tbb::parallel_for(nw::graph::cyclic_neighbor_range(nodes, num_bins), [&](auto&
            // i) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            // for (auto&& j = i.begin(); j != i.end(); ++j)
            {
              /// auto&& [hyperN, hyperN_ngh] = *j;
              // all neighbors of hyperN are hyperedges
              for (auto&& [hyperE] : nodes[hyperN]) {
                if (hyperedgedegrees[hyperE] < s) continue;
                for (auto&& [anotherhyperE] : nodes[hyperN]) {
                  // so we check compid of each hyperedge
                  // travese upper triangluar with lhs > rhs
                  // avoid self edge with lhs == rhs
                  if (hyperE >= anotherhyperE) continue;
                  // filter edges deg(e) < s
                  if (hyperedgedegrees[anotherhyperE] < s) continue; 
                  //size_t index = hyperE * M + anotherhyperE;
                  //if (0 == index % 9)
                  //std::cout << hyperE << "-" << anotherhyperE << "=" << index << std::endl;
                  // short circuit
                  //auto tmp = soverlap[index];
                  if (soverlap[hyperE][anotherhyperE] < s)
                    ++soverlap[hyperE][anotherhyperE];
                    //soverlap[index].fetch_add(1, std::memory_order_relaxed);
                    //std::atomic_fetch_add(&soverlap[index], 1);
                  //nw::graph::fetch_add(soverlap[hyperE][anotherhyperE], 1);
                    //nw::graph::fetch_add(soverlap[index], 1);
                    //soverlap[index].fetch_add(1, std::memory_order_relaxed); 
                }  // each neighbor of hyperN
              }    // each neighbor of hyperE
            }
          });
      //}, tbb::auto_partitioner());
      std::cout << "?????" << std::endl;
      for (vertex_id_t i = 0; i < M; ++i) {
        for (vertex_id_t j = 0; j < M; ++j) {
          if  (i == j) continue;
          //size_t index = i * M + j;
          
          //if (s <= soverlap[index]) 
          if (s <= soverlap[i][j]) 
          {
              //std::cout << soverlap[index] << " " << index << " !" << i <<"-" << j << std::endl;
            two_graphs[0].push_back(
                std::make_pair<vertex_id_t, vertex_id_t>(
                    std::forward<vertex_id_t>(i),
                    std::forward<vertex_id_t>(j)));
          }
        }
      }
    }
    for(size_t i = 0; i < M; ++i)
     delete soverlap[i];
    delete soverlap;
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}



template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph_spgemm_kij_cyclic_v2(ExecutionPolicy&& ep, HyperEdge& edges, HyperNode& nodes, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  size_t M = edges.size();
  size_t N = nodes.size();
  std::vector<std::map<size_t, size_t>> soverlap(M, std::map<size_t, size_t>());
  //std::atomic<int>* soverlap = new std::atomic<int>[M * M];
  //std::atomic<int> * a1 = new std::atomic<int>[ n ];
  //std::atomic<size_t> soverlap[M * M];
  // Create
  //std::vector<int> soverlap(M * M, 0);
  //std::vector<std::vector<int>> soverlap(M, std::vector<int>(M, 0));
  std::cout << "!!!!!" << std::endl;
    //std::vector<std::atomic<int>> soverlap(M * M);
//std::cout << "!!!" << soverlap.size() << std::endl;
    // Initialize.
 //for (size_t i = 0; i < M*M; ++i)
    //soverlap[i].store(0, std::memory_order_relaxed);
  //for(auto& e : soverlap)
   //     e.store(0, std::memory_order_relaxed);
  using linegraph_t = std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t>>>;
  linegraph_t two_graphs(num_bins);
  if (1 == s) {
    {
      nw::util::life_timer _(__func__);
      tbb::parallel_for(nw::graph::cyclic_neighbor_range(edges, num_bins), [&](auto& i) {
        int worker_index = tbb::this_task_arena::current_thread_index();
        std::vector<bool> visitedE(M, false);
        for (auto&& j = i.begin(); j != i.end(); ++j) {
        auto&& [hyperE, hyperE_ngh] = *j;
        std::fill(visitedE.begin(), visitedE.end(), false);
        //all neighbors of hyperedges are hypernode
        for (auto &&[hyperN] : hyperE_ngh) {
          for (auto &&[anotherhyperE] : nodes[hyperN]) {
            if (hyperE >= anotherhyperE) continue;
            if (visitedE[anotherhyperE]) continue; else visitedE[anotherhyperE] = true;  
            two_graphs[worker_index].push_back(std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE)));
          }
        }
        }
      }, tbb::auto_partitioner());
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
    //create an array of line graphs for each thread
    {
      nw::util::life_timer _(__func__);
      std::for_each(
          ep, nw::graph::counting_iterator<vertex_id_t>(0ul),
            nw::graph::counting_iterator<vertex_id_t>(N), [&](auto hyperN) {
            // tbb::parallel_for(nw::graph::cyclic_neighbor_range(nodes, num_bins), [&](auto&
            // i) {
            int worker_index = tbb::this_task_arena::current_thread_index();
            // for (auto&& j = i.begin(); j != i.end(); ++j)
            {
              /// auto&& [hyperN, hyperN_ngh] = *j;
              // all neighbors of hyperN are hyperedges
              std::vector<vertex_id_t> frontier;
              for (auto&& [hyperE] : nodes[hyperN]) {
                if (hyperedgedegrees[hyperE] >= s)
                  frontier.push_back(hyperE);
              }    // each neighbor of hyperE
              size_t n = frontier.size();
              if (1 >= n) return;
              for (size_t i = 0; i < n - 1; ++i) {
                auto hyperE = frontier[i];
                for (size_t j = i + 1; j < n; ++j) {
                  auto anotherhyperE = frontier[j];
                  if (soverlap[hyperE][anotherhyperE] < s)
                    ++soverlap[hyperE][anotherhyperE];
                    //soverlap[index].fetch_add(1, std::memory_order_relaxed);
                    //std::atomic_fetch_add(&soverlap[index], 1);
                  //nw::graph::fetch_add(soverlap[hyperE][anotherhyperE], 1);
                    //nw::graph::fetch_add(soverlap[index], 1);
                    //soverlap[index].fetch_add(1, std::memory_order_relaxed); 
                }  // each neighbor of hyperN
              }
            }
          });
      //}, tbb::auto_partitioner());
      std::cout << "?????" << std::endl;
      for (vertex_id_t i = 0; i < M; ++i) {
        for (vertex_id_t j = 0; j < M; ++j) {
          //if  (i == j) continue;
          //size_t index = i * M + j;
          
          //if (s <= soverlap[index]) 
          if (i != j && s <= soverlap[i][j]) 
          {
              //std::cout << soverlap[index] << " " << index << " !" << i <<"-" << j << std::endl;
            two_graphs[0].push_back(
                std::make_pair<vertex_id_t, vertex_id_t>(
                    std::forward<vertex_id_t>(i),
                    std::forward<vertex_id_t>(j)));
          }
        }
      }
    }
    //delete soverlap;
    return create_edgelist_with_squeeze(two_graphs);
  }//else
}




}//namespace hypergraph
}//namespace nw