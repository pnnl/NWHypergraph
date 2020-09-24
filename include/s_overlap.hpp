//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2020
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#pragma once
#include <util/timer.hpp>
#include <util/intersection_size.hpp>
#include <util/parallel_for.hpp>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_hash_map.h>
#include <util/AtomicBitVector.hpp>
#include <mutex>
#include "tbb/task_scheduler_init.h"

namespace nw {
namespace hypergraph {
/// Basic helper used for all of the inner set intersections.
///
/// This wraps `std::set_intersection` to produce the size of the set rather
/// than the set itself, and also handles the fact that our iterator value types
/// are tuples where we only care about the first element for ordering.
///
/// @tparam           A The type of the first iterator.
/// @tparam           B The type of the second iterator.
/// @tparam           C The type of the third iterator.
/// @tparam           D The type of the fourth iterator.
/// @tparam ExecutionPolicy The type of the parallel execution policy.
///
/// @param            i The beginning of the first range.
/// @param           ie The end of the first range.
/// @param            j The beginning of the second range.
/// @param           je The end of the second range.
/// @param           ep The parallel execution policy.
///
/// @returns            The size of the intersected set.

/*
*
we can pass in s  and check s against n  there when we increment n and as soon we see n==s we can immediately return without doing the full intersection.
and return true/false .
We don't have to use std::set_intersection, instead write our own as this one above and make the return type to be bool
*/
template <class A, class B, class C, class D>
bool is_intersection_size_s(A i, B&& ie, C j, D&& je, size_t s = 1) {
  // Custom comparator because we know our iterator operator* produces tuples
  // and we only care about the first value.
  static constexpr auto lt = [](auto&& x, auto&& y) { return std::get<0>(x) < std::get<0>(y); };

  // Use our own trivial loop for the intersection size when the execution
  // policy is sequential, otherwise rely on std::set_intersection.
  //
  // @todo We really don't need set intersection. You'd hope that it would be
  //       efficient with the output counter, but it just isn't. Parallelizing
  //       the intersection size seems non-trivial though.
    std::size_t n = 0;
    while (i != ie && j != je) {
      if (lt(*i, *j)) {
        ++i;
      } else if (lt(*j, *i)) {
        ++j;
      } else {
        ++n;
        if (n >= s) return true;
        ++i;
        ++j;
      }
    }
    return false;
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();
  size_t counter = 0;
  for (size_t i = 0; i < e_nbs.size(); ++i) {
    for (size_t j = i + 1; j < e_nbs.size(); ++j) {
      ++counter;
      size_t count = nw::graph::intersection_size(e_nbs[i], e_nbs[j]);    //       if intersection_size(n_nbs(i), n_nbs(j)) >= s
      if (count >= s) {
        two_graph.push_back(i, j);    //         add (i,j) to new edge_list
        std::cout << i << " " << j << std::endl;
        std::cout << i << " :";
        std::for_each(e_nbs.begin()[i].begin(), e_nbs.begin()[i].end(), [&](auto &&x) {
          auto hyperN = std::get<0>(x);
          std::cout << hyperN << " ";
        });
        std::cout << std::endl;
        std::cout << j << " :";
        std::for_each(e_nbs.begin()[j].begin(), e_nbs.begin()[j].end(), [&](auto &&x) {
          auto hyperN = std::get<0>(x);
          std::cout << hyperN << " ";
        });
        std::cout << std::endl;
      }
    }
  }
  std::cout << counter << " intersections performed" << std::endl;
  two_graph.close_for_push_back();
  return two_graph;
}

/*
* Serial efficient version
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv2(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, std::vector<index_t>& hyperedgedegrees, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();
  size_t counter = 0;
  if (1 == s) {
    for (size_t hyperE = 0; hyperE < e_nbs.size(); ++hyperE) { //O(n)
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
        auto hyperN = std::get<0>(x);
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
          //so we check compid of each hyperedge
          two_graph.push_back(hyperE, std::get<0>(y)); //         add (i,j) to new edge_list
        });
      });
    }
  }
  else {
    //O(n*average degree of hyperedges*average degree of hypernodes*average degree of hyperedges)
    for (size_t hyperE = 0; hyperE < e_nbs.size(); ++hyperE) { //O(n)
      if (hyperedgedegrees[hyperE] < s) continue;
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
        auto hyperN = std::get<0>(x);
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
          //so we check compid of each hyperedge
          auto anotherhyperE = std::get<0>(y);
          if (hyperE >= anotherhyperE) return;
          if (hyperedgedegrees[anotherhyperE] < s) return;
          ++counter;
          //O(average degree of hyperedges)
          if (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]))
            two_graph.push_back(hyperE, anotherhyperE); //         add (i,j) to new edge_list
        });
      });
    }
  }//else
  std::cout << counter << " intersections performed" << std::endl;

  two_graph.close_for_push_back();
  return two_graph;
}

/*
* Parallel efficient version
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv3(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, std::vector<index_t>& hyperedgedegrees, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();
  size_t counter = 0;
  if (1 == s) {
    //or use a parallel for loop
    std::for_each(tbb::counting_iterator<vertex_id_t>(0), tbb::counting_iterator<vertex_id_t>(e_nbs.size()), [&](auto hyperE) {
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
        auto hyperN = std::get<0>(x);
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
          //so we check compid of each hyperedge
          two_graph.push_back(hyperE, std::get<0>(y)); //         add (i,j) to new edge_list
        });
      });
    });
  }
  else {
    std::vector<vertex_id_t> frontier;
    //O(n*average degree of hyperedges*average degree of hypernodes*average degree of hyperedges)
    std::for_each(tbb::counting_iterator<vertex_id_t>(0), tbb::counting_iterator<vertex_id_t>(e_nbs.size()), [&](auto hyperE) {
      if (hyperedgedegrees[hyperE] >= s) frontier.push_back(hyperE);
    });
    std::for_each(frontier.begin(), frontier.end(), [&](auto&& hyperE) {
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
        auto hyperN = std::get<0>(x);
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
          //so we check compid of each hyperedge
          auto anotherhyperE = std::get<0>(y);
          if (hyperE >= anotherhyperE) return;
          if (hyperedgedegrees[anotherhyperE] < s) return;
          ++counter;
          //O(average degree of hyperedges)
          if (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]))
            two_graph.push_back(hyperE, anotherhyperE); 
        });
      });
    });
  }//else
  std::cout << counter << " intersections performed" << std::endl;

  two_graph.close_for_push_back();
  return two_graph;
}


template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv4(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, std::vector<index_t>& hyperedgedegrees, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();
  size_t counter = 0;
  if (1 == s) {
    //or use a parallel for loop
    std::for_each(tbb::counting_iterator<vertex_id_t>(0), tbb::counting_iterator<vertex_id_t>(e_nbs.size()), [&](auto hyperE) {
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
        auto hyperN = std::get<0>(x);
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
          //so we check compid of each hyperedge
          two_graph.push_back(hyperE, std::get<0>(y)); //         add (i,j) to new edge_list
        });
      });
    });
  }
  else {
    std::vector<vertex_id_t> frontier;
    //O(n*average degree of hyperedges*average degree of hypernodes*average degree of hyperedges)
    std::for_each(tbb::counting_iterator<vertex_id_t>(0), tbb::counting_iterator<vertex_id_t>(e_nbs.size()), [&](auto hyperE) {
      if (hyperedgedegrees[hyperE] >= s) frontier.push_back(hyperE);
    });
    std::vector<vertex_id_t>::iterator it = frontier.begin();
    std::for_each(it, frontier.end(), [&](auto&& hyperE) {
      auto begin = it + 1;
      std::for_each(begin, frontier.end(), [&](auto &&anotherhyperE) { //O(average degree of hyperedges)
          if (hyperE >= anotherhyperE) return;
          if (hyperedgedegrees[anotherhyperE] < s) return;
          ++counter;
          //O(average degree of hyperedges)
          if (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]))
            two_graph.push_back(hyperE, anotherhyperE); 
      });
    });
  }//else
  std::cout << counter << " intersections performed" << std::endl;

  two_graph.close_for_push_back();
  return two_graph;
}

/*
* Parallel efficient version
*/
template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv5(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  nw::util::life_timer _(__func__);

  size_t M = e_nbs.size();
  size_t N = n_nbs.size();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    size_t nintersections = 0, nedges = 0;
    std::vector<nw::graph::edge_list<edge_directedness>> two_graphs(num_bins, nw::graph::edge_list<edge_directedness>(0));
    std::for_each(ep, tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].open_for_push_back();
    });
      
    nw::graph::AtomicBitVector   visitedE(M);
    //or use a parallel for loop
    std::for_each(tbb::counting_iterator<vertex_id_t>(0), tbb::counting_iterator<vertex_id_t>(M), [&](auto hyperE) {
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
        auto hyperN = std::get<0>(x);
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
          //so we check compid of each hyperedge
          auto anotherhyperE = std::get<0>(y);
          if (hyperE >= anotherhyperE) return;
          if (1 == visitedE.atomic_get(anotherhyperE)) return; else visitedE.atomic_set(anotherhyperE);
          nedges++;
          two_graphs[hyperE % num_bins].push_back(hyperE, anotherhyperE);
        });
      });
      visitedE.clear();
    });
    std::cout << nintersections << " intersections performed, " 
    << nedges << " edges added" << std::endl;
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].close_for_push_back();
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //find the line graph
    /*
    std::vector<std::vector<vertex_id_t>> frontier(num_bins);

    nw::graph::parallel_for(tbb::blocked_range(0ul, M), [&](auto&& hyperE) {
      if (s <= hyperedgedegrees[hyperE]) frontier[hyperE % num_bins].push_back(hyperE);
    });
    */
    std::vector<vertex_id_t> frontier;
    //O(n*average degree of hyperedges*average degree of hypernodes*average degree of hyperedges)
    std::for_each(tbb::counting_iterator<vertex_id_t>(0), tbb::counting_iterator<vertex_id_t>(M), [&](auto hyperE) {
      if (hyperedgedegrees[hyperE] >= s) frontier.push_back(hyperE);
    });
    
    nw::graph::edge_list<edge_directedness> two_graph(0);
    two_graph.open_for_push_back();
    nw::graph::AtomicBitVector   visitedE(M);
    using label_map_t = tbb::concurrent_hash_map<vertex_id_t, vertex_id_t>;
    label_map_t relabel_map;
    std::atomic<vertex_id_t> index = 0;
    auto relabel = [&](vertex_id_t &old) {
      vertex_id_t result;
      label_map_t::const_accessor read;
      if (false == relabel_map.find(read, old)) {
        //if old has not been relabeled
        result = index++;
        label_map_t::accessor write;
        relabel_map.insert(write, std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(old), std::forward<vertex_id_t>(result)));
      }
      else
        result = read->second;
      return result;
    };

    std::mutex mtx;
    std::atomic<size_t> nvisited = 0, nbelowS = 0, nduplicates = 0, nintersections = 0, nedges = 0;
    std::for_each(ep, frontier.begin(), frontier.end(), [&](auto&& hyperE) {
      std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
        auto hyperN = std::get<0>(x);
        std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
          //so we check compid of each hyperedge
          auto anotherhyperE = std::get<0>(y);
          ++nvisited;
          //travese upper triangluar
          if (hyperE >= anotherhyperE) return;
          ++nbelowS;
          //filter edges deg(e) < s
          if (hyperedgedegrees[anotherhyperE] < s) return;
          ++nduplicates;
          //avoid duplicate intersections
          if (1 == visitedE.atomic_get(anotherhyperE)) return; else visitedE.atomic_set(anotherhyperE);
          ++nintersections;
          //O(average degree of hyperedges)
          if (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE])) {
            ++nedges;
            vertex_id_t newx = relabel(hyperE), newy = relabel(anotherhyperE);
            //std::cout << hyperE << " " << anotherhyperE << " into " << newx << " " << newy << std::endl;
            mtx.lock();
            two_graph.push_back(newx, newy);
            mtx.unlock();
          }
        });
      });
      visitedE.clear();
    });
    std::cout << nvisited << " visits, "
    << nbelowS << " hyperedges has less S nodes, "
    << nduplicates << " duplicate intersections avoid,"
    << nintersections << " intersections performed, " 
    << nedges << " edges added" << std::endl;


    two_graph.close_for_push_back();
    return two_graph;
    /*
    //after finding the line graph, relabel the line graph
    //such that the id is consecutive
    nw::graph::edge_list<edge_directedness> relabel_graph(0);
    relabel_graph.open_for_push_back();
    std::for_each(two_graph.begin(), two_graph.end(), [&](auto&& elt) {
      auto&& [x, y] = elt;
      if (0 == relabel_map.count(x)) {
        //if x has not been relabeled
        relabel_map[x] = index;
        index++;
      }
      auto newx = relabel_map[x];
      if (0 == relabel_map.count(y)) {
        //if y has not been relabeled
        relabel_map[y] = index;
        index++;
      }
      auto newy = relabel_map[y];
      //std::cout << x << " " << y << " into " << newx << " " << newy << std::endl;
      relabel_graph.push_back(newx, newy);
    });

    relabel_graph.close_for_push_back();
    return relabel_graph;
    */
  }//else
}

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv6(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  nw::util::life_timer _(__func__);

  size_t M = e_nbs.size();
  size_t N = n_nbs.size();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    //avoid intersection when s=1
    std::atomic<size_t> nedges = 0;
    std::vector<nw::graph::edge_list<edge_directedness>> two_graphs(num_bins, nw::graph::edge_list<edge_directedness>(0));
    std::for_each(ep, tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].open_for_push_back();
    });
      
    nw::graph::AtomicBitVector   visitedE(M);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            //if (hyperE >= anotherhyperE) return;
            if (1 == visitedE.atomic_get(anotherhyperE)) return; else visitedE.atomic_set(anotherhyperE);
            ++nedges;
            two_graphs[worker_index].push_back(hyperE, anotherhyperE);
          });
        });
        visitedE.clear();
      } //for
    }, tbb::auto_partitioner());
    std::cout << nedges << " edges added" << std::endl;
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].close_for_push_back();
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }
  else {
    //when s > 1
    using label_map_t = tbb::concurrent_hash_map<vertex_id_t, vertex_id_t>;
    label_map_t relabel_map;
    std::atomic<vertex_id_t> index = 0;
    auto relabel = [&](vertex_id_t &old) {
      vertex_id_t result;
      label_map_t::const_accessor read;
      if (false == relabel_map.find(read, old)) {
        //if old has not been relabeled
        result = index++;
        label_map_t::accessor write;
        relabel_map.insert(write, std::make_pair<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(old), std::forward<vertex_id_t>(result)));
      }
      else
        result = read->second;
      return result;
    };
    //create an array of line graphs for each thread
    std::vector<nw::graph::edge_list<edge_directedness>> two_graphs(num_bins, nw::graph::edge_list<edge_directedness>(0));
    std::for_each(ep, tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].open_for_push_back();
    });
    
    nw::graph::AtomicBitVector   visitedE(M);
    std::atomic<size_t> nvisited = 0, nbelowS = 0, nduplicates = 0, nintersections = 0, nedges = 0;
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        if (hyperedgedegrees[hyperE] < s) continue;
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            ++nvisited;
            //travese upper triangluar
            if (hyperE >= anotherhyperE) return;
            ++nbelowS;
            //filter edges deg(e) < s
            if (hyperedgedegrees[anotherhyperE] < s) return;
            ++nduplicates;
            //avoid duplicate intersections
            if (1 == visitedE.atomic_get(anotherhyperE)) return; else {
              visitedE.atomic_set(anotherhyperE);
            }
            ++nintersections;
            //O(average degree of hyperedges)
            //if (s <= nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE])) {
            if (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) {
              ++nedges;
              vertex_id_t newx = relabel(hyperE), newy = relabel(anotherhyperE);
              //std::cout << hyperE << " " << anotherhyperE << " into " << newx << " " << newy << std::endl;
              two_graphs[worker_index].push_back(newx, newy);
            }
          });
        });
        visitedE.clear();
      } //for
    }, tbb::auto_partitioner());
    std::cout << nvisited << " visits, "
    << nbelowS << " hyperedges has less S nodes, "
    << nduplicates << " duplicate intersections avoid,"
    << nintersections << " intersections performed, " 
    << nedges << " edges added" << std::endl;
    
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].close_for_push_back();
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto&& e){
        result.push_back(e);
      });
    });
    result.close_for_push_back();

    return result;
  }//else
}


template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv7(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, 
std::vector<index_t>& hyperedgedegrees, size_t s = 1, int num_bins = 32) {
  nw::util::life_timer _(__func__);

  size_t M = e_nbs.size();
  size_t N = n_nbs.size();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();

  if (1 == s) {
    //avoid intersection when s=1
    std::atomic<size_t> nedges = 0;
    std::vector<nw::graph::edge_list<edge_directedness>> two_graphs(num_bins, nw::graph::edge_list<edge_directedness>(0));
    std::for_each(ep, tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].open_for_push_back();
    });
      
    nw::graph::AtomicBitVector   visitedE(M);
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      for (auto hyperE = r.begin(), e = r.end(); hyperE < e; ++hyperE) {
        //all neighbors of hyperedges are hypernode
        std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto &&x) { //O(average degree of hyperedges)
          auto hyperN = std::get<0>(x);
          std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto &&y) { //O(average degree of hypernodes)
            //so we check compid of each hyperedge
            auto anotherhyperE = std::get<0>(y);
            //if (hyperE >= anotherhyperE) return;
            if (1 == visitedE.atomic_get(anotherhyperE)) return; else visitedE.atomic_set(anotherhyperE);
            ++nedges;
            two_graphs[worker_index].push_back(hyperE, anotherhyperE);
          });
        });
        visitedE.clear();
      } //for
    }, tbb::auto_partitioner());
    std::cout << nedges << " edges added" << std::endl;
    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i){
      two_graphs[i].close_for_push_back();
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
    std::vector<std::vector<std::pair<vertex_id_t, vertex_id_t>>> two_graphs(num_bins);
    std::atomic<size_t> nvisited = 0, nbelowS = 0, nduplicates = 0, nintersections = 0, nedges = 0;
    tbb::parallel_for(tbb::blocked_range<vertex_id_t>(0, M), [&](tbb::blocked_range<vertex_id_t>& r) {
      int worker_index = tbb::task_arena::current_thread_index();
      //nw::graph::AtomicBitVector visitedE(M);
      std::vector<bool> visitedE(M, false);
      //std::cout << worker_index << ": " <<  r.begin() << " " <<  r.end() << std::endl;
      for (auto hyperE = r.begin(), e = r.end(); hyperE != e; ++hyperE) {
        //visitedE.clear();
        std::fill(visitedE.begin(), visitedE.end(), false);
        if (hyperedgedegrees[hyperE] < s) continue;
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
            if (hyperedgedegrees[anotherhyperE] < s) return;

            ++nduplicates;
            //avoid duplicate intersections
            if (1 == visitedE[anotherhyperE]) return; else {
              visitedE[anotherhyperE] = 1;
            }
            
            ++nintersections;
            //O(average degree of hyperedges)
            if (is_intersection_size_s(e_nbs[hyperE].begin(), e_nbs[hyperE].end(),
            e_nbs[anotherhyperE].begin(), e_nbs[anotherhyperE].end(), s)) {
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

    nw::graph::edge_list<edge_directedness> result(0);
    result.open_for_push_back();
    //do this in serial
    vertex_id_t index = 0;
    std::unordered_map<vertex_id_t, vertex_id_t> relabel_map;
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto &&e) {
        auto &&[x, y] = e;
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
        result.push_back(newx, newy);
      });
    });
    result.close_for_push_back();

    return result;
  }//else
}



}//namespace hypergraph
}//namespace nw