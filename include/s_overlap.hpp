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

namespace nw {
namespace hypergraph {

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

}//namespace hypergraph
}//namespace nw