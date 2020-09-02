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

namespace nw {
namespace hypergraph {

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graph(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();
  uint64_t counter = 0;
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

template<directedness edge_directedness = undirected, class ExecutionPolicy, class HyperEdge, class HyperNode>
auto to_two_graphv2(ExecutionPolicy&& ep, HyperEdge& e_nbs, HyperNode& n_nbs, std::vector<index_t>& hyperedgedegrees, size_t s = 1) {
  nw::util::life_timer _(__func__);
  nw::graph::edge_list<edge_directedness> two_graph(0);
  two_graph.open_for_push_back();

  auto edges = e_nbs.begin();
  auto nodes = n_nbs.begin();
  uint64_t counter = 0;
  //O(n*average degree of hyperedges*average degree of hypernodes*average degree of hyperedges)
  for (size_t i = 0; i < e_nbs.size(); ++i) { //O(n)
    auto hyperE = i;
    if (hyperedgedegrees[hyperE] < s) continue;
    std::for_each(edges[hyperE].begin(), edges[hyperE].end(), [&](auto&& x) { //O(average degree of hyperedges)
      auto hyperN = std::get<0>(x);
      std::for_each(nodes[hyperN].begin(), nodes[hyperN].end(), [&](auto&& y) { //O(average degree of hypernodes)
        //so we check compid of each hyperedge
        auto anotherhyperE = std::get<0>(y);
        if (hyperE > anotherhyperE) return;
        ++counter;
        //O(average degree of hyperedges)
        size_t count = nw::graph::intersection_size(e_nbs[hyperE], e_nbs[anotherhyperE]);    //       if intersection_size(n_nbs(i), n_nbs(j)) >= s
        if (count >= s) {
          two_graph.push_back(hyperE, anotherhyperE);    //         add (i,j) to new edge_list
        }
      });
    });
  }
  std::cout << counter << " intersections performed" << std::endl;

  two_graph.close_for_push_back();
  return two_graph;
}

}//namespace hypergraph
}//namespace nw