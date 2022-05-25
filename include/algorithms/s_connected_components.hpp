/**
 * @file s_connected_components.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */


#pragma once

#include <nwgraph/algorithms/connected_components.hpp>
#include <nwgraph/experimental/algorithms/connected_components.hpp>

namespace nw {
namespace hypergraph {

/*
* soverlap cc using label propagation cc
*/
template<class ExecutionPolicy, class HyperNode, class SGraph>
auto linegraph_ccv1(ExecutionPolicy&& ep, HyperNode& hypernodes, SGraph& s_adj) {
  nw::util::life_timer _(__func__);
  auto E = ccv1(s_adj);
  return E;
  /*
  //if no component found, then return an empty pair
  if (E.empty()) return std::tuple(E, E);
  size_t nhypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> N(nhypernodes);
  //for each hypernode, find N[i]
  std::for_each(ep, tbb::counting_iterator<vertex_id_t>(0ul), tbb::counting_iterator<vertex_id_t>(nhypernodes), [&](auto hyperN) {
    auto hyperE = std::get<0>(*hypernodes[hyperN].begin());
    N[hyperN] = E[hyperE];
  });
  return std::tuple(N, E);
  */
}
/*
* soverlap cc using Afforest
*/
template<class ExecutionPolicy, class HyperNode, class SGraph>
auto linegraph_Afforest(ExecutionPolicy&& ep, HyperNode& hypernodes, SGraph& s_adj) {
  nw::util::life_timer _(__func__);
  nw::graph::adjacency<1> s_adj_trans(0);
  auto E = Afforest(s_adj, s_adj_trans);
  return E;
  /*
  //if no component found, then return an empty pair
  if (E.empty()) return std::tuple(E, E);
  size_t nhypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> N(nhypernodes);
  //for each hypernode, find N[i]
  std::for_each(ep, tbb::counting_iterator<vertex_id_t>(0ul), tbb::counting_iterator<vertex_id_t>(nhypernodes), [&](auto hyperN) {
    auto hyperE = std::get<0>(*hypernodes[hyperN].begin());
    N[hyperN] = E[hyperE];
  });
  return std::tuple(N, E);
  */
}
/*
* soverlap cc using label propagation cc
*/
template<class ExecutionPolicy, class SGraph>
auto linegraph_lpcc(ExecutionPolicy&& ep, SGraph& s_adj) {
  nw::util::life_timer _(__func__);
  auto E = lpcc(ep, s_adj);
  return E;
  /*
  //if no component found, then return an empty pair
  if (E.empty()) return std::tuple(E, E);
  size_t nhypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> N(nhypernodes);
  //for each hypernode, find N[i]
  std::for_each(ep, tbb::counting_iterator<vertex_id_t>(0ul), tbb::counting_iterator<vertex_id_t>(nhypernodes), [&](auto hyperN) {
    auto hyperE = std::get<0>(*hypernodes[hyperN].begin());
    N[hyperN] = E[hyperE];
  });
  return std::tuple(N, E);
  */
}
/*
* Inefficient version of relabeling
* Copy the original graph while relabeling hyperedges or hypernodes
*/
template<class ExecutionPolicy, class HyperGraph>
auto to_relabel_graph(ExecutionPolicy&& ep, HyperGraph& aos_a) {
  nw::util::life_timer _(__func__);
  auto n_nbs = num_vertices(aos_a, 1);
  auto e_nbs = num_vertices(aos_a, 0);
  edge_list<nw::graph::directedness::undirected> relabel_graph(0);
  relabel_graph.open_for_push_back();
  //we relabel the smaller set (either hypernode or hyperedge)
  if (n_nbs < e_nbs) {
    std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
      auto&& [edge, node] = elt;
      node = node + e_nbs;
      relabel_graph.push_back(edge, node);    //         add (i,j) to new edge_list
    });
  }
  else {
    std::for_each(aos_a.begin(), aos_a.end(), [&](auto&& elt) {
      auto&& [edge, node] = elt;
      edge = edge + n_nbs;
      relabel_graph.push_back(edge, node);    //         add (i,j) to new edge_list
    });
  }

  relabel_graph.close_for_push_back();
  return relabel_graph;
}

template<class ExecutionPolicy, class HyperGraph, class vertex_id_t = vertex_id_t<HyperGraph>>
auto relabelHyperCC(ExecutionPolicy&& ep, HyperGraph& aos_a) {
  nw::util::life_timer _(__func__);
  auto relabel_g = to_relabel_graph(ep, aos_a);    // 1) find 2-graph corresponding to s-overlapped hyper edges
  //relabel_g. template lexical_sort_by<0>();
  //relabel_g.uniq();
  //relabel_g.remove_self_loops();

  //relabel_g.stream_edges();

  auto s_adj  = adjacency<0>(relabel_g);        // 2) convert new edge_list to new_adjacency
  //auto s_trans_adj = adjacency<1>(relabel_g);
  //s_adj.stream_indices();

  auto labeling   = //Afforest(s_adj, s_trans_adj); 
  ccv1(s_adj);//

  auto n_nbs = num_vertices(aos_a, 1);
  auto e_nbs = num_vertices(aos_a, 0);
  std::vector<vertex_id_t> N, E;
  if (n_nbs < e_nbs) {
    E.assign(labeling.begin(), labeling.begin() + e_nbs);
    N.assign(labeling.begin() + e_nbs, labeling.end());
  }
  else {
    N.assign(labeling.begin(), labeling.begin() + n_nbs);
    E.assign(labeling.begin() + n_nbs, labeling.end());
  }
  return std::tuple(N, E);
}

}//namespace hypergraph
}//namespace nw
