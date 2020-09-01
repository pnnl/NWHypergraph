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

#include <algorithms/connected_components.hpp>

namespace nw {
namespace hypergraph {

template<class ExecutionPolicy, class HyperNode, class SGraph>
auto base_two(ExecutionPolicy&& ep, HyperNode& hypernodes, SGraph& s_adj) {
  auto E = ccv1(s_adj);
  auto nhypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> N(nhypernodes);
  //for each hypernode, find N[i]
  for (size_t hyperN = 0; hyperN < nhypernodes; ++hyperN) {
    //get the id of the first neighbor (hyperedge) of hypernode hyperN
    auto hyperE = std::get<0>(*hypernodes[hyperN].begin());
    N[hyperN] = E[hyperE];
  }
  return std::tuple(N, E);
}

template<class ExecutionPolicy, class HyperGraph>
auto to_relabel_graph(ExecutionPolicy&& ep, HyperGraph& aos_a) {

  auto n_nbs = aos_a.max()[1] + 1;
  auto e_nbs = aos_a.max()[0] + 1;
  edge_list<undirected> relabel_graph(0);
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

template<class ExecutionPolicy, class HyperGraph>
auto relabelHyperCC(ExecutionPolicy&& ep, HyperGraph& aos_a) {
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

  auto n_nbs = aos_a.max()[1] + 1;
  auto e_nbs = aos_a.max()[0] + 1;
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
