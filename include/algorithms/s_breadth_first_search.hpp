/**
 * @file s_breadth_first_search.hpp
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

#include <nwgraph/algorithms/bfs.hpp>

namespace nw {
namespace hypergraph {

template<class ExecutionPolicy, class HyperNode, class OutGraph, class InGraph, typename vertex_id_t = vertex_id_t<HyperNode>>
auto S_BFS_v0(ExecutionPolicy&& exec, vertex_id_t source_hyperedge, HyperNode& hypernodes, 
OutGraph& s_adj, InGraph& s_tran_adj, int num_bins = 32, int alpha = 15, int beta = 18) {
  std::vector<vertex_id_t> parentE = bfs_v11(s_adj, s_tran_adj, source_hyperedge, num_bins, alpha, beta);//Afforest(s_adj, s_trans_adj);                 // 3) run whatever on new_adjacency
  auto nhypernodes = hypernodes.max() + 1;
  std::vector<vertex_id_t> parentN(nhypernodes);
  //for each hypernode hyperN, find parentN[hyperN] with a bottomup round
  for (size_t hyperN = 0; hyperN < nhypernodes; ++hyperN) {
    std::for_each(hypernodes.begin()[hyperN].begin(), hypernodes.begin()[hyperN].end(), [&](auto&& x) {
        //so we check compid of each hyperedge
        auto hyperE = std::get<0>(x);
        parentN[hyperN] = hyperE;
        return;
    });
  }

  return std::tuple(parentN, parentE);
}

}//namespace hypergraph
}//namespace nw