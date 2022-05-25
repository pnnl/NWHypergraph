/**
 * @file toplexes_eg.cpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#include <nwgraph/adaptors/neighbor_range.hpp>
#include <algorithm>
#include <nwgraph/edge_list.hpp>
#include <nwgraph/adjacency.hpp>

using namespace nw::graph;


int main(int argc, char* argv[]) {
  size_t n_vtx = 5;

  bi_edge_list<directedness::directed, int> A_list(n_vtx);
  A_list.open_for_push_back();
  A_list.push_back(0, 1, 1);
  A_list.push_back(1, 2, 2);
  A_list.push_back(2, 3, 3);
  A_list.push_back(3, 4, 4);
  A_list.close_for_push_back();

  {
    biadjacency<0, int> A(A_list);
    biadjacency<1, int> AT(A_list);

    std::size_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(A)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    std::cout << std::endl;

    std::size_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(AT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    std::cout << std::endl;
  }

  {
    biadjacency<0, int> B(A_list);
    biadjacency<1, int> BT(A_list);

    std::size_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(B)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    std::cout << std::endl;

    std::size_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(BT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    std::cout << std::endl;
  }

  return 0;
}