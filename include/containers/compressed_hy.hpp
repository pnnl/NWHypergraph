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
#include <compressed.hpp>
#include <util/timer.hpp>
namespace nw {
namespace hypergraph {
/*
* To relabel bi-adjacency of a hypergraph by degree in-place, there are multiple steps:
* 1. based on relabel edge or node, we relabel by permute the id of edges or nodes
* 2. we get the permutation of step 1 and relabel the neighbor list of the other part of the bi-adjacency
* 3. sort the neighbor list of the relabeled neighbor list in step 2
**/
template<class... Attributes>
void relabel_by_degree(adjacency<0, Attributes...>& edges, adjacency<1, Attributes...>& nodes, 
int idx, std::string direction) {
  nw::util::life_timer _("relabel_by_degree");
  const int jdx = (idx + 1) % 2;
  if (0 == idx) {
    auto &&iperm = edges.permute_by_degree(direction, std::execution::par_unseq);
    assert(iperm.size() == edges.size());
    nodes.relabel_to_be_indexed(iperm);
    std::cout << "relabeling edge adjacency by degree..." << std::endl;
  }
  else {
    auto &&iperm = nodes.permute_by_degree(direction, std::execution::par_unseq);
    assert(iperm.size() == nodes.size());
    edges.relabel_to_be_indexed(iperm);
    std::cout << "relabeling nodes adjacency by degree..." << std::endl;
  }
}

} //namespace hypergraph
} //namespace nw