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
#include <containers/compressed.hpp>
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
auto relabel_by_degree(adjacency<0, Attributes...>& edges, adjacency<1, Attributes...>& nodes, 
int idx, std::string direction) {
  nw::util::life_timer _("relabel_by_degree");
  if (0 == idx) {
    auto &&perm = edges.permute_by_degree(direction, std::execution::par_unseq);
    assert(perm.size() == edges.size());
    nodes.relabel_to_be_indexed(perm, std::execution::par_unseq);
    std::cout << "relabeling edge adjacency by degree..." << std::endl;
    return perm;
  }
  else {
    auto &&perm = nodes.permute_by_degree(direction, std::execution::par_unseq);
    assert(perm.size() == nodes.size());
    edges.relabel_to_be_indexed(perm, std::execution::par_unseq);
    std::cout << "relabeling nodes adjacency by degree..." << std::endl;
    return perm;
  }
}

} //namespace hypergraph
} //namespace nw