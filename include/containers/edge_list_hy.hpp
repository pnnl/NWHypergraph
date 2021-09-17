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
#include <containers/edge_list.hpp>
#include <util/timer.hpp>
namespace nw {
namespace hypergraph {

template<int idx, directedness edge_directedness = directed, class Vector = std::vector<int>>
auto relabel_by_degree(nw::graph::edge_list<edge_directedness>& aos_a, std::string direction = "ascending", Vector&& degree = std::vector<int>(0)) {
  nw::util::life_timer _(__func__);
  int status = -4;
  aos_a.prv.push_back(nw::graph::demangle(typeid(*&aos_a).name(), nullptr, nullptr, &status) + "::" + __func__,
            "index " + std::to_string(idx) + " " + direction);

  std::vector<vertex_id_t> perm = (0 == degree.size()) ? aos_a.template perm_by_degree<idx>(direction) : aos_a.template perm_by_degree<idx>(degree, direction);

  aos_a.template relabel<idx>(perm);
  
  return perm;
}


} //namespace hypergraph
} //namespace nw