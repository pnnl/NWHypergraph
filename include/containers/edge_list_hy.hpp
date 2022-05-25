/**
 * @file edge_list_hy.hpp
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
#include <nwgraph/build.hpp>
#include <nwgraph/edge_list.hpp>
#include <nwgraph/util/timer.hpp>
namespace nw {
namespace hypergraph {

template <int idx,
          nw::graph::edge_list_graph edge_list_t,
          class Vector = std::vector<int>>
requires(false == nw::graph::is_unipartite<typename edge_list_t::bipartite_graph_base>::value)
auto relabel_by_degree(edge_list_t& aos_a,
                       std::string direction = "ascending",
                       Vector&& degree = std::vector<int>(0)) {
  nw::util::life_timer _(__func__);
  /*
  int status = -4;
  aos_a.prv.push_back(
      nw::graph::demangle(typeid(*&aos_a).name(), nullptr, nullptr, &status) +
          "::" + __func__,
      "index " + std::to_string(idx) + " " + direction);
  */
  auto&& perm =
      (0 == degree.size())
          ? perm_by_degree<idx>(aos_a, direction)
          : perm_by_degree<idx>(aos_a, degree, direction);

  relabel<idx>(aos_a, perm);

  return perm;
}

} //namespace hypergraph
} //namespace nw