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
#include "edge_list.hpp"

namespace nw {
namespace hypergraph {

template<int idx = 0, class Vector = std::vector<int>>
void relabel_by_degree_bipartite(nw::graph::edge_list<nw::graph::directed>& aos_a, std::string direction = "ascending", Vector&& degree = std::vector<int>(0)) {
    int status = -4;
    aos_a.prv.push_back(nw::graph::demangle(typeid(*&aos_a).name(), nullptr, nullptr, &status) + "::" + __func__,
                  "index " + std::to_string(idx) + " " + direction);

    std::vector<vertex_id_t> perm = degree.size() == 0 ? aos_a.perm_by_degree<idx>(direction) : aos_a.perm_by_degree<idx>(degree, direction);

    std::vector<vertex_id_t> iperm(perm.size());


    tbb::parallel_for(tbb::blocked_range(0ul, iperm.size()), [&](auto&& r) {
      for (auto i = r.begin(), e = r.end(); i != e; ++i) {
        iperm[perm[i]] = i;
      }
    });

/*
    //std::for_each(std::execution::par_unseq, base::begin(), base::end(),[&](auto&& x) {
    //    std::get<idx>(x) = iperm[std::get<idx>(x)];
    //});
    */
  }

} //namespace hypergraph
} //namespace nw