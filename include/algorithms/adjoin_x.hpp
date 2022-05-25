/**
 * @file adjoin_x.hpp
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

#include <nwgraph/util/timer.hpp>
#include <tuple>
#include <vector>
#include <nwgraph/adaptors/vertex_range.hpp>

namespace nw {
namespace hypergraph {

template<class T>
std::tuple<std::vector<T>, std::vector<T>> splitLabeling(std::vector<T>& labeling, const size_t num_realedges, const size_t num_realnodes) {
  std::vector<T> E(num_realedges), N(num_realnodes);
  if (num_realnodes < num_realedges) {
    nw::util::life_timer _("unrelabeling");
    //std::move, vector. assign, std::copy all take O(N)
    auto it = std::next(labeling.begin(), num_realedges);
    //instead of copy(vector.assign), use std::move
    std::move(labeling.begin(), it, std::back_insert_iterator(E));
    //E.assign(labeling.begin(), labeling.begin() + num_realedges);
    std::move(it, labeling.end(), std::back_insert_iterator(N));
    //N.assign(labeling.begin() + num_realedges, labeling.end());
  }
  else {
    nw::util::life_timer _("unrelabeling");
    auto it = std::next(labeling.begin(), num_realnodes);
    std::move(labeling.begin(), it, std::back_insert_iterator(N));
    //N.assign(labeling.begin(), labeling.begin() + num_realnodes);
    std::move(it, labeling.end(), std::back_insert_iterator(E));
    //E.assign(labeling.begin() + num_realnodes, labeling.end());
  }
  return std::tuple(N, E);
}

template<class ExecutionPolicy, class T>
std::tuple<std::vector<T>, std::vector<T>> splitLabeling(ExecutionPolicy& ep, std::vector<T>& labeling, const size_t num_realedges, const size_t num_realnodes) {
  std::vector<T> E(num_realedges), N(num_realnodes);
  if (num_realnodes < num_realedges) {
    nw::util::life_timer _("unrelabeling");
    //E.assign(labeling.begin(), labeling.begin() + num_realedges);
    std::for_each(ep, nw::graph::counting_iterator(0ul), nw::graph::counting_iterator(num_realedges), [&](auto i) {
      E[i] = labeling[i];
    });
    //N.assign(labeling.begin() + num_realedges, labeling.end());
    std::for_each(ep, nw::graph::counting_iterator(0ul), nw::graph::counting_iterator(num_realnodes), [&](auto i) {
      N[i] = labeling[i + num_realedges];
    }); 
  }
  else {
    nw::util::life_timer _("unrelabeling");
    std::for_each(ep, nw::graph::counting_iterator(0ul), nw::graph::counting_iterator(num_realnodes), [&](auto i) {
      N[i] = labeling[i];
    }); 
    //E.assign(labeling.begin() + num_realnodes, labeling.end());
    std::for_each(ep, nw::graph::counting_iterator(0ul), nw::graph::counting_iterator(num_realedges), [&](auto i) {
      E[i] = labeling[i + num_realnodes];
    });
  }
  return std::tuple(N, E);
}

/*
* X is the algorithm X. T is the type of the elements in the vector function f returns.
* Parameters is parameter list of the algorithm X, which is a parameter pack that
* contains zero or more function argurments.
*/
template<class X, class T, class... Args>
std::tuple<std::vector<T>, std::vector<T>> 
relabel_x(const size_t num_realedges, const size_t num_realnodes, X f, Args&&... args) {
  std::vector<T>&& labeling = f(std::forward<Args>(args)...);
  return splitLabeling<T>(labeling, num_realedges, num_realnodes);
}

template<class ExecutionPolicy, class X, class T, class... Args>
std::tuple<std::vector<T>, std::vector<T>> 
relabel_x_parallel(ExecutionPolicy& ep, const size_t num_realedges, const size_t num_realnodes, X f, Args&&... args) {
  std::vector<T>&& labeling = f(std::forward<Args>(args)...);
  return splitLabeling<ExecutionPolicy, T>(ep, labeling, num_realedges, num_realnodes);
}

}//namespace hypergraph
}//namespace nw