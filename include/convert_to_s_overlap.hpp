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
#include <iostream>
#include <vector>
#include <tuple>
#include <edge_list.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "s_overlap.hpp"
using namespace nw::graph;
namespace py = pybind11;

namespace nw {
namespace hypergraph {

template<class T>
std::tuple<py::array_t<T>, py::array_t<T>, py::array_t<T>> 
convert_to_s_overlap(py::array_t<T, py::array::c_style | py::array::forcecast> &x, 
py::array_t<T, py::array::c_style | py::array::forcecast> &y, 
py::array_t<T, py::array::c_style | py::array::forcecast> &data, size_t s = 1) {
    //sanitize check
    auto rx = x.template mutable_unchecked<1>();
    auto ry = y.template mutable_unchecked<1>();
    auto rdata = data.template mutable_unchecked<1>();
    //rx(0) = 1;
    size_t n_x = x.shape(0);
    size_t n_y = y.shape(0);
    //assume there are n_x pairs
    assert(n_x == n_y);
    size_t n_data = data.shape(0);

    //create edge list 
    if (0 == n_data) {
        std::cout << "unweighted hypergraph" << std::endl;
        nw::graph::edge_list<nw::graph::directed> aos_a(0); 
        aos_a.open_for_push_back();
        for (size_t i = 0; i < n_x; ++i) {
            aos_a.push_back(rx(i), ry(i));
        }
        aos_a.close_for_push_back();

        nw::graph::adjacency<0> hyperedges(aos_a);
        nw::graph::adjacency<1> hypernodes(aos_a);
        std::vector<nw::graph::index_t> hyperedge_degrees = aos_a.degrees<0>();

        int num_bins = 32;
        nw::graph::edge_list<nw::graph::undirected> &&linegraph = 
        to_two_graph_efficient_parallel_clean<nw::graph::undirected>(std::execution::par_unseq, hyperedges, hypernodes, hyperedge_degrees, s, num_bins);
        //where when an empty edge list is passed in, an adjacency still have two elements
        size_t n = linegraph.size();
        std::cout << "linegraph has " << n << " edges" << std::endl;
        py::array_t<T> newx({1}), newy({1});
        auto rnewx = newx.template mutable_unchecked<1>();
        auto rnewy = newy.template mutable_unchecked<1>();
        if (0 == n) return std::make_tuple(newx, newy, py::array_t<T>());
        size_t i = 0;
        for (auto&&[u, v] : linegraph) {
            rnewx(i) = u;
            rnewy(i) = v;
            ++i;
        }
        return std::make_tuple(newx, newy, py::array_t<T>());
    }
    else {
        std::cout << "weighted hypergraph" << std::endl;
        nw::graph::edge_list<nw::graph::directed, T> aos_a(0); 
        aos_a.open_for_push_back();
        for (size_t i = 0; i < n_x; ++i) {
            aos_a.push_back(rx(i), ry(i), rdata(i));
        }
        aos_a.close_for_push_back();

        nw::graph::adjacency<0, T> hyperedges(aos_a);
        nw::graph::adjacency<1, T> hypernodes(aos_a);
        std::vector<nw::graph::index_t> hyperedge_degrees = aos_a.template degrees<0>();
    
        int num_bins = 32;
        nw::graph::edge_list<nw::graph::undirected, T> &&linegraph = 
        to_two_graph_weighted_efficient_parallel_clean<nw::graph::undirected, T>(std::execution::par_unseq, hyperedges, hypernodes, hyperedge_degrees, s, num_bins);
        //where when an empty edge list is passed in, an adjacency still have two elements
        size_t n = linegraph.size();
        std::cout << "linegraph has " << n << " edges" << std::endl;
        py::array_t<T, py::array::c_style | py::array::forcecast> newx({1}), newy({1}), newdata({1});
        auto rnewx = newx.template mutable_unchecked<1>();
        auto rnewy = newy.template mutable_unchecked<1>();
        auto rnewdata = newdata.template mutable_unchecked<1>();
        if (0 == n) return std::make_tuple(newx, newy, newdata);
        size_t i = 0;
        for (auto&&[u, v, weight] : linegraph) {
            std::cout << i << ": " << u << " " << v << " " << weight << std::endl;
            rnewx(i) = u;
            rnewy(i) = v;
            rnewdata(i) = weight;
            ++i;
        }
        return std::make_tuple(newx, newy, newdata);
    }
}

}//namespace hypergraph
}//namespace nw