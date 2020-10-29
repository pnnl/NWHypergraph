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

template<class Index_t, class Data_t>
std::tuple<py::array_t<Index_t>, py::array_t<Index_t>, py::array_t<Data_t>, 
py::array_t<Index_t>, py::array_t<Index_t>, py::array_t<Data_t>> 
convert_to_s_overlap(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y, 
py::array_t<Data_t, py::array::c_style | py::array::forcecast> &data, 
size_t s = 1) {
    //sanitize check
    auto rx = x.template mutable_unchecked<1>();
    auto ry = y.template mutable_unchecked<1>();
    auto rdata = data.template mutable_unchecked<1>();
    //rx(0) = 1;
    size_t n_x = x.shape(0);
    size_t n_y = y.shape(0);
    //assume there are n_x pairs
    std::cout << "n_col=" << n_x << " n_row=" << n_y << std::endl;
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
        auto&& two_graphs = 
        to_two_graph_efficient_parallel_clean_without_sequeeze<nw::graph::undirected>(std::execution::par_unseq, hyperedges, hypernodes, hyperedge_degrees, s, num_bins);
        if (1 == s) {
            nw::graph::edge_list<nw::graph::undirected> &&linegraph = create_edgelist_without_squeeze<nw::graph::undirected>(two_graphs);
            //where when an empty edge list is passed in, an adjacency still have two elements
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n}); 
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    ++i;
                }
                return std::make_tuple(newx, newy, py::array_t<Data_t>(0), 
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));
        }
        else {
            nw::graph::edge_list<nw::graph::undirected> &&linegraph = create_edgelist_with_squeeze<nw::graph::undirected>(two_graphs);
            nw::graph::edge_list<nw::graph::undirected> &&raw_linegraph = create_edgelist_without_squeeze<nw::graph::undirected>(two_graphs);           
            //where when an empty edge list is passed in, an adjacency still have two elements
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n}), oldx({n}), oldy({n}); 
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            auto roldx = oldx.template mutable_unchecked<1>();
            auto roldy = oldy.template mutable_unchecked<1>();
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    ++i;
                }
                i = 0;
                for (auto &&[u, v] : raw_linegraph) {
                    roldx(i) = u;
                    roldy(i) = v;
                    ++i;
                }
                return std::make_tuple(newx, newy, py::array_t<Data_t>(0), oldx, oldy, py::array_t<Data_t>(0));
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));
        }//else 1 == s
    }
    else {
        std::cout << "weighted hypergraph" << std::endl;

        nw::graph::edge_list<nw::graph::directed, Data_t> aos_a(0); 
        aos_a.open_for_push_back();
        for (size_t i = 0; i < n_x; ++i) {
            aos_a.push_back(rx(i), ry(i), rdata(i));
        }
        aos_a.close_for_push_back();

        nw::graph::adjacency<0, Data_t> hyperedges(aos_a);
        nw::graph::adjacency<1, Data_t> hypernodes(aos_a);
        std::vector<nw::graph::index_t> hyperedge_degrees = aos_a.template degrees<0>();
    
        int num_bins = 32;
        auto&& two_graphs = 
        to_two_graph_weighted_efficient_parallel_clean_without_squeeze<nw::graph::undirected, Data_t>(std::execution::par_unseq, hyperedges, hypernodes, hyperedge_degrees, s, num_bins);
        if (1 == s) {
            nw::graph::edge_list<nw::graph::undirected, Data_t> &&linegraph = create_edgelist_without_squeeze<nw::graph::undirected, Data_t>(two_graphs);
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            //create results
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n});
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            py::array_t<Data_t, py::array::c_style> newdata({n}); 
            auto rnewdata = newdata.template mutable_unchecked<1>(); 
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v, weight] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    rnewdata(i) = weight;
                    ++i;
                }
                return std::make_tuple(newx, newy, newdata, 
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Index_t>(0));
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));   
        }
        else {
            nw::graph::edge_list<nw::graph::undirected, Data_t> &&linegraph = create_edgelist_with_squeeze<nw::graph::undirected, Data_t>(two_graphs);
            nw::graph::edge_list<nw::graph::undirected, Data_t> &&raw_linegraph = create_edgelist_without_squeeze<nw::graph::undirected, Data_t>(two_graphs);           
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n}), oldx({n}), oldy({n}); 
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            auto roldx = oldx.template mutable_unchecked<1>();
            auto roldy = oldy.template mutable_unchecked<1>();
            py::array_t<Data_t, py::array::c_style> newdata({n}), olddata({n}); 
            auto rnewdata = newdata.template mutable_unchecked<1>(); 
            auto rolddata = olddata.template mutable_unchecked<1>(); 
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v, weight] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    rnewdata(i) = weight;
                    ++i;
                }
                i = 0;
                for (auto &&[u, v, weight] : raw_linegraph) {
                    roldx(i) = u;
                    roldy(i) = v;
                    rolddata(i) = weight;
                    ++i;
                }
                return std::make_tuple(newx, newy, newdata, oldx, oldy, olddata);                  
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));         
        }//else 1 == s
    }//else 0 == n_data
}

}//namespace hypergraph
}//namespace nw