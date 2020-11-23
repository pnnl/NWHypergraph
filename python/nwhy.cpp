//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2020
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//


#include <pybind11/pybind11.h>
#include <pybind11/stl.h> //auto copy between stl containers and python data structures
#include <pybind11/stl_bind.h>
#include <pybind11/cast.h>
#include "convert_to_s_overlap.hpp"
#include <pybind11/stl_bind.h>
#include <vector>

namespace py = pybind11;
using namespace nw::hypergraph;

using Index_t = int;
using Data_t = int;
//PYBIND11_MAKE_OPAQUE(py::array_t<T, py::array::c_style | py::array::forcecast>);


PYBIND11_MODULE(nwhy, m) {
    m.doc() = "NWhy pybind11 module plugin"; // optional module docstring

    py::class_<NWHypergraph<Index_t, Data_t>> hypergraph_class(m, "NWHypergraph");
    hypergraph_class

    .def(py::init<>([](py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y) {
        return new NWHypergraph<Index_t, Data_t>(x, y);
    }))
    .def(py::init<>([](py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y,
    py::array_t<Data_t, py::array::c_style | py::array::forcecast> &data) {
        return new NWHypergraph<Index_t, Data_t>(x, y, data);
    }))
    .def("s_linegraph", &NWHypergraph<Index_t, Data_t>::s_linegraph, "A function which converts a hypergraph to its s line graph",
    py::arg("s") = 1, py::arg("edges") = true)
    //s_connected_component
    .def("s_connected_component", py::overload_cast<Slinegraph<Index_t, Data_t> &, bool>(&NWHypergraph<Index_t, Data_t>::s_connected_component),
     "A function which finds the connected components for its s line graph",
    py::arg("linegraph"), py::arg("return_singleton") = false)
    .def("s_connected_component", py::overload_cast<int, bool, bool>(&NWHypergraph<Index_t, Data_t>::s_connected_component), 
    "A function which finds the connected components for its s line graph",
    py::arg("s") = 1, py::arg("edges") = true, py::arg("return_singleton") = false)
    //s_distance
    .def("s_distance", py::overload_cast<Slinegraph<Index_t, Data_t> &, Index_t, Index_t>(&NWHypergraph<Index_t, Data_t>::s_distance),
     "A function which computes the distance from src to dest in its s line graph",
    py::arg("linegraph"), py::arg("src"), py::arg("dest"))
    .def("s_distance", py::overload_cast<Index_t, Index_t, int, bool>(&NWHypergraph<Index_t, Data_t>::s_distance), 
    "A function which computes the distance from src to dest in its s line graph",
    py::arg("src"), py::arg("dest"), py::arg("s") = 1, py::arg("edges") = true)
    //s_neighbor
    .def("s_neighbor", py::overload_cast<Slinegraph<Index_t, Data_t> &, Index_t>(&NWHypergraph<Index_t, Data_t>::s_neighbor),
     "A function which finds the neighbors for vertex v of its s line graph",
    py::arg("linegraph"), py::arg("v"))
    .def("s_neighbor", py::overload_cast<Index_t, int, bool>(&NWHypergraph<Index_t, Data_t>::s_neighbor), 
    "A function which finds the neighbors for vertex v of its s line graph",
    py::arg("v"), py::arg("s") = 1, py::arg("edges") = true);

    py::class_<Slinegraph<Index_t, Data_t>> slinegraph_class(m, "Slinegraph");
    slinegraph_class
    .def(py::init<>([](NWHypergraph<Index_t, Data_t>& g, int s, bool edges) {
        return new Slinegraph<Index_t, Data_t>(g, s, edges);
    }), "Init function", py::arg("g"), py::arg("s") = 1, py::arg("edges") = true)
    .def("_s", &Slinegraph<Index_t, Data_t>::getS)
    .def("s_connected_component", &Slinegraph<Index_t, Data_t>::s_connected_component,
     "A function which finds the connected components for its s line graph",
    py::arg("return_singleton") = false)
    .def("s_distance", &Slinegraph<Index_t, Data_t>::s_distance,
    "A function to compute the distance from src to dest",
    py::arg("src"), py::arg("dest"))
    .def("s_neighbor", &Slinegraph<Index_t, Data_t>::s_neighbor,
    "A function to get s neighbors of a vertex",
    py::arg("v"));

    //register version information in a module as below
    py::object version = py::cast("0.0.3");
    m.attr("_version") = version;

    //define function, its argument list, and with default argument for s
    m.def("convert_to_s_overlap", &convert_to_s_overlap<Index_t, Data_t>, "A function which converts a hypergraph to its s line graph",
    py::arg("x"), py::arg("y"), py::arg("data"), py::arg("s") = 1);
    //m.def("connected_component", &nw::graph::cc);
}