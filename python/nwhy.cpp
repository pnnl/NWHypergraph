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
/*
* About weight, it is complicated.
* It has multiple different scenarios.
* Scenario 1: unweighted multi-hypergraph collapse into weighted simple hypergraph
* Scenario 2: weighted simple hypergraph
* Scenario 3: weighted simple hypergraph into unweighted slinegraph
* Scenario 4: unweighted simple hypergraph into weighted slinegraph
*/
using Data_t = int;

PYBIND11_MODULE(nwhy, m) {
    m.doc() = "NWhy pybind11 module plugin"; // optional module docstring

    //define NWHypergraph python object
    py::class_<NWHypergraph<Index_t, Data_t>> hypergraph_class(m, "NWHypergraph");
    hypergraph_class
    .def_readonly("row", &NWHypergraph<Index_t, Data_t>::row_)
    .def_readonly("col", &NWHypergraph<Index_t, Data_t>::col_)
    .def_readonly("data", &NWHypergraph<Index_t, Data_t>::data_)
    //constuctor
    .def(py::init<>([](py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y) {
        return new NWHypergraph<Index_t, Data_t>(x, y);
    }))
    .def(py::init<>([](py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y,
    py::array_t<Data_t, py::array::c_style | py::array::forcecast> &data,
    bool collapse = false) {
        return new NWHypergraph<Index_t, Data_t>(x, y, data, collapse);
    }),
    py::arg("x"), py::arg("y"), py::arg("data"), py::arg("collapse") = false)
    //stats for hypergraph
    .def("edge_size_dist", &NWHypergraph<Index_t, Data_t>::edge_size_dist,
    "A function to get the edge size distribution of the hypergraph")
    .def("node_size_dist", &NWHypergraph<Index_t, Data_t>::node_size_dist,
    "A function to get the node size distribution of the hypergraph")
    .def("edge_incidence", &NWHypergraph<Index_t, Data_t>::edge_incidence,
    "A function to get the incident nodes of an edge in the hypergraph", py::arg("edge"))
    .def("node_incidence", &NWHypergraph<Index_t, Data_t>::node_incidence,
    "A function to get the incident edges of a node in the hypergraph", py::arg("node"))
    .def("degree", py::overload_cast<Index_t, size_t>(&NWHypergraph<Index_t, Data_t>::degree),
    "A function to get the degree of a node in the hypergraph", py::arg("node"), py::arg("size") = 1)
    .def("size", &NWHypergraph<Index_t, Data_t>::size,
    "A function to get the number of nodes that belong to an edge in the hypergraph", py::arg("edge"))
    .def("dim", &NWHypergraph<Index_t, Data_t>::dim,
    "A function to get the number of nodes that belong to an edge minus 1 in the hypergraph", py::arg("edge"))
    .def("number_of_nodes", &NWHypergraph<Index_t, Data_t>::number_of_nodes,
    "A function to get the number of nodes in the hypergraph")
    .def("order", &NWHypergraph<Index_t, Data_t>::order,
    "A function to get the number of nodes in the hypergraph")  
    .def("number_of_edges", &NWHypergraph<Index_t, Data_t>::number_of_edges,
    "A function to get the number of edges in the hypergraph")

    //return singletons 
    // a singleton is an edge of size 1 with a node of degree 1
    .def("singletons", &NWHypergraph<Index_t, Data_t>::singletons,
    "A function to get the list of singletons in the hypergraph")
    
    // return toplexes
    // a toplex is an edge where its nodes are not included in another edge
    .def("toplexes", &NWHypergraph<Index_t, Data_t>::toplexes,
    "A function to get the list of edges of toplexes in the hypergraph")

    //create slinegraph from nwhypergraph
    .def("s_linegraph", &NWHypergraph<Index_t, Data_t>::s_linegraph, 
    "A function which converts a hypergraph to its s line graph; if edges is true, then it is an edge linegraph",
    py::arg("s") = 1, py::arg("edges") = true)
    //s_connected_component
    .def("s_connected_components", py::overload_cast<Slinegraph<Index_t, Data_t> &, bool>(&NWHypergraph<Index_t, Data_t>::s_connected_components),
     "A function which finds the connected components for its s line graph",
    py::arg("linegraph"), py::arg("return_singleton") = false)
    //s_distance
    .def("distance", py::overload_cast<Slinegraph<Index_t, Data_t> &, Index_t, Index_t>(&NWHypergraph<Index_t, Data_t>::distance),
     "A function which computes the distance from src to dest in its s line graph",
    py::arg("linegraph"), py::arg("src"), py::arg("dest"))
    //s_neighbor
    .def("neighbors", py::overload_cast<Slinegraph<Index_t, Data_t> &, Index_t>(&NWHypergraph<Index_t, Data_t>::neighbors),
     "A function which finds the neighbors for vertex v of its s line graph",
    py::arg("linegraph"), py::arg("v"))
    .def("s_degree", py::overload_cast<Slinegraph<Index_t, Data_t> &, Index_t>(&NWHypergraph<Index_t, Data_t>::s_degree),
     "A function which finds the degree for vertex v of its s line graph",
    py::arg("linegraph"), py::arg("v"))
    ;

    //define Slinegraph python object
    py::class_<Slinegraph<Index_t, Data_t>> slinegraph_class(m, "Slinegraph");
    slinegraph_class
    .def_readonly("row", &Slinegraph<Index_t, Data_t>::row_)
    .def_readonly("col", &Slinegraph<Index_t, Data_t>::col_)
    .def_readonly("data", &Slinegraph<Index_t, Data_t>::data_)
    .def(py::init<>([](NWHypergraph<Index_t, Data_t>& g, int s, bool edges) {
        return new Slinegraph<Index_t, Data_t>(g, s, edges);
    }), "Init function", py::arg("g"), py::arg("s") = 1, py::arg("edges") = true)
    .def_readonly("s", &Slinegraph<Index_t, Data_t>::s_)
    .def("s_connected_components", &Slinegraph<Index_t, Data_t>::s_connected_components,
     "A function which finds the connected components for its s line graph", py::arg("return_singleton") = false)
    .def("is_s_connected", &Slinegraph<Index_t, Data_t>::is_s_connected,
    "A function which tests whether its s line graph is connected or not")
    .def("s_distance", &Slinegraph<Index_t, Data_t>::s_distance,
    "A function to compute the distance from src to dest", py::arg("src"), py::arg("dest"))
    .def("s_diameter", &Slinegraph<Index_t, Data_t>::s_diameter,
    "A function to compute the diameter of the s line graph")
    .def("s_path", &Slinegraph<Index_t, Data_t>::s_path,
    "A function which finds one shortest path from src to dest (could have multiple, only find one)", py::arg("src"), py::arg("dest")) 
    .def("s_betweenness_centrality", &Slinegraph<Index_t, Data_t>::s_betweenness_centrality,
    "A function which computes the betweenness centrality from src to all the other vertices", 
    py::arg("v"), py::arg("normalized") = true)
    .def("s_closeness_centrality", &Slinegraph<Index_t, Data_t>::s_closeness_centrality,
    "A function which computes the closeness centrality of v", py::arg("v"))
    .def("s_harmonic_closeness_centrality", &Slinegraph<Index_t, Data_t>::s_harmonic_closeness_centrality,
    "A function which computes the harmonic closeness centrality of v", py::arg("v"))
    .def("s_eccentricity", &Slinegraph<Index_t, Data_t>::s_eccentricity,
    "A function which computes the eccentrality of a vertex", py::arg("v"))        
    .def("s_neighbors", &Slinegraph<Index_t, Data_t>::s_neighbors,
    "A function to get neighbors of a vertex", py::arg("v"))
    .def("s_degree", &Slinegraph<Index_t, Data_t>::s_degree,
    "A function to get the degree of a vertex in the slinegraph", py::arg("v"));

    //register version information in a module as below
    py::object version = py::cast("0.0.3");
    m.attr("_version") = version;

    //define function, its argument list, and with default argument for s
    m.def("convert_to_s_overlap", &convert_to_s_overlap<Index_t, Data_t>, "A function which converts a hypergraph to its s line graph",
    py::arg("x"), py::arg("y"), py::arg("data"), py::arg("s") = 1);
}