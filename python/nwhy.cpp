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

    //register version information in a module as below
    py::object version = py::cast("0.0.2");
    m.attr("_version") = version;

    //define function, its argument list, and with default argument for s
    m.def("convert_to_s_overlap", &convert_to_s_overlap<Index_t, Data_t>, "A function which converts a hypergraph to its s line graph",
    py::arg("x"), py::arg("y"), py::arg("data"), py::arg("s") = 1);
}