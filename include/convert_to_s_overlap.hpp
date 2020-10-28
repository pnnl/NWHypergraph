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
#include "s_overlap.hpp"
using namespace nw::graph;


namespace nw {
namespace hypergraph {

std::tuple<std::vector<int>,std::vector<int>, 
std::vector<int>> convert_to_s_overlap(std::vector<int> x,std::vector<int> y, 
std::vector<int> data = std::vector<int>(0), size_t s = 1) {
    //TODO store a_coo
    size_t n_x = x.size();
    size_t n_y = y.size();
    size_t n_data = data.size();
    if (0 == n_data)
        std::cout << "unweighted hypergraph" << std::endl;
    return std::tuple(x, y, data);
}

}//namespace hypergraph
}//namespace nw