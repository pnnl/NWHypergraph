//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2021-2021
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#pragma once
#include <nwgraph/edge_list.hpp>
#include <execution>
#include <vector>
#include <random>
#include <math.h>

namespace nw {
namespace hypergraph {

/*
* Return a random undirected bipartite graph (in edge_list) from two given degree sequences.
* The bipartite graph is composed of two partitions. 
* By default, contiguous_id_space is false, set A has nodes 0 to
* (deg_seqa.size() - 1) and set B has nodes 0 to (deg_seqb.size() - 1).
* When contiguous_id_space is true, set A has nodes 0 to (deg_seqa.size() - 1)
* and set B has nodes deg_seqa.size() to (deg_seqa.size() + deg_seqb.size() - 1).
* Nodes from set A are connected to nodes in set B by choosing
* randomly from the possible free stubs, one in A and one in B.
* Note directed graph is not supported in configuration model.
* 
*/
template<class T>
nw::graph::bi_edge_list<nw::graph::directedness::directed> 
configuration_model(std::vector<T>& deg_seqa, std::vector<T>& deg_seqb, bool contiguous_id_space = false) {
    //validate degree sequences such that their summation are equivalent
    nw::graph::bi_edge_list<nw::graph::directedness::directed> el;
    T suma = std::reduce(std::execution::par_unseq, deg_seqa.cbegin(), deg_seqa.cend());
    T sumb = std::reduce(std::execution::par_unseq, deg_seqb.cbegin(), deg_seqb.cend());
    if (suma != sumb) {
        std::cout << "Invalid degree sequences, sum(deg_seqa)!=sum(deg_seqb)" << std::endl;
        return el;
    }
    auto populate_stubs = [&](auto& degree_sequence) -> std::vector<T> {
        std::vector<T> stubs(suma);
        std::size_t cursor = 0;
        for (std::size_t v = 0, e = degree_sequence.size(); v < e; ++v) {
            auto degree_v = degree_sequence[v];
            auto begin = stubs.begin() + cursor;
            cursor += degree_v;
            auto end = begin + cursor;
            std::fill(std::execution::par_unseq, begin, end, v);
        }
        return stubs;
    };
    auto astubs = populate_stubs(deg_seqa);
    auto bstubs = populate_stubs(deg_seqb);

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(astubs.begin(), astubs.end(), g);
    std::shuffle(bstubs.begin(), bstubs.end(), g);

    std::copy(astubs.begin(), astubs.end(), std::ostream_iterator<int>(std::cout, " "));


    el.open_for_push_back();
    for(std::size_t i = 0, e = astubs.size(); i < e; ++i) {
        el.push_back(astubs[i], bstubs[i]);
    }
    el.close_for_push_back();
    return el;
}

}//namespace hypergraph
}//namespace nw
