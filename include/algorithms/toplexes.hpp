//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2021
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#pragma once
#include <util/timer.hpp>
#include <util/intersection_size.hpp>
#include <map>
#include <vector>
#include <unordered_set>


namespace nw {
namespace hypergraph {
/*
* TODO why not use std::set_intersection
* or std::include
*/
template<typename A>
bool is_subset(A&& lhs, A&& rhs) {
    std::map<vertex_id_t, size_t> frequency;
    std::for_each(lhs.begin(), lhs.end(), [&](auto &&x) {
        auto v = std::get<0>(x);
        ++frequency[v];
    });
    bool res = true;
    std::for_each(rhs.begin(), rhs.end(), [&](auto &&x) {
        auto v = std::get<0>(x);
        if (0 < frequency[v])
            --frequency[v];
        else
        {
            res = false;
            return;
        }
    });
    return res;
}
template<typename GraphE>
auto toplexes_serial_v0(GraphE& edges) {
    nw::util::life_timer _(__func__);
    size_t M = edges.size();
    //create an empty toplex set and an empty old toplex set
    std::vector<vertex_id_t> tops;
    for (vertex_id_t e = 0; e < M; ++e) {
            //for each edge, assume it is a toplex by default
        bool flag = true;
            //make a copy of the toplex set as old_tops
        auto old_tops(tops);
            //tops.clear();
            //clear the toplex set
            //TODO Could be parallized:
            //1) if the flag can be set atomically
            //2) once the flag is set, then every thread will exist
            //3) tops.erase has to be within a critical section
        for (size_t i = 0, end = old_tops.size(); i < end; ++i) {
            vertex_id_t top = old_tops[i];
            if (e == top)
                continue; 
            if (is_subset(edges[top], edges[e])) {
                //if e is a subset of top, then e is not a toplex
                flag = false;
                break;
            } else if (is_subset(edges[e], edges[top]))
                tops.erase(tops.begin() + i);
        }//for old_tops
        if (flag)
            tops.push_back(e);
    }//for each e
    return tops;
}
template<typename GraphE>
auto toplexes_serial_v1(GraphE& edges) {
    nw::util::life_timer _(__func__);
    size_t M = edges.size();

    std::unordered_set<vertex_id_t> tops;
    for (vertex_id_t e = 0; e < M; ++e) {
        bool flag = true;
        auto old_tops(tops);
        for (auto& top : old_tops) {
            //It is necessary to differeniate e from top
            if (e == top)
                continue;
            if (is_subset(edges[top], edges[e]))
                tops.erase(top);
            else if (is_subset(edges[e], edges[top])) {
                flag = false;
                break;
            }
        }
        if (flag)
            tops.insert(e);
    }
    return tops;
}

template<typename GraphE>
auto toplexes_serial_v2(GraphE& edges) {
    nw::util::life_timer _(__func__);
    size_t M = edges.size();
    std::vector<vertex_id_t> tops;
    //std::map<Index_t, bool> tops;
    for (vertex_id_t e = 0; e < M; ++e) {
        bool flag = true;
        std::vector<vertex_id_t> old_tops(tops);
        for (size_t i = 0, end = old_tops.size(); i < end; ++i) {
            vertex_id_t top = old_tops[i];
            //TODO is this necessary?
            if (e == top)
                continue;
            auto lhs = edges[e].size();
            auto rhs = edges[top].size();
            auto s = nw::graph::intersection_size(edges[e], edges[top]);
            if (s == rhs)
                tops.erase(tops.begin() + i);
            else if (s == lhs)
            {
                flag = false;
                break;
            }
        }
        
        if (flag)
            tops.push_back(e);
    }

    return tops;
}

/*
* TODO
*/
template<typename GraphE>
auto toplexes_serial_v4(GraphE& edges) {
    nw::util::life_timer _(__func__);
    size_t M = edges.size();
    std::vector<vertex_id_t> tops(M);
    for (vertex_id_t e = 0; e < M; ++e) {
        bool flag = true;
        for (vertex_id_t another_e = 0; another_e < M; ++another_e) {
            if (is_subset(edges[e], edges[another_e]))
                ;
            else if (is_subset(edges[another_e], edges[e])) {
                flag = false;
                break;
            }
        }
        if (flag)
            tops.push_back(e);
    }
    return tops;
}

}//namespace hypergraph
}//namespace nw