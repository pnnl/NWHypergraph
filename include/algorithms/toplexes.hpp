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
    //O(m+n), where m is the size of lhs, n is the size of rhs
    auto issubset = []<class A>(A &&lhs, A &&rhs) {
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
    };
    std::vector<vertex_id_t> tops;
    for (vertex_id_t e = 0; e < M; ++e)
    {
        bool flag = true;
        std::vector<vertex_id_t> old_tops(tops);
        for (size_t i = 0, end = old_tops.size(); i < end; ++i)
        {
            vertex_id_t top = old_tops[i];
            //TODO is this necessary?
            if (e == top)
                continue;
            if (issubset(edges[top], edges[e]))
                tops.erase(tops.begin() + i);
            else if (issubset(edges[e], edges[top]))
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

template<typename GraphE>
auto toplexes_serial_v1(GraphE& edges) {
    nw::util::life_timer _(__func__);
    size_t M = edges.size();
    std::vector<vertex_id_t> tops;
    //std::map<Index_t, bool> tops;
    for (vertex_id_t e = 0; e < M; ++e)
    {
        bool flag = true;
        std::vector<vertex_id_t> old_tops(tops);
        for (size_t i = 0, end = old_tops.size(); i < end; ++i)
        {
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

}//namespace hypergraph
}//namespace nw