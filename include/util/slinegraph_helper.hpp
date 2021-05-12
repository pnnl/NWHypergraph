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
#include <vector>
#include <unordered_map>

#include <util/types.hpp>

using namespace nw::graph;
namespace nw {
namespace hypergraph {

/*
* Squeeze the edge lists such that the ids are consecutive in the new edge list
*/
template<directedness edge_directedness = nw::graph::undirected, typename... T>
auto squeeze_edgelist(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T...>>> &two_graphs) {
    nw::util::life_timer _(__func__);
    nw::graph::edge_list<edge_directedness, T...> result(0);
    //result.open_for_push_back();
    //do this in serial
    vertex_id_t index = 0;
    std::unordered_map<vertex_id_t, vertex_id_t> relabel_map;
    int num_bins = two_graphs.size();
    std::for_each(tbb::counting_iterator<int>(0),
                  tbb::counting_iterator<int>(num_bins), [&](auto i) {
                    std::for_each(
                        two_graphs[i].begin(), two_graphs[i].end(),
                        [&](auto &&elt) {
                          std::apply([&](vertex_id_t x, vertex_id_t y, T... w) {
                            if (relabel_map.end() == relabel_map.find(x)) {
                              // if x has not been relabeled
                              relabel_map[x] = index;
                              ++index;
                            }
                            auto newx = relabel_map[x];
                            if (relabel_map.end() == relabel_map.find(y)) {
                              // if y has not been relabeled
                              relabel_map[y] = index;
                              ++index;
                            }
                            auto newy = relabel_map[y];
                            // std::cout << x << " " << y << " into " << newx <<
                            // " " << newy << std::endl;
                            result.push_back(newx, newy, w...);
                          }, elt);
                        });
                  });
    result.close_for_push_back(false);

    return result;
}

/*
* Do NOT squeeze the edge lists and combine them in the new edge list.
*/
template<directedness edge_directedness = nw::graph::undirected, class... T>
auto create_edgelist_without_squeeze(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T...>>> &two_graphs) {
    nw::util::life_timer _(__func__);
    nw::graph::edge_list<edge_directedness, T...> result(0);
    //result.open_for_push_back();
    //do this in serial
    int num_bins = two_graphs.size();
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto &&e) {
        result.push_back(e);
      });
    });
    result.close_for_push_back(false);

    return result;
}

/*
* Squeeze the edge lists such that the ids are consecutive in the new edge list.
*/
template<directedness edge_directedness = nw::graph::undirected, class... T>
auto create_edgelist_with_squeeze(std::vector<std::vector<std::tuple<vertex_id_t, vertex_id_t, T...>>> &two_graphs) {
    nw::util::life_timer _(__func__);
    nw::graph::edge_list<edge_directedness, T...> result(0);
    //result.open_for_push_back();
    //do this in serial
    vertex_id_t index = 0;
    std::unordered_map<vertex_id_t, vertex_id_t> relabel_map;
    int num_bins = two_graphs.size();
    std::for_each(tbb::counting_iterator<int>(0), tbb::counting_iterator<int>(num_bins), [&](auto i) {
      std::for_each(two_graphs[i].begin(), two_graphs[i].end(), [&](auto &&elt) {
        std::apply([&](vertex_id_t x, vertex_id_t y, T... w) {
          if (relabel_map.end() == relabel_map.find(x)) {
          //if x has not been relabeled
          relabel_map[x] = index;
          ++index;
          }
          auto newx = relabel_map[x];
          if (relabel_map.end() == relabel_map.find(y)) {
          //if y has not been relabeled
            relabel_map[y] = index;
            ++index;
          }
          auto newy = relabel_map[y];
          result.push_back(newx, newy, w...);
        }, elt);
      });
    });
    result.close_for_push_back(false);

    return result;
}

}//namespace hypergraph
}//namespace nw