
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

#include <optional>
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <edge_list.hpp>
#include <util/intersection_size.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include "algorithms/slinegraph_map.hpp"
#include "algorithms/slinegraph_efficient.hpp"
#include "algorithms/toplexes.hpp"
#include <algorithms/connected_components.hpp>
#include <algorithms/delta_stepping.hpp>
#include <algorithms/betweenness_centrality.hpp>

using namespace nw::graph;
namespace py = pybind11;

namespace nw {
namespace hypergraph {

template<class Index_t, typename... Attributes> class NWHypergraph;

template<class Index_t, typename... Attributes> 
class Slinegraph {
private:
    bool edges_;
    Index_t offset_; //offset between min id to 0
    //s-linegraph generated from original hypergraph
    nw::graph::adjacency<0, Attributes...> g_;
    nw::graph::adjacency<1, Attributes...> g_t_;
public:
    py::array_t<Index_t, py::array::c_style> row_;
    py::array_t<Index_t, py::array::c_style> col_;
    py::array_t<Attributes..., py::array::c_style> data_;
    int s_;
private:
    /*
    * Note N may or may not equal to linegraph.size()
    * */
    void populate_adjacency(nw::graph::edge_list<nw::graph::undirected, Attributes...>& linegraph,
    std::size_t N) {
        size_t nedges = linegraph.size();
        std::cout << "linegraph size: " << nedges << std::endl;
        //instantiate empty adjacency
        if (0 == nedges) {
            //create an adjacency where each list is empty
            g_ = nw::graph::adjacency<0, Attributes...>(N);
            g_t_ = nw::graph::adjacency<1, Attributes...>(N);
        }
        //fill adjacency with edges
        else {
            g_ = nw::graph::adjacency<0, Attributes...>(N, linegraph);
            g_t_ = nw::graph::adjacency<1, Attributes...>(N, linegraph);
        }
    }
    void populate_py_array(nw::graph::edge_list<nw::graph::undirected, Attributes...>& linegraph) {
        pybind11::ssize_t n = linegraph.size();
        row_ = py::array_t<Index_t, py::array::c_style>(n);
        col_ = py::array_t<Index_t, py::array::c_style>(n);
        data_ = py::array_t<Attributes..., py::array::c_style>(n);
        auto rrow_ = row_.template mutable_unchecked<1>();
        auto rcol_ = col_.template mutable_unchecked<1>(); 
        auto rdata_ = data_.template mutable_unchecked<1>();
        size_t i = 0;
        for (auto &&[u, v, weight] : linegraph) {
            rrow_(i) = u;
            rcol_(i) = v;
            rdata_(i) = weight;
            ++i;
        }
    }
public:
    //constructor
    Slinegraph(NWHypergraph<Index_t, Attributes...>&g, int s = 1, bool edges = true) : 
    edges_(edges), 
    offset_(0),
    s_(s) {
        //extract linegraph from its neighbor counts
        if (edges_) {
            nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph = g.populate_edge_linegraph(s);
            populate_adjacency(linegraph, g.edges_.size());
            populate_py_array(linegraph);
        }
        else{
            nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph = 
            g.populate_node_linegraph(s);         
            populate_adjacency(linegraph, g.nodes_.size());
            populate_py_array(linegraph);
        }
    }
    Slinegraph(
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y,
    py::array_t<Attributes..., py::array::c_style | py::array::forcecast> &data,
    int s = 1, bool edges = true) : 
    edges_(edges), 
    offset_(0),
    row_(x), 
    col_(y), 
    data_(data), 
    s_(s) {   
        nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph(0);
        linegraph.open_for_push_back();
        //sanitize check
        auto rx = x.template mutable_unchecked<1>();
        auto ry = y.template mutable_unchecked<1>();
        auto rdata = data.template mutable_unchecked<1>();
        //rx(0) = 1;
        size_t n_x = x.shape(0);
        size_t n_y = y.shape(0);
        size_t n_data = data.shape(0);
        Index_t m = 0;
        if (0 == n_data) {
            //when there is no weight passed in, but request a weighted hypergraph
            //we fake a weight with value 1
            for (size_t i = 0; i < n_x; ++i) {
                linegraph.push_back({rx(i), ry(i), 1});
                m = std::max(m, rx(i));
                m = std::max(m, ry(i));
            }
        }
        else {
            for (size_t i = 0; i < n_x; ++i) {
                linegraph.push_back({rx(i), ry(i), rdata(i)});
                m = std::max(m, rx(i));
                m = std::max(m, ry(i));
            }
        }
        //offset_ = linegraph.close_for_push_back_with_shift();
        linegraph.close_for_push_back(false);
        g_ = nw::graph::adjacency<0, Attributes...>(linegraph);
        g_t_ = nw::graph::adjacency<1, Attributes...>(linegraph);
    }
    /*
    * Get all the singletons in the s line graph.
    * */ 
    py::list get_singletons() {
        auto E = nw::graph::ccv1(g_);
        // This returns the subgraph of each component.
        std::map<Index_t, py::set> comps;
        for (size_t i = 0, e = E.size(); i < e; ++i) {
            Index_t label = E[i];
            Index_t original_id = i + offset_;
            if (comps.find(label) != comps.end())
                comps[label].add(original_id);
            else {
                py::set s;
                s.add(original_id);
                comps[label] = s;
            }
        }
        if (0 != offset_) {
            //if the ids are shifted by a number, 
            //then singletons from 0 to n - 1 must be added to list
            py::list l;
            for (Index_t i = 0; i < offset_; ++i) {
                py::set s;
                s.add(i);
                l.append(s);
            }
            for (auto it = comps.begin(); it != comps.end(); ++it) {
                if (1 == it->second.size())
                    l.append(it->second);  
            }
            return l;
        }
        else {
            //if ids are not shifted, do nothing
            py::list l;
            for (auto it = comps.begin(); it != comps.end(); ++it) {
                if (1 == it->second.size())
                    l.append(it->second);  
            }
            return l;
        }
    }
    /*
    * Return a list of sets which
    * each set is a component.
    */
    py::list s_connected_components() {
        auto E = nw::graph::ccv1(g_);
        
        // This returns the subgraph of each component.
        std::map<Index_t, py::set> comps;
        for (size_t i = 0, e = E.size(); i < e; ++i) {
            Index_t label = E[i];
            Index_t original_id = i + offset_;
            if (comps.find(label) != comps.end())
                comps[label].add(original_id);
            else {
                py::set s;
                s.add(original_id);
                comps[label] = s;
            }
        }
        py::list l;
        for (auto it = comps.begin(); it != comps.end(); ++it) {
            //auto s = it->second;
            //py::ssize_t n = it->second.size();
            if (1 < it->second.size()) 
                l.append(it->second); //only return non-singletons
        }
        return l;
    }
    /*
    * Check whether slinegraph is s_connected
    * Check the neighborhood size of each vertex
    * As long as there is one non-empty neighborhood, then slinegraph is s_connected
    */
    bool is_s_connected() {
        //if every vertex in g_, if one has neighbor, then g_ is s_connected
        auto res = nw::graph::parallel_for(tbb::blocked_range<Index_t>(0, g_.size()), [&](auto& i) {            
            return g_[i].size();
        }, std::plus{}, 0ul);
        return (0 < res);
    }
    /*
    * Compute the distance from src to dest
    * return -1 if unreachable
    */
    Index_t s_distance(Index_t src, Index_t dest) {
        using distance_t = std::uint64_t;
        size_t delta = 1;
        auto dist = nw::graph::delta_stepping_v12<distance_t>(g_t_, src, delta);
        Index_t dist_of_dest = dist[dest];
        if (dist_of_dest == std::numeric_limits<Index_t>::max() - 1)
            return -1;
        return dist_of_dest;
    }
    /*
    * Compute the distance from src to dest
    * return 0 if every vertex is a singleton
    */
    auto s_diameter() {
        using distance_t = std::uint64_t;
        std::size_t delta = 1;
        std::size_t n = g_.size();
        std::vector<distance_t> bc(n);
        tbb::parallel_for(tbb::blocked_range<Index_t>(0, n), [&](tbb::blocked_range<Index_t>& r) {
            for (auto i = r.begin(), e = r.end(); i != e; ++i) {
                auto dist = nw::graph::delta_stepping_v12<distance_t>(g_t_, i, delta);
                auto tmp = std::max_element(dist.begin(), dist.end());
                bc[i] = (*tmp).load(std::memory_order_relaxed);
            }
        }, tbb::auto_partitioner());
        distance_t result = *std::max_element(bc.begin(), bc.end());
        if (result == std::numeric_limits<Index_t>::max())
            return static_cast<distance_t>(0);
        return result;
    }
    /*
    * Compute the shortest path from src to dest.
    * Note an empty list will be returned if unreachable from src to dest.
    * Current impl is in seral. More improvement needs be done (parallelize it).
    * 
    * The main idea of this implementation is:
    * Find the shortest path from src to a node u using BFS
    * Find the shortest path from dest to a node v using BFS
    * when u and v are identical, the shortest path from src to dest can be formed
    * as [src, u/v, dest]
    * To identify u and v are the same, use two bitmaps check whether either has been visited
    * by src or by dest respectively
    */
    py::list s_path(Index_t src, Index_t dest) {
        //1. create two queue to memorize the paths from src and the paths from dest
        std::queue<std::tuple<Index_t, std::vector<Index_t>>> srcQ, destQ; //node, path from node
        //2. add src and dest to its queue respectively
        srcQ.push(std::make_tuple(src, std::vector<Index_t>{src}));
        destQ.push(std::make_tuple(dest, std::vector<Index_t>{dest}));

        std::size_t n = g_.size();
        //3. create a bitmap to memorize which node has been visited
        nw::graph::AtomicBitVector<Index_t>   srcVisited(n), destVisited(n);
        srcVisited.atomic_set(src);
        destVisited.atomic_set(dest);

        py::list l;
        //4. continue if both srcQ and destQ are not empty
        while((!srcQ.empty() && !destQ.empty())) {
            auto&& [u, uPath] = srcQ.front();
            srcQ.pop();
            if (0 == destVisited.atomic_get(u)) {
                //if u has been visited by dest
                //bingo we found a path
                for(auto& v : uPath) {
                    l.append(v);
                }
                while (!destQ.empty()) {
                    auto& [v, vPath] = destQ.front();
                    destQ.pop();
                    if (v == u) {
                        //we found the path
                        //by adding reversely from the end of the path to v
                        for (size_t i = vPath.size(); i >= 0; --i) {
                            l.append(vPath[i]);
                        }
                        break;
                    }
                }
                break;
            }
            if (0 == srcVisited.atomic_get(u)) {
                //if cur_node has not been visited by src, 
                //we keep extend the path by visited all the neighbors of cur_node
                std::for_each(g_[u].begin(), g_[u].end(), [&](auto&& i) {
                    auto v = std::get<0>(i);
                    if (0 == srcVisited.atomic_get(v)) {
                        //if this neighbor has not been visited
                        if (0 == srcVisited.atomic_set(v)) {
                            uPath.push_back(v);
                            srcQ.push({v, uPath});
                        }
                    }
                    else {
                        //ignore the visited neighbors
                    }
                });
            }
            auto&& [cur_node, cur_path] = destQ.front();
            destQ.pop();
            if (0 == srcVisited.atomic_get(cur_node)) {
                //if cur_node has been visited by src
                //bingo we found a path
                while (!srcQ.empty()) {
                    auto& [node, path] = srcQ.front();
                    srcQ.pop();
                    if (node == cur_node) {
                        //we found the path from src to cur_node
                        for (size_t i = 0, e = path.size(); i < e; ++i) {
                            l.append(path[i]);
                        }
                        break;
                    }
                }
                for(auto& v : cur_path) {
                    l.append(v);
                }
                break;
            }
            if (0 == destVisited.get(cur_node)) {
                //if cur_node has not been visited by dest, 
                //we keep extend the path by visited all the neighbors of cur_node
                std::for_each(g_[cur_node].begin(), g_[cur_node].end(), [&](auto&& i) {
                    auto v = std::get<0>(i);
                    if (0 == destVisited.atomic_get(v)) {
                        //if this neighbor has not been visited
                        if (0 == destVisited.atomic_set(v)) {
                            cur_path.push_back(v);
                            destQ.push({v, cur_path});
                        }
                    }
                    else {
                        //ignore the visited neighbors
                    }
                });

            }

        }//while
        
        return l;
    }
    /*
    * Compute the betweenness centrality of vertex v to any other vertices
    * TODO compute cc then compute bc on each cc
    */
    py::list s_betweenness_centrality(bool normalized = true) {
        using score_t=float;
        using accum_t=double;
        std::vector<score_t> bc = nw::graph::betweenness_brandes<decltype(g_), score_t, accum_t>(g_, false);
        std::size_t n = bc.size();
        float scale = 1.0;
        if (normalized) {
            if (2 < n)
                scale /= ((n - 1) * (n - 2));
        }
        py::list l = py::list(n);
        if (!normalized) {
            for (std::size_t i = 0; i < n; ++i) {
                l[i] = bc[i];
            }
        }
        else {
            for (std::size_t i = 0; i < n; ++i) {
                l[i] = bc[i] * scale;
            }
        }
        return l;
    }
    /*
    * Closeness centrality of a node `v` is the reciprocal of the
    * average shortest path distance to `v` over all `n-1` reachable nodes.
    */
    py::list s_closeness_centrality(std::optional<Index_t> v = {}) {
        using distance_t = std::uint64_t;
        std::size_t delta = 1;
        Index_t n = g_.size();
        if (v.has_value()) {
            py::list l;
            //validate input
            if (v.value() >= n)
                return l;
            auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v.value(), delta);
            distance_t sum = 0ul;
            std::size_t ncomp = 0;
            for (auto& d : dist) {
                distance_t tmp = d.load(std::memory_order_relaxed);
                if (tmp != std::numeric_limits<distance_t>::max() && 0 != tmp) {
                    ++ncomp;
                    sum += tmp;
                }
            } 
            l.append<float>(1.0 * ncomp / sum);
            return l;
        }
        else {
            py::list l = py::list(n);
            //for each vertex v in the graph
            for (Index_t v = 0; v < n; ++v) {
                //get the distances from v to any other vertices
                auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v, delta);
                distance_t sum = 0ul;
                std::size_t ncomp = 0;
                for (auto &d : dist) {
                    distance_t tmp = d.load(std::memory_order_relaxed);
                    if (tmp != std::numeric_limits<distance_t>::max() && 0 != tmp) {
                        ++ncomp;
                        sum += tmp;
                    }
                }
                float res = 1.0 * ncomp / sum;
                l[v] = res;
            }
            return l;
        }
    }
    /*
    * 
    *  
    */
    py::list s_harmonic_closeness_centrality(std::optional<Index_t> v = {}) {
        using distance_t = std::uint64_t;
        std::size_t delta = 1;
        Index_t n = g_.size();
        if (v.has_value()) {
            py::list l;
            //validate input
            if (v >= n)
                return l;
            auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v.value(), delta);
            float sum = 0;
            for (auto& d : dist) {
                distance_t tmp = d.load(std::memory_order_relaxed);
                if (tmp != std::numeric_limits<distance_t>::max() && 0 != tmp)
                    sum +=  1.0 / tmp;
            }  
            l.append(1.0 * sum);
            return l;
        }
        else {
            py::list l = py::list(n);
            //for each vertex v in the graph
            for (Index_t v = 0; v < n; ++v) {
                //get the distances from v to any other vertices
                auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v, delta);
                float sum = 0;
                for (auto &d : dist) {
                    distance_t tmp = d.load(std::memory_order_relaxed);
                    if (tmp != std::numeric_limits<distance_t>::max() && 0 != tmp)
                        sum += 1.0 / tmp;
                }
                l[v] = 1.0 * sum;
            }
            return l;
        }
    }
    /*
    * Compute the eccentricity of vertex v
    * return 0 if unreachable
    */
    py::list s_eccentricity(std::optional<Index_t> v = {}) {
        using distance_t = std::uint64_t;
        std::size_t delta = 1;
        Index_t n = g_.size();
        if (v.has_value()) {
            py::list l;
            //validate input
            if (v >= n)
                return l;

            auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v.value(), delta);
            auto tmp = std::max_element(dist.begin(), dist.end());
            distance_t result = (*tmp).load(std::memory_order_relaxed);
            if (result == std::numeric_limits<distance_t>::max())
                l.append(std::numeric_limits<distance_t>::infinity());
            else
                l.append(result);
            return l;
        }
        else {
            py::list l = py::list(n);
            //for each vertex v in the graph
            for (Index_t v = 0; v < n; ++v) {
                //get the distances from v to any other vertices
                auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v, delta);
                auto tmp = std::max_element(dist.begin(), dist.end());
                distance_t result = (*tmp).load(std::memory_order_relaxed);
                if (result == std::numeric_limits<distance_t>::max())
                    l[v] = std::numeric_limits<distance_t>::infinity();
                else
                    l[v] = result;
            }     
            return l;
        }
    }
    /*
    * Get the neighbors of a vertex in the slinegraph
    */
    py::list s_neighbors(Index_t v) {
        py::list l;
        for (auto u = g_[v].begin(); u != g_[v].end(); ++u) {
            l.append(std::get<0>(*u));
        }
        return l;
    }
    py::ssize_t s_neighborhood_size(Index_t v) { return g_[v].size(); }
    py::ssize_t s_degree(Index_t v) { return g_[v].size(); }
    int getS() const { return s_; }
    bool isEdgeOverlap() const { return edges_; }
}; //slinegraph class


}//namespace hypergraph
}//namespace nw