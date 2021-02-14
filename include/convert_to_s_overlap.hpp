//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2020-2021
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
#include <algorithms/connected_components.hpp>
#include <algorithms/delta_stepping.hpp>
#include <algorithms/betweenness_centrality.hpp>
using namespace nw::graph;
namespace py = pybind11;

namespace nw {
namespace hypergraph {

//forward-declaration
template<class Index_t, typename... Attributes> class Slinegraph;

template<class Index_t, typename... Attributes>
class NWHypergraph {
private:
    //hypergraph in bi-adjacency
    nw::graph::adjacency<0, Attributes...> edges_;
    Index_t max_edge_;
    nw::graph::adjacency<1, Attributes...> nodes_;
    Index_t max_node_;

    //for each hyperedge, map stores its hyperedge neighbor list
    //which the key is the adjacent hyperedge id
    //the value is the number of overlap hypernodes
    std::vector<std::map<size_t, size_t>> edge_neighbor_count_;
    //for each hypernode, map stores its hypernode neighbor list
    //which the key is the adjacent hypernode
    //the value is the number of overlap hyperedges
    std::vector<std::map<size_t, size_t>> node_neighbor_count_;   
public:
    py::array_t<Index_t, py::array::c_style> row_;
    py::array_t<Index_t, py::array::c_style> col_;
    py::array_t<Attributes..., py::array::c_style> data_;
    friend class Slinegraph<Index_t, Attributes...>;
public:
    void populate_neighbor_count(bool edges = true) {
        if (edges){
            edge_neighbor_count_ = to_two_graph_count_neighbors_cyclic(edges_, nodes_);
        }
        else {
            node_neighbor_count_ = to_two_graph_count_neighbors_cyclic(nodes_, edges_);
        }
    }
public:
    //constructor
    NWHypergraph() {}
    NWHypergraph(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y) : row_(x), col_(y) {
        nw::graph::edge_list<nw::graph::directed, Attributes...> g(0);
        //sanitize check
        auto rx = x.template mutable_unchecked<1>();
        auto ry = y.template mutable_unchecked<1>();
        //rx(0) = 1;
        size_t n_x = x.shape(0);
        size_t n_y = y.shape(0);
        for (size_t i = 0; i < n_x; ++i) {
            g.push_back({rx(i), ry(i), 0});
        }
        edges_ = nw::graph::adjacency<0, Attributes...>(g);
        max_edge_ = edges_.size();
        nodes_ = nw::graph::adjacency<1, Attributes...>(g);
        max_node_ = nodes_.size();
    }
    NWHypergraph(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
    py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y,
    py::array_t<Attributes..., py::array::c_style | py::array::forcecast> &data,
    bool collapse = false) : row_(x), col_(y), data_(data) {
        nw::graph::edge_list<nw::graph::directed, Attributes...> g(0);
        //sanitize check
        auto rx = x.template mutable_unchecked<1>();
        auto ry = y.template mutable_unchecked<1>();
        auto rdata = data.template mutable_unchecked<1>();
        //rx(0) = 1;
        size_t n_x = x.shape(0);
        size_t n_y = y.shape(0);
        size_t n_data = data.shape(0);
        if (0 == n_data) {
            //when there is no weight passed in, but request a weighted hypergraph
            //we fake a weight
            for (size_t i = 0; i < n_x; ++i) 
                g.push_back({rx(i), ry(i), 0});
        }
        else {
            for (size_t i = 0; i < n_x; ++i) 
                g.push_back({rx(i), ry(i), rdata(i)});
        }
        
        if(collapse) {
            //Remove duplicate edges
            std::cout << "before collapse: " << g.size() << " edges" << std::endl;
            g.template lexical_sort_by<0>();
            g.uniq();
            std::cout << "after collapse: " << g.size() << " edges" << std::endl;
        }
        edges_ = nw::graph::adjacency<0, Attributes...>(g);
        max_edge_ = edges_.size();
        nodes_ = nw::graph::adjacency<1, Attributes...>(g);
        max_node_ = nodes_.size();
    }
    std::vector<std::map<size_t, size_t>> get_edge_neighbor_counts() const { return edge_neighbor_count_; }
    std::vector<std::map<size_t, size_t>> get_node_neighbor_counts() const { return node_neighbor_count_; }

    Slinegraph<Index_t, Attributes...> s_linegraph(int s = 1, bool edges = true) {
        Slinegraph<Index_t, Attributes...> slineg(*this, s, edges);
        return slineg;
    }
    /*
    * Find the connected components for a slinegraph
    */
    py::list s_connected_components(Slinegraph<Index_t, Attributes...>& linegraph, bool return_singleton = false) {
        return linegraph.s_connected_components(return_singleton);
    }
    /*
    * Find the distance from the src to dest in its slinegraph
    */
    Index_t distance(Slinegraph<Index_t, Attributes...>& linegraph, Index_t src, Index_t dest) {
        return linegraph.s_distance(src, dest);
    }
    /*
    * Find neighbors of a vertex in its slinegraph
    */
    py::list neighbors(Slinegraph<Index_t, Attributes...>& linegraph, Index_t v) {
        return linegraph.s_neighbors(v);
    } 
    py::ssize_t s_degree(Slinegraph<Index_t, Attributes...>& linegraph, Index_t v) {
        return linegraph.s_degree(v);
    }
    /*
    * Get the edge size distributation of the hypergraph
    */
    py::list edge_size_dist() const {
        auto degs = edges_.degrees();
        py::list l = py::list(max_edge_);
        tbb::parallel_for(tbb::blocked_range<Index_t>(0, max_edge_), [&](tbb::blocked_range<Index_t>& r) {
            for (auto i = r.begin(), e = r.end(); i != e; ++i) {
                l[i] = degs[i];
            }
        }, tbb::auto_partitioner());
        return l;
    }
    /*
    * Get the node size distributation of the hypergraph
    */
    py::list node_size_dist() const {
        auto degs = nodes_.degrees();
        py::list l = py::list(max_node_);
        tbb::parallel_for(tbb::blocked_range<Index_t>(0, max_node_), [&](tbb::blocked_range<Index_t>& r) {
            for (auto i = r.begin(), e = r.end(); i != e; ++i) {
                l[i] = degs[i];
            }
        }, tbb::auto_partitioner());
        return l;
    }
    /*
    * Get the incident nodes of an edge in the hypergraph
    */
    py::list edge_incidence(Index_t edge) {
        py::list l;
        /*
        std::for_each(edges_[e].begin(), edges_[e].end(), [&](auto &&x) {
            auto hyperN = std::get<0>(x);
            l.append(hyperN);
        });
        */
        for (auto &&[hyperN, w] : edges_[edge]) {
            l.append(hyperN);
        }
        return l;
    }
    /*
    * Get the incident edges of a node in the hypergraph
    */
    py::list node_incidence(Index_t node) {
        py::list l;
        std::for_each(nodes_[node].begin(), nodes_[node].end(), [&](auto &&x) {
            auto hyperE = std::get<0>(x);
            l.append(hyperE);
        });
        return l;
    }
    /*
    * Get the degree of a node in the hypergraph.
    * If min_size and/or max_size is specified, then either/both are used to filter.
    * */
    Index_t degree(Index_t node, std::size_t min_size = 1, std::optional<std::size_t> max_size = {}) {
        Index_t deg = -1;
        // validate the input
        if (node >= max_node_)
            return deg;
        if (!max_size.has_value()) {
            //if no max_size is specified
            if (1 == min_size)
                return nodes_[node].size();
            deg = 0;
            std::for_each(nodes_[node].begin(), nodes_[node].end(), [&](auto &&x) {
                auto hyperE = std::get<0>(x);
                //filter the edges whose size is no less than min_size
                if (min_size <= edges_[hyperE].size())
                    ++deg;
            });
        }
        else {
            //if max_size is specified
            //and max_size is also 1
            if (1 == min_size && 1 == max_size)
                return nodes_[node].size();
            deg = 0;
            std::for_each(nodes_[node].begin(), nodes_[node].end(), [&](auto &&x) {
                auto hyperE = std::get<0>(x);
                //filter the edges whose size is no less than min_size
                //and is no greater than max_size
                if (min_size <= edges_[hyperE].size() && max_size >= edges_[hyperE].size())
                    ++deg;
            });
        }
        return deg;
    }
    py::ssize_t number_of_nodes() const { return max_node_; }
    py::ssize_t order() const { return max_node_; }
    py::ssize_t size(Index_t edge, std::size_t min_degree = 1, std::optional<std::size_t> max_degree = {}) {
        Index_t size = -1;
        // validate the input
        if (edge >= max_edge_)
            return size;
        if (!max_degree.has_value()) {
            //if no max_size is specified
            if (1 == min_degree)
                return edges_[edge].size();
            size = 0;
            std::for_each(edges_[edge].begin(), edges_[edge].end(), [&](auto &&x) {
                auto hyperE = std::get<0>(x);
                //filter the edges whose size is no less than min_degree
                if (min_degree <= edges_[hyperE].size())
                    ++size;
            });
        }
        else {
            //if max_degree is specified
            //and max_degree is also 1
            if (1 == min_degree && 1 == max_degree)
                return edges_[edge].size();
            size = 0;
            std::for_each(edges_[edge].begin(), edges_[edge].end(), [&](auto &&x) {
                auto hyperN = std::get<0>(x);
                //filter the edges whose size is no less than min_degree
                //and is no greater than max_degree
                if (min_degree <= nodes_[hyperN].size() && max_degree >= nodes_[hyperN].size())
                    ++size;
            });
        }
        return size;
    }
    py::ssize_t dim(Index_t edge) {
        if (edge >= max_edge_)
            return -1;
        else
            return edges_[edge].size() - 1;
    }
    py::ssize_t number_of_edges() const { return max_edge_; }
    // a singleton is an edge of size 1
    py::list singletons() {
        py::list l;
        for(Index_t e = 0; e < max_edge_; ++e) {
            if (1 == edges_[e].size()) {
                std::for_each(edges_[e].begin(), edges_[e].end(), [&](auto &&x) {
                    auto hyperN = std::get<0>(x);
                    if (1 == nodes_[hyperN].size())
                        l.append(e);
                });
            }
        }
        return l;
    }
    /*
    * Verified. using intersection_size instead set_intersection
    * */
    py::list toplexes_v2() {
        std::vector<Index_t> tops;
        //std::map<Index_t, bool> tops;
        for (Index_t e = 0; e < max_edge_; ++e) {
            bool flag = true;
            std::vector<Index_t> old_tops(tops);
            for (size_t i = 0, end = old_tops.size(); i < end; ++i) {
                Index_t top = old_tops[i];
                //TODO is this necessary?
                if (e == top)
                    continue;
                auto lhs = edges_[e].size();
                auto rhs = edges_[top].size();
                auto s = nw::graph::intersection_size(edges_[e], edges_[top]);
                if (s == rhs)
                    tops.erase(tops.begin() + i);
                else if (s == lhs) {
                    flag = false;
                    break;
                }
            }
            if (flag)
                tops.push_back(e);
        }

        auto n = tops.size();
        py::list l = py::list(n);
        tbb::parallel_for(
            tbb::blocked_range<Index_t>(0, n), [&](tbb::blocked_range<Index_t> &r) {
                for (auto i = r.begin(), e = r.end(); i != e; ++i)
                    l[i] = tops[i];
            },
            tbb::auto_partitioner());
        return l;
    }
    /*
    * Verified. using unordered_set
    * */
    py::list toplexes_v1() {
        //check rhs is a subset of lhs
        //create a freq table for all the elements of lhs
        //traverse rhs and search for each element of rhs in the freq table
        //if element is found , then decrease the frequency, otherwise, return false
        //if all elements are found, return true
        //O(m+n), where m is the size of lhs, n is the size of rhs
        auto issubset = []<class A>(A&& lhs, A&& rhs) {
            std::map<Index_t, size_t> frequency;
            std::for_each (lhs.begin(), lhs.end(), [&](auto&& x) {
                auto v = std::get<0>(x);
                ++frequency[v];
            });
            bool res = true;
            std::for_each (rhs.begin(), rhs.end(), [&](auto&& x) {
                auto v = std::get<0>(x);
                if (0 < frequency[v])
                    --frequency[v];
                else {
                    res = false;
                    return;
                }
            }); 
            return res;
        };

        //create an empty toplex set and an empty old toplex set
        std::unordered_set<Index_t> tops;
        for (Index_t e = 0; e < max_edge_; ++e) {
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
            for (auto& top : old_tops) {
                if (e == top)
                    continue; 
               if (issubset(edges_[top], edges_[e])) {
                    //if e is a subset of top, then e is not a toplex
                    flag = false;
                    break;
                } else if (issubset(edges_[e], edges_[top]))
                    tops.erase(top);
            }//for old_tops
            if (flag)
                tops.insert(e);
        }//for each e

        py::list l;
        for (auto& e : tops) {
            l.append(e);
        }
        return l;
    }
    py::list toplexes() {
        //check rhs is a subset of lhs
        //create a freq table for all the elements of lhs
        //traverse rhs and search for each element of rhs in the freq table
        //if element is found , then decrease the frequency, otherwise, return false
        //if all elements are found, return true
        //O(m+n), where m is the size of lhs, n is the size of rhs
        auto issubset = []<class A>(A&& lhs, A&& rhs) {
            std::map<Index_t, size_t> frequency;
            std::for_each (lhs.begin(), lhs.end(), [&](auto&& x) {
                auto v = std::get<0>(x);
                ++frequency[v];
            });
            bool res = true;
            std::for_each (rhs.begin(), rhs.end(), [&](auto&& x) {
                auto v = std::get<0>(x);
                if (0 < frequency[v])
                    --frequency[v];
                else {
                    res = false;
                    return;
                }
            }); 
            return res;
        };

        //create an empty toplex set and an empty old toplex set
        std::vector<Index_t> tops;
        for (Index_t e = 0; e < max_edge_; ++e) {
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
                Index_t top = old_tops[i];
                if (e == top)
                    continue; 
               if (issubset(edges_[top], edges_[e])) {
                    //if e is a subset of top, then e is not a toplex
                    flag = false;
                    break;
                } else if (issubset(edges_[e], edges_[top]))
                    tops.erase(tops.begin() + i);
            }//for old_tops
            if (flag)
                tops.push_back(e);
        }//for each e

        auto n = tops.size();
        py::list l = py::list(n);
        tbb::parallel_for(
            tbb::blocked_range<Index_t>(0, n), [&](tbb::blocked_range<Index_t> &r) {
                for (auto i = r.begin(), e = r.end(); i != e; ++i)
                    l[i] = tops[i];
            }, tbb::auto_partitioner());
        return l;
    }
protected:
    void populate_edge_linegraph (nw::graph::edge_list<nw::graph::undirected, Attributes...>&& linegraph, size_t s) {
        //if neighbors have not been counted
        if (0 == edge_neighbor_count_.size()) {
            edge_neighbor_count_ = to_two_graph_count_neighbors_cyclic(edges_, nodes_);
        }
        populate_linegraph_from_neighbor_map<decltype(edges_), Attributes...>(edges_, edge_neighbor_count_, std::move(linegraph), s, false);
    }
    auto populate_node_linegraph (nw::graph::edge_list<nw::graph::undirected, Attributes...>&& linegraph, size_t s) {
        //if neighbors have not been counted
        if (0 == node_neighbor_count_.size()) {
            node_neighbor_count_ = to_two_graph_count_neighbors_cyclic(nodes_, edges_);
        }
        populate_linegraph_from_neighbor_map<decltype(nodes_), Attributes...>(nodes_, node_neighbor_count_, std::move(linegraph), s, false);
    }
}; //class NWhypergraph

template<class Index_t, typename... Attributes> 
class Slinegraph {
private:
    NWHypergraph<Index_t, Attributes...>& hyperg_;
    bool edges_;

    //s-linegraph generated from original hypergraph
    nw::graph::adjacency<0, Attributes...> g_;
    nw::graph::adjacency<1, Attributes...> g_t_;
public:
    py::array_t<Index_t, py::array::c_style> row_;
    py::array_t<Index_t, py::array::c_style> col_;
    py::array_t<Attributes..., py::array::c_style> data_;
    int s_;
private:
    void populate_adjacency(nw::graph::edge_list<nw::graph::undirected, Attributes...>& linegraph) {
        size_t nedges = linegraph.size();
        std::cout << "linegraph size: " << nedges << std::endl;

        //instantiate empty adjacency
        size_t m = edges_ ? hyperg_.edges_.size() : hyperg_.nodes_.size();
        if (0 == nedges) {
            //create an adjacency where each list is empty
            g_ = nw::graph::adjacency<0, Attributes...>(m);
            g_t_ = nw::graph::adjacency<1, Attributes...>(m);
        }
        //fill adjacency with edges
        else {
            g_ = nw::graph::adjacency<0, Attributes...>(m, linegraph);
            g_t_ = nw::graph::adjacency<1, Attributes...>(m, linegraph);
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
    Slinegraph(NWHypergraph<Index_t, Attributes...>&g, int s = 1, bool edges = true) : hyperg_(g), edges_(edges), s_(s) {
        //extract linegraph from its neighbor counts
        if (edges_) {
            nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph;
            hyperg_.populate_edge_linegraph(std::move(linegraph), s);
            populate_adjacency(linegraph);
            populate_py_array(linegraph);
        }
        else{
            nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph;
            hyperg_.populate_node_linegraph(std::move(linegraph), s);         
            populate_adjacency(linegraph);
            populate_py_array(linegraph);
        }
    }
    /*
    * Return a list of sets which
    * each set is a component.
    */
    py::list s_connected_components(bool return_singleton = false) {
        auto E = nw::graph::ccv1(g_);
        // This returns the subgraph of each component.
        std::map<Index_t, py::set> comps;
        for (size_t i = 0, e = E.size(); i < e; ++i) {
            auto label = E[i];
            comps[label].add(i);
        }
        py::list l;
        for (auto it = comps.begin(); it != comps.end(); ++it) {
            auto s = it->second;
            py::ssize_t n = it->second.size();
            if (true == return_singleton) {
                l.append(s);
            }
            else {
                if (1 < n) {
                    l.append(s);
                }
            }
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
        tbb::parallel_for(tbb::blocked_range<Index_t>(0, n), [&](tbb::blocked_range<Index_t>& r) {
            for (auto i = r.begin(), e = r.end(); i != e; ++i) {
                if (!normalized)
                    l[i] = bc[i] / 2;
                else
                    l[i] = scale * bc[i];
            }
        }, tbb::auto_partitioner());
        return l;
    }
    /*
    * Closeness centrality of a node `v` is the reciprocal of the
    * average shortest path distance to `v` over all `n-1` reachable nodes.
    */
    py::list s_closeness_centrality(std::optional<Index_t> v = {}) {
        using distance_t = std::uint64_t;
        std::size_t delta = 1;
        std::size_t n = g_.size();
        if (v.has_value()) {
            py::list l;
            //validate input
            if (v.value() >= (Index_t)g_.size())
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
            tbb::parallel_for(tbb::blocked_range<Index_t>(0, n), [&](auto& r) {
                for (Index_t v = r.begin(); v != r.end(); ++v) {
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
            });
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
        std::size_t n = g_.size();
        if (v.has_value()) {
            py::list l;
            //validate input
            if (v >= (Index_t)g_.size())
                return l;
            auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v.value(), delta);
            float sum = 0;
            for (auto& d : dist) {
                distance_t tmp = d.load(std::memory_order_relaxed);
                if (tmp != std::numeric_limits<distance_t>::max() && 0 != tmp)
                    sum +=  1.0 / tmp;
            }  
            l.append(1.0 * sum / (n - 1));
            return l;
        }
        else {
            py::list l = py::list(n);
            //for each vertex v in the graph
            tbb::parallel_for(tbb::blocked_range<Index_t>(0, n), [&](auto& r) {
                for (Index_t v = r.begin(); v != r.end(); ++v) {
                    //get the distances from v to any other vertices
                    auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v, delta);
                    float sum = 0;
                    for (auto &d : dist) {
                        distance_t tmp = d.load(std::memory_order_relaxed);
                        if (tmp != std::numeric_limits<distance_t>::max() && 0 != tmp)
                            sum += 1.0 / tmp;
                    }
                    l[v] = 1.0 * sum / (n - 1);
                }
            });
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
        std::size_t n = g_.size();
        if (v.has_value()) {
            py::list l;
            //validate input
            if (v >= (Index_t)g_.size())
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
            tbb::parallel_for(tbb::blocked_range<Index_t>(0, n), [&](auto& r) {
                for (Index_t v = r.begin(); v != r.end(); ++v) {
                    //get the distances from v to any other vertices
                    auto dist = nw::graph::delta_stepping_v12<distance_t>(g_, v, delta);
                    auto tmp = std::max_element(dist.begin(), dist.end());
                    distance_t result = (*tmp).load(std::memory_order_relaxed);
                    if (result == std::numeric_limits<distance_t>::max())
                        l[v] = std::numeric_limits<distance_t>::infinity();
                    else
                        l[v] = result;
                }
            });            
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
    py::ssize_t s_degree(Index_t v) { return g_[v].size(); }
    int getS() const { return s_; }
    bool isEdgeOverlap() const { return edges_; }
}; //slinegraph class


}//namespace hypergraph
}//namespace nw