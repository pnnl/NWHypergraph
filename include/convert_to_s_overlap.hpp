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
#include <map>
#include <edge_list.hpp>
#include <util/intersection_size.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include "s_overlap.hpp"
#include <algorithms/connected_components.hpp>
#include <algorithms/delta_stepping.hpp>
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
        std::cout << "counting neighbors\n";
        if (edges){
            edge_neighbor_count_ = to_two_graph_count_neighbors_parallel(edges_, nodes_);
        }
        else {
            node_neighbor_count_ = to_two_graph_count_neighbors_parallel(nodes_, edges_);
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
        g.open_for_push_back();
        for (size_t i = 0; i < n_x; ++i) {
            g.push_back({rx(i), ry(i), 0});
        }
        g.close_for_push_back();
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
        g.open_for_push_back();
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
        g.close_for_push_back();
        
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

    nw::graph::edge_list<nw::graph::undirected, Attributes...> populate_edge_linegraph (size_t s) {
        nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph;
        size_t M = edges_.size();
        linegraph.open_for_push_back();
        //n is the number of hyperedges, m is the number of hypernodes
        //time complexity of counting neighbors is same as the efficient: O(n*deg(edges)*deg(nodes)*deg(edges))
        //time complexity of extract slinegraph from the neighbor counts: O(n*deg(edges)) -> worst is O(n^2)
        //space complexity: O(n*total_deg(H)) -> worst is O(n^2)
        //total_deg(H)=sum of deg(each hyperedge)
        //total_deg(H) >> n?
        for (size_t hyperE = 0; hyperE < M; ++hyperE) {
            for (auto &&[anotherhyperE, val] : edge_neighbor_count_[hyperE]) {
                if (val >= s) {
                    //std::cout << hyperE << "-" << anotherhyperE << std::endl;
                    linegraph.push_back(std::make_tuple<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperE), std::forward<vertex_id_t>(anotherhyperE), val));
                }
            }
        }
        linegraph.close_for_push_back();
        return linegraph;
    }
    nw::graph::edge_list<nw::graph::undirected, Attributes...> populate_node_linegraph (size_t s) {
        nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph;
        size_t N = nodes_.size();
        linegraph.open_for_push_back();
        for (size_t hyperN = 0; hyperN < N; ++hyperN) {
            for (auto &&[anotherhyperN, val] : node_neighbor_count_[hyperN]) {
                if (val >= s) 
                    linegraph.push_back(std::make_tuple<vertex_id_t, vertex_id_t>(std::forward<vertex_id_t>(hyperN), std::forward<vertex_id_t>(anotherhyperN), val));
            }
        }
        linegraph.close_for_push_back();
        return linegraph;
    }
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
    * */
    Index_t degree(Index_t node, size_t size = 1) {
        if (node >= max_node_)
            return -1;
        else
            return nodes_[node].size();
    }
    py::ssize_t number_of_nodes() const { return max_node_; }
    py::ssize_t order() const { return max_node_; }
    py::ssize_t size(Index_t edge) {
        if (edge >= max_edge_)
            return -1;
        else
            return edges_[edge].size();
    }
    py::ssize_t dim(Index_t edge) {
        if (edge >= max_edge_)
            return -1;
        else
            return edges_[edge].size() - 1;
    }
    py::ssize_t number_of_edges() const { return max_edge_; }
    // a singleton is an edge of size 1 with a node of degree 1
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
        def toplexes(self,name=None,collapse=False,use_reps=False,return_counts=True):
        """
        Returns a :term:`simple hypergraph` corresponding to self.
        Warning
        -------
        Collapsing a hypergraph can take a long time. It may be preferable to collapse the graph first and
        pickle it, then apply the toplexes method separately.
        Parameters
        ----------
        name: str, optional, default: None
        collapse: boolean, optional, default: False
            Should the hypergraph be collapsed? This would preserve a link between duplicate maximal sets.
            If False then only one of these sets will be used and uniqueness will be up to sets of equal size.
        use_reps: boolean, optional, default: False
            If collapse=True then each toplex will be named by a representative of the set of
            equivalent edges, default is False (see collapse_edges).
        return_counts: boolean, optional, default: True
            If collapse=True then each toplex will be named by a tuple of the representative
            of the set of equivalent edges and their count
        """
        if collapse:
            if len(self.edges) > 20:  ### TODO: Determine how big is too big.
                warnings.warn('Collapsing a hypergraph can take a long time. It may be preferable to collapse the graph first and pickle it then apply the toplex method separately.')
            temp = self.collapse_edges(use_reps=use_reps,return_counts=return_counts)
        else:
            temp = self
        thdict = dict()
        for e in temp.edges:
            thdict[e] = temp.edges[e].uidset
        tops = dict()
        for e in temp.edges:
            flag = True
            old_tops = dict(tops)
            for top in old_tops:
                if thdict[e].issubset(thdict[top]):
                    flag = False
                    break
                elif set(thdict[top]).issubset(thdict[e]):
                    del tops[top]
            if flag:
                tops.update({e : thdict[e]})
        return Hypergraph(tops,name)

        */
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
    py::list toplexes() {
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
        std::vector<Index_t> tops;
        for (Index_t e = 0; e < max_edge_; ++e) {
            bool flag = true;
            std::vector<Index_t> old_tops(tops);
            for (size_t i = 0, end = old_tops.size(); i < end; ++i) {
                Index_t top = old_tops[i];
                //TODO is this necessary?
                if (e == top)
                    continue;
                if (issubset(edges_[top], edges_[e]))
                    tops.erase(tops.begin() + i);
                else if (issubset(edges_[e], edges_[top])) {
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
};

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
        nw::graph::edge_list<nw::graph::undirected, Attributes...> linegraph;
        if (edges_) {
            //if neighbors have not been counted
            if (0 == hyperg_.edge_neighbor_count_.size()) {
                hyperg_.populate_neighbor_count(edges_);
            }
            linegraph = hyperg_.populate_edge_linegraph(s);
        }
        else{
            if (0 == hyperg_.node_neighbor_count_.size()) {
                hyperg_.populate_neighbor_count(edges_);
            }
            linegraph = hyperg_.populate_node_linegraph(s);         
        }
        populate_adjacency(linegraph);
        populate_py_array(linegraph);
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
        std::cout << std::endl;
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
    * Compute the distance from src to dest
    * TODO -1 if unreachable for now
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


/*
* Deprecated
*/
template<class Index_t, class Data_t>
std::tuple<py::array_t<Index_t>, py::array_t<Index_t>, py::array_t<Data_t>, 
py::array_t<Index_t>, py::array_t<Index_t>, py::array_t<Data_t>> 
convert_to_s_overlap(py::array_t<Index_t, py::array::c_style | py::array::forcecast> &x, 
py::array_t<Index_t, py::array::c_style | py::array::forcecast> &y, 
py::array_t<Data_t, py::array::c_style | py::array::forcecast> &data, 
size_t s = 1) {
    //sanitize check
    auto rx = x.template mutable_unchecked<1>();
    auto ry = y.template mutable_unchecked<1>();
    auto rdata = data.template mutable_unchecked<1>();
    //rx(0) = 1;
    size_t n_x = x.shape(0);
    size_t n_y = y.shape(0);
    //assume there are n_x pairs
    size_t n_data = data.shape(0);
    //create edge list 
    if (0 == n_data) {
        std::cout << "reading unweighted hypergraph" << std::endl;
        nw::graph::edge_list<nw::graph::directed> aos_a(0); 
        aos_a.open_for_push_back();
        for (size_t i = 0; i < n_x; ++i) {
            aos_a.push_back(rx(i), ry(i));
        }
        aos_a.close_for_push_back();

        nw::graph::adjacency<0> hyperedges(aos_a);
        nw::graph::adjacency<1> hypernodes(aos_a);
        std::vector<nw::graph::index_t> hyperedge_degrees = aos_a.degrees<0>();

        int num_bins = 32;
        auto&& two_graphs = 
        to_two_graph_efficient_parallel_clean_without_sequeeze<nw::graph::undirected>(std::execution::par_unseq, hyperedges, hypernodes, hyperedge_degrees, s, num_bins);
        if (1 == s) {
            nw::graph::edge_list<nw::graph::undirected> &&linegraph = create_edgelist_without_squeeze<nw::graph::undirected>(two_graphs);
            //where when an empty edge list is passed in, an adjacency still have two elements
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n}); 
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    ++i;
                }
                return std::make_tuple(newx, newy, py::array_t<Data_t>(0), 
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));
        }
        else {
            nw::graph::edge_list<nw::graph::undirected> &&linegraph = create_edgelist_with_squeeze<nw::graph::undirected>(two_graphs);
            nw::graph::edge_list<nw::graph::undirected> &&raw_linegraph = create_edgelist_without_squeeze<nw::graph::undirected>(two_graphs);           
            //where when an empty edge list is passed in, an adjacency still have two elements
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n}), oldx({n}), oldy({n}); 
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            auto roldx = oldx.template mutable_unchecked<1>();
            auto roldy = oldy.template mutable_unchecked<1>();
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    ++i;
                }
                i = 0;
                for (auto &&[u, v] : raw_linegraph) {
                    roldx(i) = u;
                    roldy(i) = v;
                    ++i;
                }
                return std::make_tuple(newx, newy, py::array_t<Data_t>(0), oldx, oldy, py::array_t<Data_t>(0));
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));
        }//else 1 == s
    }
    else {
        std::cout << "reading weighted hypergraph" << std::endl;

        nw::graph::edge_list<nw::graph::directed, Data_t> aos_a(0); 
        aos_a.open_for_push_back();
        for (size_t i = 0; i < n_x; ++i) {
            aos_a.push_back(rx(i), ry(i), rdata(i));
        }
        aos_a.close_for_push_back();

        nw::graph::adjacency<0, Data_t> hyperedges(aos_a);
        nw::graph::adjacency<1, Data_t> hypernodes(aos_a);
        std::vector<nw::graph::index_t> hyperedge_degrees = aos_a.template degrees<0>();
    
        int num_bins = 32;
        auto&& two_graphs = 
        to_two_graph_weighted_efficient_parallel_clean_without_squeeze<nw::graph::undirected, Data_t>(std::execution::par_unseq, hyperedges, hypernodes, hyperedge_degrees, s, num_bins);
        if (1 == s) {
            nw::graph::edge_list<nw::graph::undirected, Data_t> &&linegraph = create_edgelist_without_squeeze<nw::graph::undirected, Data_t>(two_graphs);
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            //create results
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n});
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            py::array_t<Data_t, py::array::c_style> newdata({n}); 
            auto rnewdata = newdata.template mutable_unchecked<1>(); 
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v, weight] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    rnewdata(i) = weight;
                    ++i;
                }
                return std::make_tuple(newx, newy, newdata, 
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Index_t>(0));
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));   
        }
        else {
            nw::graph::edge_list<nw::graph::undirected, Data_t> &&linegraph = create_edgelist_with_squeeze<nw::graph::undirected, Data_t>(two_graphs);
            nw::graph::edge_list<nw::graph::undirected, Data_t> &&raw_linegraph = create_edgelist_without_squeeze<nw::graph::undirected, Data_t>(two_graphs);           
            pybind11::ssize_t n = linegraph.size();
            std::cout << "linegraph has " << n << " edges" << std::endl;
            std::cout << "rawlinegraph has " << raw_linegraph.size() << " edges" << std::endl;
            py::array_t<Index_t, py::array::c_style> newx({n}), newy({n}), oldx({n}), oldy({n}); 
            auto rnewx = newx.template mutable_unchecked<1>();
            auto rnewy = newy.template mutable_unchecked<1>(); 
            auto roldx = oldx.template mutable_unchecked<1>();
            auto roldy = oldy.template mutable_unchecked<1>();
            py::array_t<Data_t, py::array::c_style> newdata({n}), olddata({n}); 
            auto rnewdata = newdata.template mutable_unchecked<1>(); 
            auto rolddata = olddata.template mutable_unchecked<1>(); 
            if (0 != n) {
                size_t i = 0;
                for (auto &&[u, v, weight] : linegraph) {
                    rnewx(i) = u;
                    rnewy(i) = v;
                    rnewdata(i) = weight;
                    ++i;
                }
                i = 0;
                for (auto &&[u, v, weight] : raw_linegraph) {
                    roldx(i) = u;
                    roldy(i) = v;
                    rolddata(i) = weight;
                    ++i;
                }
                return std::make_tuple(newx, newy, newdata, oldx, oldy, olddata);                  
            }
            else
                return std::make_tuple(py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0),
                py::array_t<Index_t>(0), py::array_t<Index_t>(0), py::array_t<Data_t>(0));         
        }//else 1 == s
    }//else 0 == n_data
}

}//namespace hypergraph
}//namespace nw