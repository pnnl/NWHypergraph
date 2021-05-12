//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2020
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0
// International License https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#pragma once
#include <util/timer.hpp>

#include "algorithms/slinegraph_efficient.hpp"
#include "algorithms/slinegraph_map.hpp"
#include "algorithms/slinegraph_naive.hpp"

using namespace nw::graph;

namespace nw {
namespace hypergraph {

auto twograph_reader(int loader_version, bool verbose, std::bitset<8> &features,
                     adjacency<0> &edges, adjacency<1> &nodes,
                     std::vector<nw::graph::index_t> &edgedegrees, size_t s = 1,
                     int num_bins = 32) {
  switch (loader_version) {
    case 0: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_efficient_parallel_portal<undirected>(
              verbose, features, std::execution::par_unseq, edges, nodes,
              edgedegrees, s, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max= " << s_adj.max() << std::endl;
      return s_adj;
    }
    case 1: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_efficient_parallel_cyclic_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max= " << s_adj.max() << std::endl;
      return s_adj;
    }
    case 2: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_naive_parallel_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, s, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max= " << s_adj.max() << std::endl;
      return s_adj;
    }
    case 3: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_map_blocked_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max= " << s_adj.max() << std::endl;
      return s_adj;
    }
    case 4: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_map_cyclic_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max= " << s_adj.max() << std::endl;
      return s_adj;
    }
    case 5: {
      std::vector<std::map<size_t, size_t>> neighbor_count =
          to_two_graph_count_neighbors_blocked(edges, nodes);
      nw::graph::edge_list<undirected> &&linegraph =
          populate_linegraph_from_neighbor_map<undirected>(neighbor_count, s);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max= " << s_adj.max() << std::endl;
      return s_adj;
    }
    case 6: {
      std::vector<std::map<size_t, size_t>> neighbor_count =
          to_two_graph_count_neighbors_cyclic(edges, nodes);
      nw::graph::edge_list<undirected> &&linegraph =
          populate_linegraph_from_neighbor_map<undirected>(neighbor_count, s);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max= " << s_adj.max() << std::endl;
      return s_adj;
    }
    default: {
      std::cerr << "unknown soverlap computation loader" << std::endl;
      return nw::graph::adjacency<0>(0);
    }
  }
};

}  // namespace hypergraph
}  // namespace nw
