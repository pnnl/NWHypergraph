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
#include "algorithms/slinegraph_efficient.hpp"
#include "algorithms/slinegraph_map.hpp"
#include "algorithms/slinegraph_naive.hpp"
#include "algorithms/slinegraph_adjoin.hpp"
#include "algorithms/slinegraph_kij.hpp"

namespace nw {
namespace hypergraph {

/*
* Below are different versions of sline graph computation algorithms.
*/
enum soverlap_version {
  Efficient_Blocked = 0,
  Efficient_Cyclic = 1,
  Naive = 2,
  Map_Blocked = 3,
  Map_Cyclic = 4,
  Ensemble_Blocked = 5,
  Ensemble_Cyclic = 6,
  Map_Frontier_Blocked = 7,
  Map_Frontier_Cyclic = 8,
  HashMap_Frontier_Blocked = 9,
  HashMap_Frontier_Cyclic = 10,
  Efficient_Frontier_Blocked = 11,
  Efficient_Frontier_Cyclic = 12,
  HashMap_Blocked = 13,
  HashMap_Cyclic = 14,
  Vector_Blocked = 15,
  Vector_Cyclic = 16,
  Static_HashMap_Blocked = 17,
  Static_HashMap_Cyclic = 18,
  Efficient_Blocked_Size = 19,
  Map_Blocked_Size = 20,
  HashMap_Blocked_Size = 21,
  Vector_Blocked_Size = 22,
  Static_HashMap_Blocked_Size = 23,
  SPGEMM_KIJ = 99, //TODO not ready
  Unknown
};

soverlap_version operator++(soverlap_version& value, int) {
  soverlap_version current = value;
  if (Unknown < current + 1) value = Efficient_Blocked;
  else value = static_cast<soverlap_version>(value + 1);

  return current;
}

std::ostream& operator<<(std::ostream& out, const soverlap_version value){
    static std::map<soverlap_version, std::string> strings;
    if (0 == strings.size()){
#define INSERT_ELEMENT(p) strings[p] = #p
        INSERT_ELEMENT(Efficient_Blocked);     
        INSERT_ELEMENT(Efficient_Cyclic);     
        INSERT_ELEMENT(Naive);     
        INSERT_ELEMENT(Map_Blocked);     
        INSERT_ELEMENT(Map_Cyclic);     
        INSERT_ELEMENT(Ensemble_Blocked);      
        INSERT_ELEMENT(Ensemble_Cyclic);     
        INSERT_ELEMENT(Map_Frontier_Blocked);     
        INSERT_ELEMENT(Map_Frontier_Cyclic);  
        INSERT_ELEMENT(HashMap_Frontier_Blocked);      
        INSERT_ELEMENT(HashMap_Frontier_Cyclic);     
        INSERT_ELEMENT(Efficient_Frontier_Blocked);     
        INSERT_ELEMENT(Efficient_Frontier_Cyclic); 
        INSERT_ELEMENT(HashMap_Blocked);
        INSERT_ELEMENT(HashMap_Cyclic);  
        INSERT_ELEMENT(Vector_Blocked);  
        INSERT_ELEMENT(Vector_Cyclic);
        INSERT_ELEMENT(Static_HashMap_Blocked);
        INSERT_ELEMENT(Static_HashMap_Cyclic); 
        INSERT_ELEMENT(Efficient_Blocked_Size);  
        INSERT_ELEMENT(Map_Blocked_Size);
        INSERT_ELEMENT(HashMap_Blocked_Size);
        INSERT_ELEMENT(Vector_Blocked_Size); 
        INSERT_ELEMENT(Static_HashMap_Blocked_Size); 
        INSERT_ELEMENT(Unknown);           
#undef INSERT_ELEMENT
    }   

    return out << strings[value];
}
/**
* This function prints the versions of the sline graph algorithms. 
*
*/
void print_soverlap_version() {
  for (int i = Efficient_Blocked; i < Unknown; ++i) {
    std::cout << i << ")" << static_cast<soverlap_version>(i) << " ";
  }
  std::cout << std::endl;
}
/**
* This function computes the s-line graph of a hypergraph and returns the edge list of the s-line graph. 
* Only unweighted hypergraph is supported for now.
*
* @param[in] version s-line graph algorithm version id
* @param[in] verbose flag to control verbose stats collection
* @param[in] features a fixed-size sequence of 8 bits to control which heuristics are enabled in the efficient version
* @param[in] edges adjacency for hyperedges
* @param[in] nodes adjacency for hypernodes
* @param[in] iperm permutation array of hyperedges/nodes IDs when adjoin is enabled
* @param[in] nrealedges the number of real hyperedges when adjoin is enabled
* @param[in] nrealnodes the number of real nodes when adjoin is enabled
* @param[in] s the number of the overlapping vertices between hyperedge pairs
* @param[in] num_threads the number of threads to run
* @param[in] num_bins the number of bins to divide the workload
* @returns the edge list of the s-line graph
*
*/
auto twograph_reader(int version, bool verbose, std::bitset<8> &features,
                     nw::graph::adjacency<0> &edges, nw::graph::adjacency<1> &nodes,
                     std::vector<nw::graph::index_t> &edgedegrees,
                     std::vector<vertex_id_t>& iperm,
                     size_t nrealedges, size_t nrealnodes, 
                     size_t s,
                     int num_threads, 
                     int num_bins = 32) {
  switch (version) {
    case Efficient_Blocked: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_efficient_blocked_portal<undirected>(
              verbose, features, std::execution::par_unseq, edges, nodes,
              edgedegrees, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Efficient_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_efficient_parallel_cyclic_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Naive: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_naive_parallel_portal<undirected>(
              verbose, edges, nodes, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Map_Blocked: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_map_blocked_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Map_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_map_cyclic_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Ensemble_Blocked: {
      std::vector<std::unordered_map<size_t, size_t>>&& neighbor_count =
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
    case Ensemble_Cyclic: {
      std::vector<std::unordered_map<size_t, size_t>>&& neighbor_count =
          to_two_graph_count_neighbors_cyclic(edges, nodes);
      nw::graph::edge_list<undirected> &&linegraph =
          populate_linegraph_from_neighbor_map<undirected>(neighbor_count, s);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Map_Frontier_Blocked: {
      nw::graph::edge_list<undirected>&& linegraph =
          to_two_graph_map_frontier_blocked_portal<undirected>(
              std::execution::par_unseq, edges, nodes, edgedegrees,
              iperm,
              nrealedges, nrealnodes, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have
      // two elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Map_Frontier_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_map_frontier_cyclic_portal<undirected>(
              std::execution::par_unseq, edges, nodes, edgedegrees,
              iperm,
              nrealedges, nrealnodes, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have
      // two elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case HashMap_Frontier_Blocked: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_hashmap_frontier_blocked_portal<undirected>(
              std::execution::par_unseq, edges, nodes, edgedegrees,
              iperm,
              nrealedges, nrealnodes, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have
      // two elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case HashMap_Frontier_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_hashmap_frontier_cyclic_portal<undirected>(
              std::execution::par_unseq, edges, nodes, edgedegrees,
              iperm,
              nrealedges, nrealnodes, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have
      // two elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Efficient_Frontier_Blocked: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_efficient_frontier_blocked_portal<undirected>(
              std::execution::par_unseq, edges, nodes, edgedegrees,
              iperm,
              nrealedges, nrealnodes, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have
      // two elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Efficient_Frontier_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_efficient_frontier_cyclic_portal<undirected>(
              std::execution::par_unseq, edges, nodes, edgedegrees,
              iperm,
              nrealedges, nrealnodes, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have
      // two elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case HashMap_Blocked: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_hashmap_blocked_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case HashMap_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_hashmap_cyclic_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Vector_Blocked: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_vector_blocked_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Vector_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_vector_cyclic_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Static_HashMap_Blocked: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_static_hashmap_blocked_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Static_HashMap_Cyclic: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_static_hashmap_cyclic_portal<undirected>(
              verbose, std::execution::par_unseq, edges, nodes, edgedegrees, s,
              num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Efficient_Blocked_Size: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_efficient_blocked_vary_size<undirected>(
              std::execution::par_unseq, edges, nodes,
              edgedegrees, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Map_Blocked_Size: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_map_blocked_vary_size<undirected>(
              std::execution::par_unseq, edges, nodes,
              edgedegrees, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case HashMap_Blocked_Size: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_hashmap_blocked_vary_size<undirected>(
              std::execution::par_unseq, edges, nodes,
              edgedegrees, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }  
    case Vector_Blocked_Size: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_vector_blocked_vary_size<undirected>(
              std::execution::par_unseq, edges, nodes,
              edgedegrees, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    } 
    case Static_HashMap_Blocked_Size: {
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_static_hashmap_blocked_vary_size<undirected>(
              std::execution::par_unseq, edges, nodes,
              edgedegrees, s, num_threads, num_bins);
      // where when an empty edge list is passed in, an adjacency still have two
      // elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    } 
    case SPGEMM_KIJ: {
      std::cout << "graph edges = " << nodes.size() << std::endl;
      nw::graph::edge_list<undirected> &&linegraph =
          to_two_graph_spgemm_kij_cyclic_v2<undirected>(
              std::execution::seq, edges, nodes, edgedegrees,
              s, num_bins);
      // where when an empty edge list is passed in, an adjacency still have
      // two elements
      if (0 == linegraph.size()) return nw::graph::adjacency<0>(0, 0);
      nw::graph::adjacency<0> s_adj(linegraph);
      std::cout << "line graph edges = " << linegraph.size()
                << ", adjacency size = " << s_adj.size()
                << ", max = " << s_adj.max() << std::endl;
      return s_adj;
    }
    case Unknown:
    default: {
      std::cerr << "unknown soverlap computation loader" << std::endl;
      return nw::graph::adjacency<0>(0);
    }
  }
};

}  // namespace hypergraph
}  // namespace nw
