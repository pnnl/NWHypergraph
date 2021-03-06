/**
 * @file hypergraph_io.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#pragma once
#include <iostream>
#include <fstream>

#include <nwgraph/edge_list.hpp>
#include <nwgraph/adjacency.hpp>
using namespace nw::graph;

namespace nw {
namespace hypergraph {

std::string AdjHypergraphHeader = "AdjacencyHypergraph";
std::string WghAdjHypergraphHeader = "WeightedAdjacencyHypergraph";
/* The detail of the file format is at
* https://github.com/jshun/ligra/tree/ligra-h#input-format-for-ligra-h-applications
*/

/* 
* fill biadjacency.
* return adjacency<0> and adjacency<1>.
*/
template <std::unsigned_integral vertex_id_t>
auto adj_hypergraph_fill(std::istream& inputStream) {
  //100, 972 100 972
  size_t n0, m0, n1, m1;
  inputStream >> n0;
  inputStream >> m0;
  inputStream >> n1;
  inputStream >> m1;

  vertex_id_t tmp;
  std::vector<vertex_id_t> v0(n0 + 1);
  for (size_t i = 0; i < n0; ++i) {
    inputStream >> tmp;
    v0[i] = tmp;
  }
  std::vector<vertex_id_t> e0(m0);
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    e0[i] = tmp;
  }
  v0[n0] = m0;

  std::vector<vertex_id_t> v1(n1 + 1);
  for (size_t i = 0; i < n1; ++i) {
    inputStream >> tmp;
    v1[i] = tmp;
  }
  std::vector<vertex_id_t> e1(m1);
  for (size_t i = 0; i < m1; ++i) {
    inputStream >> tmp;
    e1[i] = tmp;
  }
  v1[n1] = m1;

  //in adj_hypergraph, <0> is hypernodes, <1> is hyperedges
  biadjacency<1> N(n0, m0);
  N.move(std::move(v0), std::move(e0));
  biadjacency<0> E(n1, m1);
  E.move(std::move(v1), std::move(e1));
  return std::tuple(E, N);
}

/* 
* Fill weighted biadjacency.
* Force weights to be vertex_id_t type.
* return adjacency<0, vertex_id_t> and adjacency<1, vertex_id_t>.
*/
template <std::unsigned_integral vertex_id_t, typename T>
auto weighted_adj_hypergraph_fill(std::istream& inputStream) {
  //100, 972, 100, 972
  //the number of weights are the same as the number of half edges (offsets)
  size_t n0, m0, n1, m1;
  inputStream >> n0;
  inputStream >> m0;
  inputStream >> n1;
  inputStream >> m1;

  vertex_id_t tmp;
  std::vector<vertex_id_t> v0(n0 + 1);
  for (size_t i = 0; i < n0; ++i) {
    inputStream >> tmp;
    v0[i] = tmp;
  }
  std::vector<vertex_id_t> e0(m0);
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    e0[i] = tmp;
  }
  v0[n0] = m0;
  std::vector<T> w0(m0);
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    w0[i] = tmp;
  }

  std::vector<vertex_id_t> v1(n1 + 1);
  for (size_t i = 0; i < n1; ++i) {
    inputStream >> tmp;
    v1[i] = tmp;
  }
  std::vector<vertex_id_t> e1(m1);
  for (size_t i = 0; i < m1; ++i) {
    inputStream >> tmp;
    e1[i] = tmp;
  }
  v1[n1] = m1;
  std::vector<T> w1(m1);
  for (size_t i = 0; i < m1; ++i) {
    inputStream >> tmp;
    w1[i] = tmp;
  }

  //in adj_hypergraph, <0> is hypernodes, <1> is hyperedges
  biadjacency<1, T> N(n0, m0);
  N.move(std::move(v0), std::move(e0), std::move(w0));
  biadjacency<0, T> E(n1, m1);
  E.move(std::move(v1), std::move(e1), std::move(w1));
  return std::tuple(E, N);
}

/*
* combine biadjacency into adjacency. Return g and its transpose
* Transform a bipartite graph into a general graph
*/
template <std::unsigned_integral vertex_id_t>
auto adj_hypergraph_pair_fill_and_adjoin(std::istream& inputStream,
const size_t n0, const size_t m0, const size_t n1, const size_t m1) {
  vertex_id_t tmp;
  std::string buffer;
  size_t N = n0 + n1;
  size_t M = m0 + m1;
  std::vector<vertex_id_t> v0(N + 1), v1(N + 1);
  std::vector<vertex_id_t> e0(M), e1(M);
  for (size_t i = 0; i < n0; ++i)
  {
    inputStream >> tmp;
    v0[i] = tmp;
    v1[i] = tmp + m1;
  }
  for (size_t i = 0; i < m0; ++i)
  {
    inputStream >> tmp;
    //increment each neighbor of hypernodes
    //by n0 (num of hypernodes)
    e0[i] = tmp + n0;
    e1[i] = tmp;
  }
  for (size_t i = n0, e = N; i < e; ++i)
  {
    inputStream >> tmp;
    //increment each index of hyperedges
    //by m0 (the last index of first adjacency)
    v0[i] = tmp + m0;
    v1[i] = tmp;
  }
  for (size_t i = m0, e = M; i < e; ++i)
  {
    inputStream >> tmp;
    e0[i] = tmp;
    e1[i] = tmp + n1;
  }
  v0[N] = M;
  v1[N] = M;
  //create use move constructor
  // adjacency<0> is hyperedges, adjacency<1> is hypernodes
  return std::tuple(adjacency<0>(std::move(v1), std::move(e1)),
  adjacency<1>(std::move(v0), std::move(e0)));
}
/*
* combine biadjacency into adjacency. Return g and its transpose
* Transform a bipartite graph into a general graph
*/
template <std::unsigned_integral vertex_id_t>
auto adj_hypergraph_pair_fill_and_adjoin_v1(std::istream& inputStream,
const size_t n0, const size_t m0, const size_t n1, const size_t m1) {
  vertex_id_t tmp;
  size_t N = n0 + n1;
  size_t M = m0 + m1;
  //from 0 to n0 - 1 are the hypernode csr indices
  //from next 0 to m0 - 1 are the hypernode incident lists
  //from next 0 to n1 - 1 are the hyperedge csr indices
  //from next 0 to m1 - 1 are the hyperedge incident lists
  //v0-e0 is graph
  //v1-e1 is transpose
  std::vector<vertex_id_t> v0(N + 1), v1(N + 1);
  std::vector<vertex_id_t> e0(M), e1(M);
  for (size_t i = 0, j = n1; i < n0; ++i, ++j) {
    inputStream >> tmp;
    v0[i] = tmp;
    //increment each index of hypernodes
    //by m0 (the last index of first adjacency)
    v1[j] = tmp + m0;
  }
  for (size_t i = 0, j = m1; i < m0; ++i, ++j) {
    inputStream >> tmp;
    //increment each neighbor of hypernodes
    //by n0 (num of hypernodes)
    e0[i] = tmp + n0;
    e1[j] = tmp;
  }
  for (size_t i = n0, j = 0, e = N; i < e; ++i, ++j) {
    inputStream >> tmp;
    //increment each index of hyperedges
    //by m0 (the last index of first adjacency)
    v0[i] = tmp + m0;
    v1[j] = tmp;
  }
  for (size_t i = m0, j = 0, e = M; i < e; ++i, ++j) {
    inputStream >> tmp;
    e0[i] = tmp;
    //increment each neighbor of hypernodes
    //by n1 (num of hyperedges)
    e1[j] = tmp + n0;
  }
  v0[N] = M;
  v1[N] = M;
  //create use move constructor
  // adjacency<0> is hyperedges, adjacency<1> is hypernodes
  return std::tuple(adjacency<0>(std::move(v1), std::move(e1)),
  adjacency<1>(std::move(v0), std::move(e0)));
}

/*
* combine biadjacency into an adjacency
* Transform a bipartite graph into a general graph
*/
template <std::unsigned_integral vertex_id_t>
auto adjacency_fill_adjoin(std::istream& inputStream,
const size_t n0, const size_t m0, const size_t n1, const size_t m1) {
  vertex_id_t tmp;
  size_t N = n0 + n1;
  size_t M = m0 + m1;
  std::vector<vertex_id_t> v0(N + 1);
  std::vector<vertex_id_t> e0(M);
  if (n0 >= n1) {
    //if the first adjacency (hypernode) is larger the the second adjacency (hyperedge)
    //we relabel the latter (increment each hyperedge by n0) and leave indices of first untouched
    //but we need to increment the neighbors of first by n0
    //here is what v0 contains: v0={0, ..., n0-1, n0,..., n0+n1-2, n0+n1-1, n0+n1}
    //where v0[n0+n1]=m0+m1, which is populated right before adjacency construction
    //here is what e0 contains: e0={0, ..., m0-1, m0,..., m0+m1-1}
    for (size_t i = 0; i < n0; ++i) {
      //read in hypernode indices
      inputStream >> tmp;
      v0[i] = tmp;
    }
    for (size_t i = 0; i < m0; ++i) {
      inputStream >> tmp;
      //increment each neighbor of hypernodes
      //by n0 (num of hypernodes)
      e0[i] = tmp + n0;
    }
    for (size_t i = n0, e = N; i < e; ++i) {
      inputStream >> tmp;
      //increment each index of hyperedges
      //by m0 (the last index of first adjacency)
      v0[i] = tmp + m0;
    }
    for (size_t i = m0, e = M; i < e; ++i) {
      inputStream >> tmp;
      e0[i] = tmp;
    }
  }
  else {
    //if the first adjacency is smaller the the second adjacency
    //we relabel the first
    //here is what v0 contains: v0={0, ..., n1-1, n1,..., n1+n0-1, n1+n0}
    //where v0[n0+n1]=m0+m1, which is populated later before return
    //here is what e0 contains: e0={0, ..., m1-1, m1,..., m1+m0-1}
    for (size_t i = n1, e = N; i < e; ++i) {
      inputStream >> tmp;
      //increment the first index by
      //the last offset of the second adjacency
      v0[i] = tmp + m1;
    }
    for (size_t i = m1, e = M; i < e; ++i) {
      inputStream >> tmp;
      //the offset of the first index remains the same
      e0[i] = tmp;
    }
    for (size_t i = 0, e = n1; i < e; ++i) {
      inputStream >> tmp;
      //the second index remains the same
      v0[i] = tmp;
    }
    for (size_t i = 0, e = m1; i < e; ++i) {
      inputStream >> tmp;
      //increment each neighbor of the second adjacency
      //by n1 the num of indices of first adjacency
      e0[i] = tmp + n1;
    }
  }
  v0[N] = M;
  return std::tuple(v0, e0);
}

/* 
* Fill edgelist.
* Convention is that column 0 stores hyperedges,
* column 1 stores hypernodes.
*/
template<class Edgelist>
void adjacency_fill(std::istream& inputStream, Edgelist& A,
const size_t nnodes, const size_t m0, const size_t nedges, const size_t m1) {
  using vertex_id_t = typename graph_traits<Edgelist>::vertex_id_type;

  A.open_for_push_back();
  vertex_id_t tmp;
  //skip the node-edge adjacency
  std::vector<vertex_id_t> v0(nnodes + 1);
  for (size_t i = 0; i < nnodes; ++i) {
    inputStream >> tmp;
    v0[i] = tmp;
  }
  std::vector<vertex_id_t> e0(m0);
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    e0[i] = tmp;
  }
  v0[nnodes] = m0;
  /*
  for (vertex_id_t v = 0; v < nnodes; ++v) {
    for (size_t j = v0[v]; j < v0[v + 1]; ++j) {
      A.push_back(e0[j], v);
    }
  }
*/
  //populate the edge-node adjacency
  std::vector<vertex_id_t> v1(nedges + 1);
  for (size_t i = 0; i < nedges; ++i) {
    inputStream >> tmp;
    v1[i] = tmp;
  }
  std::vector<vertex_id_t> e1(m1);
  for (size_t i = 0; i < m1; ++i) {
    inputStream >> tmp;
    e1[i] = tmp;
  }
  v1[nedges] = m1;
  for (vertex_id_t e = 0; e < nedges; ++e) {
    for (size_t j = v1[e]; j < v1[e + 1]; ++j) {
      A.push_back(e, e1[j]);
    }
  }
  A.close_for_push_back();
}
/* 
* Fill edgelist.
* Convention is that column 0 stores hyperedges,
* column 1 stores hypernodes.
*/
template<class Edgelist, class T>
void weighted_adjacency_fill(std::istream& inputStream, Edgelist& A,
const size_t n0, const size_t m0, const size_t n1, const size_t m1) {
  using vertex_id_t = typename graph_traits<Edgelist>::vertex_id_type;

  A.open_for_push_back();
  vertex_id_t tmp;
  //skip the node-edge adjacency indices
  for (size_t i = 0; i < n0; ++i) {
    inputStream >> tmp;
  }
  //skip the node-edge adjacency offsets
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
  }
  //skip the node-edge adjacency weights
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
  }

  //populate the edge-node adjacency indices
  std::vector<vertex_id_t> v1(n1 + 1);
  for (size_t i = 0; i < n1; ++i) {
    inputStream >> tmp;
    v1[i] = tmp;
  }
  //populate the edge-node adjacency offsets
  std::vector<vertex_id_t> e1(m1);
  for (size_t i = 0; i < m1; ++i) {
    inputStream >> tmp;
    e1[i] = tmp;
  }
  v1[n1] = m1;
  T weight;
  //read the edge-node adjacency weights and insert an edge pair
  for (vertex_id_t e = 0; e < n1; ++e) {
    for (size_t j = v1[e]; j < v1[e + 1]; ++j) {
      inputStream >> weight;
      A.push_back(e, e1[e], weight);
    }
  }
  A.close_for_push_back();
}

/* 
* Fill edgelist, where each edge is from hyperedge to hypernode.
* Convention is that column 0 are hyperedges,
* column 1 are hypernodes.
*/
template<class Edgelist>
void edgelist_fill_adjoin(std::istream& inputStream, Edgelist& A,
const size_t nnodes, const size_t m0, const size_t nedges, const size_t m1) {
  using vertex_id_t = typename graph_traits<Edgelist>::vertex_id_type;

  A.open_for_push_back();
  vertex_id_t tmp;
  //skip the node-edge adjacency
  std::vector<vertex_id_t> v0(nnodes + 1);
  for (size_t i = 0; i < nnodes; ++i) {
    inputStream >> tmp;
    v0[i] = tmp;
  }
  std::vector<vertex_id_t> e0(m0);
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    e0[i] = tmp;
  }
  v0[nnodes] = m0;
  /*
  if (nnodes >= nedges) {
    for (vertex_id_t v = 0; v < nnodes; ++v) {
      for (size_t j = v0[v]; j < v0[v + 1]; ++j) {
        A.push_back(e0[j] + nnodes, v);
      }
    }
  } else {
    for (vertex_id_t v = 0; v < nnodes; ++v) {
      for (size_t j = v0[v]; j < v0[v + 1]; ++j) {
        A.push_back(e0[j], v + nedges);
      }
    }
  }
*/
  //populate the edge-node adjacency
  std::vector<vertex_id_t> v1(nedges + 1);
  for (size_t i = 0; i < nedges; ++i) {
    inputStream >> tmp;
    v1[i] = tmp;
  }
  std::vector<vertex_id_t> e1(m1);
  for (size_t i = 0; i < m1; ++i) {
    inputStream >> tmp;
    e1[i] = tmp;
  }
  v1[nedges] = m1;
  if (nnodes >= nedges) {
    // shift hyperedge id by nnodes
    for (vertex_id_t e = 0; e < nedges; ++e) {
      for (size_t j = v1[e]; j < v1[e + 1]; ++j) {
        A.push_back(e + nnodes, e1[j]);
      }
    }
  } else {
    // shift hypernode id by nedges
    for (vertex_id_t e = 0; e < nedges; ++e) {
      for (size_t j = v1[e]; j < v1[e + 1]; ++j) {
        A.push_back(e, e1[j] + nedges);
      }
    }
  }
  A.close_for_push_back();
}

template <std::unsigned_integral vertex_id_t>
auto read_and_adjoin_adj_hypergraph_pair(const std::string& filename, size_t& nreal_edges, size_t& nreal_nodes) {
  std::ifstream inputStream(filename, std::ifstream::in);
  if (!inputStream.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  std::string buffer;
  inputStream >> buffer;

  if (AdjHypergraphHeader.c_str() != buffer) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  //100, 972 100 972
  size_t m0, m1;
  inputStream >> nreal_nodes;
  inputStream >> m0;
  inputStream >> nreal_edges;
  inputStream >> m1;
  //std::cout << nreal_nodes << " " << m0 << " " << nreal_edges << " " << m1 << std::endl;

  auto&& [v0, e0] = adjacency_fill_adjoin<vertex_id_t>(inputStream, nreal_nodes, m0, nreal_edges, m1);
  return std::tuple(adjacency<0>(std::move(v0), std::move(e0)), adjacency<1>(std::move(v0), std::move(e0))); 
}

template <std::unsigned_integral vertex_id_t>
auto read_and_adjoin_adj_hypergraph(const std::string& filename, size_t& nreal_edges, size_t& nreal_nodes) {
  std::ifstream inputStream(filename, std::ifstream::in);
  if (!inputStream.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  std::string buffer;
  inputStream >> buffer;

  if (AdjHypergraphHeader.c_str() != buffer) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  //100, 972 100 972
  size_t m0, m1;
  inputStream >> nreal_nodes;
  inputStream >> m0;
  inputStream >> nreal_edges;
  inputStream >> m1;
  //std::cout << nreal_nodes << " " << m0 << " " << nreal_edges << " " << m1 << std::endl;
  auto&& [v0, e0] = adjacency_fill_adjoin<vertex_id_t>(inputStream, nreal_nodes, m0, nreal_edges, m1);
  return adjacency<0>(std::move(v0), std::move(e0));
}

/*
 * Read unweighted adjacency hypergraph as bi-adjacency
 **/
template <std::unsigned_integral vertex_id_t, typename... Attributes>
auto read_adj_hypergraph(const std::string& filename) {
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  std::string header;
  file >> header;

  if (AdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  return adj_hypergraph_fill<vertex_id_t>(file);
}

/*
 * Read the adjacency graph as an edgelist.
 **/
template<nw::graph::directedness sym, typename... Attributes>
nw::graph::bi_edge_list<sym, Attributes...> read_adjacency(const std::string& filename) {
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  std::string header;
  file >> header;

  if (AdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported adjacency format" << std::endl;
    throw;
  }
  size_t nnodes, m0, nedges, m1;
  file >> nnodes;
  file >> m0;
  file >> nedges;
  file >> m1;

  nw::graph::bi_edge_list<sym, Attributes...> A(nedges, nnodes);
  A.reserve(m1);
  adjacency_fill(file, A, nnodes, m0, nedges, m1);

  return A;
}

/*
 * Read the adjacency graph as an edgelist.
 **/
template<nw::graph::directedness sym, typename... Attributes>
nw::graph::bi_edge_list<sym, Attributes...> read_weighted_adjacency(const std::string& filename) {
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  std::string header;
  file >> header;

  if (WghAdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported adjacency format" << std::endl;
    throw;
  }
  size_t nnodes, m0, nedges, m1;
  file >> nnodes;
  file >> m0;
  file >> nedges;
  file >> m1;

  nw::graph::bi_edge_list<sym, Attributes...> A(nedges, nnodes);
  A.reserve(m1);
  (weighted_adjacency_fill<nw::graph::bi_edge_list<sym, Attributes...>, Attributes>(file, A, nnodes, m0, nedges, m1), ...);

  return A;
}
/*
 * Read the adjacency hypergraph into adjoin graph in the format of an edgelist.
 **/
template<nw::graph::directedness sym, typename... Attributes>
nw::graph::edge_list<sym, Attributes...> read_adjacency_adjoin(const std::string& filename, size_t& nreal_edges, size_t& nreal_nodes) {
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  std::string header;
  file >> header;

  if (AdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported adjacency format" << std::endl;
    throw;
  }

  size_t m0, m1;
  file >> nreal_nodes;
  file >> m0;
  file >> nreal_edges;
  file >> m1;

  nw::graph::edge_list<sym, Attributes...> A(nreal_edges + nreal_nodes);
  A.reserve(m1);
  edgelist_fill_adjoin<nw::graph::edge_list<sym, Attributes...>>(file, A, nreal_nodes, m0, nreal_edges, m1);

  return A;
}

template <std::unsigned_integral vertex_id_t, typename... Attributes>
auto read_weighted_adj_hypergraph(const std::string& filename) {
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  std::string header;
  file >> header;

  if (WghAdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  return (weighted_adj_hypergraph_fill<vertex_id_t, Attributes>(file), ...);
}

auto write_adj_hypergraph(const std::string& filename, biadjacency<0>& E, biadjacency<1>& N) {
  using vertex_id_t = typename graph_traits<decltype(E)>::vertex_id_type; 
  nw::util::life_timer _(__func__);
  std::ofstream file(filename, std::ios::out | std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    return false;
  }
  std::vector<vertex_id_t> v0(N.get_indices());
  v0.pop_back();
  std::vector<vertex_id_t> e0(std::get<0>(N.get_to_be_indexed()));
  std::vector<vertex_id_t> v1(E.get_indices());
  v1.pop_back();
  std::vector<vertex_id_t> e1(std::get<0>(E.get_to_be_indexed()));
  //write header at line 1
  file << AdjHypergraphHeader.c_str() << std::endl;
  size_t n0 = v0.size();
  size_t m0 = v0[n0];
  size_t n1 = v1.size();
  size_t m1 = v1[n1];
  std::cout << n0 << " " << m0 << " " << n1 << " " << m1 << std::endl;
  //write the stats of biadjacency followed by header
  file << n0 << std::endl;
  file << m0 << std::endl;
  file << n1 << std::endl;
  file << m1 << std::endl;


  for (size_t i = 0; i < n0; ++i) {
    file << v0[i] << std::endl;
  }
  
  for (size_t i = 0; i < m0; ++i) {
    file << e0[i] << std::endl;
  }

  for (size_t i = 0; i < n1; ++i) {
    file << v1[i] << std::endl;
  }
  
  for (size_t i = 0; i < m1; ++i) {
    file << e1[i] << std::endl;
  }
  return true;
}
}//namespace hypergraph
}//namespace nw