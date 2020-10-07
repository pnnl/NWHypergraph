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
#include <fstream>

#include <edge_list.hpp>
#include <compressed.hpp>
#include <util/timer.hpp>
using namespace nw::graph;

namespace nw {
namespace hypergraph {

std::string AdjHypergraphHeader = "AdjacencyHypergraph";
std::string WghAdjHypergraphHeader = "WeightedAdjacencyHypergraph";


void adj_hypergraph_fill(std::istream& inputStream, adjacency<0>& E, adjacency<1>& N,
size_t n0, size_t m0, size_t n1, size_t m1) {
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
  v0[n0 + 1] = m0;
  N.move(std::move(v0), std::move(e0));
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
  v1[n1 + 1] = m1;
  E.move(std::move(v1), std::move(e1));
}

auto adj_hypergraph_fill_and_relabel(std::istream& inputStream,
size_t n0, size_t m0, size_t n1, size_t m1) {
  vertex_id_t tmp;
  std::vector<vertex_id_t> v0(n0 + n1 + 1);
  std::vector<vertex_id_t> e0(m0 + m1);
  if (n0 > n1) {
    //if the first adjacency is larger the the second adjacency
    //we relabel the latter
  for (size_t i = 0; i < n0; ++i) {
    inputStream >> tmp;
    v0[i] = tmp;
  }
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    e0[i] = tmp + n0;
  }
  for (size_t i = n0, e = n0 + n1; i < e; ++i) {
    inputStream >> tmp;
    v0[i] = tmp + m0;
  }
  for (size_t i = 0, e = m0 + m1; i < e; ++i) {
    inputStream >> tmp;
    e0[i] = tmp;
  }
  }
  else {
    //if the first adjacency is smaller the the second adjacency
    //we relabel the first
  for (size_t i = 0; i < n0; ++i) {
    inputStream >> tmp;
    v0[i] = tmp + m1;
  }
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    e0[i] = tmp;
  }
  for (size_t i = n0, e = n0 + n1; i < e; ++i) {
    inputStream >> tmp;
    v0[i] = tmp;
  }
  for (size_t i = 0, e = m0 + m1; i < e; ++i) {
    inputStream >> tmp;
    e0[i] = tmp + n1;
  }
  }
  v0[n0 + n1 + 1] = m0 + m1;
  return adjacency<0>(std::move(v0), std::move(e0));
}

auto read_and_relabel_adj_hypergraph(std::istream& inputStream, size_t& nreal_edges, size_t& nreal_nodes) {
  std::string header;
  inputStream >> header;

  if (AdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  //100, 972 100 972
  size_t m0, m1;
  inputStream >> nreal_nodes;
  inputStream >> m0;
  inputStream >> nreal_edges;
  inputStream >> m1;

  return adj_hypergraph_fill_and_relabel(inputStream, nreal_nodes, m0, nreal_edges, m1);
}

auto read_and_relabel_adj_hypergraph(const std::string& filename, size_t& nreal_edges, size_t& nreal_nodes) {
  nw::util::life_timer _(__func__);
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  return read_and_relabel_adj_hypergraph(file, nreal_edges, nreal_nodes);
}

template <typename... Attributes>
auto read_adj_hypergraph(std::istream& inputStream) {
  std::string header;
  inputStream >> header;

  if (AdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  //100, 972 100 972
  size_t n0, m0, n1, m1;
  inputStream >> n0;
  inputStream >> m0;
  inputStream >> n1;
  inputStream >> m1;

  //in adj_hypergraph, <0> is hypernodes, <1> is hyperedges
  adjacency<0> E(n1, m1);
  adjacency<1> N(n0, m0);
  adj_hypergraph_fill(inputStream, E, N, n0, m0, n1, m1);
  return std::tuple(E, N);
}

template <typename... Attributes>
auto read_adj_hypergraph(const std::string& filename) {
  nw::util::life_timer _(__func__);
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  return read_adj_hypergraph<Attributes...>(file);
}

template <typename... Attributes>
auto read_weighted_adj_hypergraph(std::istream& inputStream) {
  std::string header;
  inputStream >> header;

  if (WghAdjHypergraphHeader.c_str() != header) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  size_t n0, m0, w0, n1, m1, w1;
  inputStream >> n0;
  inputStream >> m0;
  inputStream >> w0;
  inputStream >> n1;
  inputStream >> m1;
  inputStream >> w1;
  std::cout << n0 << " " << m0 << " " << w0 << std::endl;
  std::cout << n1 << " " << m1 << " " << w1 << std::endl;
  // assert(n0 == n1);

  //in adj_hypergraph, <0> is hypernodes, <1> is hyperedges
  adjacency<0> E(n1, m1);
  adjacency<1> N(n0, m0);
  //TODO adj_hypergraph_fill(inputStream, E, N, n0, m0, n1, m1);
  return std::tuple(E, N);
}

template <typename... Attributes>
auto read_weighted_adj_hypergraph(const std::string& filename) {
  std::ifstream file(filename, std::ifstream::in);
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  return read_weighted_adj_hypergraph<Attributes...>(file);
}


}//namespace hypergraph
}//namespace nw