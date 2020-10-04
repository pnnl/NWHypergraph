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

using namespace nw::graph;

namespace nw {
namespace hypergraph {

std::string AdjHypergraphHeader = "AdjacencyHypergraph";
std::string WghAdjHypergraphHeader = "WeightedAdjacencyHypergraph";


void adj_hypergraph_fill(std::istream& inputStream, adjacency<0>& E, adjacency<1>& N,
size_t n0, size_t m0, size_t n1, size_t m1) {
  vertex_id_t tmp;
  std::vector<vertex_id_t> v0(n0);
  for (size_t i = 0; i < n0; ++i) {
    inputStream >> tmp;
    v0[i] = tmp;
  }
  std::vector<vertex_id_t> e0(m0);
  for (size_t i = 0; i < m0; ++i) {
    inputStream >> tmp;
    e0[i] = tmp;
  }
  //N.copy(v0, e0);
  N.move(std::move(v0), std::move(e0));
  std::vector<vertex_id_t> v1(n1);
  for (size_t i = 0; i < n1; ++i) {
    inputStream >> tmp;
    v1[i] = tmp;
  }
  std::vector<vertex_id_t> e1(m1);
  for (size_t i = 0; i < m1; ++i) {
    inputStream >> tmp;
    e1[i] = tmp;
  }
  //E.copy(v1, e1);
  E.move(std::move(v1), std::move(e1));
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