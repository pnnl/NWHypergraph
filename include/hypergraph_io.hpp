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

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <fcntl.h>
#include <iostream>
#include <sstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <edge_list.hpp>
using namespace nw::graph;

namespace nw {
namespace hypergraph {

std::string AdjHypergraphHeader = "AdjacencyHypergraph";
std::string WghAdjHypergraphHeader = "WeightedAdjacencyHypergraph";

template <directedness sym, typename... Attributes>
edge_list<sym, Attributes...> read_adj_hypergraph(const std::string& filename) {
  std::ifstream file(filename);

  edge_list<sym, Attributes...> A = read_adj_hypergraph<sym, Attributes...>(filename);
  A.set_origin(filename);
  
  return A;
}

template <typename T>
void adj_hypergraph_fill(std::istream& inputStream, edge_list<undirected, T>& A,
size_t n0, size_t m0, size_t n1, size_t m1) {
  A.open_for_push_back();
  std::vector<size_t> v0(n0);
  for (size_t i = 0; i < n0; ++i) {
    std::string buffer;
    size_t      tmp;
    std::getline(inputStream, buffer);
    std::stringstream(buffer) >> tmp;
    v0[i] = tmp;
  }
  std::vector<size_t> v1(n1);
  for (size_t i = 0; i < n1; ++i) {
    std::string buffer;
    size_t      tmp;
    std::getline(inputStream, buffer);
    std::stringstream(buffer) >> tmp;
    v1[i] = tmp;
  }



    
    //A.push_back(d0, d1);


  A.close_for_push_back();
}

template <directedness sym, typename... Attributes>
edge_list<sym, Attributes...> read_adj_hypergraph(std::istream& inputStream) {
  std::string header;
  inputStream >> header;

  if (header != AdjHypergraphHeader) {
    std::cerr << "Unsupported format" << std::endl;
    throw;
  }
  size_t n0, m0, n1, m1;
  std::string string_input;
  std::stringstream(string_input) >> n0;
  std::stringstream(string_input) >> m0;
  std::stringstream(string_input) >> n1;
  std::stringstream(string_input) >> m1;
  std::cout << n0 << " " << m0 << std::endl;
  std::cout << n1 << " " << m1 << std::endl;
  // assert(n0 == n1);

  edge_list<sym, Attributes...> A(0);
  adj_hypergraph_fill(inputStream, A, n0, m0, n1, m1);
  return A;
}

template <directedness sym, typename... Attributes>
edge_list<sym, Attributes...> read_weighted_adj_hypergraph(const std::string& filename) {


  edge_list<sym, Attributes...> A(0);
  
  return A;
}


}//namespace hypergraph
}//namespace nw