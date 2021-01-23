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
#include <util/timer.hpp>
#include <mmio.hpp>

void mm_fill_relabeling(std::istream& inputStream, nw::graph::edge_list<nw::graph::directed>& A, size_t nedges, size_t nnodes, 
size_t nNonzeros, bool file_symmetry, bool pattern) {

  A.reserve(nNonzeros);
  A.open_for_push_back();

  for (size_t i = 0; i < nNonzeros; ++i) {
    size_t d0, d1;
    double d2;

    if (pattern) {
      inputStream >> d0 >> d1;
    } else {
      inputStream >> d0 >> d1 >> d2;
    }
    if (nedges > nnodes)
      d1 += nedges;
    else
      d0 += nnodes;
    A.push_back(d0, d1);
  }
  A.close_for_push_back();
}

template<typename T>
void mm_fill_relabeling(std::istream& inputStream, nw::graph::edge_list<nw::graph::directed, T>& A, size_t nedges, size_t nnodes, 
size_t nNonzeros, bool file_symmetry, bool pattern) {

  A.reserve((file_symmetry ? 2 : 1) * nNonzeros);
  A.open_for_push_back();
  for (size_t i = 0; i < nNonzeros; ++i) {
    std::string buffer;
    size_t      d0, d1;
    T           v(1.0);

    std::getline(inputStream, buffer);
    if (pattern) {
      std::stringstream(buffer) >> d0 >> d1;
    } else {
      std::stringstream(buffer) >> d0 >> d1 >> v;
    }
    if (nedges > nnodes)
      d1 += nedges;
    else
      d0 += nnodes;
    A.push_back(d0, d1, v);

    if (file_symmetry && (d0 != d1)) {
      A.push_back(d1, d0, v);
    }
  }
  A.close_for_push_back();
}

void mm_fill_relabeling(std::istream& inputStream, nw::graph::edge_list<nw::graph::undirected>& A, size_t nedges, size_t nnodes, 
size_t nNonzeros, bool file_symmetry, bool pattern) {

  A.reserve(nNonzeros);
  A.open_for_push_back();

  for (size_t i = 0; i < nNonzeros; ++i) {
    size_t d0, d1;
    double d2;

    if (pattern) {
      inputStream >> d0 >> d1;
    } else {
      inputStream >> d0 >> d1 >> d2;
    }
    if (nedges > nnodes)
      d1 += nedges;
    else
      d0 += nnodes;
    A.push_back(d0, d1);
  }
  A.close_for_push_back();
}

template<typename T>
void mm_fill_relabeling(std::istream& inputStream, nw::graph::edge_list<nw::graph::undirected, T>& A, size_t nedges, size_t nnodes, 
size_t nNonzeros, bool file_symmetry, bool pattern) {
  assert(file_symmetry);
  A.reserve(nNonzeros);
  A.open_for_push_back();
  for (size_t i = 0; i < nNonzeros; ++i) {
    std::string buffer;
    size_t      d0, d1;
    T           v(1.0);
    std::getline(inputStream, buffer);
    if (pattern) {
      std::stringstream(buffer) >> d0 >> d1;
    } else {
      std::stringstream(buffer) >> d0 >> d1 >> v;
    }
    if (nedges > nnodes)
      d1 += nedges;
    else
      d0 += nnodes;
    A.push_back(d0, d1, v);
  }
  A.close_for_push_back();
}

//loader for mmio
template<nw::graph::directedness sym, typename... Attributes>
nw::graph::edge_list<sym, Attributes...> read_mm_relabeling(const std::string& filename, size_t& numRealEdges, size_t& numRealNodes) {
  nw::util::life_timer _(__func__);
  std::ifstream inputStream(filename);
  std::string              string_input;
  bool                     file_symmetry = false;
  std::vector<std::string> header(5);

  // %%MatrixMarket matrix coordinate integer symmetric
  std::getline(inputStream, string_input);
  std::stringstream h(string_input);
  for (auto& s : header)
    h >> s;

  if (header[0] != "%%MatrixMarket") {
    std::cerr << "Not Matrix Market format" << std::endl;
    return nw::graph::edge_list<sym, Attributes...>(0);
  }
  if (header[4] == "symmetric") {
    file_symmetry = true;
  } else if (header[4] == "general") {
    file_symmetry = false;
  } else {
    std::cerr << "Bad format (symmetry): " << header[4] << std::endl;
    throw;
  }

  while (std::getline(inputStream, string_input)) {
    if (string_input[0] != '%') break;
  }
  size_t nNonzeros;
  std::stringstream(string_input) >> numRealEdges >> numRealNodes >> nNonzeros;

  nw::graph::edge_list<sym, Attributes...> A(numRealEdges);
  mm_fill_relabeling(inputStream, A, numRealEdges, numRealNodes, nNonzeros, file_symmetry, (header[3] == "pattern"));
  A.set_origin(filename);
  
  return A;
}

template <size_t w_idx, int idx, typename... Attributes>
void hy_adjacency_stream(std::ofstream& outputStream, nw::graph::adjacency<idx, Attributes...>& A, 
size_t size_col0, size_t size_col1,
const std::string& file_symmetry, std::string& w_type) {
  outputStream << "%%MatrixMarket matrix coordinate " << w_type << " " << file_symmetry << "\n%%\n";

  outputStream << size_col0 << " " << size_col1 << " "
               << std::accumulate(A.begin(), A.end(), 0, [&](int a, auto b) { return a + (int)(b.end() - b.begin()); })
               << std::endl;

  for (auto first = A.begin(); first != A.end(); ++first) {
    for (auto v = (*first).begin(); v != (*first).end(); ++v) {
      outputStream << first - A.begin() + (1) << " " << std::get<0>(*v) + (1);
      if (w_idx != 0) outputStream << " " << std::get<w_idx>(*v);
      outputStream << std::endl;
    }
  }
}

template <size_t w_idx = 0, typename idxtype = void, int idx, typename... Attributes>
void write_mm_hy(const std::string& filename, nw::graph::adjacency<idx, Attributes...>& A, size_t size_col0, size_t size_col1,
const std::string& file_symmetry = "general") {
  /*if (file_symmetry == "symmetric" && sym == directedness::directed) {
    std::cerr << "cannot save directed matrix as symmetric matrix market" << std::endl;
  }*/

  std::string w_type = "pattern";
  if (std::numeric_limits<idxtype>::is_integer)
    w_type = "integer";
  else if (std::is_floating_point<idxtype>::value)
    w_type = "real";

  std::ofstream outputStream(filename);
  hy_adjacency_stream<w_idx>(outputStream, A, size_col0, size_col1, file_symmetry, w_type);
}
