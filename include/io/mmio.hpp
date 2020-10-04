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

template<nw::graph::directedness sym, typename... Attributes>
nw::graph::edge_list<sym, Attributes...> read_mm_relabeling(std::istream& inputStream, size_t& numRealEdges, size_t& numRealNodes) {
  std::string              string_input;
  bool                     file_symmetry = false;
  std::vector<std::string> header(5);

  // %%MatrixMarket matrix coordinate integer symmetric
  std::getline(inputStream, string_input);
  std::stringstream h(string_input);
  for (auto& s : header)
    h >> s;

  if (header[0] != "%%MatrixMarket") {
    std::cerr << "Unsupported format" << std::endl;
    throw;
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

  return A;
}

//loader for mmio
template<nw::graph::directedness sym, typename... Attributes>
nw::graph::edge_list<sym, Attributes...> read_mm_relabeling(const std::string& filename, size_t& numRealEdges, size_t& numRealNodes) {
  std::ifstream inputFile(filename);
  std::string type;
  inputFile >> type;

  if ("%%MatrixMarket" == type) {
    std::cout << "Reading matrix market input " << filename << " (slow)" << std::endl;
    nw::util::life_timer _("read and relable mm");
    nw::graph::edge_list<sym, Attributes...> A = read_mm_relabeling<sym, Attributes...>(inputFile, numRealEdges, numRealNodes);
    A.set_origin(filename);

    return A;
  }
  else {
    //std::cerr << "Did not recognize graph input file " << file << "\n";
    return nw::graph::edge_list<sym, Attributes...>(0);
  }
}
