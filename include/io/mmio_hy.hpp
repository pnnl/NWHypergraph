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
#include <nwgraph/util/timer.hpp>
#include <nwgraph/io/mmio.hpp>

namespace nw {
namespace hypergraph {

template<class EdgeList>
void mm_fill_adjoin(std::istream& inputStream, 
EdgeList& A, 
size_t nedges, size_t nnodes, 
size_t nNonzeros, bool file_symmetry, bool pattern) {
  //A.reserve(nNonzeros);
  A.reserve((file_symmetry ? 2 : 1) * nNonzeros);
  A.open_for_push_back();
  if (nedges > nnodes) {
    if (pattern) {
      for (size_t i = 0; i < nNonzeros; ++i) {
        size_t d0, d1;
        inputStream >> d0 >> d1;
        
        A.push_back(d0 - 1, d1 + nedges - 1);
        if (file_symmetry && (d0 != d1)) {
          A.push_back(d1 + nedges - 1, d0 - 1);
        }
      }
    }
    else {
      for (size_t i = 0; i < nNonzeros; ++i) {
        size_t d0, d1;
        double d2;
        inputStream >> d0 >> d1 >> d2;

        A.push_back(d0 - 1, d1 + nedges - 1);
        if (file_symmetry && (d0 != d1)) {
          A.push_back(d1 + nedges - 1, d0 - 1);
        }
      }
    }
  }
  else {
    if (pattern) {
      for (size_t i = 0; i < nNonzeros; ++i) {
        size_t d0, d1;
        inputStream >> d0 >> d1;
        A.push_back(d0 + nnodes - 1, d1 - 1);
        if (file_symmetry && (d0 != d1)) {
          A.push_back(d1 - 1, d0 + nnodes - 1);
        }
      }
    }
    else {
      for (size_t i = 0; i < nNonzeros; ++i) {
        size_t d0, d1;
        double d2;
        inputStream >> d0 >> d1 >> d2;

        A.push_back(d0 + nnodes - 1, d1 - 1);
        if (file_symmetry && (d0 != d1)) {
          A.push_back(d1 - 1, d0 + nnodes - 1);
        }
      }
    }
  }
  A.close_for_push_back();
}

template<typename EdgeList, typename T>
void mm_fill_adjoin(std::istream& inputStream, 
EdgeList& A, 
size_t nedges, size_t nnodes, 
size_t nNonzeros, bool file_symmetry, bool pattern) {

  A.reserve((file_symmetry ? 2 : 1) * nNonzeros);
  A.open_for_push_back();
  if (nedges > nnodes) {
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

      d1 += nedges;

      A.push_back(d0 - 1, d1 - 1, v);

      if (file_symmetry && (d0 != d1)) {
        A.push_back(d1 - 1, d0 - 1, v);
      }
    }
  }
  else {
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

      d0 += nnodes;

      A.push_back(d0 - 1, d1 - 1, v);

      if (file_symmetry && (d0 != d1)) {
        A.push_back(d1 - 1, d0 - 1, v);
      }
    }    
  }
  A.close_for_push_back();
}

template <directedness sym, typename... Attributes>
nw::graph::bi_edge_list<sym, Attributes...> read_mm(const std::string& filename) {
  std::ifstream            inputStream(filename);
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
  size_t n0, n1, nNonzeros;
  std::stringstream(string_input) >> n0 >> n1 >> nNonzeros;

  // bipartite edge list
  if (file_symmetry) {
    std::cerr << "Can not populate bipartite graph with symmetric matrix"
              << std::endl;
    throw;
  }
  nw::graph::bi_edge_list<sym, Attributes...> A(n0, n1);
  mm_fill(inputStream, A, nNonzeros, file_symmetry, (header[3] == "pattern"));

  return A;
}

//loader for mmio_adjoin
template<nw::graph::directedness sym, typename... Attributes>
nw::graph::edge_list<sym, Attributes...> read_mm_adjoin(const std::string& filename, size_t& numRealEdges, size_t& numRealNodes) {
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

  nw::graph::edge_list<sym, Attributes...> A(numRealEdges + numRealNodes);
  mm_fill_adjoin<nw::graph::edge_list<sym, Attributes...>, Attributes...>(inputStream, A, numRealEdges, numRealNodes, nNonzeros, file_symmetry, (header[3] == "pattern"));
  
  return A;
}

template <size_t w_idx, int idx, typename... Attributes>
void hy_adjacency_stream(std::ofstream& outputStream, nw::graph::biadjacency<idx, Attributes...>& A, 
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
void write_mm_hy(const std::string& filename, nw::graph::biadjacency<idx, Attributes...>& A, size_t size_col0, size_t size_col1,
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

}//namespace hypergraph
}//namespace nw