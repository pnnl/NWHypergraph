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

#include <nwgraph/edge_list.hpp>
#include <nwgraph/adjacency.hpp>
#include <nwgraph/util/timer.hpp>
using namespace nw::graph;

namespace nw {
namespace hypergraph {

/*
* Read the csv file row by row, split the numbers by comma and insert to edge list
**/
template <class EdgeList>
void csv_fill(std::istream& inputStream, EdgeList& A) {
  using vertex_id_t = typename graph_traits<EdgeList>::vertex_id_type;
  std::string line, word, temp;
  vertex_id_t e = 0;
  A.open_for_push_back();
  while (std::getline(inputStream, line)) {
    // read an entire row and
    // store it in a string variable 'line'

    // used for breaking words
    std::stringstream s(line);

    // read every column data of a row and
    // store it in a string variable, 'word'
    while (std::getline(s, word, ',')) {
      // add all the column data
      // of a row to a vector
      try {
        // to validate whether we are handling a standard csv file
        vertex_id_t v = std::stoi(word) - 1;
        A.push_back(e, v);
      } catch (const std::invalid_argument& ia) {
        //if it is not a valid csv file (there is a header in it),
        //std::stoi will throw an invalid_argument exception
        //std::cerr << "Invalid argument: " << ia.what() << std::endl;
        std::cerr << "Not a csv file" << std::endl;
        return;
      }
    }
    ++e;
  }
  A.close_for_push_back();
}

/*
* Convert (adjoin) a biparite edge list into a unipartite edge list
**/
template <class BiEdgeList, class EdgeList>
void populate_adjoin_edge_list(BiEdgeList& A, EdgeList& B,
                               const size_t nedges, const size_t nnodes) {
  B.open_for_push_back();
  if (nedges > nnodes) {
    for (auto&& [d0, d1] : A) {
      B.push_back(d0, d1 + nedges);
    }
  } else {
    for (auto&& [d0, d1] : A) {
      B.push_back(d0 + nnodes, d1);
    }
  }
  B.close_for_push_back();
}

/*
* Read the csv file as an edge list.
* Assuming csv file is index-1 based.
**/
template <directedness sym, typename... Attributes>
nw::graph::bi_edge_list<sym, Attributes...> read_csv(
    const std::string& filename) {
  std::ifstream file(filename, std::ifstream::in);
  // test file open
  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    throw;
  }
  nw::graph::bi_edge_list<sym, Attributes...> A(0, 0);
  // Read the data from the file as edge_list
  csv_fill(file, A);
  file.close();
  return A;
}

/*
* Read the csv file as an edge list.
* Assuming csv file is index-1 based.
**/
template <directedness sym, typename... Attributes>
nw::graph::edge_list<sym, Attributes...> read_csv_adjoin(const std::string& filename, size_t& numRealEdges, size_t& numRealNodes) {
    std::ifstream file(filename, std::ifstream::in);
    // test file open
    if (!file.is_open()) {
      std::cerr << "Can not open file: " << filename << std::endl;
      throw;
    }
    nw::graph::bi_edge_list<sym, Attributes...> A;
    // Read the data from the file as edge_list
    csv_fill(file, A);
    file.close();
    //adjoin the edge list
    numRealEdges = A.num_vertices()[0];
    numRealNodes = A.num_vertices()[1];
    nw::graph::edge_list<sym, Attributes...> B(numRealEdges + numRealNodes);
    B.reserve(A.size());
    populate_adjoin_edge_list(A, B, numRealEdges, numRealNodes);
    return B;
}

}//namespace hypergraph
}//namespace nw