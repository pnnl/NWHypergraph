//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2021
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#pragma once
#include <nwgraph/edge_list.hpp>
#include <nwgraph/io/mmio.hpp>

#include "hypergraph_io.hpp"
#include "csv_io.hpp"
#include "mmio_hy.hpp"

using namespace nw::graph;

namespace nw {
namespace hypergraph {

/*
 * This loader loads matrix market, adjacency graph/hypergraph or csv format into edge list.
 **/
template <class... Attributes>
nw::graph::bi_edge_list<nw::graph::directedness::directed, Attributes...> load_graph(std::string& file) {
  std::ifstream in(file);
  std::string type;
  in >> type;

  if (type == "NW") {
    nw::util::life_timer _("deserialize");
    nw::graph::bi_edge_list<nw::graph::directedness::directed, Attributes...> aos_a(0, 0);
    aos_a.deserialize(file);
    return aos_a;
  }
  else if (type == "%%MatrixMarket") {
    std::cout << "Reading matrix market input " << file << " (slow)" << std::endl;
    //using el_t = nw::graph::bi_edge_list<Directedness, Attributes...>;
    nw::util::life_timer _("read mm");
    return nw::hypergraph::read_mm<nw::graph::directedness::directed, Attributes...>(file);
  }
  else if (type == AdjHypergraphHeader.c_str()) {
    std::cout << "Reading adjacency input " << file << " (slow)" << std::endl;
    nw::util::life_timer _("read adjacency");
    return read_adjacency<nw::graph::directedness::directed, Attributes...>(file);
  }
  else {
    std::cout << "Reading CSV input " << file << " (slow)" << std::endl;
    nw::util::life_timer _("read csv");
    return read_csv<nw::graph::directedness::directed, Attributes...>(file);
  }
}

template <nw::graph::directedness Directedness, class... Attributes>
nw::graph::edge_list<Directedness, Attributes...> load_graph(std::string& file) {
  std::ifstream in(file);
  std::string   type;
  in >> type;

  if (type == "NW") {
    nw::util::life_timer                   _("deserialize");
    nw::graph::edge_list<Directedness, Attributes...> aos_a(0);
    aos_a.deserialize(file);
    return aos_a;
  } else if (type == "%%MatrixMarket") {
    std::cout << "Reading matrix market input " << file << " (slow)\n";
    nw::util::life_timer _("read mm");
    return nw::graph::read_mm<Directedness, Attributes...>(file);
  } else {
    std::cerr << "Did not recognize graph input file " << file << "\n";
    exit(1);
  }
}

/*
 * This loader loads weighted matrix market and adjacency graph/hypergraph into edge list.
 **/
template <class... Attributes>
nw::graph::bi_edge_list<nw::graph::directedness::directed, Attributes...> load_weighted_graph(std::string file) {
  std::ifstream in(file);
  std::string type;
  in >> type;

  if (type == "NW") {
    nw::util::life_timer _("deserialize");
    nw::graph::bi_edge_list<nw::graph::directedness::directed, Attributes...> aos_a(0, 0);
    aos_a.deserialize(file);
    return aos_a;
  }
  else if (type == "%%MatrixMarket") {
    std::cout << "Reading matrix market input " << file << " (slow)" << std::endl;
    //using el_t = nw::graph::bi_edge_list<Directedness, Attributes...>;
    nw::util::life_timer _("read mm");
    return nw::hypergraph::read_mm<nw::graph::directedness::directed, Attributes...>(file);
  }
  else if (type == WghAdjHypergraphHeader.c_str()) {
    std::cout << "Reading adjacency input " << file << " (slow)" << std::endl;
    nw::util::life_timer _("read adjacency");
    return read_weighted_adjacency<nw::graph::directedness::directed, Attributes...>(file);
  }
  else {
    std::cout << "Unknown input format: " << file << std::endl;
    return nw::graph::bi_edge_list<nw::graph::directedness::directed, Attributes...>(0, 0);
  }
}

/*
 * This loader loads matrix market, adjacency graph/hypergraph or csv format into adjoin 
 * graph (in the format of edge list).
 **/
template <directedness Directedness = nw::graph::directedness::directed, class... Attributes>
nw::graph::edge_list<Directedness, Attributes...> load_adjoin_graph(std::string file, size_t& numRealEdges, size_t& numRealNodes) {
  std::ifstream in(file);
  std::string type;
  in >> type;

  if (type == "%%MatrixMarket") {
    std::cout << "Reading matrix market input " << file << " (slow)" << std::endl;
    nw::util::life_timer _("read mm adjoin");
    return read_mm_adjoin<Directedness, Attributes...>(file, numRealEdges, numRealNodes);
  }
  else if (type == AdjHypergraphHeader.c_str() || type == WghAdjHypergraphHeader.c_str()) {
    //return nw::graph::edge_list<Directedness, Attributes...>(0);   
    std::cout << "Reading adjacency input " << file << " (slow)" << std::endl;
    nw::util::life_timer _("read adjacency adjoin");
    return read_adjacency_adjoin<Directedness, Attributes...>(file, numRealEdges, numRealNodes); 
  }
  else {
    std::cout << "Reading CSV input " << file << " (slow)" << std::endl;
    nw::util::life_timer _("read csv adjoin");
    return read_csv_adjoin<Directedness, Attributes...>(file, numRealEdges, numRealNodes);
  }
}

/*
 * This loader loads adjacency graph/hypergraph format into bi-adjacency.
 **/
template<std::unsigned_integral vertex_id_t, class... Attributes>
std::tuple<nw::graph::biadjacency<0, Attributes...>, nw::graph::biadjacency<1, Attributes...>> 
load_adjacency(std::string file) {
  nw::util::life_timer _("read adjacency");
  std::ifstream in(file);
  std::string type;
  in >> type;
  if (type == AdjHypergraphHeader.c_str()) {
    std::cout << "Reading adjacency input " << file << " (slow)" << std::endl;
    return read_adj_hypergraph<vertex_id_t>(file);
  }
  else {
    std::cerr << "Did not recognize graph input file " << file << std::endl;;
    return std::make_tuple(nw::graph::biadjacency<0, Attributes...>(), nw::graph::biadjacency<1, Attributes...>());
  }
}

/*
 * This loader loads adjacency graph/hypergraph format into bi-adjacency.
 **/
template<std::unsigned_integral vertex_id_t, class... Attributes>
std::tuple<nw::graph::biadjacency<0, Attributes...>, nw::graph::biadjacency<1, Attributes...>> 
load_weighted_adjacency(std::string file) {
  nw::util::life_timer _("read adjacency");
  std::ifstream in(file);
  std::string type;
  in >> type;
  if (type == "WeightedAdjacencyHypergraph") {
    std::cout << "Reading weighted adjacency input " << file << " (slow)" << std::endl;
    return read_weighted_adj_hypergraph<vertex_id_t, Attributes...>(file);
  }
  else {
    std::cerr << "Did not recognize graph input file " << file << std::endl;;
    return std::make_tuple(nw::graph::biadjacency<0, Attributes...>(), nw::graph::biadjacency<1, Attributes...>());
  }
}

}//namespace hypergraph
}//namespace nw