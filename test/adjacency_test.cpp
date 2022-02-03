//
// This file is part of NWHypergraph
// (c) Pacific Northwest National Laboratory 2018-2021
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include <catch2/catch.hpp>
#include <nwgraph/edge_list.hpp>
#include <nwgraph/adjacency.hpp>
#include <nwgraph/adaptors/neighbor_range.hpp>
#include <algorithm>

using namespace nw::graph;


TEST_CASE("bi-adjacency", "[bi-adjacency]") {

  SECTION("bi-edgelist constructor with push_back") {
     /// For a bipartite graph, the edge list format must be directed
    bi_edge_list<directedness::directed> A_list;
    A_list.open_for_push_back();
    A_list.push_back(0, 1);
    A_list.push_back(1, 2);
    A_list.push_back(2, 3);
    A_list.push_back(3, 4);
    A_list.close_for_push_back();
    A_list.stream_stats();
    A_list.stream_edges();  
    REQUIRE(4 == A_list.num_vertices()[0]);
    REQUIRE(5 == A_list.num_vertices()[1]);
    REQUIRE(4 == num_vertices(A_list));
    REQUIRE(5 == num_vertices(A_list, 1));   
  }

  SECTION("bi-adjacency constructor passing a unweighted edge list") {
   
    /// Test for bipartite graph using adjacency
    /// For a bipartite graph, the edge list format must be directed
    bi_edge_list<directedness::directed> A_list;
    A_list.open_for_push_back();
    A_list.push_back(0, 1);
    A_list.push_back(1, 2);
    A_list.push_back(2, 3);
    A_list.push_back(3, 4);
    A_list.close_for_push_back();

    /// pass in edge list without any constraint
    /// assume the edge list to be bipartite
    biadjacency<0> A(A_list);
    A.stream_stats();
    A.stream_indices();
    biadjacency<1> AT(A_list);
    AT.stream_stats();
    AT.stream_indices();
    for (auto&& [u, neighborhood] : neighbor_range(A)) {
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == A.num_vertices()[0]);
    REQUIRE(5 == A.num_vertices()[1]);
    REQUIRE(4 == num_vertices(A));
    REQUIRE(5 == num_vertices(A, 1));

    for (auto&& [u, neighborhood] : neighbor_range(AT)) {
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(5 == AT.num_vertices()[0]);
    REQUIRE(4 == AT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(AT));
    REQUIRE(4 == num_vertices(AT, 1));
  }

  SECTION("bi-adjacency constructor passing a weighted edge list") {
    /// Test for bipartite graph using adjacency
    /// For a bipartite graph, the edge list format must be directed
    bi_edge_list<directedness::directed, int> A_list;
    A_list.open_for_push_back();
    A_list.push_back(0, 1, 1);
    A_list.push_back(1, 2, 2);
    A_list.push_back(2, 3, 3);
    A_list.push_back(3, 4, 4);
    A_list.close_for_push_back();
    /// pass in edge list without any constraint
    /// assume the edge list to be bipartite
    biadjacency<0, int> A(A_list);
    biadjacency<1, int> AT(A_list);

    for (auto&& [u, neighborhood] : neighbor_range(A)) {
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == A.num_vertices()[0]);
    REQUIRE(5 == A.num_vertices()[1]);
    REQUIRE(4 == num_vertices(A));
    REQUIRE(5 == num_vertices(A, 1));

    for (auto&& [u, neighborhood] : neighbor_range(AT)) {
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(5 == AT.num_vertices()[0]);
    REQUIRE(4 == AT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(AT));
    REQUIRE(4 == num_vertices(AT, 1));
  }

  SECTION("bi-adjacency constructor passing partition size and a weighted edge list") {
    /// Test for bipartite graph using adjacency
    /// For a bipartite graph, the edge list format must be directed
    bi_edge_list<directedness::directed, int> A_list;
    A_list.open_for_push_back();
    A_list.push_back(0, 1, 1);
    A_list.push_back(1, 2, 2);
    A_list.push_back(2, 3, 3);
    A_list.push_back(3, 4, 4);
    A_list.close_for_push_back();
    biadjacency<0, int> B(A_list);
    biadjacency<1, int> BT(A_list);

    for (auto&& [u, neighborhood] : neighbor_range(B)) {
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == B.num_vertices()[0]);
    REQUIRE(5 == B.num_vertices()[1]);
    REQUIRE(4 == num_vertices(B));
    REQUIRE(5 == num_vertices(B, 1));

    for (auto&& [u, neighborhood] : neighbor_range(BT)) {
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(5 == BT.num_vertices()[0]);
    REQUIRE(4 == BT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(BT));
    REQUIRE(4 == num_vertices(BT, 1));
  }
  SECTION("bi-adjacency copy constructor passing a unweighted csr (indices and offsets)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    //IMPORTANT: remember the number of vertices is the indices.size() - 1
    // the last element in indices indicates the last offest of the edges
    biadjacency<0> C(v1.size() - 1, v0, e0);
    biadjacency<1> CT(v0.size() - 1, v1, e1);

    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }

  SECTION("bi-adjacency copy constructor passing a weighted csr (indices and offsets and weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    
    biadjacency<0, double> C(v1.size() - 1, v0, e0, w0);
    biadjacency<1, double> CT(v0.size() - 1, v1, e1, w1);

    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }
  SECTION("bi-adjacency copy constructor passing a weighted csr (indices and offsets and two weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<uint32_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
    
    biadjacency<0, double, float> C(v1.size() - 1, v0, e0);
    biadjacency<1, double, float> CT(v0.size() - 1, v1, e1);

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));


    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }
  
  SECTION("bi-adjacency move constructor passing a unweighted csr (indices and offsets)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    size_t N0 = v0.size() - 1, N1 = v1.size() - 1;
    biadjacency<0> C(N1, std::move(v0), std::move(e0));
    //NOTE: Do NOT use like follow, the move constructor above will erase v0
    //v0.size() becomes 0 after the above move constructor
    //biadjacency<1> CT(v0.size() - 1, std::move(v1), std::move(e1));
    biadjacency<1> CT(N0, std::move(v1), std::move(e1));

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }

  SECTION("bi-adjacency move constructor passing a weighted csr (indices and offsets and weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    size_t N0 = v0.size() - 1, N1 = v1.size() - 1;    
    biadjacency<0, double> C(N1, std::move(v0), std::move(e0), std::move(w0));
    biadjacency<1, double> CT(N0, std::move(v1), std::move(e1), std::move(w1));

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }
  SECTION("bi-adjacency move constructor passing a weighted csr (indices and offsets and two weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<uint32_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
        
    size_t N0 = v0.size() - 1, N1 = v1.size() - 1;
    biadjacency<0, double, float> C(N1, std::move(v0), std::move(e0));
    biadjacency<1, double, float> CT(N0, std::move(v1), std::move(e1));

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));


    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }
  
  SECTION("bi-adjacency default constructor followed by copying a unweighted csr (indices and offsets)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    biadjacency<0> C(v0.size() - 1, v1.size() - 1, e0.size());
    C.copy(v0, e0);
    biadjacency<1> CT(v1.size() - 1, v0.size() - 1, e1.size());
    CT.copy(v1, e1);

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }
  
  SECTION("bi-adjacency default constructor followed by copying a weighted csr (indices and offsets and weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    biadjacency<0, double> C(v0.size() - 1, v1.size() - 1, e0.size());
    C.copy(v0, e0, w0);
    biadjacency<1, double> CT(v1.size() - 1, v0.size() - 1, e1.size());
    CT.copy(v1, e1, w1);

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }

  SECTION("bi-adjacency default constructor followed by copying a weighted csr (indices and offsets and two weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<uint32_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
    
    biadjacency<0, double, float> C(v0.size() - 1, v1.size() - 1, std::get<0>(e0).size());
    C.copy(v0, e0);
    biadjacency<1, double, float> CT(v1.size() - 1, v0.size() - 1, std::get<0>(e1).size());
    CT.copy(v1, e1);

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));


    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }

  SECTION("bi-adjacency default constructor followed by moving a unweighted csr (indices and offsets)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    biadjacency<0> C(v0.size() - 1, v1.size() - 1, e0.size());
    C.move(std::move(v0), std::move(e0));
    biadjacency<1> CT(v1.size() - 1, v0.size() - 1, e1.size());
    CT.move(std::move(v1), std::move(e1));

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }
  
  SECTION("bi-adjacency default constructor followed by moving a weighted csr (indices and offsets and weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<uint32_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    biadjacency<0, double> C(v0.size() - 1, v1.size() - 1, e0.size());
    C.move(std::move(v0), std::move(e0), std::move(w0));
    biadjacency<1, double> CT(v1.size() - 1, v0.size() - 1, e1.size());
    CT.move(std::move(v1), std::move(e1), std::move(w1));

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));

    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }

  SECTION("bi-adjacency default constructor followed by moving a weighted csr (indices and offsets and two weights)") {
    std::vector<uint32_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<uint32_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
    
    size_t N0 = v0.size() - 1, N1 = v1.size() - 1;
    biadjacency<0, double, float> C(N0, N1, std::get<0>(e0).size());
    C.copy(std::move(v0), std::move(e0));
    biadjacency<1, double, float> CT(N1, N0, std::get<0>(e1).size());
    CT.copy(std::move(v1), std::move(e1));

    uint32_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max - u < 0 ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(4 == C.num_vertices()[0]);
    REQUIRE(5 == C.num_vertices()[1]);
    REQUIRE(4 == num_vertices(C));
    REQUIRE(5 == num_vertices(C, 1));


    uint32_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max - u < 0? u : col1_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(5 == CT.num_vertices()[0]);
    REQUIRE(4 == CT.num_vertices()[1]);
    REQUIRE(5 == num_vertices(CT));
    REQUIRE(4 == num_vertices(CT, 1));
  }
}