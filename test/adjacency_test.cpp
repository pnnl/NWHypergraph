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
#include <containers/edge_list.hpp>
#include <adaptors/neighbor_range.hpp>
#include <algorithm>

using namespace nw::graph;


TEST_CASE("bi-adjacency", "[bi-adjacency]") {

  SECTION("bi-adjacency constructor passing a unweighted edge list") {
   
    /// Test for bipartite graph using adjacency
    /// For a bipartite graph, the edge list format must be directed
    edge_list<directedness::directed> A_list;
    A_list.open_for_push_back();
    A_list.push_back(0, 1);
    A_list.push_back(1, 2);
    A_list.push_back(2, 3);
    A_list.push_back(3, 4);
    A_list.close_for_push_back();
    /// pass in edge list without any constraint
    /// assume the edge list to be bipartite
    adjacency<0> A(A_list);
    adjacency<1> AT(A_list);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(A)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(AT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == col1_max);
  }

  SECTION("bi-adjacency constructor passing a weighted edge list") {
    /// Test for bipartite graph using adjacency
    /// For a bipartite graph, the edge list format must be directed
    edge_list<directedness::directed, int> A_list;
    A_list.open_for_push_back();
    A_list.push_back(0, 1, 1);
    A_list.push_back(1, 2, 2);
    A_list.push_back(2, 3, 3);
    A_list.push_back(3, 4, 4);
    A_list.close_for_push_back();
    /// pass in edge list without any constraint
    /// assume the edge list to be bipartite
    adjacency<0, int> A(A_list);
    adjacency<1, int> AT(A_list);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(A)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(AT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }

  SECTION("bi-adjacency constructor passing partition size and a weighted edge list") {
    /// Test for bipartite graph using adjacency
    /// For a bipartite graph, the edge list format must be directed
    edge_list<directedness::directed, int> A_list;
    A_list.open_for_push_back();
    A_list.push_back(0, 1, 1);
    A_list.push_back(1, 2, 2);
    A_list.push_back(2, 3, 3);
    A_list.push_back(3, 4, 4);
    A_list.close_for_push_back();
    adjacency<0, int> B(A_list.max()[0] + 1, A_list);
    adjacency<1, int> BT(A_list.max()[1] + 1, A_list);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(B)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(BT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
  SECTION("bi-adjacency copy constructor passing a unweighted csr (indices and offsets)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    adjacency<0> C(v0, e0);
    adjacency<1> CT(v1, e1);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == col1_max);
  }

  SECTION("bi-adjacency copy constructor passing a weighted csr (indices and offsets and weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    
    adjacency<0, double> C(v0, e0, w0);
    adjacency<1, double> CT(v1, e1, w1);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
  SECTION("bi-adjacency copy constructor passing a weighted csr (indices and offsets and two weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<vertex_id_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
    
    adjacency<0, double, float> C(v0, e0);
    adjacency<1, double, float> CT(v1, e1);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(3 == col0_max);


    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
  
  SECTION("bi-adjacency move constructor passing a unweighted csr (indices and offsets)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    adjacency<0> C(std::move(v0), std::move(e0));
    adjacency<1> CT(std::move(v1), std::move(e1));

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == col1_max);
  }

  SECTION("bi-adjacency move constructor passing a weighted csr (indices and offsets and weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    
    adjacency<0, double> C(std::move(v0), std::move(e0), std::move(w0));
    adjacency<1, double> CT(std::move(v1), std::move(e1), std::move(w1));

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
  SECTION("bi-adjacency move constructor passing a weighted csr (indices and offsets and two weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<vertex_id_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
    
    adjacency<0, double, float> C(std::move(v0), std::move(e0));
    adjacency<1, double, float> CT(std::move(v1), std::move(e1));

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(3 == col0_max);


    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
  
  SECTION("bi-adjacency default constructor followed by copying a unweighted csr (indices and offsets)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    adjacency<0> C(v0.size() - 1, e0.size());
    C.copy(v0, e0);
    adjacency<1> CT(v1.size() - 1, e1.size());
    CT.copy(v1, e1);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
  
  SECTION("bi-adjacency default constructor followed by copying a weighted csr (indices and offsets and weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    adjacency<0, double> C(v0.size() - 1, e0.size());
    C.copy(v0, e0, w0);
    adjacency<1, double> CT(v1.size() - 1, e1.size());
    CT.copy(v1, e1, w1);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }

  SECTION("bi-adjacency default constructor followed by copying a weighted csr (indices and offsets and two weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<vertex_id_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
    
    adjacency<0, double, float> C(v0.size() - 1, std::get<0>(e0).size());
    C.copy(v0, e0);
    adjacency<1, double, float> CT(v1.size() - 1, std::get<0>(e1).size());
    CT.copy(v1, e1);

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(3 == col0_max);


    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }

  SECTION("bi-adjacency default constructor followed by moving a unweighted csr (indices and offsets)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    adjacency<0> C(v0.size() - 1, e0.size());
    C.move(std::move(v0), std::move(e0));
    adjacency<1> CT(v1.size() - 1, e1.size());
    CT.move(std::move(v1), std::move(e1));

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v] : neighborhood)
        std::cout << "edge " << u << " to " << v << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
  
  SECTION("bi-adjacency default constructor followed by moving a weighted csr (indices and offsets and weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::vector<vertex_id_t> e0{1, 2, 3, 4}, e1{0, 1, 2, 3};
    std::vector<double> w0{1, 2, 3, 4}, w1{1, 2, 3, 4};
    adjacency<0, double> C(v0.size() - 1, e0.size());
    C.move(std::move(v0), std::move(e0), std::move(w0));
    adjacency<1, double> CT(v1.size() - 1, e1.size());
    CT.move(std::move(v1), std::move(e1), std::move(w1));

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(3 == col0_max);

    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }

  SECTION("bi-adjacency default constructor followed by moving a weighted csr (indices and offsets and two weights)") {
    std::vector<vertex_id_t> v0{0, 1, 2, 3, 4}, v1{0, 0, 1, 2, 3, 4};
    std::tuple<std::vector<vertex_id_t>, std::vector<double>, std::vector<float>> e0{
      {1, 2, 3, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}
    }, e1 {
      {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}
    };
    
    adjacency<0, double, float> C(v0.size() - 1, std::get<0>(e0).size());
    C.copy(std::move(v0), std::move(e0));
    adjacency<1, double, float> CT(v1.size() - 1, std::get<0>(e1).size());
    CT.copy(std::move(v1), std::move(e1));

    vertex_id_t col0_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(C)) {
      col0_max = col0_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood) {
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
        REQUIRE(w1 == w2);
      }
    }
    REQUIRE(3 == col0_max);


    vertex_id_t col1_max = 0;
    for (auto&& [u, neighborhood] : neighbor_range(CT)) {
      col1_max = col1_max < u ? u : col0_max;
      for (auto&& [v, w1, w2] : neighborhood)
        std::cout << "edge " << u << " to " << v << " has weight " << w1
                  << std::endl;
    }
    REQUIRE(4 == col1_max);
  }
}