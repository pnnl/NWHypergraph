#pragma once

#include <edge_list.hpp>
#include <edge_range.hpp>
#include <mmio.hpp>
#include <util/timer.hpp>
#include <util/traits.hpp>

#include <tbb/global_control.h>
#include <tbb/parallel_for.h>

#include <iomanip>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <bitset>

#include "io/hypergraph_io.hpp"

using namespace nw::graph;

namespace nw::hypergraph {
namespace tools {

#if defined (EXECUTION_POLICY)
constexpr inline bool WITH_TBB = true;
#else
constexpr inline bool WITH_TBB = false;
#endif

auto set_n_threads(long n) {
  if constexpr (WITH_TBB) {
    return tbb::global_control(tbb::global_control::max_allowed_parallelism, n);
  }
  else {
    return 0;
  }
}

long get_n_threads() {
  if constexpr (WITH_TBB) {
    return tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
  }
  else {
    return 1;
  }
}

std::vector<long> parse_n_threads(const std::vector<std::string>& args) {
  std::vector<long> threads;
  if constexpr (WITH_TBB) {
    if (args.size() == 0) {
      threads.push_back(tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism));
    }
    else {
      for (auto&& n : args) {
        threads.push_back(std::stol(n));
      }
    }
  }
  else {
    threads.push_back(1);
  }
  //if no thread number is given by user, set 1 thread as default
  if (threads.empty())
    threads.push_back(1);
  return threads;
}

std::vector<long> parse_ids(const std::vector<std::string>& args) {
  std::vector<long> ids;
  for (auto&& n : args) {
    ids.push_back(std::stol(n));
  }
  return ids;
}
/*
* parse the option string as a bitset
*/
template<size_t N = 8>
std::bitset<N> parse_features(const std::vector<std::string>& args) {
  std::bitset<N> ids;
  for (auto&& n : args) {
    auto i = std::stoul(n);
    if (i > N) {
      std::cerr << "Wrong feature specified" << std::endl;
      exit(1);
    }
    ids[i] = true;
  }
  return ids;
}

template <directedness Directedness, class... Attributes>
nw::graph::edge_list<Directedness, Attributes...> load_graph(std::string file) {
  std::ifstream in(file);
  std::string type;
  in >> type;

  if (type == "BGL17") {
    nw::util::life_timer _("deserialize");
    nw::graph::edge_list<Directedness, Attributes...> aos_a(0);
    aos_a.deserialize(file);
    return aos_a;
  }
  else if (type == "%%MatrixMarket") {
    std::cout << "Reading matrix market input " << file << " (slow)" << std::endl;
    nw::util::life_timer _("read mm");
    return nw::graph::read_mm<Directedness, Attributes...>(file);
  }
  else {
    //std::cerr << "Did not recognize graph input file " << file << "\n";
    return nw::graph::edge_list<Directedness, Attributes...>(0);
  }
}

template<class... Attributes>
std::tuple<nw::graph::adjacency<0, Attributes...>, nw::graph::adjacency<1, Attributes...>> 
load_adjacency(std::string file) {
  std::ifstream in(file);
  std::string type;
  in >> type;
  if (type == "AdjacencyHypergraph") {
    std::cout << "Reading adjacency hypergraph input " << file << " (slow)" << std::endl;
    return nw::hypergraph::read_adj_hypergraph<Attributes...>(file);
  }
  else if (type == "WeightedAdjacencyHypergraph") {
    std::cout << "Reading weighted adjacency hypergraph input " << file << " (slow)" << std::endl;
    return nw::hypergraph::read_weighted_adj_hypergraph<Attributes...>(file);
  }
  else {
    std::cerr << "Did not recognize graph input file " << file << std::endl;;
    exit(1);
  }
}

template <int Adj, directedness Directedness, class... Attributes>
nw::graph::adjacency<Adj, Attributes...> build_adjacency(nw::graph::edge_list<Directedness, Attributes...>& graph) {
  nw::util::life_timer _("build adjacency");
  return { graph };
}

template <class Graph>
auto build_degrees(Graph&& graph)
{
  nw::util::life_timer _("degrees");
  std::vector<nw::graph::vertex_id_t> degrees(graph.size());
  tbb::parallel_for(edge_range(graph), [&](auto&& edges) {
    for (auto&& [i, j] : edges) {
      __atomic_fetch_add(&degrees[j], 1, __ATOMIC_ACQ_REL);
    }
  });
  return degrees;
}

template <class Graph>
auto build_random_sources(Graph&& graph, size_t n, long seed)
{
  using Id = typename nw::graph::vertex_id<std::decay_t<Graph>>::type;

  auto sources = std::vector<nw::graph::vertex_id_t>(n);
  auto degrees = build_degrees(graph);
  auto     gen = std::mt19937(seed);
  auto     dis = std::uniform_int_distribution<Id>(0, graph.max());

  for (auto& id : sources) {
    for (id = dis(gen); degrees[id] == 0; id = dis(gen)) {}
  }
  return sources;
}

/// Load a set of vertices from a file.
///
/// This will load a set of vertices from the passed `file` and verify that we
/// have the expected number `n`.
template <class Graph>
auto load_sources_from_file(Graph&&, std::string file, size_t n = 0) {
  std::vector sources = read_mm_vector<nw::graph::vertex_id_t>(file);
  if (n && sources.size() != n) {
    std::cerr << file << " contains " << sources.size() << " sources, however options require " << n << "\n";
    exit(1);
  }
  return sources;
}

/// Helper to time an operation.
template <class Op>
auto time_op(Op&& op) {
  if constexpr (std::is_void_v<decltype(op())>) {
    auto start = std::chrono::high_resolution_clock::now();
    op();
    std::chrono::duration<double> end = std::chrono::high_resolution_clock::now() - start;
    return std::tuple{ end.count() };
  }
  else {
    auto start = std::chrono::high_resolution_clock::now();
    auto e = op();
    std::chrono::duration<double> end = std::chrono::high_resolution_clock::now() - start;
    return std::tuple { end.count(), std::move(e) };
  }
}

template <class Op, class Check>
auto time_op_verify(Op&& op, Check&& check) {
  auto&& [t, e] = time_op(std::forward<Op>(op));
  return std::tuple(t, check(std::move(e)));
}


}//bench
}//nw::hypergraph
