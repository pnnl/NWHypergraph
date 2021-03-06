/**
 * @file common.hpp
 *
 * @copyright SPDX-FileCopyrightText: 2022 Battelle Memorial Institute
 * @copyright SPDX-FileCopyrightText: 2022 University of Washington
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * 
 * Author: Xu Tony Liu
 * 
 */

#pragma once

#include <nwgraph/edge_list.hpp>
#include <nwgraph/adaptors/edge_range.hpp>
#include <nwgraph/io/mmio.hpp>
#include <nwgraph/util/timer.hpp>
#include <nwgraph/util/traits.hpp>

#include <tbb/global_control.h>
#include <tbb/parallel_for.h>

#include <variant>
#include <iomanip>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>
#include <bitset>

#include "io/loader.hpp"
#include "containers/compressed_hy.hpp"
#include "containers/edge_list_hy.hpp"

using namespace nw::graph;

namespace nw::hypergraph {
namespace bench {

constexpr inline bool WITH_TBB = true;

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


/*
 * This graph reader loads the graph into edge list then convert to bi-adjacency.
 **/
template<std::unsigned_integral vertex_id_t, 
typename... Attributes>
auto graph_reader(std::string file) {
  auto aos_a = load_graph<Attributes...>(file);
  if (0 == aos_a.size()) {
    auto&&[hyperedges, hypernodes] = load_adjacency<vertex_id_t, Attributes...>(file);
    std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
    return std::tuple(hyperedges, hypernodes, std::vector<vertex_id_t>());
  }
  nw::graph::biadjacency<0, Attributes...> hyperedges(aos_a);
  nw::graph::biadjacency<1, Attributes...> hypernodes(aos_a);

  std::cout << "num_hyperedges = " << hyperedges.size()
            << " num_hypernodes = " << hypernodes.size() << std::endl;
  return std::tuple(hyperedges, hypernodes, std::vector<vertex_id_t>());
}

/*
 * This graph reader loads weighted graph into edge list then convert to bi-adjacency.
 **/
template<std::unsigned_integral vertex_id_t, 
typename... Attributes>
auto weighted_graph_reader(std::string file) {
  auto aos_a = load_weighted_graph<Attributes...>(file);
  if (0 == aos_a.size()) {
    auto&&[hyperedges, hypernodes] = load_weighted_adjacency<vertex_id_t, Attributes...>(file);
    std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
    return std::tuple(hyperedges, hypernodes, std::vector<vertex_id_t>());
  }
  nw::graph::biadjacency<0, Attributes...> hyperedges(aos_a);
  nw::graph::biadjacency<1, Attributes...> hypernodes(aos_a);

  std::cout << "num_hyperedges = " << hyperedges.size()
            << " num_hypernodes = " << hypernodes.size() << std::endl;
  return std::tuple(hyperedges, hypernodes, std::vector<vertex_id_t>());
}

/*
 * This graph reader loads the graph into edge list then convert to bi-adjacency.
 * And can relabel the ids of either hyperedges or hypernodes by degree
 * in the direction of ascending (default) or descending.
 **/
template<std::unsigned_integral vertex_id_t, 
typename... Attributes>
auto graph_reader_relabel(std::string file, int idx, std::string direction) {
  auto aos_a = load_graph<Attributes...>(file);
  std::vector<vertex_id_t> iperm;
  if (0 == aos_a.size()) {
    auto&& [hyperedges, hypernodes] = load_adjacency<vertex_id_t, Attributes...>(file);
    // Run relabeling. This operates directly on the incoming edglist.
    if (-1 != idx) {
      auto&& iperm = nw::hypergraph::relabel_by_degree(
          hyperedges, hypernodes, idx, direction);
        std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
      return std::tuple(hyperedges, hypernodes, iperm);
    }
    std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
    return std::tuple(hyperedges, hypernodes, iperm); 
  }
  // Run relabeling. This operates directly on the incoming edglist.
  if (-1 != idx) {
    std::cout << "relabeling edge_list by degree..." << std::endl;
    if (1 == idx) {
      iperm = nw::hypergraph::relabel_by_degree<1>(aos_a, direction);
    }
    else {
      iperm = nw::hypergraph::relabel_by_degree<0>(aos_a, direction);
    }
  }
  nw::graph::biadjacency<0, Attributes...> hyperedges(aos_a);
  nw::graph::biadjacency<1, Attributes...> hypernodes(aos_a);

  std::cout << "num_hyperedges = " << hyperedges.size()
            << " num_hypernodes = " << hypernodes.size() << std::endl;
  return std::tuple(hyperedges, hypernodes, iperm);
}

/*
 * This weighted graph reader loads the graph into edge list then convert to bi-adjacency.
 * And can relabel the ids of either hyperedges or hypernodes by degree
 * in the direction of ascending (default) or descending.
 **/
template<std::unsigned_integral vertex_id_t, 
typename... Attributes>
auto weighted_graph_reader_relabel(std::string file, int idx, std::string direction) {
  auto aos_a = load_weighted_graph<Attributes...>(file);
  std::vector<vertex_id_t> iperm;
  if (0 == aos_a.size()) {
    auto&& [hyperedges, hypernodes] = load_weighted_adjacency<vertex_id_t, Attributes...>(file);
    // Run relabeling. This operates directly on the incoming edglist.
    if (-1 != idx) {
      auto&& iperm = nw::hypergraph::relabel_by_degree(
          hyperedges, hypernodes, idx, direction);
        std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
      return std::tuple(hyperedges, hypernodes, iperm);
    }
    std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
    return std::tuple(hyperedges, hypernodes, iperm); 
  }
  // Run relabeling. This operates directly on the incoming edglist.
  if (-1 != idx) {
    std::cout << "relabeling edge_list by degree..." << std::endl;
    if (1 == idx) {
      iperm = nw::hypergraph::relabel_by_degree<1>(aos_a, direction);
    }
    else {
      iperm = nw::hypergraph::relabel_by_degree<0>(aos_a, direction);
    }
  }
  nw::graph::biadjacency<0, Attributes...> hyperedges(aos_a);
  nw::graph::biadjacency<1, Attributes...> hypernodes(aos_a);

  std::cout << "num_hyperedges = " << hyperedges.size()
            << " num_hypernodes = " << hypernodes.size() << std::endl;
  return std::tuple(hyperedges, hypernodes, iperm);
}

/*
 * This graph reader can read either mtx or hygra (adjacency) file 
 * in adjoin fashion.
 * Adjoin will increment the ids of the smaller set of the hyperedges or hypernodes,
 * such that the ids of both are now contiguous.
 * And it also can relabel the ids of either hyperedges or hypernodes by degree
 * in the direction of ascending (default) or descending.
 **/
template<directedness edge_directedness = directedness::undirected, 
std::unsigned_integral vertex_id_t, 
typename... Attributes>
auto graph_reader_adjoin_and_relabel(std::string file, const int idx, std::string direction, 
size_t& nrealedges, size_t& nrealnodes) {
  auto aos_a = load_adjoin_graph<edge_directedness, Attributes...>(
      file, nrealedges, nrealnodes);
  std::vector<vertex_id_t> iperm;
  if (0 == aos_a.size()) {
    //we get adjoin graph and its transpose
    auto&& [g, gt] = read_and_adjoin_adj_hypergraph_pair<vertex_id_t>(file, nrealedges, nrealnodes);
    if (-1 != idx) {
      // Run relabeling. This operates directly on adjoin graph and its transpose.
      // Note that g and gt are symmetric,
      // therefore the permutation of g is also the permuation of its transpose
      nw::util::ms_timer t("relabel_by_degree");
      auto&& iperm =
            g.permute_by_degree(direction, std::execution::par_unseq);
      g.relabel_to_be_indexed(iperm, std::execution::par_unseq);
      auto&& iperm_t =
            gt.permute_by_degree(direction, std::execution::par_unseq);
      gt.relabel_to_be_indexed(iperm_t, std::execution::par_unseq);
      t.stop();
      std::cout << t << std::endl;

      std::cout << "num_realedges = " << nrealedges
                << " num_realnodes = " << nrealnodes << std::endl;
      if (0 == idx)
        return std::tuple(g, gt, iperm);
      else 
        return std::tuple(g, gt, iperm_t);
    }

    std::cout << "num_realedges = " << nrealedges
            << " num_realnodes = " << nrealnodes << std::endl; 
    std::cout << "num_vertices = " << g.size() << std::endl;  
    return std::tuple(g, gt, iperm);
  }
  else {
    // Run relabeling. This operates directly on the incoming edglist.
    if (-1 != idx) {
      //since the adjoin graph is symmetric, when we need to relabel it
      //we need to operate on both column 0 and column 1
      nw::util::ms_timer t("relabel_edge_list_by_degree");
      if (1 == idx) {
        iperm = nw::graph::relabel_by_degree<1>(aos_a, direction);
      }
      else {
        iperm = nw::graph::relabel_by_degree<0>(aos_a, direction);
      }
      t.stop();
      std::cout << t << std::endl;
    }

    nw::graph::adjacency<0, Attributes...> g(aos_a);
    nw::graph::adjacency<1, Attributes...> g_t(aos_a);
    std::cout << "num_realedges = " << nrealedges
            << " num_realnodes = " << nrealnodes << std::endl;  
    std::cout << "num_vertices = " << g.size() << std::endl;  
    return std::tuple(g, g_t, iperm);
  }
}

/*
 * This graph reader loads the graph into adjoin graph (edge list) then convert to bi-adjacency.
 * Adjoin will increment the ids of the smaller set of the hyperedges or hypernodes,
 * such that the ids of both are now contiguous.
 **/
template <directedness edge_directedness = directedness::directed, 
std::unsigned_integral vertex_id_t, 
typename... Attributes>
auto graph_reader_adjoin(std::string file, size_t& nrealedges,
                         size_t& nrealnodes) {
  auto aos_a = load_adjoin_graph<edge_directedness, Attributes...>(
      file, nrealedges, nrealnodes);
  if (0 == aos_a.size()) {
    //after adjoin, we have adjoin graph and its transpose
    auto&& [g, g_t] =
        read_and_adjoin_adj_hypergraph_pair<vertex_id_t>(file, nrealedges, nrealnodes);
    std::cout << "num_realedges = " << nrealedges
            << " num_realnodes = " << nrealnodes << std::endl;
    std::cout << "num_vertices = " << g.size() << std::endl;
    return std::tuple(g, g_t, std::vector<vertex_id_t>());
  }
  nw::graph::adjacency<0, Attributes...> g(aos_a);
  nw::graph::adjacency<1, Attributes...> g_t(aos_a);
  std::cout << "num_realedges = " << nrealedges
            << " num_realnodes = " << nrealnodes << std::endl;
  std::cout << "num_vertices = " << g.size() << std::endl;
  return std::tuple(g, g_t, std::vector<vertex_id_t>());
}

template <std::unsigned_integral vertex_id_t = vertex_id_t<nw::graph::biadjacency<0>>, typename... Attributes>
auto adjoin_graph_reader(std::string file, const int idx, std::string direction,
                         size_t& nrealedges, size_t& nrealnodes) {
  if (-1 != idx)
    return graph_reader_adjoin_and_relabel<directedness::undirected,
                                           vertex_id_t, Attributes...>(
        file, idx, direction, nrealedges, nrealnodes);
  else
    return graph_reader_adjoin<directedness::undirected, vertex_id_t,
                               Attributes...>(file, nrealedges, nrealnodes);
}

template<std::unsigned_integral vertex_id_t = vertex_id_t<nw::graph::biadjacency<0>>, 
typename... Attributes>
auto graph_reader(std::string file, int idx, std::string direction) {
  if (-1 != idx)
    return graph_reader_relabel<vertex_id_t, Attributes...>(file, idx, direction);
  else
    return graph_reader<vertex_id_t, Attributes...>(file);
}

template<std::unsigned_integral vertex_id_t = vertex_id_t<nw::graph::biadjacency<0>>, 
typename... Attributes>
auto weighted_graph_reader(std::string file, int idx, std::string direction) {
  if (-1 != idx)
    return weighted_graph_reader_relabel<vertex_id_t, Attributes...>(file, idx, direction);
  else
    return weighted_graph_reader<vertex_id_t, Attributes...>(file);
}

template <int Adj, directedness Directedness, class... Attributes>
nw::graph::adjacency<Adj, Attributes...> build_adjacency(nw::graph::edge_list<Directedness, Attributes...>& graph) {
  nw::util::life_timer _("build adjacency");
  return { graph };
}

template <class Graph>
auto build_random_sources(Graph&& graph, size_t n, long seed)
{
  using vertex_id_t = typename graph_traits<Graph>::vertex_id_type;  

  auto sources = std::vector<vertex_id_t>(n);
  auto degrees = nw::graph::degrees(graph);
  auto     gen = std::mt19937(seed);
  auto     dis = std::uniform_int_distribution<vertex_id_t>(0, graph.max());

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
  std::vector sources = read_mm_vector<vertex_id_t<Graph>>(file);
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

template <class... Extra>
class Times
{
  using Sample = std::tuple<double, Extra...>;
  std::map<std::tuple<std::string, long, long>, std::vector<Sample>> times_ = {};

 public:
  decltype(auto) begin() const { return times_.begin(); }
  decltype(auto)   end() const { return times_.end(); }

  template <class Op>
  auto record(std::string file, long id, long thread, Op&& op, Extra... extra) {
    return std::apply([&](auto time, auto&&... rest) {
      append(file, id, thread, time, extra...);
      return std::tuple { std::forward<decltype(rest)>(rest)... };
    }, time_op(std::forward<Op>(op)));
  }

  template <class Op, class Verify>
  void record(std::string file, long id, long thread, Op&& op, Verify&& verify, Extra... extra) {
    auto&& [time, result] = time_op(std::forward<Op>(op));
    verify(std::forward<decltype(result)>(result));
    append(file, id, thread, time, extra...);
  }

  void append(std::string file, long id, long thread, double trial, Extra... extra) {
    times_[std::tuple(file, id, thread)].emplace_back(trial, extra...);
  }

  void print(std::ostream& out) const {
    std::size_t n = 4;
    for (auto&& [config, samples] : times_) {
      n = std::max(n, std::get<0>(config).size());
    }

    out << std::setw(n + 2) << std::left << "File";
    out << std::setw(10) << std::left << "Version";
    out << std::setw(10) << std::left << "Threads";
    out << std::setw(20) << std::left << "Min";
    out << std::setw(20) << std::left << "Avg";
    out << std::setw(20) << std::left << "Max";
    out << "\n";

    for (auto&& [config, samples] : times_) {
      auto [file, id, threads] = config;
      auto [min, max, avg] = minmaxavg(samples);

      out << std::setw(n + 2) << std::left << file;
      out << std::setw(10) << std::left << id;
      out << std::setw(10) << std::left << threads;
      out << std::setw(20) << std::left << std::setprecision(6) << std::fixed << min;
      out << std::setw(20) << std::left << std::setprecision(6) << std::fixed << avg;
      out << std::setw(20) << std::left << std::setprecision(6) << std::fixed << max;
      out << "\n";
    }
  }

 private:
  static auto average(const std::vector<Sample>& times) {
    double total = 0.0;
    for (auto&& sample : times) {
      total += std::get<0>(sample);
    }
    return total / times.size();
  }

  static auto minmax(const std::vector<Sample>& times) {
    return std::apply([](auto... minmax) {
      return std::tuple(std::get<0>(*minmax)...);
    }, std::minmax_element(times.begin(), times.end(), [](auto&& a, auto&& b) {
      return std::get<0>(a) < std::get<0>(b);
    }));
  }

  static auto minmaxavg(const std::vector<Sample>& times) {
    return std::apply([&](auto... minmax) {
      return std::tuple(minmax..., average(times));
    }, minmax(times));
  }
};

template <class... Extra>
class Times_WithS
{
  using Sample = std::tuple<double, Extra...>;
  std::map<std::tuple<std::string, long, long, long>, std::vector<Sample>> times_ = {};

 public:
  decltype(auto) begin() const { return times_.begin(); }
  decltype(auto)   end() const { return times_.end(); }

  template <class Op>
  auto record(std::string file, long id, long thread, long s, Op&& op, Extra... extra) {
    return std::apply([&](auto time, auto&&... rest) {
      append(file, id, thread, s, time, extra...);
      return std::tuple { std::forward<decltype(rest)>(rest)... };
    }, time_op(std::forward<Op>(op)));
  }

  template <class Op, class Verify>
  void record(std::string file, long id, long thread, long s, Op&& op, Verify&& verify, Extra... extra) {
    auto&& [time, result] = time_op(std::forward<Op>(op));
    verify(std::forward<decltype(result)>(result));
    append(file, id, thread, s, time, extra...);
  }

  void append(std::string file, long id, long thread, long s, double trial, Extra... extra) {
    times_[std::tuple(file, id, thread, s)].emplace_back(trial, extra...);
  }

  void print(std::ostream& out) const {
    std::size_t n = 4;
    for (auto&& [config, samples] : times_) {
      n = std::max(n, std::get<0>(config).size());
    }

    out << std::setw(n + 2) << std::left << "File";
    out << std::setw(10) << std::left << "Version";
    out << std::setw(10) << std::left << "Threads";
    out << std::setw(10) << std::left << "S";
    out << std::setw(20) << std::left << "Min";
    out << std::setw(20) << std::left << "Avg";
    out << std::setw(20) << std::left << "Max";
    out << "\n";

    for (auto&& [config, samples] : times_) {
      auto [file, id, threads, s] = config;
      auto [min, max, avg] = minmaxavg(samples);

      out << std::setw(n + 2) << std::left << file;
      out << std::setw(10) << std::left << id;
      out << std::setw(10) << std::left << threads;
      out << std::setw(10) << std::left << s;    
      out << std::setw(20) << std::left << std::setprecision(6) << std::fixed << min;
      out << std::setw(20) << std::left << std::setprecision(6) << std::fixed << avg;
      out << std::setw(20) << std::left << std::setprecision(6) << std::fixed << max;
      out << "\n";
    }
  }

 private:
  static auto average(const std::vector<Sample>& times) {
    double total = 0.0;
    for (auto&& sample : times) {
      total += std::get<0>(sample);
    }
    return total / times.size();
  }

  static auto minmax(const std::vector<Sample>& times) {
    return std::apply([](auto... minmax) {
      return std::tuple(std::get<0>(*minmax)...);
    }, std::minmax_element(times.begin(), times.end(), [](auto&& a, auto&& b) {
      return std::get<0>(a) < std::get<0>(b);
    }));
  }

  static auto minmaxavg(const std::vector<Sample>& times) {
    return std::apply([&](auto... minmax) {
      return std::tuple(minmax..., average(times));
    }, minmax(times));
  }
};



}//bench
}//nw::hypergraph
