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
#include "io/mmio.hpp"
#include "containers/compressed_hy.hpp"
#include "containers/edge_list_hy.hpp"

using namespace nw::graph;

namespace nw::hypergraph {
namespace bench {

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

/*
 * This graph reader can read either mtx or hygra (adjacency) file.
 **/
template<directedness edge_directedness = directed, typename... Attributes>
auto graph_reader(std::string file) {
  auto aos_a = load_graph<edge_directedness, Attributes...>(file);
  if (0 == aos_a.size()) {
    auto&&[hyperedges, hypernodes] = load_adjacency<Attributes...>(file);
    std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
    return std::tuple(hyperedges, hypernodes, std::vector<vertex_id_t>());
  }
  nw::graph::adjacency<0, Attributes...> hyperedges(aos_a);
  nw::graph::adjacency<1, Attributes...> hypernodes(aos_a);

  std::cout << "num_hyperedges = " << hyperedges.size()
            << " num_hypernodes = " << hypernodes.size() << std::endl;
  return std::tuple(hyperedges, hypernodes, std::vector<vertex_id_t>());
}

/*
 * This graph reader can read either mtx or hygra (adjacency) file.
 * And can relabel the ids of either hyperedges or hypernodes by degree
 * in the direction of ascending (default) or descending.
 **/
template<directedness edge_directedness = directed, typename... Attributes>
auto graph_reader_relabel(std::string file, int idx, std::string direction) {
  auto aos_a = load_graph<edge_directedness, Attributes...>(file);
  if (0 == aos_a.size()) {
    auto&& [hyperedges, hypernodes] = load_adjacency<Attributes...>(file);
    // Run relabeling. This operates directly on the incoming edglist.
    if (-1 != idx) {
      auto&& iperm = nw::hypergraph::relabel_by_degree(
          hyperedges, hypernodes, idx, direction);
        std::cout << "num_hyperedges = " << hyperedges.size()
              << " num_hypernodes = " << hypernodes.size() << std::endl;
      return std::tuple(hyperedges, hypernodes, iperm);
    }
  }
  std::vector<vertex_id_t> iperm;
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
  nw::graph::adjacency<0, Attributes...> hyperedges(aos_a);
  nw::graph::adjacency<1, Attributes...> hypernodes(aos_a);

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
template<directedness edge_directedness = directed, typename... Attributes>
auto graph_reader_adjoin_and_relabel(std::string file, const int idx, std::string direction, 
size_t& nrealedges, size_t& nrealnodes) {
  auto aos_a =
      read_mm_adjoin<edge_directedness, Attributes...>(file, nrealedges, nrealnodes);
  std::vector<vertex_id_t> iperm, perm;
  if (0 == aos_a.size()) {
    //we get adjoin graph and its transpose
    auto&& [g, gt] = read_and_adjoin_adj_hypergraph_pair(file, nrealedges, nrealnodes);
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
    return std::tuple(g, gt, iperm);
  }
  else {
    // Run relabeling. This operates directly on the incoming edglist.
    if (-1 != idx) {
      //since the adjoin graph is symmetric, when we need to relabel it
      //we need to operate on both column 0 and column 1
      std::cout << "relabeling edge_list by degree..." << std::endl;
      nw::util::ms_timer t("relabel_by_degree");
      if (1 == idx) {
        perm = aos_a.template perm_by_degree<1>(direction);
        iperm = aos_a.relabel(perm);
      }
      else {
        perm = aos_a.template perm_by_degree<0>(direction);
        iperm = aos_a.relabel(perm);
      }
      t.stop();
      std::cout << t << std::endl;
    }

    nw::graph::adjacency<0, Attributes...> g(aos_a);
    nw::graph::adjacency<1, Attributes...> g_t(aos_a);
    std::cout << "num_realedges = " << nrealedges
            << " num_realnodes = " << nrealnodes << std::endl;    
    return std::tuple(g, g_t, iperm);
  }
}

/*
 * This graph reader can read either mtx or hygra (adjacency) file 
 * in adjoin fashion.
 * Adjoin will increment the ids of the smaller set of the hyperedges or hypernodes,
 * such that the ids of both are now contiguous.
 **/
template <directedness edge_directedness = directed, typename... Attributes>
auto graph_reader_adjoin(std::string file, size_t& nrealedges,
                         size_t& nrealnodes) {
  auto aos_a = read_mm_adjoin<edge_directedness, Attributes...>(
      file, nrealedges, nrealnodes);
  if (0 == aos_a.size()) {
    //after adjoin, we have adjoin graph and its transpose
    auto&& [g, g_t] =
        read_and_adjoin_adj_hypergraph_pair(file, nrealedges, nrealnodes);
    std::cout << "num_realedges = " << nrealedges
            << " num_realnodes = " << nrealnodes << std::endl;
      std::cout << "num_hyperedges = " << g.size()
            << " num_hypernodes = " << g_t.size() << std::endl;
    return std::tuple(g, g_t, std::vector<vertex_id_t>());
  }
  nw::graph::adjacency<0, Attributes...> g(aos_a);
  nw::graph::adjacency<1, Attributes...> g_t(aos_a);
  std::cout << "num_realedges = " << nrealedges
            << " num_realnodes = " << nrealnodes << std::endl;
  std::cout << "num_hyperedges = " << g.size()
            << " num_hypernodes = " << g_t.size() << std::endl;
  return std::tuple(g, g_t, std::vector<vertex_id_t>());
}

template<directedness edge_directedness = directed, typename... Attributes>
auto graph_reader(std::string file, int idx, std::string direction, bool adjoin, 
size_t& nrealedges, size_t& nrealnodes) {
  if (adjoin) {
    //adjoin hypergraph will be forced to be symmetric (undirected)
    if (-1 != idx)
      return graph_reader_adjoin_and_relabel<undirected, Attributes...>(file, idx, direction, 
      nrealedges, nrealnodes);
    else
      return graph_reader_adjoin<undirected, Attributes...>(file, nrealedges, nrealnodes);
  }
  else {
    if (-1 != idx)
      return graph_reader_relabel<edge_directedness, Attributes...>(file, idx, direction);
    else
      return graph_reader<edge_directedness, Attributes...>(file);
  }
}

template<directedness edge_directedness = directed, typename... Attributes>
auto graph_reader(std::string file, int idx, std::string direction) {
  if (-1 != idx)
    return graph_reader_relabel<edge_directedness, Attributes...>(file, idx, direction);
  else
    return graph_reader<edge_directedness, Attributes...>(file);
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
