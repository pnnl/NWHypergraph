//
// This file is part of the Graph Standard Library (aka nw::graph aka NWGraph)
// (c) 2020 Pacific Northwest National Laboratory
//
// Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
// https://creativecommons.org/licenses/by-nc-sa/4.0/
//
// Author: Xu Tony Liu
//

#include <execution>
#include <unordered_set>
#include <edge_list.hpp>
#include <algorithms/connected_components.hpp>
#include <docopt.h>
#include "Log.hpp"
#include "common.hpp"
#include "io/mmio.hpp"
#include "algorithms/relabel_x.hpp"



using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;
using namespace nw::graph;

static constexpr const char USAGE[] =
    R"(hyccrelabel.exe: nw::graph hypergraph connected components benchmark driver.
  Usage:
      hyccrelabel.exe (-h | --help)
      hyccrelabel.exe [-f FILE...] [--version ID...] [-n NUM] [--succession STR] [--relabel] [--clean] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

  Options:
      -h, --help            show this screen
      --version ID          algorithm version to run [default: 0]
      -f FILE               input file paths (can have multiples)
      -n NUM                number of trials [default: 1]
      --relabel             relabel the graph or not
      -c, --clean           clean the graph or not
      --direction DIR       graph relabeling direction - ascending/descending [default: descending]
      --succession STR      successor/predecessor [default: successor]
      --log FILE            log times to a file
      --log-header          add a header to the log file
      -d, --debug           run in debug mode
      -v, --verify          verify results
      -V, --verbose         run in verbose mode
)";

int main(int argc, char* argv[]) {
  std::vector<std::string> strings(argv + 1, argv + argc);
  auto                     args = docopt::docopt(USAGE, strings, true);

  bool verify  = args["--verify"].asBool();
  bool verbose = args["--verbose"].asBool();
  bool debug   = args["--debug"].asBool();
  long trials  = args["-n"].asLong() ?: 1;

  std::vector<long> ids     = parse_ids(args["--version"].asStringList());
  std::vector<long> threads = parse_n_threads(args["THREADS"].asStringList());

  std::vector<std::string> files;
  for (auto&& file : args["-f"].asStringList()) {
    files.emplace_back(file);
  }

  Times<bool> times;

  // Appease clang.
  //
  // These are captured in lambdas later, and if I use structured bindings
  // they have to be listed as explicit captures (this is according to the 17
  // standard). That's a little bit noisy where it happens, so I just give
  // them real symbols here rather than the local bindings.
  for (auto&& file : files) {
    auto reader = [&](std::string file, bool verbose, size_t& nrealedges, size_t& nrealnodes) {
      //auto aos_a   = load_graph<directed>(file);
      auto aos_a   = read_mm_relabeling<nw::graph::directed>(file, nrealedges, nrealnodes);
      if (0 == aos_a.size()) {
        auto&& g = read_and_relabel_adj_hypergraph(file, nrealedges, nrealnodes);
        return std::tuple(g, nw::graph::adjacency<1>(0, 0));
      }
      // Run relabeling. This operates directly on the incoming edglist.
      if (args["--relabel"].asBool()) {
        auto degrees = aos_a.degrees();
        aos_a.relabel_by_degree<0>(args["--direction"].asString(), degrees);
      }
      // Clean up the edgelist to deal with the normal issues related to
      // undirectedness.
      if (args["--clean"].asBool()) {
        aos_a.swap_to_triangular<0>(args["--succession"].asString());
        aos_a.lexical_sort_by<0>();
        aos_a.uniq();
        aos_a.remove_self_loops();
      }

      nw::graph::adjacency<0> g(aos_a);
      nw::graph::adjacency<1> g_t(aos_a);
      if (verbose) {
        g_t.stream_stats();
        g.stream_stats();
      }
      std::cout << "num_hyperedges = " << nrealedges << " num_hypernodes = " << nrealnodes << std::endl;
      return std::tuple(g, g_t);
    };
    size_t num_realedges, num_realnodes;
    auto&& [g, g_t]     = reader(file, verbose, num_realedges, num_realnodes);

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        auto verifier = [&](auto&& result) {
          auto&& [N, E] = result;
          std::unordered_set<vertex_id_t> uni_comps(E.begin(), E.end());
          std::cout << uni_comps.size() << " components found" << std::endl;
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, std::forward<decltype(op)>(op), verifier, true); };
        using Graph = nw::graph::adjacency<0>;
        using Transpose = nw::graph::adjacency<1>;
        using ExecutionPolicy = decltype(std::execution::par_unseq);
        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { 
                auto lpf = nw::graph::ccv1<Graph>;
                using LabelPropagationF = decltype(lpf);
                return nw::hypergraph::relabel_x_parallel<ExecutionPolicy, LabelPropagationF, vertex_id_t>(std::execution::par_unseq, num_realedges, num_realnodes, lpf, g); });
              break;
            case 1:
              record([&] { 
                auto af = nw::graph::Afforest<Graph, Transpose>;
                using AfforestF = decltype(af);
                return nw::hypergraph::relabel_x_parallel<ExecutionPolicy, AfforestF, vertex_id_t>(std::execution::par_unseq, num_realedges, num_realnodes, af, g, g_t, 2); });
              break;
            case 2:
              record([&] { 
                auto lpf = nw::graph::ccv5<Graph>;
                using LabelPropagationF = decltype(lpf);
                return nw::hypergraph::relabel_x_parallel<ExecutionPolicy, LabelPropagationF, vertex_id_t>(std::execution::par_unseq, num_realedges, num_realnodes, lpf, g); });
              break;
            case 3:
              record([&] { 
                auto svf = nw::graph::compute_connected_components_v2<Graph>;
                using SVF = decltype(svf);
                return nw::hypergraph::relabel_x_parallel<ExecutionPolicy, SVF, vertex_id_t>(std::execution::par_unseq, num_realedges, num_realnodes, svf, g); });
              break;
            default:
              std::cout << "Unknown version v" << id << "\n";
          }
        }
      }
    }
  }

  times.print(std::cout);

  if (args["--log"]) {
    auto file   = args["--log"].asString();
    bool header = args["--log-header"].asBool();
    log("cc", file, times, header, "Time(s)");
  }

  return 0;
}
