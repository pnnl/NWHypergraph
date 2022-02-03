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
#include <nwgraph/edge_list.hpp>
#include <nwgraph/algorithms/connected_components.hpp>
#include <nwgraph/experimental/algorithms/connected_components.hpp>
#include <docopt.h>
#include "Log.hpp"
#include "common.hpp"
#include "io/mmio_hy.hpp"
#include "algorithms/adjoin_x.hpp"


using namespace nw::hypergraph::bench;
using namespace nw::hypergraph;
using namespace nw::graph;

static constexpr const char USAGE[] =
    R"(adjoincc.exe: nw::graph hypergraph connected components benchmark driver.
  Usage:
      adjoincc.exe (-h | --help)
      adjoincc.exe [-f FILE...] [--version ID...] [-n NUM] [--succession STR] [--relabel] [--clean] [--direction DIR] [-dvV] [--log FILE] [--log-header] [THREADS]...

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
    size_t num_realedges = 0, num_realnodes = 0;
    using vertex_id_t = uint32_t;
    auto&& [g, g_t, iperm]    = graph_reader_adjoin<directedness::undirected, vertex_id_t>(file, num_realedges, num_realnodes);
    std::cout << "num_hyperedges = " << num_realedges << " num_hypernodes = " << num_realnodes << std::endl;
    std::cout << "size of the merged adjacency = " << g.size() << std::endl;

    for (auto&& thread : threads) {
      auto _ = set_n_threads(thread);
      for (auto&& id : ids) {
        auto verifier = [&](auto&& result) {
          auto&& [N, E] = result;
          std::unordered_set<vertex_id_t> uni_comps(N.begin(), N.end());
          std::cout << uni_comps.size() << " components found" ;
          if (verbose) {
            std::cout << ":" << std::endl;
            for(auto it = uni_comps.begin(); it != uni_comps.end(); ++it)
              std::cout << *it << std::endl;
            for(size_t i = 0, e = E.size(); i < e; ++i)
              std::cout << i << " " << E[i] << std::endl;
          }
          std::cout << std::endl;
        };

        auto record = [&](auto&& op) { times.record(file, id, thread, std::forward<decltype(op)>(op), verifier, true); };
        using Graph = nw::graph::adjacency<0>;
        using Transpose = nw::graph::adjacency<1>;
        using ExecutionPolicy = decltype(std::execution::par_unseq);
        for (int j = 0, e = trials; j < e; ++j) {
          switch (id) {
            case 0:
              record([&] { 
                auto l = nw::graph::Afforest<Graph, Transpose>(g, g_t, 2);
                return splitLabeling<ExecutionPolicy, vertex_id_t>(std::execution::par_unseq, l, num_realedges, num_realnodes);
              });
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
            case 4:
              record([&] { 
                auto lpf = nw::graph::ccv1<Graph, vertex_id_t>;
                using LabelPropagationF = decltype(lpf);
                return nw::hypergraph::relabel_x< LabelPropagationF, vertex_id_t>(num_realedges, num_realnodes, lpf, g); });
              break;
            case 5:
              record([&] { 
                auto af = nw::graph::Afforest<Graph, Transpose>;
                using AfforestF = decltype(af);
                return nw::hypergraph::relabel_x<AfforestF, vertex_id_t>(num_realedges, num_realnodes, af, g, g_t, 2); });
              break;
            case 6:
              record([&] { 
                auto lpf = nw::graph::lpcc<ExecutionPolicy, Graph>;
                using LabelPropagationF = decltype(lpf);
                int num_bins = 32;
                return nw::hypergraph::relabel_x_parallel<ExecutionPolicy, LabelPropagationF, vertex_id_t>(std::execution::par_unseq, num_realedges, num_realnodes, lpf, std::execution::par_unseq, g, num_bins); });
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
