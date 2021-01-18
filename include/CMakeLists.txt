# -----------------------------------------------------------------------------
# Define the primary NW_HYPERGRAPH target.
# -----------------------------------------------------------------------------
#
# This target keeps track of all of the header sources for NW_HYPERGRAPH, and sets up all
# of the associated flags and dependencies such that we can correctly compile
# and link executables that are dependent on NW_HYPERGRAPH. Ultimately this target will
# control the install target too, once we have one.
#
# Executable targets should be defined with `target_link_libraries(exec NW_HYPERGRAPH)`
# and then they'll have the appropriate `-I`, `-D`, and `-l` stuff set
# automagically though the transitive dependency infrastructure.

# @note
#   If NW_HYPERGRAPH starts to have .cpp sources, change this from an interface library to
#   a STATIC (or DYNAMIC) library, and add those sources here. Change other
#   instances of INTERFACE to PUBLIC.
# @note
add_library(NW_HYPERGRAPH INTERFACE)

# Extend this as needed.
target_sources(NW_HYPERGRAPH
  INTERFACE
  containers/edge_list_hy.hpp
  containers/compressed_hy.hpp
  io/hypergraph_io.hpp
  io/mmio.hpp
  s_overlap.hpp
  algorithms/relabel_x.hpp
  algorithms/hyper_connected_components.hpp
  algorithms/hyper_breadth_first_search.hpp
  algorithms/s_connected_components.hpp
  algorithms/s_breadth_first_search.hpp
  algorithms/toplexes.hpp
  )

target_include_directories(NW_HYPERGRAPH INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# -----------------------------------------------------------------------------
# Check to see if our C++ compiler supports Parallelism TS. If it does than
# nothing weird will happen. If it doesn't we define the `NEED_` macros and that
# will pull in `detail/numeric.hpp`.
# -----------------------------------------------------------------------------
include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "#include<numeric>
    int main() {
        int* ptr;
        std::inclusive_scan(ptr, ptr, ptr);
    }"
  HAVE_INCLUSIVE_SCAN)

if (NOT HAVE_INCLUSIVE_SCAN)
  target_compile_definitions(NW_HYPERGRAPH INTERFACE NW_HYPERGRAPH_NEED_INCLUSIVE_SCAN)
endif ()

check_cxx_source_compiles(
  "#include<numeric>
    int main() {
        int* ptr;
        std::exclusive_scan(ptr, ptr, ptr, 0);
    }"
  HAVE_EXCLUSIVE_SCAN)

if (NOT HAVE_EXCLUSIVE_SCAN)
  target_compile_definitions(NW_HYPERGRAPH INTERFACE NW_HYPERGRAPH_NEED_EXCLUSIVE_SCAN)
endif ()

check_cxx_source_compiles(
  "#include<numeric>
    int main() {
        int* ptr;
        std::reduce(ptr, ptr, 0);
    }"
  HAVE_REDUCE)

if (NOT HAVE_REDUCE)
  target_compile_definitions(NW_HYPERGRAPH INTERFACE NW_HYPERGRAPH_NEED_REDUCE)
endif ()

# -----------------------------------------------------------------------------
# Handle requests for parallel execution
# -----------------------------------------------------------------------------
if (NW_HYPERGRAPH_EXECUTION_POLICY)
  set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(TBB REQUIRED)
  find_package(Threads REQUIRED)
  target_compile_definitions(NW_HYPERGRAPH INTERFACE EXECUTION_POLICY)
  target_link_libraries(NW_HYPERGRAPH INTERFACE TBB::tbb Threads::Threads)
  if (NW_HYPERGRAPH_USE_TBBMALLOC)
    target_link_libraries(NW_HYPERGRAPH INTERFACE TBB::tbbmalloc) 
  endif ()
endif()

# -----------------------------------------------------------------------------
# Handle Apple-specific things
# -----------------------------------------------------------------------------
if (APPLE)
  if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    target_compile_definitions(NW_HYPERGRAPH INTERFACE _Atomic=volatile)
  endif()
endif()
