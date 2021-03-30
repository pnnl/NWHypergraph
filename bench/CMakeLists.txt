# -----------------------------------------------------------------------------
# The date library from MIT provides C++20 std date functionality in C++17
# -----------------------------------------------------------------------------
find_package(Filesystem REQUIRED Final Experimental)
message(STATUS ${CXX_FILESYSTEM_IS_EXPERIMENTAL})
add_definitions(-Wno-deprecated -Wno-deprecated-declarations)

# -----------------------------------------------------------------------------
# Generate the configuration files used for logging, and -I its built path.
# -----------------------------------------------------------------------------
get_filename_component(NW_HYPERGRAPH_CXX_COMPILER ${CMAKE_CXX_COMPILER} NAME)
configure_file(config.h.in config.h)
#include_directories(${CMAKE_CURRENT_BINARY_DIR})

# -----------------------------------------------------------------------------
# Create an interface library for the benchmarks to depend on.
# -----------------------------------------------------------------------------
add_library(bench INTERFACE)
target_sources(bench INTERFACE
  common.hpp
  Log.hpp
  config.h)
target_include_directories(bench INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(bench INTERFACE NW_HYPERGRAPH docopt date::date std::filesystem)


# -----------------------------------------------------------------------------
# NW_HYPERGRAPH drivers
# -----------------------------------------------------------------------------
add_executable(hybfs.exe hyperbfs.cpp)
target_link_libraries(hybfs.exe bench)

add_executable(relabelbfs.exe relabelbfs.cpp)
target_link_libraries(relabelbfs.exe bench)

add_executable(hycc.exe hypercc.cpp)
target_link_libraries(hycc.exe bench)

add_executable(relabelcc.exe relabelcc.cpp)
target_link_libraries(relabelcc.exe bench)

add_executable(scc.exe soverlapcc.cpp)
target_link_libraries(scc.exe bench)

add_executable(sbfs.exe soverlapbfs.cpp)
target_link_libraries(sbfs.exe bench)

add_executable(ssssp.exe soverlapsssp.cpp)
target_link_libraries(ssssp.exe bench)

add_executable(sbc.exe soverlapbc.cpp)
target_link_libraries(sbc.exe bench)

add_executable(top.exe toplexes.cpp)
target_link_libraries(top.exe bench)

