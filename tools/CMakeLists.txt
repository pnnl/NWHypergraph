# -----------------------------------------------------------------------------
# Create an interface library for the benchmarks to depend on.
# -----------------------------------------------------------------------------
add_library(tools INTERFACE)
target_sources(tools INTERFACE)
target_include_directories(tools INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(tools INTERFACE NW_HYPERGRAPH docopt date::date)


# -----------------------------------------------------------------------------
# NW_HYPERGRAPH tools
# -----------------------------------------------------------------------------

add_executable(mm2adj.exe mmio_to_hyadj.cpp)
target_link_libraries(mm2adj.exe tools)

add_executable(adj2mm.exe hyadj_to_mmio.cpp)
target_link_libraries(adj2mm.exe tools)

add_executable(hystat.exe hyperstats.cpp)
target_link_libraries(hystat.exe tools)

add_executable(mm2mm.exe mmio_to_mmio.cpp)
target_link_libraries(mm2mm.exe tools)

add_executable(hypergen.exe hypergenerator.cpp)
target_link_libraries(hypergen.exe tools)

add_executable(imdb2mtx.exe IMDB_to_mtx.cpp)
target_link_libraries(imdb2mtx.exe tools xtensor)

add_executable(mtx2bin.exe mtx_to_binary.cpp)
target_link_libraries(mtx2bin.exe tools)