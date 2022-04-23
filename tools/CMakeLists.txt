# -----------------------------------------------------------------------------
# Create an interface library for the benchmarks to depend on.
# -----------------------------------------------------------------------------
add_library(tools INTERFACE)
target_sources(tools INTERFACE)
target_include_directories(tools INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(tools INTERFACE NW_HYPERGRAPH docopt date::date)

include(FetchContent)
set(XTENSOR_USE_TBB 1)

# xsimd and TBB dependencies are automatically
# searched when the following is executed

find_package(xtensor QUIET)

if (NOT xtensor_FOUND)
  message("Installed xtensor not found, using cmake config fetched version")
  FetchContent_Declare(
    xtl
    GIT_REPOSITORY https://github.com/xtensor-stack/xtl
  )

  FetchContent_GetProperties(xtl)
  if(NOT xtl_POPULATED)
    message("  Fetching xtl...")

    FetchContent_Populate(xtl)
    add_subdirectory(${xtl_SOURCE_DIR} ${xtl_BINARY_DIR})
  endif()

  FetchContent_Declare(xtensor
    GIT_REPOSITORY https://github.com/xtensor-stack/xtensor
  )
  FetchContent_GetProperties(xtensor)
  if(NOT xtensor_POPULATED)
    message("  Fetching xtensor...")
    FetchContent_Populate(xtensor)
    add_subdirectory(${xtensor_SOURCE_DIR} ${xtensor_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()
endif()

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