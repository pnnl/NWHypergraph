# -----------------------------------------------------------------------------
# The date library from MIT provides C++20 std date functionality in C++17
# -----------------------------------------------------------------------------
#find_package(Filesystem REQUIRED Final Experimental)
#message(STATUS ${CXX_FILESYSTEM_IS_EXPERIMENTAL})
#add_definitions(-Wno-deprecated -Wno-deprecated-declarations)

# -----------------------------------------------------------------------------
# Generate the configuration files used for logging, and -I its built path.
# -----------------------------------------------------------------------------
get_filename_component(NW_HYPERGRAPH_CXX_COMPILER ${CMAKE_CXX_COMPILER} NAME)
configure_file(config.h.in config.h)
#include_directories(${CMAKE_CURRENT_BINARY_DIR})

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
# Create an interface library for the benchmarks to depend on.
# -----------------------------------------------------------------------------
add_library(bench INTERFACE)
target_sources(bench INTERFACE
  common.hpp
  Log.hpp
  #config.h
  )
target_include_directories(bench INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(bench INTERFACE NW_HYPERGRAPH docopt date::date)


# -----------------------------------------------------------------------------
# NW_HYPERGRAPH drivers
# -----------------------------------------------------------------------------
#set policy to suppress warnings
cmake_policy(SET CMP0115 NEW)

add_executable(hybfs.exe hyperbfs.cpp)
target_link_libraries(hybfs.exe bench)

add_executable(adjoinbfs.exe adjoinbfs.cpp)
target_link_libraries(adjoinbfs.exe bench)

add_executable(hycc.exe hypercc.cpp)
target_link_libraries(hycc.exe bench)

add_executable(adjoincc.exe adjoincc.cpp)
target_link_libraries(adjoincc.exe bench)

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

add_executable(imdb.exe imdb.cpp)
target_link_libraries(imdb.exe bench xtensor)