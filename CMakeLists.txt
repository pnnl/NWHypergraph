##########################################################################
#	CMake Build Rules for the NW hypergraph (v0.3)                              #
##########################################################################
# Basic Usage:                                                           #
#   cmake ..								 #
#   make  								 #
#   make test                                                            #
# For more information about CMake, see http://www.cmake.org             #
##########################################################################
#
# Useful things to know 
# make VERBOSE=1
# cmake --trace
# cmake .. -DCMAKE_CXX_COMPILER=icpx
# cmake .. -DCMAKE_BUILD_TYPE=Release (or Debug)
# cmake .. -DNW_HYPERGRAPH_BUILD_TEST=ON (or OFF)
# cmake .. -DNW_HYPERGRAPH_BUILD_BENCH=ON (or OFF)
# cmake .. -DNW_HYPERGRAPH_BUILD_PYBIND=OFF (or ON)
# cmake .. -DNW_HYPERGRAPH_BUILD_TOOLS=OFF (or ON)
#
##########################################################################
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set(NW_HYPERGRAPH_VERSION 0.3)

if(APPLE)
  set(CMAKE_APPLE_SILICON_PROCESSOR arm64)
endif()

project(
  NW_HYPERGRAPH
  VERSION ${NW_HYPERGRAPH_VERSION}
  LANGUAGES CXX)

# -----------------------------------------------------------------------------
# Extend the module path so we can find our custom modules
# -----------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# -----------------------------------------------------------------------------
# Select a C++ standard
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------------------------------------------------------
# Build types and default flags
# -----------------------------------------------------------------------------
message("Architecture is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL arm64)
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-elide-constructors -fconcepts-diagnostics-depth=3 " CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG " CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Ofast -g -DNDEBUG" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG " CACHE STRING "" FORCE)
else()
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fno-elide-constructors -fconcepts-diagnostics-depth=3 " CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -DNDEBUG " CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-Ofast -g -march=native -DNDEBUG" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -march=native -DNDEBUG " CACHE STRING "" FORCE)
endif()

# Default to Release build
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (default Release)" FORCE)
endif()

# Control compiler-specific flags.
include(CompilerFlags)

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
set(NW_HYPERGRAPH_SANITIZE "" CACHE STRING "-fsanitize compiler options (e.g., =address)")
if (NW_HYPERGRAPH_SANITIZE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=${NW_HYPERGRAPH_SANITIZE}")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=${NW_HYPERGRAPH_SANITIZE}")
endif (NW_HYPERGRAPH_SANITIZE)

option(NW_HYPERGRAPH_EXECUTION_POLICY "with(out) C++ parallel execution policy" ON)
option(NW_HYPERGRAPH_BUILD_TEST "Determines whether to build tests." ON)
option(NW_HYPERGRAPH_BUILD_BENCH "Determines whether to build performance benchmarks." ON)
option(NW_HYPERGRAPH_BUILD_PYBIND "Determines whether to build python binding modules." OFF)
option(NW_HYPERGRAPH_BUILD_TOOLS "Determines whether to build tools." OFF)

add_subdirectory(include)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../NWgr/include)
include(NWGraph)

if (NW_HYPERGRAPH_BUILD_BENCH)
  add_subdirectory(bench)
  include(Date)
  include(Docopt)
  #include(RangeV3)
endif()

if (NW_HYPERGRAPH_BUILD_PYBIND)
  add_subdirectory(python)
endif()

if (NW_HYPERGRAPH_BUILD_TOOLS)
  add_subdirectory(tools)
  include(Docopt)
endif()

if (NW_HYPERGRAPH_BUILD_TEST)
  include(CTest)
  add_subdirectory(test)
endif()
