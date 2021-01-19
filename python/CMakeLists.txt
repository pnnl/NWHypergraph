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
add_library(dep INTERFACE)
target_sources(dep INTERFACE
  config.h)
target_include_directories(dep INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(dep INTERFACE NW_HYPERGRAPH)

# -----------------------------------------------------------------------------
# Config Python and Pybind11 for the benchmarks to depend on.
# -----------------------------------------------------------------------------
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" OFF)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/python/pybind11/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

# add pybind11 dependency
find_package(Python COMPONENTS Interpreter Development)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pybind11/include)

# -----------------------------------------------------------------------------
# NW_HYPERGRAPH pybind11 module to verify that we're correctly generating source lists.
# -----------------------------------------------------------------------------

#add_library(hello MODULE hello.cpp)
pybind11_add_module(nwhy nwhy.cpp)
target_link_libraries(nwhy PRIVATE dep pybind11::module)
#target_link_libraries(hello PRIVATE dep 
#pybind11::headers pybind11::module pybind11::lto pybind11::windows_extras)

#pybind11_extension(hello)
#pybind11_strip(hello)

#set_target_properties(hello PROPERTIES CXX_VISIBILITY_PRESET "hidden"
#                                       CUDA_VISIBILITY_PRESET "hidden")