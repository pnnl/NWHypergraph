# Fetch Catch2 header from github
include(Catch2)

macro (nwhy_add_exe EXENAME)
  add_executable(${EXENAME}.exe ${EXENAME}.cpp)
  target_link_libraries(${EXENAME}.exe NW_HYPERGRAPH)
endmacro (nwhy_add_exe)

macro (nwhy_add_test TESTNAME)
  add_executable(${TESTNAME}.exe ${TESTNAME}.cpp)
  target_link_libraries(${TESTNAME}.exe NW_HYPERGRAPH Catch2::Catch2 catch_main)
  target_compile_definitions(${TESTNAME}.exe PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data/")
  target_compile_definitions(${TESTNAME}.exe PRIVATE DATA_FILE="${CMAKE_SOURCE_DIR}/data/line.mmio")
  add_test(NAME ${TESTNAME}_test COMMAND ${TESTNAME}.exe)
  #  catch_discover_tests(${TESTNAME})
endmacro (nwhy_add_test)

add_library(catch_main STATIC catch_main.cpp)
target_link_libraries(catch_main Catch2::Catch2)

# Add Catch2 tests
nwhy_add_test(adjacency_test)

# Example programs (not catch2)
nwhy_add_exe(toplexes_eg)


