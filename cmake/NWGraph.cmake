# -----------------------------------------------------------------------------
# Allow our executables to use NWGraph library.
# -----------------------------------------------------------------------------
include(FetchContent)

FetchContent_Declare(
  nwgr
  GIT_REPOSITORY https://github.com/NWmath/NWgr.git
  GIT_TAG master
)

FetchContent_MakeAvailable(nwgr)

include_directories(${nwgr_SOURCE_DIR}/include)