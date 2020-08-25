# -----------------------------------------------------------------------------
# Allow our executables to use docopt.
# -----------------------------------------------------------------------------
include(ExternalProject)

ExternalProject_Add(
  date

  GIT_REPOSITORY    "git@github.com:HowardHinnant/date.git"
  GIT_TAG           master
  GIT_SHALLOW       ON

  BUILD_ALWAYS      OFF
  CMAKE_ARGS        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  INSTALL_COMMAND   ""
)
ExternalProject_Get_property(date SOURCE_DIR)
set(DATE_INCLUDE_DIR ${SOURCE_DIR}/include)
include_directories(${DATE_INCLUDE_DIR})

include(FetchContent)

FetchContent_Declare(
  date
  GIT_REPOSITORY https://github.com/HowardHinnant/date.git
  GIT_TAG master)

FetchContent_MakeAvailable(date)

target_compile_options(date INTERFACE -Wno-deprecated-declarations)