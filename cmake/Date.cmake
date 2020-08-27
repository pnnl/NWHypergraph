# -----------------------------------------------------------------------------
# Allow our executables to use docopt.
# -----------------------------------------------------------------------------
include(FetchContent)

FetchContent_Declare(
  date
  GIT_REPOSITORY https://github.com/HowardHinnant/date.git
  GIT_TAG master)

FetchContent_MakeAvailable(date)

include_directories(${date_SOURCE_DIR}/include)