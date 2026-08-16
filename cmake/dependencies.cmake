# Fetched only when NANOH3_BUILD_TESTS or NANOH3_BUILD_BENCHMARKS is ON, which
# is only when nano-h3 is the top-level project. The library itself needs none
# of this: nanoh3.hpp includes <cmath>, <cstdint> and its own .inl, nothing else.
include(FetchContent)

# H3 is C. Enabled here rather than in project(LANGUAGES ...) so that consuming
# the header never requires a C compiler: this file is only included when
# nano-h3 is the top-level project. Without it CMake fails at generate time with
# "No known features for C compiler", because the top-level scope has no C
# feature table for the test target to resolve against.
enable_language(C)

# H3 is the test oracle, not a dependency of the library. Same pin the header
# was transcribed from; changing it invalidates the bit-identity contract.
FetchContent_Declare(
  h3
  GIT_REPOSITORY https://github.com/uber/h3.git
  GIT_TAG        v4.1.0
  GIT_SHALLOW    TRUE
)
# These are h3's own option names, checked against its CMakeLists: they are
# BUILD_BENCHMARKS and friends, NOT H3_BUILD_BENCHMARKS. Setting the H3_-
# prefixed spelling does nothing at all and leaves h3 building ~30 benchmark,
# filter and fuzzer binaries we never run.
set(BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(BUILD_FUZZERS OFF CACHE BOOL "" FORCE)
set(BUILD_FILTERS OFF CACHE BOOL "" FORCE)
set(BUILD_GENERATORS OFF CACHE BOOL "" FORCE)
set(ENABLE_DOCS OFF CACHE BOOL "" FORCE)
set(ENABLE_LINTING OFF CACHE BOOL "" FORCE)
set(ENABLE_FORMAT OFF CACHE BOOL "" FORCE)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
FetchContent_MakeAvailable(h3)

# BOTH directories are required, not one for tidiness: h3api.h is generated into
# the binary dir, so the source dir alone gives
#   fatal error: h3api.h: No such file or directory
set(H3_INCLUDE_DIRS
  ${h3_SOURCE_DIR}/src/h3lib/include
  ${h3_BINARY_DIR}/src/h3lib/include
  CACHE INTERNAL "H3 include directories"
)

if(NANOH3_BUILD_TESTS)
  FetchContent_Declare(
    doctest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG        v2.4.11
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(doctest)
endif()

if(NANOH3_BUILD_BENCHMARKS)
  FetchContent_Declare(
    nanobench
    GIT_REPOSITORY https://github.com/martinus/nanobench.git
    GIT_TAG        v4.3.11
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(nanobench)
endif()
