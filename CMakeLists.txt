# 3.21, not 3.20: PROJECT_IS_TOP_LEVEL was added in 3.21. On 3.20 it expands to
# nothing, `option(X "desc" ${PROJECT_IS_TOP_LEVEL})` collapses to the legal
# two-argument form, and every option below silently defaults OFF -- a standalone
# build would configure cleanly and build zero tests.
cmake_minimum_required(VERSION 3.21)
project(nano-h3 VERSION 0.2.0 LANGUAGES CXX)

add_library(nano-h3 INTERFACE)
add_library(nano-h3::nano-h3 ALIAS nano-h3)
target_include_directories(nano-h3 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(nano-h3 INTERFACE cxx_std_17)

# OFF whenever nano-h3 is consumed as a subproject, so a FetchContent consumer
# never inherits h3/doctest/nanobench declarations from us.
option(NANOH3_BUILD_TESTS      "Build the differential test suite" ${PROJECT_IS_TOP_LEVEL})
option(NANOH3_BUILD_BENCHMARKS "Build the nanobench comparison"    ${PROJECT_IS_TOP_LEVEL})
option(NANOH3_INSTALL          "Generate install and export rules" ${PROJECT_IS_TOP_LEVEL})

if(NANOH3_BUILD_TESTS OR NANOH3_BUILD_BENCHMARKS)
  include(cmake/dependencies.cmake)
endif()

if(NANOH3_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(NANOH3_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

# Gated because add_subdirectory(... EXCLUDE_FROM_ALL) suppresses a subproject's
# install rules but FetchContent_MakeAvailable does not. Without the gate the two
# consumption paths disagree about what `cmake --install` emits, and a
# FetchContent consumer finds our headers in its own install tree.
if(NANOH3_INSTALL)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(TARGETS nano-h3 EXPORT nano-h3-targets)
  install(EXPORT nano-h3-targets
          FILE nano-h3-targets.cmake
          NAMESPACE nano-h3::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nano-h3)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/nano-h3-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/nano-h3-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nano-h3)

  # ARCH_INDEPENDENT: without it the generated file compares pointer sizes and
  # rejects valid consumers of a header-only library. SameMinorVersion, not
  # SameMajorVersion: at 0.x the latter would accept 0.2.0 for a request of
  # 0.1.0, which is exactly the silent skew the version macros exist to catch.
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/nano-h3-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
    ARCH_INDEPENDENT)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nano-h3-config.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/nano-h3-config-version.cmake
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nano-h3)
endif()
