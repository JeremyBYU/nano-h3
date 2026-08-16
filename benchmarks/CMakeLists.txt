add_executable(bench-nanoh3 bench_nanoh3.cpp)
target_link_libraries(bench-nanoh3 PRIVATE nano-h3::nano-h3 h3 nanobench)
target_include_directories(bench-nanoh3 PRIVATE ${H3_INCLUDE_DIRS})
