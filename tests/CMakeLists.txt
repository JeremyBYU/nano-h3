add_executable(nanoh3-tests
  test_main.cpp
  nanoh3_test.cpp
)
target_link_libraries(nanoh3-tests PRIVATE nano-h3::nano-h3 h3 doctest::doctest)
target_include_directories(nanoh3-tests PRIVATE ${H3_INCLUDE_DIRS})

add_test(NAME nanoh3-tests COMMAND nanoh3-tests)
