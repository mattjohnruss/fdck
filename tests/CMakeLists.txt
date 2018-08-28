find_package(Catch2 REQUIRED)

add_library(tests_main tests_main.cpp)
target_link_libraries(tests_main Catch2::Catch2)

add_executable(test_config test_config.cpp)
target_link_libraries(test_config library tests_main)

#enable_testing()

#include(CTest)
#include(Catch)
#catch_discover_tests(test_blah)
