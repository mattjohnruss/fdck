find_package(Catch2 REQUIRED)

add_library(tests_main tests_main.cpp)
target_link_libraries(tests_main Catch2::Catch2)

add_executable(test_blah test_blah.cpp)
target_link_libraries(test_blah library tests_main)

#enable_testing()

#include(CTest)
#include(Catch)
#catch_discover_tests(test_blah)
