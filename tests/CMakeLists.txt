find_package(Catch2 REQUIRED)

add_library(tests_main tests_main.cpp)
target_link_libraries(tests_main Catch2::Catch2)

file(GLOB TESTS_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test_*.cpp)

message("Scanning for tests:")

foreach(test_src ${TESTS_SRC})
    string(REPLACE ".cpp" "" test ${test_src})
    message("Found test \"${test}\"")
    add_executable(${test} ${test_src})
    target_link_libraries(${test} PUBLIC library tests_main)
endforeach()

#enable_testing()

#include(CTest)
#include(Catch)
#catch_discover_tests(test_blah)
