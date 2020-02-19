add_library(tests_main tests_main.cpp)

message(STATUS "Looking for Catch2")
if(EXISTS ${fdck_SOURCE_DIR}/vendor/catch2)
    target_include_directories(tests_main PUBLIC ${fdck_SOURCE_DIR}/vendor/catch2/single_include)
    message("Using bundled Catch2 at ${fdck_SOURCE_DIR}/vendor/catch2")
else()
    message(FATAL_ERROR "Catch2 not found in vendor/catch2")
endif()

file(GLOB TESTS_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test_*.cpp)

message(STATUS "Scanning for tests:")

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
