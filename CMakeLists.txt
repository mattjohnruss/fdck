cmake_minimum_required(VERSION 3.9)

project(fd_diff++ C CXX)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

include_directories(${fd_diff++_SOURCE_DIR}/inc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

add_subdirectory(library)
add_subdirectory(drivers)
