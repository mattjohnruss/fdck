cmake_minimum_required(VERSION 3.9)

project(fd_diff++ C CXX)

message("Build type: ${CMAKE_BUILD_TYPE}")

if(NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif(NOT CMAKE_CROSSCOMPILING)

include_directories(${fd_diff++_SOURCE_DIR}/inc)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -static-libgcc -static-libstdc++")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")

add_subdirectory(library)
add_subdirectory(drivers)
