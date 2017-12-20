cmake_minimum_required(VERSION 3.0)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(fd_diff++ C CXX)

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
#find_package(SuperLU REQUIRED)
#find_package(Umfpack REQUIRED)
#find_package(OpenMP REQUIRED)

#include_directories(${SUPERLU_INCLUDES})
include_directories(${fd_diff++_SOURCE_DIR}/inc)

add_executable(main src/main.cpp src/problem.cpp)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DNDEBUG -march=native -std=c++14 -Wall -pedantic -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -march=native -std=c++17 -Wall -pedantic -O3")

if(EXISTS ${fd_diff++_SOURCE_DIR}/vendor/eigen)
    include_directories(vendor/eigen)
else()
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(main Eigen3::Eigen)
endif()

#target_link_libraries(main Eigen3::Eigen umfpack superlu)
