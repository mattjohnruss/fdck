cmake_minimum_required(VERSION 3.9)

project(fd_diff++ C CXX)

message("Build type: ${CMAKE_BUILD_TYPE}")

if(CMAKE_CROSSCOMPILING)
    message("Cross compiling!")
endif(CMAKE_CROSSCOMPILING)

include_directories(${fd_diff++_SOURCE_DIR}/inc)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-libgcc")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -static-libgcc -static-libstdc++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# force RelWithDebInfo builds to use -O3 optimisation, which is what Release
# builds use by default
string(REPLACE "-O2" "-O3" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

if(NOT CMAKE_CROSSCOMPILING)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif(CMAKE_BUILD_TYPE STREQUAL "Release")
endif(NOT CMAKE_CROSSCOMPILING)

#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
set(CMAKE_EXE_LINKER_FLAGS "-static")

add_subdirectory(library)
add_subdirectory(drivers)
add_subdirectory(tests)
