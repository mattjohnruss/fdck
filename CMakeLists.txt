cmake_minimum_required(VERSION 3.9)

project(fd_diff++ C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Platform ID (CXX): ${CMAKE_CXX_PLATFORM_ID}")

if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --shell-file shell_minimal.html")
endif()

if(CMAKE_CROSSCOMPILING)
    message(STATUS "Cross compiling!")
endif()

include_directories(${fd_diff++_SOURCE_DIR}/inc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wextra -Wpedantic")

# force RelWithDebInfo builds to use -O3 optimisation, which is what Release
# builds use by default
string(REPLACE "-O2" "-O3" CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

if(NOT CMAKE_CROSSCOMPILING)
    #set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
    #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")

    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
        else()
            message(STATUS "Compiling with Clang - disabling LTO!")
        endif()
    endif()
endif()

set(CMAKE_EXE_LINKER_FLAGS "-static")

add_subdirectory(library)
add_subdirectory(drivers)
add_subdirectory(tests)
