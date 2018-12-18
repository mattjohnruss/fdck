file(GLOB LIBRARY_SRC src/*.cpp)
file(GLOB LIBRARY_INC inc/*.h)

add_library(library STATIC ${LIBRARY_SRC} ${LIBRARY_INC})

target_include_directories(library PUBLIC inc)

message(STATUS "Looking for Eigen")
if(EXISTS ${fd_diff++_SOURCE_DIR}/vendor/eigen)
    target_include_directories(library PUBLIC ${fd_diff++_SOURCE_DIR}/vendor/eigen)
    message("Using bundled Eigen at ${fd_diff++_SOURCE_DIR}/vendor/eigen")
else()
    message(FATAL_ERROR "Eigen not found in vendor/Eigen")
endif()

message(STATUS "Looking for spdlog")
if(EXISTS ${fd_diff++_SOURCE_DIR}/vendor/spdlog)
    target_include_directories(library PUBLIC ${fd_diff++_SOURCE_DIR}/vendor/spdlog/include)
    message("Using bundled spdlog at ${fd_diff++_SOURCE_DIR}/vendor/spdlog")
else()
    message(FATAL_ERROR "spdlog not found in vendor/spdlog")
endif()
