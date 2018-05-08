file(GLOB LIBRARY_SRC src/*.cpp)
file(GLOB LIBRARY_INC inc/*.h)

add_library(library STATIC ${LIBRARY_SRC} ${LIBRARY_INC})

target_include_directories(library PUBLIC inc)

if(EXISTS ${fd_diff++_SOURCE_DIR}/vendor/eigen)
    target_include_directories(library PUBLIC ${fd_diff++_SOURCE_DIR}/vendor/eigen)
    message("Using bundled Eigen at ${fd_diff++_SOURCE_DIR}/vendor/eigen")
else()
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    target_link_libraries(library Eigen3::Eigen)
    message("Using system Eigen")
endif()
