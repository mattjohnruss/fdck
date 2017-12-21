file(GLOB DRIVERS_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(driver_src ${DRIVERS_SRC})
    string(REPLACE ".cpp" "" driver ${driver_src})
    add_executable(${driver} ${driver_src})
    target_link_libraries(${driver} PUBLIC library)
endforeach(driver_src)
