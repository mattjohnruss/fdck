file(GLOB DRIVERS_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(driver_src ${DRIVERS_SRC})
    string(REPLACE ".cpp" "" driver ${driver_src})
    message("Found driver \"${driver}\"")
    add_executable(${driver} ${driver_src})
    set_target_properties(${driver} PROPERTIES EXCLUDE_FROM_ALL 1)
    target_link_libraries(${driver} PUBLIC library)
endforeach()
