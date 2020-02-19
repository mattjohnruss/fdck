file(GLOB DRIVERS_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

message(STATUS "Scanning for experimental drivers:")

foreach(driver_src ${DRIVERS_SRC})
    string(REPLACE ".cpp" "" driver ${driver_src})
    message("Found experimental driver \"${driver}\"")
    add_executable(${driver} ${driver_src})
    set_target_properties(${driver} PROPERTIES EXCLUDE_FROM_ALL 1)
    target_link_libraries(${driver} PUBLIC library)
endforeach()
