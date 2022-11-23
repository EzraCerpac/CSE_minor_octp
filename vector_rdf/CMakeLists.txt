cmake_minimum_required(VERSION 3.23)
project(TW3740TU C)

set(CMAKE_C_STANDARD 99)

SET(GCC_COVERAGE_COMPILE_FLAGS "-lm")
SET(GCC_COVERAGE_LINK_FLAGS    "-lm")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})
link_libraries(${GCC_COVERAGE_LINK_FLAGS})

add_executable(rdf_force rdf_force.c)