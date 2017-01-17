cmake_minimum_required(VERSION 2.8.11)
project(HYBRIDSIM)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

list( APPEND CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -fopenmp -march=native  -ffast-math -Wall -O3")

include_directories ("${PROJECT_BINARY_DIR}")

include_directories ("${PROJECT_SOURCE_DIR}/lib")
add_subdirectory (lib)

include_directories ("${PROJECT_SOURCE_DIR}/src")
add_subdirectory (src)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Particles")
add_subdirectory(lib/Particles)

include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")
add_subdirectory(lib/Hydrodynamics)

#find_package(Eigen3 REQUIRED)


