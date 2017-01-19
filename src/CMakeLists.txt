
include_directories ("${PROJECT_SOURCE_DIR}/lib/Particles")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")


add_executable(test test.cpp) 
target_link_libraries(test Random Molecule Particle MPC)

add_executable(thread_random thread_random.cpp)
target_link_libraries(thread_random Random)
