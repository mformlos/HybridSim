
include_directories ("${PROJECT_SOURCE_DIR}/lib/Particles")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")


add_executable(test test.cpp) 
target_link_libraries(test Particle MPC)
