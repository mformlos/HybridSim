
include_directories ("${PROJECT_SOURCE_DIR}/lib/Particles")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")
include_directories ("${PROJECT_SOURCE_DIR}/lib/MolecularDynamics")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")

add_executable(hybridsim_auto hybridsim_auto.cpp)
target_link_libraries(hybridsim_auto Random System MPC) 

add_executable(hybridsim hybridsim.cpp)
target_link_libraries(hybridsim Random System Molecule Particle MPC)

add_executable(moltest mol_test.cpp)
target_link_libraries(moltest Random Molecule Particle System)

add_executable(test test.cpp) 
target_link_libraries(test Random System Molecule Particle MPC)

add_executable(thread_random thread_random.cpp)
target_link_libraries(thread_random Random)
