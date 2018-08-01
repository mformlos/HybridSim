
include_directories ("${PROJECT_SOURCE_DIR}/lib/Particles")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Hydrodynamics")
include_directories ("${PROJECT_SOURCE_DIR}/lib/MolecularDynamics")
include_directories ("${PROJECT_SOURCE_DIR}/lib/Analysis")

add_executable(hybridsim_auto hybridsim_auto.cpp)
target_link_libraries(hybridsim_auto Random System MPC) 

add_executable(AFMsim AFMsim.cpp)
target_link_libraries(AFMsim Random System)

add_executable(AFMConstraintSim AFMConstraintSim.cpp)
target_link_libraries(AFMConstraintSim Random System)

add_executable(hybridsim hybridsim.cpp)
target_link_libraries(hybridsim Random System Molecule Particle MPC)

add_executable(moltest mol_test.cpp)
target_link_libraries(moltest Random Molecule Particle System)

add_executable(LangevinTest Langevin_test.cpp)
target_link_libraries(LangevinTest Random Molecule Particle System)

add_executable(test test.cpp) 
target_link_libraries(test Random System Molecule Particle MPC)

add_executable(thread_random thread_random.cpp)
target_link_libraries(thread_random Random)

add_executable(copy_test copy_test.cpp) 
target_link_libraries(copy_test Particle)

add_executable(evaltest evaltest.cpp) 

add_executable(msd msd.cpp)
target_link_libraries(msd Molecule)

add_executable(rgyr_corr rgyr_corr.cpp)

add_executable(rgyr_corr_multiple rgyr_corr_multiple.cpp)

add_executable(rgyr_corr_multiple_stepfile rgyr_corr_multiple_stepfile.cpp)

add_executable(evec_corr evec_corr.cpp)

add_executable(stress_tensor stress_tensor.cpp)
target_link_libraries(stress_tensor Molecule)

