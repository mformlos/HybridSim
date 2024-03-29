#include <sys/time.h>
#include <omp.h>
#include "Particle.h"
#include "Molecule.h"
#include "MPC.h"
#include "Velocity_Profile.h"
#include "System.h"

int main() {
    Particle test_part; 
    unsigned Steps = 500; 
  
    int tid, proc_num; 
    unsigned Lx = 100, Ly = 100, Lz = 100; 
    
    System sys_test(Lx,Ly,Lz, 0.0); 
    
    sys_test.addMolecules("/home/formanek/HYBRIDSIM/input/SCNP-0-chain"); 
    sys_test.addLinks("/home/formanek/HYBRIDSIM/input/SCNP-0-bonds"); 
    sys_test.initializePositions("/home/formanek/HYBRIDSIM/input/SCNP-0-config"); 
    sys_test.initializeVelocitiesRandom(1.0); 
    
    
    test_part.Velocity(0) = 1.0; 
    
    MPC test_mpc(Lx, Ly, Lz, 10, 1.0, 0.0);
    test_mpc.initializeRandom();   
    std::cout << "Temperature after initialization: " << test_mpc.virtualTemperature() << std::endl;
    Molecule mol(1); 
    
    mol.Monomers.front().Position << Rand::real_uniform()*Lx, Rand::real_uniform()*Ly, Rand::real_uniform()*Lz; 
    //mol.Monomers.front().Velocity << Rand::real_uniform() - 0.5, Rand::real_uniform() - 0.5, Rand::real_uniform() - 0.5; 
    
    unsigned N_s {sys_test.NumberOfParticles()};
    std::cout << "number of solute particles: " << N_s << std::endl; 
    test_mpc.initializeSoluteVector(N_s); 
    
    std::cout << "part pos: " << mol.Monomers.front().Position.transpose() << std::endl; 
    std::cout << "part vel: " << mol.Monomers.front().Velocity.transpose() << std::endl; 
    VelocityProfile vel_prof{0.2}; 
    
    ofstream temp_file{}; 
    temp_file.open("temperature.dat", ios::out | ios::trunc);     
    timeval start, end; 

    gettimeofday(&start, NULL);

    /*#pragma omp parallel 
    {
    for (unsigned n = 0; n < Steps; n++) {
        #pragma omp single 
        {
            //if (n%5==0) test_mpc.sortVector(); 
            test_mpc.stream(0.1); 
            try {
                test_mpc.sort(); 
            }
            catch(...) {
                std::cout << "exception caught"; 
                //return EXIT_FAILURE; 
            }
        }
        #pragma omp for schedule(static) 
        for (unsigned i = 0; i < test_mpc.NumberOfCells; i++) {
            Vector3d COMVel {test_mpc.CenterOfMassVelocity(i)};
            test_mpc.collide(i, COMVel); 
        }
        std::cout << "current temperature: " << test_mpc.virtualTemperature() << std::endl; 
    }
    }*/
    #pragma omp parallel private(tid)
    {
        #ifdef _OPENMP
        tid = omp_get_thread_num(); 
        if (tid == 0) {
            proc_num = omp_get_num_procs();
            std::cout << "You are using " << proc_num << " processes. " << std::endl;  
        }
        #endif
        std::cout << "Hello from process: " << tid << std::endl;
        Rand::seed(tid);
        Rand::warmup(10000); 
        
        for (unsigned n = 0; n < Steps; n++) {
            #pragma omp single 
            {
                test_mpc.shiftGrid();
                test_mpc.updateBoxShift(0.1);  
                test_mpc.getSolute(sys_test.Molecules);    
            }
        
            #pragma omp for schedule(static) 
            for (unsigned part = 0; part < test_mpc.NumberOfParticles; part++) {
                test_mpc.streamPlusCellAssignment(test_mpc.Fluid[part], 0.1); 
            }
            
            #pragma omp for schedule(static) 
            for (unsigned sol = 0; sol < test_mpc.Solute.size(); sol++) {
                test_mpc.updateSoluteCellIndex(test_mpc.Solute[sol]); 
            }
            
            
            #pragma omp single 
            {
                if (n%10==0) {
                    test_mpc.sortVector();
                    if (n > 500) {
                        test_mpc(vel_prof); 
                    }
                }
                
            } 
            
            #pragma omp sections
            {
                #pragma omp section
                test_mpc.sortOnly(); 
                
                #pragma omp section 
                for (unsigned Index = 0; Index < test_mpc.NumberOfCells; Index++) {
                    test_mpc.drawRotation(Index); 
                }   
            }

            #pragma omp for schedule(static) 
            for (unsigned Index = 0; Index < test_mpc.NumberOfCells; Index++) {
                test_mpc.calculateCOMVel(Index); 
            }
            #pragma omp for schedule(static)
            for (unsigned part = 0; part < test_mpc.NumberOfParticles; part++) {
                test_mpc.rotate(part); 
            }
            
            #pragma omp for schedule(static)
            for (unsigned sol = 0; sol < test_mpc.Solute.size(); sol++) {
                test_mpc.rotate(sol); 

            }

            //std::cout << "current temperature " << test_mpc.virtualTemperature() << std::endl; 
            #pragma omp single
            {
                test_mpc.returnSolute(sys_test.Molecules);
                if(n%50==0) {
                    double temp = test_mpc.virtualTemperature();
                    std::cout << "Step: " << n << " Temperature: " << temp << std::endl;
                    temp_file << n << " " << temp << std::endl; 
                }
            } 
        }
    }
   
    gettimeofday(&end, NULL);
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/test_mpc.NumberOfParticles/Steps << std::endl;

    ofstream fluid_profile{}; 
    fluid_profile.open("fluid_profile.dat", ios::out | ios::trunc); 
    vel_prof.print_result(fluid_profile);     
    temp_file.close(); 
    fluid_profile.close(); 
    return 0; 
}
