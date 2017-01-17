#include <sys/time.h>
#include <omp.h>
#include "Particle.h"
#include "MPC.h"

int main() {
    Particle test_part; 
    unsigned Steps = 50; 
    unsigned N = 10000; 
    
    test_part.Velocity(0) = 1.0; 
    
    MPC test_mpc(100, 100, 100, 10, 1.0);
    test_mpc.initialize_random();   
    std::cout << "Temperature after initialization: " << test_mpc.virtualTemperature() << std::endl;

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
    #pragma omp parallel 
    {
        for (unsigned n = 0; n < Steps; n++) {
            #pragma omp for schedule(static) 
            for (unsigned part = 0; part < test_mpc.NumberOfParticles; part++) {
                test_mpc.streamPlusCellAssignment(test_mpc.Fluid[part], 0.1); 
                //test_mpc.stream(test_mpc.Fluid[part], 0.1); 
            }
            #pragma omp single 
            {
                if (n%10==0) {
                    //test_mpc.updateParticleCellIndex();
                    test_mpc.sortVector(); 
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

            //std::cout << "current temperature " << test_mpc.virtualTemperature() << std::endl; 
            #pragma omp single
            {
                if(n%10==0) std::cout << "Step: " << n << " Temperature: " << test_mpc.virtualTemperature() << std::endl;
            } 
        }
    }
   
    gettimeofday(&end, NULL);
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/test_mpc.NumberOfParticles/Steps << std::endl;

    return 0; 
}
