#include <sys/time.h>
#include <omp.h>
#include "MPC.h"
#include "Velocity_Profile.h"
#include "System.h"

int main() {

    unsigned Lx = 40, Ly = 40, Lz = 40, Steps = 10000, COMwrapInterval = 1000, MPCInterval;   
    int tid, procs, maxt, inpar, dynamic, nested, nthreads;
    double StepSize = 0.01, StepSizeMPC = 0.1, Shear = 0.0;
    
    MPCInterval = (unsigned)(StepSizeMPC/StepSize); 
    std::cout << "MPC every " << MPCInterval << " steps" << std::endl; 
    
    System sys(1.2, 1.5, Lx,Ly,Lz, Shear, 0.0, false, true); 
    MPC mpc(Lx, Ly, Lz, 5, 1.0, Shear);
    VelocityProfile vel_prof{0.2}; 

    sys.addMolecules("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-chain", 5.0); 
    sys.addLinks("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-bonds"); 
    sys.initializePositions("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-config"); 
    sys.initializeVelocitiesRandom(1.0); 
    

    Vector3d boxCenter {Lx*0.5, Ly*0.5, Lz*0.5}; 
    sys.setMoleculeCOM(0, boxCenter); 
    
    
    sys.updateVerletLists(); 
    sys.calculateForces(); 
    
    mpc.initializeRandom();  
    mpc.initializeSoluteVector(sys.NumberOfParticles()); 
    
    

    timeval start, end;
    gettimeofday(&start, NULL);
    
    ofstream gyrfile {}; 
    gyrfile.open("rgyrmpc2.dat", ios::out | ios::trunc);
    FILE* pdb {}; 
    pdb = fopen("pdb_mpc.pdb", "w");
    
    std::cout << "omega: " << sys.Molecules.front().RotationFrequency().transpose() << std::endl; 
    
    #pragma omp parallel private(tid)
    {
        #ifdef _OPENMP
        tid = omp_get_thread_num(); 
        if (tid == 0) {
            printf("Thread %d getting environment info...\n", tid);

            /* Get environment information */
            procs = omp_get_num_procs();
            nthreads = omp_get_num_threads();
            maxt = omp_get_max_threads();
            inpar = omp_in_parallel();
            dynamic = omp_get_dynamic();
            nested = omp_get_nested();
    
            /* Print environment information */
            printf("Number of processors = %d\n", procs);
            printf("Number of threads = %d\n", nthreads);
            printf("Max threads = %d\n", maxt);
            printf("In parallel? = %d\n", inpar);
            printf("Dynamic threads enabled? = %d\n", dynamic);
            printf("Nested parallelism enabled? = %d\n", nested);
        }
        #endif
        Rand::seed(tid);
        Rand::warmup(10000); 
        
        for (unsigned n = 0; n < Steps; n++) {
            #pragma omp single 
            {
                mpc.updateBoxShift(StepSize); 
                sys.delrx=mpc.delrx; 
                if (n%COMwrapInterval == 0) sys.wrapMoleculesCOM();
                if (n%1000 ==0) sys.propagate(StepSize, true); 
                else sys.propagate(StepSize); 
            }
            
            if (n%MPCInterval == MPCInterval-1) {
                #pragma omp single 
                {
                    mpc.shiftGrid();
                    mpc.getSolute(sys.Molecules); 
                    for (auto& sol : mpc.Solute) {
                        mpc.updateSoluteCellIndex(sol); 
                    }
                }
                
                #pragma omp for schedule(static) 
                for (unsigned part = 0; part < mpc.NumberOfParticles; part++) {
                    mpc.streamPlusCellAssignment(mpc.Fluid[part], 0.1); 
                }
                
                #pragma omp single 
                {
                    if (n%(10*MPCInterval)==MPCInterval-1) {
                        mpc.sortVector();
                    }
                }    
                
                #pragma omp sections
                {
                    #pragma omp section
                    mpc.sortOnly(); 
                
                    #pragma omp section 
                    for (unsigned Index = 0; Index < mpc.NumberOfCells; Index++) {
                        mpc.drawRotation(Index); 
                    }   
                }
                
                #pragma omp for schedule(static) 
                for (unsigned Index = 0; Index < mpc.NumberOfCells; Index++) {
                    mpc.calculateCOMVel(Index); 
                }
                #pragma omp for schedule(static)
                for (unsigned part = 0; part < mpc.NumberOfParticles; part++) {
                    mpc.rotate(part); 
                }
            
                #pragma omp for schedule(static)
                for (unsigned sol = 0; sol < mpc.Solute.size(); sol++) {
                    mpc.rotateSolute(sol); 
                }
                
                #pragma omp single
                {
                    mpc.returnSolute(sys.Molecules);
                } 
            }
            #pragma omp single 
            {    
                if (n%1000 ==0) {
                    sys.printStatistics(gyrfile, n*0.01);
                    double temp = mpc.virtualTemperature();
                    
                    mpc(vel_prof); 
                    sys.printPDB(pdb, n);
                    std::cout << n << " " << temp << std::endl;
                }
            }
        }
    } //end of parallel section 
    
    gettimeofday(&end, NULL);
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/mpc.NumberOfParticles/Steps << std::endl;
    ofstream vel_prof_file{}; 
    vel_prof_file.open("vel_prof2.dat", ios::out | ios::trunc); 
    vel_prof.print_result(vel_prof_file);
    vel_prof_file.close();      
    
    
}
