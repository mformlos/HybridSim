
#include <sys/time.h>
#include "System.h"
#include "HelperFunctions.h"


int main() {
    System sys_test(200,200,200, 0.0, 0.0, false, false); 
    std::vector<unsigned long long> OutputSteps; 
    std::vector<unsigned long long>::iterator OutputStepsIt{};
    unsigned long long TotalSteps {100000000};
    
    
    sys_test.addMolecules("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-PRECURSOR-chain", 1.0); 
    //sys_test.addLinks("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-bonds"); 
    sys_test.initializePositions("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-PRECURSOR-config"); 
    sys_test.initializeVelocitiesRandom(1.0); 
    initializeStepVector(OutputSteps, "/home/formanek/HYBRIDSIM/input/filetimeoutput"); 
    OutputStepsIt = OutputSteps.begin();    
    
    std::cout << "radius of gyration: " << sys_test.Molecules[0].radiusOfGyration() << std::endl; 
    

    FILE* pdb {}; 
    pdb = fopen("pdb_equil_langevin.pdb", "w");
    
    Vector3d NewCOM (100., 100., 100.); 
    sys_test.setMoleculeCOM(0, NewCOM);
     
   
    
      
    sys_test.printPDB(pdb, 0); 
    
    sys_test.calculateForcesBrute(true); 
    std::cout << sys_test.PotentialEnergy() << std::endl; 
    ofstream gyr{"./results/NH-SCNP-0-langevin-stats-2-mass.dat"}; 
    ofstream ext{"./results/ext.dat"};
    ofstream comstat("./results/adsorptionCOM.dat"); 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    bool anchored {false}; 
    bool force_set {false}; 
    Vector3d force(0.0, 0.0, 0.005);
    
    for (unsigned long long i = 0; i < TotalSteps; i++) {
        if (i == *OutputStepsIt) {
            sys_test.propagateLangevin(0.01, 1., 0.05, true); 
            //sys_test.propagate(0.001, true); 
            sys_test.printStatistics(gyr, i*0.01);
            sys_test.printForceExtension(ext, i*0.01, 2);
            sys_test.printPDB(pdb,i);
            comstat << 0.01*i << " " << sys_test.Molecules[0].centerOfMassPosition().transpose() << std::endl; 
            std::cout << 0.01*i << std::endl; 
            OutputStepsIt++;
        }        
        else sys_test.propagateLangevin(0.01, 1., 0.05); //sys_test.propagate(0.001); //

        if (!anchored && i > 15000000000) {
            if (sys_test.Molecules[0].Monomers[0].Position(2) <= 1.05 && sys_test.Molecules[0].Monomers[0].Position(2) >= 0.9 ) {
                Vector3d pos {sys_test.Molecules[0].Monomers[0].Position};
                pos(2) = 0.0; 
                sys_test.setAnchor(0, pos); 
                Vector3d relPos {pos - sys_test.Molecules[0].Monomers[0].Position}; 
                double radius2 {relPos.squaredNorm()};
                std::cout << "relPos: " << relPos.transpose() << "\nradius2: " << radius2 << std::endl;
                std::cout << "Epot in FENE for Anchor: " << FENE_Potential(radius2) << std::endl;  
                std::cout << "Force in FENE for Anchor: " << FENE_Force(radius2) << std::endl;
                sys_test.SurfaceEnergy = 1.5; 
                sys_test.calculateForcesBrute(false);
                anchored = true; 
                std::cout << "Chain begin anchored at time " << i*0.01 << std::endl; 
            }
        }
        
        if (anchored && !force_set && i > 1000000) {
            std::cout << "starting force application" << std::endl; 
            sys_test.setDrive(force, 199); 
            force_set = true; 
            sys_test.calculateForcesBrute(false);
        }
        
        if (force_set && i%100000==0) {
            force(2) += 0.01; 
            sys_test.changeDrive(force, 199);
            sys_test.calculateForcesBrute(false);
        } 
        
        
        
    }
    
    fclose(pdb); 
    
    gettimeofday(&end, NULL); 
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per step: " << realTime/TotalSteps << std::endl;
    
    return 0; 
}
