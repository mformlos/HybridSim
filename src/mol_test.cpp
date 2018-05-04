#include "System.h"
#include "HelperFunctions.h"

int main() {
    System sys_test(40,40,40, 0.0); 
    std::vector<unsigned long long> OutputSteps; 
    std::vector<unsigned long long>::iterator OutputStepsIt{};
    
    sys_test.addMolecules("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-chain", 5.0); 
    sys_test.addLinks("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-bonds"); 
    sys_test.initializePositions("/home/formanek/HYBRIDSIM/input/SCNPs/SCNP-0-config"); 
    sys_test.initializeVelocitiesRandom(1.0); 
    initializeStepVector(OutputSteps, "/home/formanek/HYBRIDSIM/input/teststeps"); 
    OutputStepsIt = OutputSteps.begin();    
    unsigned monocount {0}, bondcount {0};  
    for (auto& mol : sys_test.Molecules) {
        //std::cout << "new molecule" << std::endl; 
        for (auto& mono : mol.Monomers) {
            monocount++; 
            for (auto& bond : mono.Bonds) {
                bondcount++; 
                //std::cout << &mono << " " << bond << std::endl;   
            }
        }
    }
    std::cout << "we have " << monocount << " monomers and " << bondcount << " bonds in total" << std::endl; 
    
    std::vector<Vector3d> COMPos;
    for (auto& mol : sys_test.Molecules) {
        COMPos.push_back(mol.centerOfMassPosition()); 
    }
    
    for (unsigned i = 0; i < COMPos.size(); i++) {
        std::cout << i << " " << COMPos[i].transpose() << std::endl; 
        for (unsigned j = i+1; j < COMPos.size(); j++) {
            double distance {(COMPos[i] - COMPos[j]).norm()}; 
            if (distance < 100) std::cout << i << " " << j << " " << distance << std::endl;     
        } 
    }
    
    double rgyr_mean {}; 
    unsigned mol_count {}; 
    for (auto& mol : sys_test.Molecules) {
        double rgyr {mol.radiusOfGyration()}; 
        rgyr_mean += rgyr; 
        ++mol_count; 
        std::cout << rgyr << std::endl; 
    }
    rgyr_mean /= mol_count; 
    std::cout << "mean radius of gyration: " << rgyr_mean << std::endl; 
    
    Particle test_part{}; 
    test_part.Position(0) = 24.0; 
    test_part.Position(1) = 13.0; 
    test_part.Position(2) = 37.0; 
    FILE* pdb {}; 
    pdb = fopen("pdb_equil_brute.pdb", "w");
    Vector3d vec {-COMPos.front()}; 
    vec[0] += sys_test.BoxSize[0]*0.5; 
    vec[1] += sys_test.BoxSize[1]*0.5; 
    vec[2] += sys_test.BoxSize[2]*0.5; 
    sys_test.Molecules.front().translate(vec);
    sys_test.printPDB(pdb, 0);
    
    /*vec[0] = sys_test.BoxSize[0]*0.5+5.0; 
    vec[1] = sys_test.BoxSize[1]*0.5+3.0;
    vec[2] = sys_test.BoxSize[2]*0.5+12.0;  
    sys_test.Molecules.front().translate(vec);
    */
    std::cout << "center of mass out of box: " << sys_test.Molecules.front().centerOfMassPosition().transpose() << std::endl;
    std::cout << "distance to test particle: " << relative(sys_test.Molecules.front().Monomers.back(), test_part, sys_test.BoxSize, 0.5).transpose() << std::endl; 
    sys_test.printPDB(pdb, 0); 
    
    wrapCOM(sys_test.Molecules.front(), sys_test.BoxSize, 0.0, 0.5); 
    sys_test.printPDB(pdb, 0); 
    
   
    
    std::cout << "new center of mass: " << sys_test.Molecules.front().centerOfMassPosition().transpose() << std::endl; 
    
    std::cout << "distance to test particle: " << relative(sys_test.Molecules.front().Monomers.back(), test_part, sys_test.BoxSize, 0.5).transpose() << std::endl; 
    //sys_test.updateVerletLists(); 
    sys_test.calculateForcesBrute(true); 
    std::cout << sys_test.PotentialEnergy() << std::endl; 
    ofstream gyr{"./results/NH-SCNP-0-stats.dat"}; 

    
    for (unsigned long long i = 0; i < 10000000; i++) {
        if (i == *OutputStepsIt) {
            sys_test.propagate(0.001, true); 
            sys_test.printStatistics(gyr, i*0.001);
            std::cout << 0.001*i << std::endl; 
            OutputStepsIt++;
        }
        else sys_test.propagate(0.001); 
    }
    
    fclose(pdb); 
    
    
    return 0; 
}
