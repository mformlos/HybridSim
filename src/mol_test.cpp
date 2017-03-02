#include "System.h"

int main() {
    System sys_test(50,50,50); 
    
    sys_test.addMolecules("/home/formanek/HYBRIDSIM/input/SCNP-0-chain", 10.0); 
    sys_test.addLinks("/home/formanek/HYBRIDSIM/input/SCNP-0-bonds"); 
    sys_test.initializePositions("/home/formanek/HYBRIDSIM/input/SCNP-0-config"); 
    sys_test.initializeVelocitiesRandom(1.0); 
    
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
    
    Vector3d vec {-COMPos.front()}; 
    vec[0] += sys_test.BoxSize[0]*0.5; 
    vec[1] += sys_test.BoxSize[1]*0.5; 
    vec[2] += sys_test.BoxSize[2]*0.5; 
    sys_test.Molecules.front().translate(vec);
    
    std::cout << "new center of mass: " << sys_test.Molecules.front().centerOfMassPosition().transpose() << std::endl; 
    
    std::cout << "normal distance: " << (sys_test.Molecules.front().Monomers.front().Position - sys_test.Molecules.front().Monomers.back().Position).transpose(); 
    std::cout << " , with boundaries: " << relative(sys_test.Molecules.front().Monomers.back(), sys_test.Molecules.front().Monomers.front(), sys_test.BoxSize, 0.0).transpose() << std::endl; 
    sys_test.updateVerletLists(); 
    sys_test.calculateForces(); 
    ofstream gyr{"rgyrequil2.dat"}; 
    
    for (int i = 0; i < 10000; i++) {
        if (!(i%100)) {
            sys_test.propagate(0.01, true); 
            double rgyr = sys_test.Molecules.front().radiusOfGyration();
            Vector3d COM = sys_test.Molecules.front().centerOfMassPosition();
            std::cout << i << " " << rgyr << " " << sys_test.Molecules.front().Epot << " COM: " << COM.transpose() << std::endl;
            gyr << i << " " << rgyr << " " << sys_test.Molecules.front().Epot << " COM: " << COM.transpose() << std::endl; 
        } 
        else sys_test.propagate(0.01); 
    }
    
    
    
    
    return 0; 
}
