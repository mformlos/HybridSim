#include "System.h"

int main() {
    System sys_test(10,10,10); 
    
    sys_test.addMolecules("/home/formanek/HYBRIDSIM/input/filebondunlinked"); 
    sys_test.addLinks("/home/formanek/HYBRIDSIM/input/bondlinked"); 
    
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
    
    
    
    return 0; 
}
