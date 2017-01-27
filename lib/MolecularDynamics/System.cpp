#include "System.h"

System::System(unsigned Lx, unsigned Ly, unsigned Lz) : 
    Cutoff {1.5}, 
    VerletRadius {2.0}, 
    VerletRadiusSq {4.0} {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        Cells[0] = unsigned(BoxSize[0]/VerletRadius); 
        Cells[1] = unsigned(BoxSize[1]/VerletRadius); 
        Cells[2] = unsigned(BoxSize[2]/VerletRadius);     
        CellSideLength[0] = BoxSize[0]/(double)Cells[0]; 
        CellSideLength[1] = BoxSize[1]/(double)Cells[1];
        CellSideLength[2] = BoxSize[2]/(double)Cells[2];
        CellList = std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>>(Cells[0], std::vector<std::vector<std::forward_list<MDParticle*>>>(Cells[1], std::vector<std::forward_list<MDParticle*>>(Cells[2], std::forward_list<MDParticle*>()))); 
    }
    

System::System(double aCutoff, double aVerletRadius, unsigned Lx, unsigned Ly, unsigned Lz) : 
    Cutoff {aCutoff}, 
    VerletRadius {aVerletRadius}, 
    VerletRadiusSq {aVerletRadius*aVerletRadius} {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        Cells[0] = unsigned(BoxSize[0]/VerletRadius); 
        Cells[1] = unsigned(BoxSize[1]/VerletRadius); 
        Cells[2] = unsigned(BoxSize[2]/VerletRadius);     
        CellSideLength[0] = BoxSize[0]/(double)Cells[0]; 
        CellSideLength[1] = BoxSize[1]/(double)Cells[1];
        CellSideLength[2] = BoxSize[2]/(double)Cells[2];
        CellList = std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>>(Cells[0], std::vector<std::vector<std::forward_list<MDParticle*>>>(Cells[1], std::vector<std::forward_list<MDParticle*>>(Cells[2], std::forward_list<MDParticle*>()))); 
    }
        

void System::updateVerletLists() {
    //std::cout << "rebuilding Verlet list... " << std::endl;  
    for (auto& sheet : CellList) {
        for (auto& row : sheet) {
            for (auto& list : row) {
                list.clear(); 
            }    
        }
    }
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.VerletList.clear(); 
        }
    }
    
    std::array<int, 3> CellNumber{}; 
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            Vector3d imPos{image(mono, BoxSize, 0.0)}; 
            for (unsigned i = 0; i < 3; i++) {
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            mono.VerletPosition = mono.Position; 
            //std::cout << "CellNumber: " << CellNumber[0] << " " << CellNumber[1] << " " << CellNumber[2] << " Position: " << mono.Position.transpose() << " image position: " << imPos.transpose()<< std::endl;
            //CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono); 
            try {
                CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono); 
            }
            catch (std::exception& e) {
                std::cout << "The following exception occurred : " << e.what() << std::endl;
                std::cout << "for a particle with position "  << mono.Position.transpose() << std::endl;
                throw std::out_of_range("OOR");  
                return; 
            }
        }
    }
    int l, m, n; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            for (unsigned i = 0; i < 3; i++) {
                Vector3d imPos{image(mono, BoxSize, 0.0)}; 
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            for (int i = CellNumber[0]-1; i < CellNumber[0]+2; i++) {
                l = i - floor((double)i/Cells[0])*Cells[0]; 
                for (int j = CellNumber[1]-1; j < CellNumber[1]+2; j++) {
                    m = j - floor((double)j/Cells[1])*Cells[1];
                    for (int k = CellNumber[2]-1; k < CellNumber[2]+2; k++) {
                        n = k - floor((double)k/Cells[2])*Cells[2]; 
                        //std::cout << l << " " << m << " " << n << std::endl; 
                        for (auto& other : CellList[l][m][n]) {
                            if (other == &mono) continue; 
                            Vector3d relPos {relative(mono, *other, BoxSize, 0.0)}; 
                            //TODO: boundary conditions 
                            double distance {relPos.squaredNorm()}; 
                            if (distance <= VerletRadiusSq) {
                                mono.VerletList.push_front(other); 
                            }
                        }
                    }  
                }
            }    
        }
    }
}    

void System::checkVerletLists() {
    Vector3d displacement {Vector3d::Zero()}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            displacement = mono.Position - mono.VerletPosition; 
            //TODO: boundary conditions; 
            if (displacement.norm() > (VerletRadius - Cutoff)*0.5) {
                updateVerletLists(); 
                return; 
            }
        }
    }
}

void System::calculateForces() {
    double force_abs {}; 
    Vector3d force {}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) mono.Force = Vector3d::Zero();     
    }
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            for (auto& other : mono.VerletList) {
                Vector3d relPos {relative(mono, *other, BoxSize, 0.0)};  
                double radius2 {relPos.squaredNorm()};
                force_abs = RLJ_Force(radius2); 
                force = relPos*force_abs; 
                mono.Force -= force;  
            }
            
            for (auto& bonded : mono.Bonds) {
                Vector3d relPos {relative(mono, *bonded, BoxSize, 0.0)}; 
                double radius2 {relPos.squaredNorm()}; 
                force_abs = FENE_Force(radius2); 
                force = relPos*force_abs; 
                mono.Force -= force; 
                bonded -> Force += force; 
            } 
        }
    }
}

void System::addMolecules(std::string filename) {
    std::ifstream file {filename}; 
    std::string line; 
    unsigned current_mol, mol, bond1, bond2, monos, numberOfLines;   
    current_mol = 1; 
    monos = 0; 
    if (file.is_open()) {
        std::cout << "file " << filename << " successfully opened" << std::endl; 
        file >> numberOfLines; 
        while (file >> mol >> bond1 >> bond2) {
            if (mol != current_mol) {
                Molecules.push_back(Molecule(monos, 10.)); 
                Molecules[mol-2].setChainBonds(); 
                current_mol = mol;               
            }
            monos = bond2; 
        }
    }
    Molecules.push_back(Molecule(monos, 10.)); 
    Molecules[mol-1].setChainBonds(); 
    current_mol = mol;     
    std::cout << "initialized " << Molecules.size() << " molecules" << std::endl;
}

void System::addLinks(std::string filename) {
    std::ifstream file {filename};
    unsigned mol, bond1, bond2, numberOfLines, count;  
    if (file.is_open()) {
        file >> numberOfLines;
        count = 0;
        while (file >> mol >> bond1 >> bond2) {
            Molecules[mol-1].setLink(bond1-1, bond2-1);
            count++;  
        }    
    } 
    std::cout << "set " << count << " out of " << numberOfLines << " links" << std::endl; 
}

void System::initializePositions(std::string filename) {
    std::ifstream file {filename};
    double x, y, z;  
    unsigned count {}; 
    if (file.is_open()) {  
        for (auto& mol : Molecules) {
            for (auto& mono : mol.Monomers) {
                if (file >> x >> y >> z) {
                    mono.Position(0) = x; 
                    mono.Position(1) = y; 
                    mono.Position(2) = z; 
                    count++; 
                }
                else {
                    std::cout << "only " << count << " monomers were initialized" << std::endl; 
                    return; 
                }
            }
        }
    }
    std::cout << "all " << count << " monomers were initialized" << std::endl; 
    if (file >> x >> y >> z) std::cout << "...but there is more data..." << std::endl; 
}

void System::initializeVelocitiesRandom(double Temperature) {
    Vector3d COMVel {Vector3d::Zero()}; 
    unsigned totalMonomers {}; 
    double EKin {}; 
    for (auto& mol : Molecules) {
        totalMonomers += mol.NumberOfMonomers; 
        for (auto& mono : mol.Monomers) {
            for (unsigned i = 0; i < 3; i++) {
                mono.Velocity(i) = Rand::real_uniform() - 0.5; 
            }
            COMVel += mono.Velocity; 
        }
    }
    COMVel /= totalMonomers; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity -= COMVel; 
            EKin += mono.Velocity.squaredNorm(); 
        }
    }    
    double scaling = sqrt(3.*totalMonomers*Temperature/EKin);
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity *= scaling; 
        }
    }  
}

void System::propagate(double dt) {
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5; 
            mono.Position += mono.Velocity*dt; 
            //TODO: boundary 
        }
    }
    checkVerletLists(); 
    calculateForces(); 
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}




    
