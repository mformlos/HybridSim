#include "System.h"

System::System(unsigned Lx, unsigned Ly, unsigned Lz, double gamma) : 
    Cutoff {1.5}, 
    VerletRadius {2.0}, 
    VerletRadiusSq {4.0},
    delrx {0.0}, 
    Shear {gamma} {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        Cells[0] = unsigned(BoxSize[0]/(sqrt(2)*VerletRadius)); 
        Cells[1] = unsigned(BoxSize[1]/(sqrt(2)*VerletRadius)); 
        Cells[2] = unsigned(BoxSize[2]/(sqrt(2)*VerletRadius));     
        CellSideLength[0] = BoxSize[0]/(double)Cells[0]; 
        CellSideLength[1] = BoxSize[1]/(double)Cells[1];
        CellSideLength[2] = BoxSize[2]/(double)Cells[2];
        CellList = std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>>(Cells[0], std::vector<std::vector<std::forward_list<MDParticle*>>>(Cells[1], std::vector<std::forward_list<MDParticle*>>(Cells[2], std::forward_list<MDParticle*>()))); 
    }
    

System::System(double aCutoff, double aVerletRadius, unsigned Lx, unsigned Ly, unsigned Lz, double gamma) : 
    Cutoff {aCutoff}, 
    VerletRadius {aVerletRadius}, 
    VerletRadiusSq {aVerletRadius*aVerletRadius},
    delrx {0.0}, 
    Shear {gamma} {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        Cells[0] = unsigned(BoxSize[0]/(sqrt(2)*VerletRadius)); 
        Cells[1] = unsigned(BoxSize[1]/(sqrt(2)*VerletRadius)); 
        Cells[2] = unsigned(BoxSize[2]/(sqrt(2)*VerletRadius));      
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
            Vector3d imPos{image(mono, BoxSize, delrx)}; 
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
                throw CellAllocationException(mono, CellNumber);  
            }
        }
    }
    int l, m, n; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            Vector3d imPos{image(mono, BoxSize, delrx)}; 
            for (unsigned i = 0; i < 3; i++) { 
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            /*if (mono.Position(1) > BoxSize[1]-CellSideLength[1] && delrx > CellSideLength[1]) {
                            std::cout << delrx << " " << CellNumber[0] << " " << CellNumber[1] << " " << CellNumber[2] << endl;;  
            }*/
            for (int i = CellNumber[0]-2; i < CellNumber[0]+3; i++) {
                
                for (int j = CellNumber[1]-2; j < CellNumber[1]+3; j++) {
                    
                    for (int k = CellNumber[2]-2; k < CellNumber[2]+3; k++) {
                        l = i - floor((double)j/Cells[1])*(int)(delrx/CellSideLength[0]); 
                        l -= floor((double)l/Cells[0])*Cells[0]; 
                        m = j - floor((double)j/Cells[1])*Cells[1];
                        n = k - floor((double)k/Cells[2])*Cells[2]; 
                        /*if (mono.Position(1) > BoxSize[1]-CellSideLength[1] && delrx > CellSideLength[1]) {
                            std::cout << l << " " << m << " " << n << endl;;  
                        }*/
                        for (auto& other : CellList[l][m][n]) {
                            if (other == &mono) continue; 
                            Vector3d relPos {relative(mono, *other, BoxSize, delrx)}; 
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
            if (displacement.norm() > (VerletRadius - Cutoff)*0.5) {
                updateVerletLists(); 
                return; 
            }
        }
    }
}

void System::calculateForces(bool calcEpot) {
    //unsigned count_bonds {0}; 
    double force_abs {}; 
    Vector3d force {}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) mono.Force = Vector3d::Zero(); 
        if (calcEpot) mol.Epot = 0.0;    
    }
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            for (auto& other : mono.VerletList) {
                Vector3d relPos {relative(mono, *other, BoxSize, delrx)};  
                double radius2 {relPos.squaredNorm()};
                if (calcEpot) {
                    mol.Epot += 0.5*RLJ_Potential(radius2);

                }
                force_abs = RLJ_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(RLJException(mono.Identifier, mono.Position, mono.Velocity, other -> Identifier, other -> Position, other -> Velocity, force_abs)); 
                }
                force = relPos*force_abs; 
                mono.Force -= force;  
            }
            
            for (auto& bonded : mono.Bonds) {
                //count_bonds++; 
                Vector3d relPos {relative(mono, *bonded, BoxSize, delrx)}; 
                double radius2 {relPos.squaredNorm()}; 
                if (calcEpot) {
                    //double fene {FENE_Potential(radius2)}; 
                    //if (radius2 > 2.25) std::cout << "fene: " << fene << " bond radius: " << radius2 << " at mono: " << mono_count <<  std::endl; 
                    mol.Epot += FENE_Potential(radius2);
                    
                }    
                force_abs = FENE_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(FENEException(mono.Identifier, bonded->Identifier, force_abs)); 
                }
                force = relPos*force_abs; 
                mono.Force -= force; 
                bonded -> Force += force; 
            } 
        }
    }
}

void System::calculateForcesBrute(bool calcEpot) {
    double force_abs {}; 
    Vector3d force {}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) mono.Force = Vector3d::Zero(); 
        if (calcEpot) mol.Epot = 0.0;    
    }
    for (auto& mol : Molecules) {
        for (unsigned i = 0; i < mol.NumberOfMonomers; i++) {
            for (unsigned j = i+1; j < mol.NumberOfMonomers; j++) {            
                Vector3d relPos {relative(mol.Monomers[i], mol.Monomers[j], BoxSize, delrx)};
                //Vector3d relPos {mol.Monomers[j].Position - mol.Monomers[i].Position}; 
                double radius2 {relPos.squaredNorm()}; 
                if (calcEpot) {
                    mol.Epot += RLJ_Potential(radius2); 
                }
                force_abs = RLJ_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    MDParticle* mono {&mol.Monomers[i]}; 
                    MDParticle* other {&mol.Monomers[j]};  
                    throw(RLJException(mono -> Identifier, mono -> Position, mono -> Velocity, other -> Identifier, other -> Position, other -> Velocity, force_abs)); 
                }
                //std::cout << i << " " << j << " " << force_abs << std::endl;
                force = relPos*force_abs; 
                mol.Monomers[i].Force -= force; 
                mol.Monomers[j].Force += force;     
            }
            for (auto& bonded : mol.Monomers[i].Bonds) {
                Vector3d relPos {relative(mol.Monomers[i], *bonded, BoxSize, delrx)}; 
                //Vector3d relPos {bonded->Position - mol.Monomers[i].Position};
                double radius2 {relPos.squaredNorm()}; 
                if (calcEpot) {
                    mol.Epot += FENE_Potential(radius2); 
                }
                force_abs = FENE_Force(radius2); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(FENEException(mol.Monomers[i].Identifier, bonded->Identifier, force_abs)); 
                }
                force = relPos*force_abs; 
                mol.Monomers[i].Force -= force; 
                bonded -> Force += force; 
            }        
        }
    }
}
            




bool System::addMolecules(std::string filename, double mass) {
    std::ifstream file(filename, ios::in);
    if (!file.is_open()) return false; 
    std::string line; 
    unsigned current_mol, mol, bond1, bond2, monos, numberOfLines, mono_count;   
    current_mol = 1; 
    monos = 0; 
    mono_count = 0; 
    if (file.is_open()) {
        std::cout << "file " << filename << " successfully opened" << std::endl; 
        file >> numberOfLines; 
        if (numberOfLines == 0) {
            std::cout << "no molecules in the system" << std::endl; 
            return true; 
        }
        while (file >> mol >> bond1 >> bond2) {
            if (mol != current_mol) {
                mono_count += monos;
                Molecules.push_back(Molecule(monos, mass, mono_count)); 
                Molecules[mol-2].setChainBonds(); 
                current_mol = mol;               
            }
            monos = bond2; 
        }
    }
    Molecules.push_back(Molecule(monos, mass)); 
    Molecules[mol-1].setChainBonds(); 
    current_mol = mol;     
    std::cout << "initialized " << Molecules.size() << " molecules" << std::endl;
    
    /*for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            std::cout << mono.Identifier << " "; 
        }
        std::cout << std::endl; 
    }*/
    
    return true; 
}

bool System::addLinks(std::string filename) {
    std::ifstream file {filename};
    if (!file.is_open()) return false; 
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
    return true; 
}

bool System::initializePositions(std::string filename) {
    std::ifstream file {filename};
    if (!file.is_open()) return false; 
    double x, y, z;  
    unsigned count {0}; 
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
                    return false; 
                }
            }
        }
    }
    std::cout << "all " << count << " monomers were initialized" << std::endl; 
    if (file >> x >> y >> z) std::cout << "...but there is more data..." << std::endl; 
    return true; 
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
        }
        mol.removeAngularMomentum(); 
    }
    for (auto& mol : Molecules) {
        COMVel = mol.centerOfMassVelocity(); 
        for (auto& mono : mol.Monomers) {
            mono.Velocity -= COMVel; 
            EKin += mono.Mass*mono.Velocity.squaredNorm(); 
        }
    }    
    double scaling = sqrt(3.*totalMonomers*Temperature/EKin);
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity *= scaling; 
        }
    }  
}


void System::setMoleculeCOM(unsigned molIndex, Vector3d newCOM) {
    if (molIndex >= Molecules.size()) {
        std::cout << "Molecule number " << molIndex << " does not exist." << std::endl; 
        return; 
    }
    Vector3d currentCOM {Molecules[molIndex].centerOfMassPosition()}; 
    newCOM -= currentCOM; 
    for (auto& mono : Molecules[molIndex].Monomers) {
        mono.Position += newCOM; 
        if (mono.Position(1) >= BoxSize[1]) { 
            std::cout << mono.Position(1) << " COM to big"; 
            exit(0); 
        }
    }
}

void System::centerMolecule(unsigned molIndex) {
    Vector3d BoxCenter(BoxSize[0]*0.5, BoxSize[1]*0.5, BoxSize[2]*0.5); 
    setMoleculeCOM(molIndex, BoxCenter); 
}

void System::wrapMoleculesCOM() {
    for (auto& mol : Molecules) {
        wrapCOM(mol, BoxSize, Shear, delrx); 
    }
}
void System::propagate(double dt, bool calcEpot) {
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5; 
            mono.Position += mono.Velocity*dt; 
            //TODO: boundary 
        }
    }
    //checkVerletLists(); 
    //calculateForces(calcEpot); 
    calculateForcesBrute(calcEpot); 
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}



unsigned System::NumberOfParticles() {
    unsigned count {}; 
    for (auto& mol : Molecules) {
        count += mol.NumberOfMonomers; 
    }
    return count; 
}

unsigned System::NumberOfMolecules() {return Molecules.size();}

double System::KineticEnergy() {
    double Ekin = 0.0; 
    for (auto& mol : Molecules) {
        Ekin += mol.KineticEnergy(); 
    }
    return Ekin; 
}

double System::PotentialEnergy() {
    double Epot {0.0}; 
    for (auto& mol : Molecules) {
        Epot += mol.PotentialEnergy(); 
    }
    return Epot; 
}

std::tuple<double, Matrix3d> System::GyrationTensor() {
    double rgyr {}; 
    Matrix3d gyrTensor {Matrix3d::Zero()}; 
    std::tuple<double, Matrix3d> gyrMol {}; 
    for (auto& mol : Molecules) {
        gyrMol = mol.GyrationTensor(); 
        gyrTensor += std::get<1>(gyrMol); 
        rgyr += std::get<0>(gyrMol); 
    }
    gyrTensor /= Molecules.size(); 
    rgyr /= Molecules.size(); 
    return std::make_tuple(rgyr, gyrTensor);    
}

Vector3d System::RotationFrequency() {
    Vector3d omega {Vector3d::Zero()}; 
    for (auto& mol : Molecules) {
        omega += mol.RotationFrequency(); 
    }
    omega /= Molecules.size(); 
    return omega; 
}

void System::printPDB(FILE* pdb, int step, bool velocs) {
    int mol_count{0}; 
    fprintf(pdb, "MODEL     %d \n", step);
    for (auto& mol : Molecules) {
		mol_count++;
		int mono_count{0}; 
		for (auto& mono : mol.Monomers) {
		    if (velocs) {
		        fprintf(pdb, "ATOM %6d  C   GLY    %2d     %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n", mono_count+1, mol_count, mono.Position(0), mono.Position(1), mono.Position(2), mono.Velocity(0), mono.Velocity(1), mono.Velocity(2));
		    }
			else {
			    fprintf(pdb, "ATOM %6d  C   GLY    %2d     %7.3f %7.3f %7.3f \n", mono_count+1, mol_count, mono.Position(0), mono.Position(1), mono.Position(2));
		    }
		    mono_count++; 
		}	
        fprintf(pdb, "TER \n");
    }
    fprintf(pdb, "ENDMDL \n");   
    fflush(pdb);       
}

void System::printStatistics(std::ofstream& os, double time) {
    os.precision(6); 
    double Ekin {KineticEnergy()}, Epot {PotentialEnergy()}; 
    Vector3d Omega {RotationFrequency()}; 
    std::tuple<double, Matrix3d> GyrTuple {GyrationTensor()};
    Matrix3d GyrTensor {std::get<1>(GyrTuple)}; 
    os.precision(6); 
    os.width(14); 
    os << time << " "; 
    os.width(14); 
    os << Ekin << " ";
    os.width(14); 
    os << Epot << " "; 
    for (unsigned i = 0; i < 3; i++) {
        os.width(14);
        os << Omega(i) << " "; 
    }
    os.width(14);
    os << std::get<0>(GyrTuple); 
    for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = 0; j < 3; j++) {
            os.width(14); 
            os << GyrTensor(i,j) << " ";
        } 
    }
    os << std::endl;  
}


    
