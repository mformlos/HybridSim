#include "System.h"

System::System(unsigned Lx, unsigned Ly, unsigned Lz, double gamma, double energy, bool AdsorptionOn, bool PBCon) : 
    Cutoff {1.2}, 
    VerletRadius {1.5}, 
    VerletRadiusSq {2.25},
    delrx {0.0}, 
    Shear {gamma},
    SurfaceEnergy {energy},
    Epot{0.0},
    Adsorption {AdsorptionOn},
    PBC{PBCon}
     {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        /*Cells[0] = unsigned(BoxSize[0]/(sqrt(2)*VerletRadius)); 
        Cells[1] = unsigned(BoxSize[1]/(sqrt(2)*VerletRadius)); 
        Cells[2] = unsigned(BoxSize[2]/(sqrt(2)*VerletRadius));*/   
        Cells[0] = (unsigned)(BoxSize[0]/Cutoff); 
        Cells[1] = (unsigned)(BoxSize[1]/Cutoff); 
        Cells[2] = (unsigned)(BoxSize[2]/Cutoff);    
        CellSideLength[0] = BoxSize[0]/(double)Cells[0]; 
        CellSideLength[1] = BoxSize[1]/(double)Cells[1];
        CellSideLength[2] = BoxSize[2]/(double)Cells[2];
        CellList = std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>>(Cells[0], std::vector<std::vector<std::forward_list<MDParticle*>>>(Cells[1], std::vector<std::forward_list<MDParticle*>>(Cells[2], std::forward_list<MDParticle*>()))); 
    }
    

System::System(double aCutoff, double aVerletRadius, unsigned Lx, unsigned Ly, unsigned Lz, double gamma, double energy, bool AdsorptionOn, bool PBCon) : 
    Cutoff {aCutoff}, 
    VerletRadius {aVerletRadius}, 
    VerletRadiusSq {aVerletRadius*aVerletRadius},
    delrx {0.0}, 
    Shear {gamma},
    SurfaceEnergy {energy},
    Epot{0.0},
    Adsorption {AdsorptionOn},
    PBC {PBCon} {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly; 
        BoxSize[2] = Lz;
        /*Cells[0] = unsigned(BoxSize[0]/(sqrt(2)*VerletRadius)); 
        Cells[1] = unsigned(BoxSize[1]/(sqrt(2)*VerletRadius)); 
        Cells[2] = unsigned(BoxSize[2]/(sqrt(2)*VerletRadius));*/
        Cells[0] = unsigned(BoxSize[0]/Cutoff); 
        Cells[1] = unsigned(BoxSize[1]/Cutoff); 
        Cells[2] = unsigned(BoxSize[2]/Cutoff);     
        CellSideLength[0] = BoxSize[0]/(double)Cells[0]; 
        CellSideLength[1] = BoxSize[1]/(double)Cells[1];
        CellSideLength[2] = BoxSize[2]/(double)Cells[2];
        CellList = std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>>(Cells[0], std::vector<std::vector<std::forward_list<MDParticle*>>>(Cells[1], std::vector<std::forward_list<MDParticle*>>(Cells[2], std::forward_list<MDParticle*>()))); 
    }
        
bool System::setNeighbourDirections(std::string filename) {
    std::ifstream file (filename, std::ios::in); 
    if (!file.is_open()) {
        return false; 
    }
    int x {}, y {}, z{}; 
    unsigned n {0};
    while(file >> x >> y >> z) {
        std::cout << x << " "<< y << " " << z << std::endl; 
        NeighbourDirections[n][0] = x; 
        NeighbourDirections[n][1] = y;
        NeighbourDirections[n][2] = z;
        n++; 
    }
    if (n != 16) {
        std::cout << n << " instead of 16 neighbour directions supplied!" << std::endl; 
        return false; 
    }
    return true; 
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
            Vector3d imPos; 
            if (PBC) imPos = image(mono, BoxSize, delrx);
            else imPos = mono.Position; 
            for (unsigned i = 0; i < 3; i++) {
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            mono.VerletPosition = mono.Position; 
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
            Vector3d imPos; 
            if (PBC) imPos = image(mono, BoxSize, delrx);
            else imPos = mono.Position; 
            for (unsigned i = 0; i < 3; i++) { 
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            for (int i = CellNumber[0]-2; i < CellNumber[0]+3; i++) {
                
                for (int j = CellNumber[1]-2; j < CellNumber[1]+3; j++) {
                    
                    for (int k = CellNumber[2]-2; k < CellNumber[2]+3; k++) {
                        l = i - floor((double)j/Cells[1])*(int)(delrx/CellSideLength[0]); 
                        l -= floor((double)l/Cells[0])*Cells[0]; 
                        m = j - floor((double)j/Cells[1])*Cells[1];
                        n = k - floor((double)k/Cells[2])*Cells[2]; 
                        for (auto& other : CellList[l][m][n]) {
                            if (other == &mono) continue; 
                            Vector3d relPos; 
                            if (PBC) relPos = relative(mono, *other, BoxSize, delrx);
                            else relPos = other -> Position - mono.Position; 
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

void System::updateCellLists() {
    std::array<int, 3> CellNumber{}; 
    for (auto& sheet : CellList) {
        for (auto& row : sheet) {
            for (auto& list : row) {
                list.clear(); 
            }    
        }
    }
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            Vector3d imPos; 
            if (PBC) imPos = image(mono, BoxSize, delrx);
            else imPos = mono.Position; 
            for (unsigned i = 0; i < 3; i++) {
                CellNumber[i] = (int)(imPos(i)/CellSideLength[i]); 
            }
            try {
                CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mono); 
            }
            catch (std::exception& e) {
                throw CellAllocationException(mono, CellNumber);  
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

void System::calculateForcesVerlet(bool calcEpot) {
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
                Vector3d relPos;
                if (PBC) relPos = relative(mono, *other, BoxSize, delrx); 
                else relPos = other -> Position - mono.Position;  
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
                Vector3d relPos;
                if (PBC) relPos = relative(mono, *bonded, BoxSize, delrx); 
                else relPos = bonded -> Position - mono.Position;  
                double radius2 {relPos.squaredNorm()}; 
                if (calcEpot) {
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
            
            if (Adsorption) {
                double radius2 {pow(mono.Position(2),2)}; 
                if (calcEpot) mol.Epot += LJSurf_Potential(radius2, SurfaceEnergy); 
                force_abs = LJSurf_Force(radius2, SurfaceEnergy); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(LJSurfException(mono.Identifier, mono.Position, mono.Velocity, force_abs)); 
                } 
                mono.Force(2) += mono.Position(2)*force_abs; 
            } 
        }
    }
    
    for (auto& drivepair : Driven) {
        drivepair.DrivenParticle -> Force += drivepair.Force; 
    }
    
    for (auto& anch : Anchors) {
        Vector3d relPos;
        if (PBC) relPos = relative(anch.AnchoredParticle -> Position, anch.AnchorPoint, BoxSize, delrx);
        else relPos = anch.AnchorPoint - anch.AnchoredParticle -> Position;  
        double radius2 {relPos.squaredNorm()}; 
        if (calcEpot) { 
            Molecules[0].Epot += FENE_Potential(radius2); 
        }    
        force_abs = FENE_Force(radius2); 
            if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                throw(FENEException(anch.AnchoredParticle -> Identifier, -1, force_abs)); 
            }
        force = relPos*force_abs; 
        anch.AnchoredParticle -> Force -= force;   
    }    
}

void System::calculateForcesCellList(bool calcEpot) {
    int xbox_displacement{(unsigned)(delrx/CellSideLength[0])};
    double force_abs {}; 
    Vector3d force {}; 
    Vector3d relPos;
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) mono.Force = Vector3d::Zero(); 
        if (calcEpot) Epot = 0.0;    
    }
    //loop over all cells 
    for (int i = 0; i < Cells[0]; i++) {
        for (int j = 0; j < Cells[1]; j++) {
            for (int k = 0; k < Cells[2]; k++) {
                //std::cout << "Current Cell: " << i << " " << j << " " << k << std::endl; 
                ///loop over all monomers in this cell
                
                for (auto first_it = CellList[i][j][k].begin(); first_it != CellList[i][j][k].end(); first_it++) {
                    MDParticle* first = *first_it;
                    //std::vector<int> Neighbours {}; 
                    
                    //std::cout << "first: " << first -> Identifier << std::endl; 
                    /// loop over all further monomers in this cell
                    if (first_it !=  CellList[i][j][k].end()) {
                        for (auto second_it = std::next(first_it,1); second_it != CellList[i][j][k].end(); second_it++){
                            MDParticle* second = *second_it;  
                            //std::cout << "second: " << second -> Identifier << std::endl; 
                            if (PBC) relPos = relative(*first, *second, BoxSize, delrx); 
                            else relPos = second -> Position - first -> Position;  
                            double radius2 {relPos.squaredNorm()};
                            if (calcEpot) {
                                Epot += RLJ_Potential(radius2);

                            }
                            force_abs = RLJ_Force(radius2); 
                            if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                                throw(RLJException(first -> Identifier, first -> Position, first -> Velocity, second -> Identifier, second -> Position, second -> Velocity, force_abs)); 
                            }
                            //if (force_abs > 0) Neighbours.push_back(second -> Identifier); 
                            force = relPos*force_abs; 
                            first -> Force -= force;  
                            second -> Force += force; 
                        }
                    }    
                    ////loop over all neighbors 
                    int l{}, m{}, n{}; 
                    /// all except the top layer boxes 
                    if (j != Cells[1] - 1) {  
                        for (unsigned dir = 0; dir < 13; dir++) {
                            
                            l = i + NeighbourDirections[dir][0];
                            m = j + NeighbourDirections[dir][1];
                            n = k + NeighbourDirections[dir][2];
                            l -= floor((double)l/Cells[0])*Cells[0]; 
                            m -= floor((double)m/Cells[1])*Cells[1];
                            n -= floor((double)n/Cells[2])*Cells[2];
                            //std::cout << "Middle Cell: " << l << " " << m << " " << n << std::endl; 
                            //l -= floor((double)j/Cells[1]) *(int)(delrx/CellSideLength[0]); 
                            /// loop over all monomers in this neighbouring cell
                            for (auto& second : CellList[l][m][n]) {
                                //std::cout << "second: " << second -> Identifier << std::endl; 
                                if (PBC) relPos = relative(*first, *second, BoxSize, delrx); 
                                else relPos = second -> Position - first -> Position;  
                                double radius2 {relPos.squaredNorm()};
                                if (calcEpot) {
                                    Epot += RLJ_Potential(radius2);

                                }
                                force_abs = RLJ_Force(radius2); 
                                //std::cout << force_abs << std::endl;    
                                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                                    throw(RLJException(first -> Identifier, first -> Position, first -> Velocity, second -> Identifier, second -> Position, second -> Velocity, force_abs)); 
                                }
                                //if (force_abs > 0) Neighbours.push_back(second -> Identifier); 
                                force = relPos*force_abs; 
                                first -> Force -= force;  
                                second -> Force += force; 
                            }
                        }
                    }
                    /// for the top layer boxes
                    else { 
                        //std::cout << "Current Cell: " << i << " " << j << " " << k << " delrx: " << delrx << " displ: " << xbox_displacement << std::endl;
                       
                        for (unsigned dir = 0; dir < 16; dir++) {
                            
                            l = i + NeighbourDirections[dir][0];
                            m = j + NeighbourDirections[dir][1];
                            n = k + NeighbourDirections[dir][2];
                            
                            if (NeighbourDirections[dir][1] == 1) {
                                l -= xbox_displacement;     
                            }
                            
                            l -= floor((double)l/Cells[0])*Cells[0]; 
                            m -= floor((double)m/Cells[1])*Cells[1];
                            n -= floor((double)n/Cells[2])*Cells[2];
                            
                            //std::cout << "Top Cell: " << l << " " << m << " " << n << std::endl; 
                            //l -= floor((double)j/Cells[1]) *(int)(delrx/CellSideLength[0]); 
                            /// loop over all monomers in this neighbouring cell
                            for (auto& second : CellList[l][m][n]) {
                                //std::cout << "second: " << second -> Identifier << std::endl;
                                if (PBC) relPos = relative(*first, *second, BoxSize, delrx); 
                                else relPos = second -> Position - first -> Position;  
                                double radius2 {relPos.squaredNorm()};
                                if (calcEpot) {
                                    Epot += RLJ_Potential(radius2);

                                }
                                force_abs = RLJ_Force(radius2); 
                                //std::cout << force_abs << std::endl; 
                                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                                    throw(RLJException(first -> Identifier, first -> Position, first -> Velocity, second -> Identifier, second -> Position, second -> Velocity, force_abs)); 
                                }
                                //if (force_abs > 0) Neighbours.push_back(second -> Identifier); 
                                force = relPos*force_abs; 
                                first -> Force -= force;  
                                second -> Force += force; 
                            }
                        }    
                    
                    }
                    /// calculation of FENE etc. 
                    for (auto& bonded : first -> Bonds) {
                        Vector3d relPos;
                        if (PBC) relPos = relative(*first, *bonded, BoxSize, delrx); 
                        else relPos = bonded -> Position - first -> Position;  
                        double radius2 {relPos.squaredNorm()}; 
                        if (calcEpot) {
                            Epot += FENE_Potential(radius2);
                        }    
                        force_abs = FENE_Force(radius2); 
                        if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                            throw(FENEException(first -> Identifier, bonded->Identifier, force_abs)); 
                        }
                        force = relPos*force_abs; 
                        first -> Force -= force; 
                        bonded -> Force += force; 
                    }
                    /*if (first -> Identifier == 18 || first -> Identifier == 19 ) {
                        std::cout << "Particle " << first -> Identifier << "'s position: " << first -> Position << std::endl; 
                        std::cout << "cell: " << i << " " << j << " " << k << std::endl; 
                        std::cout << "neighbours: " << std::endl; 
                        for(auto& neigh : Neighbours) std::cout << neigh << " "; 
                        std::cout << std::endl; 
                    }*/
                }   
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
    /// loop over all molecules
    for (unsigned molm = 0; molm < Molecules.size(); molm++) {
        Molecule& mol {Molecules[molm]};  
        /// intermolecular interactions
        for (unsigned i = 0; i < mol.NumberOfMonomers; i++) {
            for (unsigned j = i+1; j < mol.NumberOfMonomers; j++) {
                Vector3d relPos;             
                if (PBC) relPos = relative(mol.Monomers[i], mol.Monomers[j], BoxSize, delrx);
                else relPos = mol.Monomers[j].Position - mol.Monomers[i].Position; 
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
                Vector3d relPos;
                if (PBC) relPos = relative(mol.Monomers[i], *bonded, BoxSize, delrx); 
                else relPos = bonded->Position - mol.Monomers[i].Position;
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
            if (Adsorption) {
                double radius2 {pow(mol.Monomers[i].Position(2),2)}; 
                if (calcEpot) mol.Epot += LJSurf_Potential(radius2, SurfaceEnergy); 
                force_abs = LJSurf_Force(radius2, SurfaceEnergy); 
                if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                    throw(LJSurfException(mol.Monomers[i].Identifier, mol.Monomers[i].Position, mol.Monomers[i].Velocity, force_abs)); 
                } 
                mol.Monomers[i].Force(2) += mol.Monomers[i].Position(2)*force_abs; 
            }         
        }
        /// loop over other molecules with higher index
        for (unsigned moln = molm + 1; moln < Molecules.size(); moln++) {
            Molecule& other {Molecules[moln]}; 
            for (unsigned i = 0; i < mol.NumberOfMonomers; i++) {
                for (unsigned j = 0; j < other.NumberOfMonomers; j++) { 
                    Vector3d relPos;             
                    if (PBC) relPos = relative(mol.Monomers[i], other.Monomers[j], BoxSize, delrx);
                    else relPos = other.Monomers[j].Position - mol.Monomers[i].Position; 
                    double radius2 {relPos.squaredNorm()}; 
                    if (calcEpot) {
                        double pot {0.5*RLJ_Potential(radius2)};
                        mol.Epot += pot; 
                        other.Epot += pot;
                    }
                    force_abs = RLJ_Force(radius2); 
                    if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
                        MDParticle* mono {&mol.Monomers[i]}; 
                        MDParticle* othermono {&other.Monomers[j]};  
                        throw(RLJException(mono -> Identifier, mono -> Position, mono -> Velocity, othermono -> Identifier, othermono -> Position, othermono -> Velocity, force_abs)); 
                    }
                    //std::cout << i << " " << j << " " << force_abs << std::endl;
                    force = relPos*force_abs; 
                    mol.Monomers[i].Force -= force; 
                    other.Monomers[j].Force += force;   
                }
            }
        }
    }
    
    for (auto& drivepair : Driven) {
        drivepair.DrivenParticle -> Force += drivepair.Force; 
    }
    
    for (auto& anch : Anchors) {
        Vector3d relPos;
        if (PBC) relPos = relative(anch.AnchoredParticle -> Position, anch.AnchorPoint, BoxSize, delrx);
        else relPos = anch.AnchorPoint - anch.AnchoredParticle -> Position;   
        double radius2 {relPos.squaredNorm()};
        if (calcEpot) { 
            Molecules[0].Epot += FENE_Potential(radius2); 
        }    
        force_abs = FENE_Force(radius2); 
        if (fabs(force_abs) > 1e4 || std::isinf(force_abs) || std::isnan(force_abs)) {
            throw(FENEException(anch.AnchoredParticle -> Identifier, -1, force_abs)); 
        }
        force = relPos*force_abs; 
        anch.AnchoredParticle -> Force -= force;   
    } 
}


bool System::addMolecules(std::string filename, double mass) {
    std::ifstream file(filename, ios::in);
    if (!file.is_open()) return false;
    std::string line;
    unsigned current_mol{0}, mol, monos, mono_start{0};

    if (file.is_open()) {
        std::cout << "file " << filename << " successfully opened" << std::endl;
        while (file >> mol >> monos) {
        	if (mol != current_mol) {
        		std::cout << "order of molecules is wrong" << std::endl;
        		return false;
        	}
        	Molecules.push_back(Molecule(monos, mass, mono_start));
        	Molecules[mol].setChainBonds();
            mono_start += monos;
            current_mol++;
        }
    }
    std::cout << "initialized " << Molecules.size() << " molecules with a total of " << mono_start << " monomers." << std::endl;

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

bool System::initializePositionsPDB(std::string filename) {
    std::ifstream file {filename};
    if (!file.is_open()) return false;
    std::string dump, line{}; 
    double x, y, z, vx, vy, vz;
    unsigned count {0}, Npart {NumberOfParticles()}; 
    unsigned mol{0}, mono{0}; 
    if (file.is_open()) {
        while(getline(file, line)) {
            std::istringstream iss(line); 
            if (iss >> dump >> dump >> dump >> dump >> dump >> x >> y >> z >> vx  >> vy >> vz)        
            { 
			    if (count >= Npart) {
			        std::cout << "all " << count << " monomer positions were initialized, but there is more data" << std::endl; 
			        return true; 
			    }
			    Molecules[mol].Monomers[mono].Position(0) = x;
			    Molecules[mol].Monomers[mono].Position(1) = y;
			    Molecules[mol].Monomers[mono].Position(2) = z;
			    Molecules[mol].Monomers[mono].Velocity(0) = vx;
			    Molecules[mol].Monomers[mono].Velocity(1) = vy;
			    Molecules[mol].Monomers[mono].Velocity(2) = vz;
			    mono++; 
			    if (mono >= Molecules[mol].NumberOfMonomers) {
			        mono = 0; 
			        mol++;    
			    }
			    count++;
			}
        }
    }    
       
    if (count != Npart) {
       	std::cout << "only " << count << " monomers were initialized" << std::endl;
        return false;
    }
    //std::cout << "all " << count << " monomer positions were initialized" << std::endl;
    file.close();
    return true;
}


bool System::initializeVelocities(std::string filename) {
    std::ifstream file {filename};
    if (!file.is_open()) return false; 
    double x, y, z;  
    unsigned count {0}; 
    if (file.is_open()) {  
        for (auto& mol : Molecules) {
            for (auto& mono : mol.Monomers) {
                if (file >> x >> y >> z) {
                    mono.Velocity(0) = x; 
                    mono.Velocity(1) = y; 
                    mono.Velocity(2) = z; 
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
    if (!PBC) {
        std::cout << "periodic boundary conditions not activated!" << std::endl; 
        return;
    }
    for (auto& mol : Molecules) {
        wrapCOM(mol, BoxSize, Shear, delrx); 
    }
}

bool System::setDrive(Vector3d aForce, unsigned PartNumber) {
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            if (mono.Identifier == PartNumber) {
                Driven.push_back(Drive(aForce, &mono));  
                return true; 
            }
        }
    }
    std::cout << "Particle identifier does not exist!" << std::endl; 
    return false; 
}

bool System::changeDrive(Vector3d aForce, unsigned PartNumber) {
    for (auto& drivepair : Driven) {
        if (drivepair.DrivenParticle -> Identifier == PartNumber) {
            drivepair.Force = aForce; 
            return true; 
        }
    }
    std::cout << "Particle " << PartNumber << " was not set to be driven!" << std::endl;
    return false; 
}

bool System::setConstraint(double aConstraint, unsigned PartNumber) {
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            if (mono.Identifier == PartNumber) {
                Constraints.push_back(Constraint(aConstraint, &mono));  
                return true; 
            }
        }
    }
    std::cout << "Particle identifier does not exist!" << std::endl; 
    return false; 
}

bool System::changeConstraint(double aConstraint, unsigned PartNumber) {
    for (auto& constraintpair : Constraints) {
        if (constraintpair.ConstrainedParticle -> Identifier == PartNumber) {
            constraintpair.ConstraintPoint = aConstraint; 
            return true; 
        }
    }
    std::cout << "Particle " << PartNumber << " was not set to be constraint!" << std::endl;
    return false; 
}


/*void System::setSMDParticles(std::vector<unsigned> PartNumbers) {
    for (auto& number : PartNumbers) {
        if (number > Molecules[0].NumberOfMonomers) {
            std::cout << "SMD only works with 1 Molecule so far" << std::endl; 
            return;
        }
        SMDParticles.push_back(&Molecules[0].Monomers[number]); 
    }
}

void System::setupSMD(Vector3d pos, Vector3d vel) {
    SMDdrivingParticle.Position = pos; 
    SMDdrivingParticle.Velocity = vel; 
}*/

void System::setAnchor(unsigned PartNumber, Vector3d AnchorPos) {
    if (PartNumber > Molecules[0].NumberOfMonomers) {
        std::cout << "Anchoring currently only works with 1 Molecule" << std::endl; 
        return;  
    }
    Anchors.push_back(Anchor(AnchorPos, &Molecules[0].Monomers[PartNumber]));    
}

void System::setAnchorAuto(unsigned PartNumber) {
    if (PartNumber > Molecules[0].NumberOfMonomers) {
        std::cout << "Anchoring currently only works with 1 Molecule" << std::endl; 
        return;  
    }
    Vector3d AnchorPos (Molecules[0].Monomers[PartNumber].Position); 
    AnchorPos(2) = 0.0; 
    Anchors.push_back(Anchor(AnchorPos, &Molecules[0].Monomers[PartNumber]));    
}

void System::propagate(double dt, bool calcEpot) {
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5; 
            mono.Position += mono.Velocity*dt; 
            //TODO: boundary 
        }
    }
    updateCellLists();
    calculateForcesCellList(calcEpot);
    //checkVerletLists(); 
    //calculateForces(calcEpot); 
    /*if (SMD) {
        SMDdrivingParticle.Position += SMDdrivingParticle.Velocity*dt;
    }*/
    //calculateForcesBrute(calcEpot); 
    
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}

void System::propagateLangevin(double dt, double Temperature, double gamma, bool calcEpot) {
    double velcexp {exp(-gamma*dt)}; 
    double velcsqrt {sqrt(2.*gamma*Temperature)}; 
    double poscexp {(1.-velcexp)/(gamma)}; 
    double poscsqrt {velcsqrt/gamma}; 
    double tau1 {(1.-exp(-gamma*dt))/gamma}; 
    double tau2 {(1.-exp(-2.*gamma*dt))/(2.*gamma)}; 
    double zc11 {sqrt(tau2)};
    double zc21 {(tau1-tau2)/zc11}; 
    double zc22 {sqrt(dt - tau1*tau1/tau2)}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            Vector3d Z1; 
            Vector3d Z2; 
            Vector3d O1 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Vector3d O2 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Z1 = zc11*O1;
            Z2 = zc21*O1 + zc22*O2; 
            Vector3d veltilde {mono.Velocity + mono.Force*dt/(2.*mono.Mass)}; 
            mono.Velocity = velcexp*veltilde + velcsqrt*Z1/sqrt(mono.Mass); 
            mono.Position += poscexp*veltilde + poscsqrt*Z2/sqrt(mono.Mass); 
        }
    }
    /*if (SMD) {
        SMDdrivingParticle.Position += SMDdrivingParticle.Velocity*dt;
    }*/
    for (auto& constraintpair : Constraints) {
        constraintpair.ConstrainedParticle -> Position(2) = constraintpair.ConstraintPoint; 
    }
    calculateForcesBrute(calcEpot); 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}

/*void System::propagateLangevin(double dt, double Temperature, double gamma, bool calcEpot) { //Angel's version with gamma -> gamma*m 

    double velcsqrt {sqrt(2.*gamma*Temperature)}; 
    
    double poscsqrt {velcsqrt/gamma}; 
    double tau1 {(1-exp(-gamma*dt))/gamma}; 
    double tau2 {(1-exp(-2.*gamma*dt))/(2.*gamma)}; 
    double zc11 {sqrt(tau2)};
    double zc21 {(tau1-tau2)/zc11}; 
    double zc22 {sqrt(dt - tau1*tau1/tau2)}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            double velcexp {exp(-gamma*dt/mono.Mass)};
            double poscexp {(1-velcexp)*mono.Mass/(gamma)};  
            Vector3d Z1; 
            Vector3d Z2; 
            Vector3d O1 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Vector3d O2 {Rand::real_normal(), Rand::real_normal(), Rand::real_normal()}; 
            Z1 = zc11*O1;
            Z2 = zc21*O1 + zc22*O2; 
            Vector3d veltilde {mono.Velocity + mono.Force*dt/(2.*mono.Mass)}; 
            mono.Velocity = velcexp*veltilde + velcsqrt*Z1/mono.Mass; 
            mono.Position += poscexp*veltilde + poscsqrt*Z2; 
        }
    }
    calculateForcesBrute(calcEpot); 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            mono.Velocity += (mono.Force/mono.Mass)*dt*0.5;  
        }
    }
}*/


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

std::vector<double> System::calculateExtension(unsigned dim) {
    std::vector<double> Extension; 
    for (auto& mol : Molecules) {
        Extension.push_back(abs(mol.Monomers.front().Position(dim) - mol.Monomers.back().Position(dim))); 
        Extension.push_back((mol.Monomers.front().Position-mol.Monomers.back().Position).norm());
    }
    return Extension; 
}

void System::printPDB(FILE* pdb, int step, bool velocs) {
    wrapMoleculesCOM();
    int mol_count{0}; 
    fprintf(pdb, "MODEL     %d \n", step);
    for (auto& mol : Molecules) {
		mol_count++;
		int mono_count{0}; 
		for (auto& mono : mol.Monomers) {
		    if (velocs) {
		        fprintf(pdb, "ATOM %6d  C   GLY   %3d     %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n", mono_count+1, mol_count, mono.Position(0), mono.Position(1), mono.Position(2), mono.Velocity(0), mono.Velocity(1), mono.Velocity(2));
		    }
			else {
			    fprintf(pdb, "ATOM %6d  C   GLY   %3d     %7.3f %7.3f %7.3f \n", mono_count+1, mol_count, mono.Position(0), mono.Position(1), mono.Position(2));
		    }
		    mono_count++; 
		}	
        fprintf(pdb, "TER \n");
    }
    fprintf(pdb, "ENDMDL \n");   
    fflush(pdb);       
}

void System::printStatistics(std::ofstream& os, double time) {
    double Ekin {KineticEnergy()}; 
    Vector3d Omega {RotationFrequency()}; 
    std::tuple<double, Matrix3d> GyrTuple {GyrationTensor()};
    Matrix3d GyrTensor {std::get<1>(GyrTuple)}; 
    os.precision(10); 
    os.width(16); 
    os << time << " "; 
    os.precision(6); 
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


void System::printForceExtension(std::ofstream& os, double time, unsigned dim) {
    std::vector<double> Extension {calculateExtension(dim)}; 
    os.precision(10); 
    os.width(16); 
    os << time << " ";
    os.precision(6);  
    os.width(14); 
    if (Driven.empty()) os << 0.0 << " "; 
    else os << Driven.front().Force(dim) << " "; 
    for (auto& ext : Extension) {
        os.width(14);
        os << ext << " "; 
    }
    os << std::endl;  
}

void System::printExtensionForce(std::ofstream& os, double time, unsigned dim) {
    os.precision(10); 
    os.width(16); 
    os << time << " ";
    os.precision(6);  
    os.width(14); 
    if (Constraints.empty()) os << "NaN ";
    else os << Constraints.front().ConstraintPoint << " ";  
    os << Constraints.front().ConstrainedParticle -> Force(dim); 
    os << std::endl;  
}

    
