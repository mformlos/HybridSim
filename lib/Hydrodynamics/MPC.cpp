#include "MPC.h"

MPC::MPC(unsigned Lx, unsigned Ly, unsigned Lz, unsigned N_c, double aTemperature, double aShear) : 
    c {cos(130.0*M_PI/180.0)}, 
    s {sin(130.0*M_PI/180.0)},
    Temperature {aTemperature},
    Shear {aShear}, 
    delrx {0.0},
    GridShift {Vector3d::Zero() }  
    {
        BoxSize[0] = Lx; 
        BoxSize[1] = Ly;
        BoxSize[2] = Lz;
        NumberOfCells = Lx*Ly*Lz; 
        NumberOfParticles = N_c*NumberOfCells; 
        CellList = std::vector<std::forward_list<MPCParticle*>>(NumberOfCells, std::forward_list<MPCParticle*>()); 
        CellData = std::vector<CellMembers>(NumberOfCells, CellMembers());  
        Fluid.reserve(NumberOfParticles); 
        for (unsigned i = 0; i < NumberOfParticles; i++) {
            Fluid.push_back(MPCParticle(1.0)); 
        }  
        
    }
    
void MPC::initializeRandom() {
    Vector3d COMVel(Vector3d::Zero()); 
    double EKin { };
    double VelScaling { };  
    
    for (auto& part : Fluid) {
        for (unsigned dim = 0; dim < 3; dim++) {
            part.Position(dim) = BoxSize[dim]*(Rand::real_uniform()); 
            part.Velocity(dim) = Rand::real_uniform()-0.5; 
        }
        COMVel += part.Velocity; 
    }
    
    COMVel /= NumberOfParticles; 
    
    for (auto& part : Fluid) {
        part.Velocity -= COMVel; 
        EKin += part.Velocity.squaredNorm(); 
    }
    
    VelScaling = sqrt(3.*NumberOfParticles*Temperature/EKin); 
    
    for(auto& part : Fluid) {
        part.Velocity *= VelScaling; 
    }
}

/*double MPC::virtualTemperature() {
    double virtTemp { };
    for (auto& part : Fluid) {
        virtTemp += part.Velocity.squaredNorm(); 
    } 
    return (virtTemp/(3.*NumberOfParticles));    
}*/

double MPC::virtualTemperature() {
    double virtTemp { }; 
    for (auto& part : Fluid) {
        int Index {part.CellIndex}; 
        virtTemp += (part.Velocity - CellData[Index].CellCOMVel).squaredNorm();  
        //virtTemp += part.Velocity.squaredNorm();  
        
    }
    return (virtTemp/(3.*(NumberOfParticles-filledCells()))); 
}

inline unsigned MPC::filledCells() {
    unsigned count {}; 
    for (auto& data : CellData) {
        if (data.CellThermo) ++count; 
    }
    return count; 
}

void MPC::shiftGrid() {
    GridShift << Rand::real_uniform() -0.5, Rand::real_uniform()-0.5, Rand::real_uniform()-0.5;
}



void MPC::streamPlusCellAssignment(MPCParticle& part, double dt) {
    part.Position += part.Velocity*dt + GridShift; 
    wrap(part, BoxSize, Shear, delrx);
    part.CellIndex = (int)part.Position(2) + (int)part.Position(1)*BoxSize[2] + (int)part.Position(0)*BoxSize[2]*BoxSize[1]; 
}

void MPC::updateSoluteCellIndex(MPCParticle& solute) {
    solute.Position += GridShift; 
    wrap(solute, BoxSize, Shear, delrx);
    solute.CellIndex = (int)solute.Position(2) + (int)solute.Position(1)*BoxSize[2] + (int)solute.Position(0)*BoxSize[2]*BoxSize[1];
}


void MPC::sortOnly() {
    for (auto& Cell : CellList) {
        Cell.clear(); 
    }
    for (auto& part : Fluid) {
        try {
            CellList.at(part.CellIndex).push_front(&part); 
        }
        catch (std::exception& e) {
            std::cout << "The following exception occurred : " << e.what() << std::endl;
            std::cout << "for a particle with position "  << part.Position.transpose() << " and cell index " << part.CellIndex << std::endl;
            throw std::out_of_range("OOR");  
            return; 
        }    
    }
    for (auto& sol : Solute) {
        try {
            CellList.at(sol.CellIndex).push_front(&sol); 
        }
        catch (std::exception& e) {
            std::cout << "The following exception occurred : " << e.what() << std::endl;
            std::cout << "for a particle with position "  << sol.Position.transpose() << " and cell index " << sol.CellIndex << std::endl;
            throw std::out_of_range("OOR");  
            return; 
        }    
    }
}


Vector3d MPC::CenterOfMassVelocity(unsigned Index) {
    Vector3d COMVel(Vector3d::Zero()); 
    double totalMass { }; 
    for (auto& part : CellList[Index]) {
        COMVel += part -> Velocity * part -> Mass; 
        totalMass += part -> Mass; 
    }
    if (totalMass > 0) COMVel /= totalMass; 
    return COMVel; 
}


void MPC::sortVector() {
    std::sort(Fluid.begin(), Fluid.end()); 
}


void MPC::calculateCOMVel(unsigned Index) {
    double totalMass { }; 
    double Ekin { };
    unsigned count { }; 
    Vector3d COMVel {Vector3d::Zero()}; 
    for (auto& part : CellList[Index]) {
        COMVel += part -> Velocity * part -> Mass; 
        Ekin += part -> Mass *  (part -> Velocity).squaredNorm(); 
        totalMass += part -> Mass; 
        ++count; 
    }
    if (count > 0) CellData[Index].CellCOMVel = COMVel / totalMass;
    else CellData[Index].CellCOMVel = Vector3d::Zero();
    if (count < 2) CellData[Index].CellThermo = false; 
    else {
        CellData[Index].CellThermo = true; 
        Ekin = 0.5*(Ekin - totalMass*CellData[Index].CellCOMVel.squaredNorm());  
        CellData[Index].CellScaling = sqrt(Rand::real_gamma(1.5*(count-1), Temperature)/Ekin);     
    }
}

void MPC::drawRotation(unsigned Index) {
    double phi { };
	double theta { };
	Vector3d RotationAxis { };
	Matrix3d RotationMatrix { };
	phi = 2.*M_PI*(Rand::real_uniform());
	theta = 2.*(Rand::real_uniform()-0.5);
	RotationAxis(0) = sqrt(1-theta*theta)*cos(phi);
	RotationAxis(1) = sqrt(1-theta*theta)*sin(phi);
	RotationAxis(2) = theta;
	RotationMatrix(0,0) = RotationAxis(0)*RotationAxis(0) + (1 - RotationAxis(0)*RotationAxis(0))*c;
	RotationMatrix(0,1) = RotationAxis(0)*RotationAxis(1)*(1 - c) - RotationAxis(2)*s;
	RotationMatrix(0,2) = RotationAxis(0)*RotationAxis(2)*(1 - c) + RotationAxis(1)*s;
	RotationMatrix(1,0) = RotationAxis(0)*RotationAxis(1)*(1 - c) + RotationAxis(2)*s;
	RotationMatrix(1,1) = RotationAxis(1)*RotationAxis(1) + (1 - RotationAxis(1)*RotationAxis(1))*c;
	RotationMatrix(1,2) = RotationAxis(1)*RotationAxis(2)*(1 - c) - RotationAxis(0)*s;
	RotationMatrix(2,0) = RotationAxis(0)*RotationAxis(2)*(1 - c) - RotationAxis(1)*s;
	RotationMatrix(2,1) = RotationAxis(1)*RotationAxis(2)*(1 - c) + RotationAxis(0)*s;
	RotationMatrix(2,2) = RotationAxis(2)*RotationAxis(2) + (1 - RotationAxis(2)*RotationAxis(2))*c;
    CellData[Index].CellRotation = RotationMatrix;     
}

void MPC::rotate(unsigned i) {
    int Index{Fluid[i].CellIndex}; 
    Fluid[i].Velocity = CellData[Index].CellCOMVel + CellData[Index].CellRotation*(Fluid[i].Velocity - CellData[Index].CellCOMVel); 
    if (CellData[Index].CellThermo) Fluid[i].Velocity = CellData[Index].CellScaling * Fluid[i].Velocity + (1-CellData[Index].CellScaling)*CellData[Index].CellCOMVel; 
    Fluid[i].Position -= GridShift;
    wrap(Fluid[i],BoxSize, Shear, delrx); 
}

void MPC::rotateSolute(unsigned i) {
    int Index{Solute[i].CellIndex}; 
    Solute[i].Velocity = CellData[Index].CellCOMVel + CellData[Index].CellRotation*(Solute[i].Velocity - CellData[Index].CellCOMVel); 
    /*if (CellData[Index].CellThermo) Solute[i].Velocity = CellData[Index].CellScaling * Solute[i].Velocity + (1-CellData[Index].CellScaling)*CellData[Index].CellCOMVel;*/ 
}

void MPC::updateBoxShift(double dt) {
    delrx += Shear*BoxSize[1]*dt; 
    delrx -= BoxSize[0]*floor(delrx/BoxSize[0]); 
}

void MPC::getSolute(const std::vector<MDParticle>& sol) {
    vector<MPCParticle>::iterator it {Solute.begin()};
    for (auto& mono : sol) {
        *it = mono; 
        wrap(*it, BoxSize, Shear, delrx);
        ++it; 
    }
}

void MPC::getSolute(const std::vector<Molecule>& Molecules) {
    vector<MPCParticle>::iterator it {Solute.begin()}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) {
            *it = mono; 
            wrap(*it, BoxSize, Shear, delrx); 
            ++it; 
        }
    }
}

void MPC::returnSolute(std::vector<Molecule>& Molecules) {
    vector<MPCParticle>::iterator it {Solute.begin()}; 
    for (auto& mol : Molecules) {
        for (auto& mono : mol.Monomers) { 
            //std::cout << "vel before: " << mono.Velocity.transpose() << " ";
            wrapVelocityBack(mono, *it, BoxSize, Shear, delrx); 
            //std::cout << "vel after: " << mono.Velocity.transpose() << std::endl; 
            ++it; 
        }
    }
}

void MPC::initializeSoluteVector(unsigned N) {
    Solute.reserve(N); 
    for (unsigned i = 0; i < N; i++) {
        Solute.push_back(MPCParticle()); 
    }
} 

void MPC::printFluid(FILE* file, double time) {
    fprintf(file, "%14.6f \n", time); 
    for (auto& part : Fluid) {
        fprintf(file, "%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f \n", part.Position(0), part.Position(1), part.Position(2), part.Velocity(0), part.Velocity(1), part.Velocity(2)); 
    }
}



//obsolete functions: 
/*void MPC::stream(MPCParticle& part, double dt) {
    part.Position += part.Velocity*dt + GridShift; 
    wrap(part); 
}*/

/*void MPC::stream(double dt) {
    for (auto& part : Fluid) {
        stream(part, dt); 
    }
}*/

/*void MPC::sort() {
    for (auto& Cell : CellList) {
        Cell.clear(); 
    }

    for (auto& part : Fluid) {
        part.CellIndex = (int)part.Position(2) + (int)part.Position(1)*BoxSize[2] + (int)part.Position(0)*BoxSize[2]*BoxSize[1]; 
        try {
            CellList.at(part.CellIndex).push_front(&part); 
        }
        catch (std::exception& e) {
            std::cout << "The following exception occurred : " << e.what() << std::endl;
            std::cout << "for a particle with position "  << part.Position.transpose() << " and cell index " << part.CellIndex << std::endl;
            throw std::out_of_range("OOR");  
            return; 
        }
    }
}*/



/*void MPC::collide(unsigned Index, Vector3d COMVel) {
    double phi { };
	double theta { };
	Vector3d RotationAxis { };
	Matrix3d RotationMatrix { };
	phi = 2.*M_PI*(Rand::real_uniform());
	theta = 2.*(Rand::real_uniform()-0.5);
	RotationAxis(0) = sqrt(1-theta*theta)*cos(phi);
	RotationAxis(1) = sqrt(1-theta*theta)*sin(phi);
	RotationAxis(2) = theta;
	RotationMatrix(0,0) = RotationAxis(0)*RotationAxis(0) + (1 - RotationAxis(0)*RotationAxis(0))*c;
	RotationMatrix(0,1) = RotationAxis(0)*RotationAxis(1)*(1 - c) - RotationAxis(2)*s;
	RotationMatrix(0,2) = RotationAxis(0)*RotationAxis(2)*(1 - c) + RotationAxis(1)*s;
	RotationMatrix(1,0) = RotationAxis(0)*RotationAxis(1)*(1 - c) + RotationAxis(2)*s;
	RotationMatrix(1,1) = RotationAxis(1)*RotationAxis(1) + (1 - RotationAxis(1)*RotationAxis(1))*c;
	RotationMatrix(1,2) = RotationAxis(1)*RotationAxis(2)*(1 - c) - RotationAxis(0)*s;
	RotationMatrix(2,0) = RotationAxis(0)*RotationAxis(2)*(1 - c) - RotationAxis(1)*s;
	RotationMatrix(2,1) = RotationAxis(1)*RotationAxis(2)*(1 - c) + RotationAxis(0)*s;
	RotationMatrix(2,2) = RotationAxis(2)*RotationAxis(2) + (1 - RotationAxis(2)*RotationAxis(2))*c;
    
    for (auto& part : CellList[Index]) {
        part -> Velocity = COMVel + RotationMatrix*(part -> Velocity - COMVel); 
    }
} */

/*inline void MPC::wrap(MPCParticle& part) {
    double cy {floor(part.Position(1)/BoxSize[1])}; 
    part.Position(0) -= BoxSize[0]*floor(part.Position(0)/BoxSize[0]); 
    part.Position(0) -= cy*delrx; 
    part.Position(0) -= BoxSize[0]*floor(part.Position(0)/BoxSize[0]); 
    part.Position(1) -= BoxSize[1]*cy; 
    part.Position(2) -= BoxSize[2]*floor(part.Position(2)/BoxSize[2]); 
    part.Velocity(0) -= cy*Shear*BoxSize[1]; 
}*/


