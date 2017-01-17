#include "MPC.h"

MPC::MPC(unsigned Lx, unsigned Ly, unsigned Lz, unsigned N_c, double aTemperature) : 
    c {cos(130.0*M_PI/180.0)}, 
    s {sin(130.0*M_PI/180.0)},
    Temperature(aTemperature) 
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
    
void MPC::initialize_random() {
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

double MPC::virtualTemperature() {
    double virtTemp { };
    for (auto& part : Fluid) {
        virtTemp += part.Velocity.squaredNorm(); 
    } 
    return (virtTemp/(3.*NumberOfParticles));    
}

void MPC::stream(MPCParticle& part, double dt) {
    part.Position += part.Velocity*dt; 
    wrap(part); 
} 

void MPC::streamPlusCellAssignment(MPCParticle& part, double dt) {
    part.Position += part.Velocity*dt; 
    wrap(part);
    part.CellIndex = (int)part.Position(0) + (int)part.Position(1)*BoxSize[0] + (int)part.Position(2)*BoxSize[0]*BoxSize[1]; 
}

void MPC::stream(double dt) {
    for (auto& part : Fluid) {
        stream(part, dt); 
    }
}

inline void MPC::wrap(MPCParticle& part) {
    for (unsigned dim = 0; dim < 3; dim++) {
        part.Position(dim) -= BoxSize[dim]*floor(part.Position(dim)/BoxSize[dim]); 
    }
}

void MPC::sort() {
    for (auto& Cell : CellList) {
        Cell.clear(); 
    }

    for (auto& part : Fluid) {
        part.CellIndex = (int)part.Position(0) + (int)part.Position(1)*BoxSize[0] + (int)part.Position(2)*BoxSize[0]*BoxSize[1]; 
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
}



void MPC::updateParticleCellIndex() {
    for (auto& part : Fluid) {
        part.CellIndex = (int)part.Position(0) + (int)part.Position(1)*BoxSize[0] + (int)part.Position(2)*BoxSize[0]*BoxSize[1];
    } 
}

void MPC::collide(unsigned Index, Vector3d COMVel) {
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
    Vector3d COMVel {Vector3d::Zero()}; 
    for (auto& part : CellList[Index]) {
        COMVel += part -> Velocity * part -> Mass; 
        totalMass += part -> Mass; 
    }
    if (totalMass > 0) CellData[Index].CellCOMVel = COMVel / totalMass;
    else CellData[Index].CellCOMVel = Vector3d::Zero();  
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
}
