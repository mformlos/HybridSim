#include "Molecule.h" 

Molecule::Molecule(unsigned N) :
    NumberOfMonomers {N}, 
    Epot { } 
    {
        Monomers.reserve(N); 
        for (unsigned i = 0; i < NumberOfMonomers; i++) {
            Monomers.push_back(MDParticle()); 
        }
    }

Molecule::Molecule(unsigned N, double Mass) :
    NumberOfMonomers {N},
    Epot { }
    {
        Monomers.reserve(N); 
        for (unsigned i = 0; i < NumberOfMonomers; i++) {
            Monomers.push_back(MDParticle(Mass)); 
        }
    }

MDParticle& Molecule::operator[](unsigned i) {return Monomers[i];}

const MDParticle& Molecule::operator[](unsigned i) const {return Monomers[i];}

void Molecule::push_back(MDParticle& part) {
    ++NumberOfMonomers; 
    Monomers.push_back(part); 
}

void Molecule::setChainBonds() {
    for (unsigned i = 0; i < NumberOfMonomers-1; i++) {
        Monomers[i].setBond(Monomers[i+1]); 
    }
}

void Molecule::setLink(unsigned first, unsigned second) {
    Monomers[first].setBond(Monomers[second]); 
}

Vector3d Molecule::centerOfMassPosition() {
    Vector3d COMPos {Vector3d::Zero()}; 
    for (auto& mono : Monomers) {
        COMPos += mono.Position; 
    }
    COMPos /= (double)NumberOfMonomers; 
    return COMPos; 
}

Vector3d Molecule::centerOfMassVelocity() {
    Vector3d COMVel {Vector3d::Zero()}; 
    for (auto& mono : Monomers) {
        COMVel += mono.Velocity; 
    }
    COMVel /= NumberOfMonomers; 
    return COMVel; 
}

void Molecule::translate(Vector3d vec) {
    for (auto& mono : Monomers) {
        mono.Position += vec; 
    }
}

void Molecule::removeAngularMomentum() {
    Vector3d COMPos {centerOfMassPosition()}; 
    Vector3d omega {RotationFrequency()}; 
    for (auto& mono : Monomers) {
        mono.Velocity += (mono.Position-COMPos).cross(omega); 
    }
}

double Molecule::KineticEnergy() {
    double Ekin {0.0}; 
    for (auto& mono : Monomers) {
        Ekin += mono.Mass*mono.Velocity.squaredNorm(); 
    }
    Ekin /= 2; 
    return Ekin; 
}

double Molecule::radiusOfGyration() {
    double rgyr {}; 
    for (unsigned i = 0; i < NumberOfMonomers; i++) {
        for (unsigned j = i+1; j < NumberOfMonomers; j++) {
            Vector3d distance {Monomers[i].Position - Monomers[j].Position}; 
            rgyr += distance.squaredNorm(); 
        }
    }
    rgyr /= pow(NumberOfMonomers,2); 
    rgyr = sqrt(rgyr); 
    return rgyr; 
}

std::tuple<double, Matrix3d> Molecule::GyrationTensor() {
    Matrix3d gyrTensor {Matrix3d::Zero()}; 
    Vector3d COMPos {centerOfMassPosition()}; 
    double rgyr { }; 
    for (auto& mono : Monomers) {
        Vector3d relPos {mono.Position-COMPos}; 
        for (unsigned alpha = 0; alpha < 3; alpha++) {
            for (unsigned beta = alpha; beta < 3; beta++) {
                gyrTensor(alpha, beta) += relPos(alpha)*relPos(beta); 
            }
        }
    }
    gyrTensor /= NumberOfMonomers; 
    gyrTensor(1,0) = gyrTensor(0,1);
    gyrTensor(2,0) = gyrTensor(0,2);
    gyrTensor(2,1) = gyrTensor(1,2); 
    rgyr = gyrTensor(0,0) + gyrTensor(1,1) + gyrTensor(2,2);
    rgyr = sqrt(rgyr); 
    return std::make_tuple(rgyr, gyrTensor);  
}

Vector3d Molecule::RotationFrequency() {
    Vector3d omega {Vector3d::Zero()}; 
    Matrix3d inertiaTensor {Matrix3d::Zero()}; 
    Vector3d angularMomentum {Vector3d::Zero()};
    Vector3d COMPos {centerOfMassPosition()}; 
    Vector3d COMVel {centerOfMassVelocity()};  
    for (auto& mono : Monomers) {   
        Vector3d relPos {mono.Position - COMPos}; 
        Vector3d relVel {mono.Velocity - COMVel}; 
        double rsqr {relPos.squaredNorm()};
        for (unsigned alpha = 0; alpha < 3; alpha++) {
            inertiaTensor(alpha, alpha) += rsqr; 
            for (unsigned beta = alpha; beta < 3; beta++) {
                inertiaTensor(alpha, beta) -= relPos(alpha)*relPos(beta);       
            }
        }
        angularMomentum += relPos.cross(relVel); 
    }
    inertiaTensor(1,0) = inertiaTensor(0,1); 
    inertiaTensor(2,0) = inertiaTensor(0,2); 
    inertiaTensor(2,1) = inertiaTensor(1,2); 
    omega = inertiaTensor.ldlt().solve(angularMomentum); 
    return omega; 
}
 

