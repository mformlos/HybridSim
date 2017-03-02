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


 

