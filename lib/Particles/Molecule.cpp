#include "Molecule.h" 

Molecule::Molecule(unsigned N) :
    NumberOfMonomers {N} 
    {
        Monomers.reserve(N); 
        for (unsigned i = 0; i < NumberOfMonomers; i++) {
            Monomers.push_back(MDParticle()); 
        }
    }

Molecule::Molecule(unsigned N, double Mass) :
    NumberOfMonomers {N} 
    {
        Monomers.reserve(N); 
        for (unsigned i = 0; i < NumberOfMonomers; i++) {
            Monomers.push_back(MDParticle(Mass)); 
        }
    }

MDParticle& Molecule::operator[](unsigned i) {return Monomers[i];}

const MDParticle& Molecule::operator[](unsigned i) const {return Monomers[i];}

void Molecule::setChainBonds() {
    for (unsigned i = 0; i < NumberOfMonomers-1; i++) {
        Monomers[i].setBond(Monomers[i+1]); 
    }
}

void Molecule::setLink(unsigned first, unsigned second) {
    Monomers[first].setBond(Monomers[second]); 
}
 

