#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "Particle.h"

class Molecule {
public: 

    unsigned NumberOfMonomers; 
    double Epot; 
    std::vector<MDParticle> Monomers; 
    
    Molecule(unsigned); //Standard Initialization of N Particles
    Molecule(unsigned, double); //Initialization of N Particles of mass M
    
    MDParticle& operator[](unsigned); //random access
    const MDParticle& operator[](unsigned) const; 
    void push_back(MDParticle&); 
    
    void setChainBonds(); 
    void setLink(unsigned, unsigned); 
    void translate(Vector3d); 
    void removeAngularMomentum(); 
    
    //Property getter functions
    Vector3d centerOfMassPosition(); 
    Vector3d centerOfMassVelocity(); 
    double KineticEnergy(); 
    double PotentialEnergy() {return Epot;}
    double radiusOfGyration(); 
    std::tuple<double, Matrix3d> GyrationTensor(); 
    Vector3d RotationFrequency(); 
    
    
}; 



#endif
