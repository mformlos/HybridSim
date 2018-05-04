#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <forward_list>
#include <../Eigen/Eigen/Dense>
#include "Potentials.h"
using namespace Eigen;


class Particle {
public: 
    //Members:
    Vector3d Position; 
    Vector3d Velocity; 
    double Mass; 
    
    //Constructor:
    Particle();  //Initialize everything to 0, mass to 1.0
    Particle(double);  //initialize only mass
    Particle(Vector3d, Vector3d, double); 
    ~Particle() = default;
    
    
}; 



class MPCParticle : public Particle {
public:
    int CellIndex;

    MPCParticle(); 
    MPCParticle(double); 
    MPCParticle(Vector3d, Vector3d, double); 
    ~MPCParticle() = default;
    bool operator < (const MPCParticle& part) const; 
};


class MDParticle : public MPCParticle {
public: 
    //Members: 
    int Identifier; 
    Vector3d Force; 
    Vector3d VerletPosition; 
    std::forward_list<MDParticle*> VerletList; 
    std::forward_list<MDParticle*> Bonds; 

    
    //Constructors:
    MDParticle(int);  //Initialize everything to 0 except mass; 
    MDParticle(int, double);
    MDParticle(int, Vector3d, Vector3d, double);  //Initialize Position and Velocity; 
    ~MDParticle() = default; 
    
    void setBond(MDParticle&);
}; 
#endif


