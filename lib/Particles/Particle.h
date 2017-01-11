#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <iostream>
#include <forward_list>
#include <../Eigen/Eigen/Dense>
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
    MPCParticle(); 
    MPCParticle(double); 
    MPCParticle(Vector3d, Vector3d, double); 
    ~MPCParticle() = default;
};


class MDParticle : public MPCParticle {
public: 
    //Members: 
    Vector3d Force; 
    Vector3d VerletPosition; 
    std::forward_list<MDParticle*> VerletList; 
    
    //Constructors:
    MDParticle();  //Initialize everything to 0 except mass; 
    MDParticle(double);
    MDParticle(Vector3d, Vector3d, double);  //Initialize Position and Velocity; 
    ~MDParticle() = default; 
}; 
#endif


