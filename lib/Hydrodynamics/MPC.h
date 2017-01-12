#ifndef LIB_MPC_H_
#define LIB_MPC_H_

#include <stdlib.h>
#include "Particle.h"
#include "Rand.h"

class MPC {
public: 
    double c; 
    double s; 
    double Temperature; 
    std::vector<MPCParticle> Fluid; 
    std::vector<std::forward_list<MPCParticle*>> CellList; 
    unsigned NumberOfCells; 
    unsigned NumberOfParticles; 
    std::array<unsigned,3> BoxSize;
    
    MPC(unsigned Lx, unsigned Ly, unsigned Lz, unsigned N_c, double); 
    
    //Initialization: 
    void initialize_random(); 
    
    
    //MPC routine 
    inline void stream(MPCParticle&, double dt); 
    void stream(double dt); 
    void sort(); 
    void collide(unsigned, Vector3d); //collide in specific box
    void collide(); //collide in all boxes
    bool sortByPosition(MPCParticle, MPCParticle); 
    void sortVector(); 
    
    
    //Boundary conditions: 
    inline void wrap(MPCParticle&);
    
    
    //Properties: 
    double virtualTemperature(); 
    Vector3d CenterOfMassVelocity(unsigned); 
    
    bool operator() (MPCParticle, MPCParticle); 


}; 

#endif
