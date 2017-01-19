#ifndef LIB_MPC_H_
#define LIB_MPC_H_

#include <stdlib.h>
#include "Particle.h"
#include "Rand.h"

struct CellMembers{
    Vector3d CellCOMVel;
    double CellScaling;  
    Matrix3d CellRotation;
    bool CellThermo;  
    
    CellMembers() : 
        CellCOMVel {Vector3d::Zero()}, 
        CellScaling {0.0}, 
        CellRotation {Matrix3d::Zero()},
        CellThermo {false} {}
};

class MPC {
public: 
    double c; 
    double s; 
    double Temperature;
    double Shear;
    double delrx;  
    Vector3d GridShift;
    std::vector<MPCParticle> Fluid; 
    std::vector<std::forward_list<MPCParticle*>> CellList; 
    std::vector<CellMembers> CellData;  
    unsigned NumberOfCells; 
    unsigned NumberOfParticles; 
    std::array<unsigned,3> BoxSize;
    
    MPC(unsigned Lx, unsigned Ly, unsigned Lz, unsigned N_c, double T, double gamma); 
    
    //Initialization: 
    void initialize_random(); 
    
    
    //MPC routine 
    void updateBoxShift(double dt);
    void stream(MPCParticle&, double dt); 
    void stream(double dt); 
    void streamPlusCellAssignment(MPCParticle&, double);
    void sort(); 
    void sortOnly(); 
    void updateParticleCellIndex(); 
    void calculateCOMVel(unsigned); 
    void drawRotation(unsigned); 
    void rotate(unsigned); //rotate a particle

    
    void collide(unsigned, Vector3d); //collide in specific box
    void collide(); //collide in all boxes

    void sortVector(); 
    
    
    //Boundary conditions: 
    inline void wrap(MPCParticle&);
    
    
    //Properties: 
    double virtualTemperature(); 
    
    Vector3d CenterOfMassVelocity(unsigned); 
     


}; 

#endif
