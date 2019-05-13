#ifndef LIB_MPC_H_
#define LIB_MPC_H_

#include <stdlib.h>
#include <fstream>
#include <list>
#include "Particle.h"
#include "Molecule.h"
#include "Rand.h"
#include "BoundaryConditions.h"

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
    unsigned Rho; 
    double STDNoHI; 
    Vector3d GridShift;
    std::vector<MPCParticle> Fluid; 
    std::vector<MPCParticle> Solute; 
    std::vector<std::forward_list<MPCParticle*>> CellList; 
    std::vector<CellMembers> CellData;  
    unsigned NumberOfCells; 
    unsigned NumberOfParticles; 
    std::array<unsigned,3> BoxSize;
    
    MPC(unsigned Lx, unsigned Ly, unsigned Lz, unsigned N_c, double T, double gamma); 
    MPC(unsigned Lx, unsigned Ly, unsigned Lz, unsigned N_c, double T, double gamma, bool HI);
    
    //Initialization: 
    void initializeRandom(); 
    void initializeProfile();
    bool initializeFile(std::string);
    
    //MPC routine 
    void updateBoxShift(double dt);
    void shiftGrid(); 
    //void stream(MPCParticle&, double dt); 
    //void stream(double dt); 
    void streamPlusCellAssignment(MPCParticle&, double);
    //void sort(); 
    void sortOnly();
    void sortOnly(std::vector<MDParticle>&);  
    void updateSoluteCellIndex(MPCParticle&);
    void updateSoluteCellIndexWithList(MPCParticle&, std::list<unsigned int>&); 
    void calculateCOMVel(unsigned); 
    void calculateCOMVelNoHI(unsigned); 
    void drawRotation(unsigned); 
    void rotate(unsigned); //rotate a particle
    void rotateSolute(unsigned); 
    void rotateSoluteNoHI(unsigned); 
    void getSolute(const std::vector<MDParticle>&); 
    void getSolute(const std::vector<Molecule>&); 
    void returnSolute(std::vector<Molecule>&); 
    void initializeSoluteVector(unsigned N); 
    void collisionNoHI(unsigned); 
    
     

    
    
    //void collide(unsigned, Vector3d); //collide in specific box
    //void collide(); //collide in all boxes

    void sortVector(); 
    
    
    //Boundary conditions: 
    //inline void wrap(MPCParticle&);
    inline bool isInBox(const Particle&); 
    
    
    //Properties: 
    double virtualTemperature(); 
    unsigned filledCells(); 
    Vector3d CenterOfMassVelocity(unsigned); 
    
    template<class UnaryFunc> 
    void operator() (UnaryFunc& ufunc) const {
        for (auto first = Fluid.cbegin(); first != Fluid.cend(); ++first) {
            ufunc(*first); 
        }        
    }
     
    
    void printFluid(FILE*, unsigned long); 


}; 

#endif
