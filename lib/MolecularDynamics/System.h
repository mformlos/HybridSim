#ifndef LIB_SYSTEM_H
#define LIB_SYSTEM_H 

#include <string>
#include <fstream>
#include "Molecule.h"
#include "Potentials.h"
#include "Rand.h"
#include "BoundaryConditions.h"
#include "../Exceptions/LibraryException.h"

class System {
public: 
    double Cutoff; 
    double VerletRadius; 
    double VerletRadiusSq; 
    double delrx; 
    double Shear;
    double SMDconstant; 
    bool SMD; 
    
    std::array<unsigned,3> BoxSize; 
    std::array<unsigned,3> Cells; 
    std::array<double,3> CellSideLength; 
    
    std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>> CellList; 
    
    std::vector<Molecule> Molecules; 
    
    std::vector<MDParticle*> SMDParticles; 
    Particle SMDdrivingParticle; 
    
    System(unsigned, unsigned, unsigned, double, bool SMDon = false, double k = 0.0); //initialize only boxsize and Shear; 
    System(double, double, unsigned, unsigned, unsigned, double, bool SMDon = false, double k = 0.0); //initialize with cutoffs
    
    void updateVerletLists(); 
    void checkVerletLists(); 
    
    void calculateForces(bool calcEpot=false); 
    void calculateForcesBrute(bool calcEpot=false); 
    
    void wrapMoleculesCOM(); 
    
    // Initialize Molecules 
    bool addMolecules(std::string, double mass = 1.0); 
    bool addLinks(std::string);  
    bool initializePositions(std::string); 
    void initializeVelocitiesRandom(double); 
    void setMoleculeCOM(unsigned, Vector3d); 
    void centerMolecule(unsigned); 
    void setSMDParticles(std::vector<unsigned>); 
    void setupSMD(Vector3d, Vector3d); 
    
    void propagate(double dt, bool calcEpot=false); 
    void propagateLangevin(double dt, double Temperature, double gamma=0.05, bool calcEpot=false); 
    
    //getter:
    
    unsigned NumberOfParticles(); 
    unsigned NumberOfMolecules(); 
    double KineticEnergy(); 
    double PotentialEnergy(); 
    std::tuple<double, Matrix3d> GyrationTensor(); 
    Vector3d RotationFrequency(); 
    
    void printPDB(FILE* pdb, int step, bool velocs = false); 
    void printStatistics(std::ofstream& os, double time); 
};

#endif
