#ifndef LIB_SYSTEM_H
#define LIB_SYSTEM_H 

#include <string>
#include <fstream>
#include "Molecule.h"
#include "Potentials.h"
#include "Rand.h"
#include "BoundaryConditions.h"
#include "../Exceptions/LibraryException.h"

struct Anchor {
    Vector3d AnchorPoint; 
    MDParticle* AnchoredParticle; 
    
    Anchor(Vector3d point, MDParticle* part) :
        AnchorPoint {point}, 
        AnchoredParticle {part} {}
};

struct Drive {
    Vector3d Force; 
    MDParticle* DrivenParticle; 
    
    Drive(Vector3d f, MDParticle* driven) : 
        Force {f}, 
        DrivenParticle {driven} {}
};


class System {
public: 
    double Cutoff; 
    double VerletRadius; 
    double VerletRadiusSq; 
    double delrx; 
    double Shear;
    double SurfaceEnergy; 
    bool Adsorption; 
    
    std::array<unsigned,3> BoxSize; 
    std::array<unsigned,3> Cells; 
    std::array<double,3> CellSideLength; 
    
    std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>> CellList; 
    
    std::vector<Molecule> Molecules; 
    
    
    std::vector<Anchor> Anchors; 
    std::vector<Drive> Driven; 
    
    
    
    System(unsigned, unsigned, unsigned, double, double SurfE = 1.0, bool AdsorptionOn = false); //initialize only boxsize and Shear; 
    System(double, double, unsigned, unsigned, unsigned, double, double SurfE = 1.0, bool AdsorptionOn = false); //initialize with cutoffs
    
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
    bool setDrive(Vector3d, unsigned); 
    bool changeDrive(Vector3d, unsigned); 
    void setAnchor(unsigned, Vector3d); 
    void setAnchorAuto(unsigned); 
    
    void propagate(double dt, bool calcEpot=false); 
    void propagateLangevin(double dt, double Temperature, double gamma=0.05, bool calcEpot=false); 
    
    //getter:
    
    unsigned NumberOfParticles(); 
    unsigned NumberOfMolecules(); 
    double KineticEnergy(); 
    double PotentialEnergy(); 
    std::tuple<double, Matrix3d> GyrationTensor(); 
    Vector3d RotationFrequency(); 
    std::vector<double> calculateExtension(unsigned); 
    
    
    void printPDB(FILE* pdb, int step, bool velocs = false); 
    void printStatistics(std::ofstream& os, double time); 
    void printForceExtension(std::ofstream& os, double time, unsigned dim); 
};


#endif
