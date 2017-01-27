#ifndef LIB_SYSTEM_H
#define LIB_SYSTEM_H 

#include <string>
#include <fstream>
#include "Molecule.h"
#include "Potentials.h"
#include "Rand.h"
#include "BoundaryConditions.h"


class System {
public: 
    double Cutoff; 
    double VerletRadius; 
    double VerletRadiusSq; 
    
    std::array<unsigned,3> BoxSize; 
    std::array<unsigned,3> Cells; 
    std::array<double,3> CellSideLength; 
    
    std::vector<std::vector<std::vector<std::forward_list<MDParticle*>>>> CellList; 
    
    std::vector<Molecule> Molecules; 
    
    System(unsigned, unsigned, unsigned); //initialize only boxsize; 
    System(double, double, unsigned, unsigned, unsigned); //initialize with cutoffs
    
    void updateVerletLists(); 
    void checkVerletLists(); 
    
    void calculateForces(); 
    
    // Initialize Molecules 
    void addMolecules(std::string); 
    void addLinks(std::string);  
    void initializePositions(std::string); 
    void initializeVelocitiesRandom(double); 
    
    void propagate(double dt); 

};

#endif
