
#include <iostream>
#include <vector>
#include "Particle.h"
#include <../Eigen/Eigen/Dense>

int main() {
    
    std::vector<MDParticle> Solute; 
    std::vector<MPCParticle> CopiedSolute; 
    
    Solute.push_back(MDParticle(0,1.)); 
    Solute.front().Position(0) = 1.0;
    Solute.front().Position(1) = 2.0;
    Solute.front().Position(2) = 3.0; 
    Solute.front().Velocity(0) = -1.;
    Solute.front().Velocity(1) = -2.;
    Solute.front().Velocity(2) = -3.;

    CopiedSolute.push_back(MPCParticle()); 
    
    std::cout << "before any copying: " << std::endl; 
    std::cout << "Solute: " << Solute.front().Position.transpose() << " " << Solute.front().Velocity.transpose() <<std::endl;
    std::cout << "Copy: " << CopiedSolute.front().Position.transpose() << " " << CopiedSolute.front().Velocity.transpose() <<std::endl;
    
    std::vector<MPCParticle>::iterator it {CopiedSolute.begin()}; 
    
    for (auto& mono : Solute) *it = mono; 
    
    std::cout << "after copying: " << std::endl; 
    std::cout << "Solute: " << Solute.front().Position.transpose() << " " << Solute.front().Velocity.transpose() <<std::endl;
    std::cout << "Copy: " << CopiedSolute.front().Position.transpose() << " " << CopiedSolute.front().Velocity.transpose() <<std::endl;
    
    it -> Velocity *= 2.; 
    
    
    std::cout << "after multiplication: " << std::endl; 
    std::cout << "Solute: " << Solute.front().Position.transpose() << " " << Solute.front().Velocity.transpose() <<std::endl;
    std::cout << "Copy: " << CopiedSolute.front().Position.transpose() << " " << CopiedSolute.front().Velocity.transpose() <<std::endl;
    
    
    
    return 0;
}
