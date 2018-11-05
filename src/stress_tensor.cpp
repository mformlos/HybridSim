#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "Molecule.h"

int main(int argc, char* argv[]) {

    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, BondFile{};
    double StartTime{}, EndTime{}, DeltaTSim{0.01};
    int StartStep{}, EndStep{}, Step{}, SamplingStep{}, N{}; 
    Molecule Mol(200); 
    Matrix3d StressAverage{Matrix3d::Zero()};
    Matrix3d StressCurrent{Matrix3d::Zero()};
    Matrix3d StressSTD{Matrix3d::Zero()}; 
    Vector3d COMPos; 
    
    if (argc != 5) {
            std::cout << "usage: ./stress_tensor DIRECTORY STARTSTEP ENDSTEP SAMPLINGSTEP" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    EndStep = std::stoi(argv[3]);
    SamplingStep = std::stoi(argv[4]); 
    
    std::cout << "Directory: " << Directory << std::endl; 
    std::cout << "StartStep: " << StartStep << "   EndStep: " << EndStep << "   SamplingStep: " << SamplingStep << std::endl; 
    
    N = ((EndStep-StartStep)/SamplingStep);
    
    StartTime = StartStep*DeltaTSim;
    EndTime = EndStep*DeltaTSim; 
    
    ConfigFileStart = Directory+"/configs/config";
    BondFile = Directory+"/bonds"; 
    
    Mol.setChainBonds();
    Mol.addLinks(BondFile, 1); 
    
    
    std::string name; 
    name = Directory+"/data/tau_xx"; 
    std::ifstream test(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_xx(name, std::ios::out | std::ios::trunc);  
    tau_xx.precision(8); 
    
    name = Directory+"/data/tau_xy"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_xy(name, std::ios::out | std::ios::trunc);  
    tau_xy.precision(8);
    
    name = Directory+"/data/tau_xz";
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    }  
    std::ofstream tau_xz(name, std::ios::out | std::ios::trunc);  
    tau_xz.precision(8);
    
    name = Directory+"/data/tau_yx"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_yx(name, std::ios::out | std::ios::trunc);  
    tau_yx.precision(8);
    
    name = Directory+"/data/tau_yy"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_yy(name, std::ios::out | std::ios::trunc);  
    tau_yy.precision(8);
    
    name = Directory+"/data/tau_yz"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_yz(name, std::ios::out | std::ios::trunc);  
    tau_yz.precision(8);
    
    name = Directory+"/data/tau_zx"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_zx(name, std::ios::out | std::ios::trunc);  
    tau_zx.precision(8);
    
    name = Directory+"/data/tau_zy"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_zy(name, std::ios::out | std::ios::trunc);  
    tau_zy.precision(8);
    
    name = Directory+"/data/tau_zz"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_zz(name, std::ios::out | std::ios::trunc);  
    tau_zz.precision(8);
    
    
    
    Step = StartStep; 
    
    while (Step <= EndStep) {
        ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
        Mol.initializePositions(ConfigFile); 
        //Mol.calculateInternalForces(); 
        Mol.calculateSpringForces(); 
        StressCurrent = Mol.StressTensor(); 
        //StressAverage += StressCurrent; 
        //StressSTD += StressCurrent.cwiseProduct(StressCurrent); 
        Step += SamplingStep; 
        tau_xx << Step << " " << StressCurrent(0,0) << std::endl;
        tau_xy << Step << " " << StressCurrent(0,1) << std::endl;
        tau_xz << Step << " " << StressCurrent(0,2) << std::endl;
        tau_yx << Step << " " << StressCurrent(1,0) << std::endl;
        tau_yy << Step << " " << StressCurrent(1,1) << std::endl;
        tau_yz << Step << " " << StressCurrent(1,2) << std::endl;
        tau_zx << Step << " " << StressCurrent(2,0) << std::endl;
        tau_zy << Step << " " << StressCurrent(2,1) << std::endl;
        tau_zz << Step << " " << StressCurrent(2,2) << std::endl;
        
        
    }
    
    tau_xx.close();
    tau_xy.close();
    tau_xz.close();
    tau_yx.close();
    tau_yy.close();
    tau_yz.close();
    tau_zx.close();
    tau_zy.close();
    tau_zz.close();
    
    //StressAverage /= N;
    //StressSTD /= N;
    
    //StressSTD -= StressAverage.cwiseProduct(StressAverage); 
    //StressSTD = StressSTD.array().sqrt(); 
    
    
    
    /*for (unsigned i = 0; i < 3; i++) {
        for (unsigned j = i; j < 3; j++) {
            OutputStream << StressAverage(i,j) << " " << StressSTD(i,j) << std::endl; 
        }
    } */  
    
    //OutputStream.close(); 
     
    return 1; 
}
