#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "Molecule.h"


int main(int argc, char* argv[]) {

    std::string Directory{}, ConfigFile{}, ConfigFileStart{};  
    double StartTime{}, EndTime{}, MaxMSDTime{}, SamplingTime{}, DeltaTSim{0.01}, DeltaMSDTime {}, CurrentDt{}; 
    int StartStep{}, EndStep{}, MaxMSDStep{}, SamplingStep{}, CurrentSamplingStep{}, CurrentDtStep{}, DeltaMSDStep{};   
    std::map<double,double> msd; 
    std::map<double,double>::iterator msd_iter; 
    std::map<double,double> msd_count; 
    std::map<double,double>::iterator msd_count_iter; 
    Molecule MolFirst(200);
    Molecule MolSecond(200);  
    Vector3d COMFirst; 
    Vector3d COMSecond; 
    
    if (argc != 6) {
            std::cout << "usage: ./msd DIRECTORY STARTSTEP ENDSTEP MAXMSDSTEP SAMPLINGSTEP DELTAMSDSTEP" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    EndStep = std::stoi(argv[3]); 
    MaxMSDStep = std::stoi(argv[4]); 
    SamplingStep = std::stoi(argv[5]); 
    DeltaMSDStep = std::stoi(argv[6]);
    
    StartTime = StartStep*DeltaTSim;
    EndTime = EndStep*DeltaTSim; 
    MaxMSDTime = MaxMSDStep*DeltaTSim; 
    SamplingTime = SamplingStep*DeltaTSim; 
    DeltaMSDTime = DeltaMSDStep*DeltaTSim; 
    
    CurrentSamplingStep = StartStep; 
    ConfigFileStart = Directory+"/configs/config";
    
    std::string MSDout = Directory+"/data/MSD"; 
    std::ifstream test (MSDout); 
    if (test.good()) {
        std::cout << "File " << MSDout << " already exists! Aborting ..." << std::endl; 
        return EXIT_FAILURE; 
    }
    std::ofstream MSDoutStream(MSDout, std::ios::out | std::ios::trunc); 
    
   
    
    
    while (CurrentSamplingStep < EndStep) {
        //std::cout << CurrentSamplingStep << std::endl; 
        CurrentDtStep = DeltaMSDStep; 
        CurrentDt = DeltaMSDTime;     
        ConfigFile = ConfigFileStart+std::to_string(CurrentSamplingStep)+".pdb";
        //std::cout << "CurrentFirstFile: " << ConfigFile << std::endl; 
        if (!(MolFirst.initializePositions(ConfigFile))) return EXIT_FAILURE; 
        COMFirst = MolFirst.centerOfMassPosition(); 
        while(CurrentDtStep <= MaxMSDStep && CurrentSamplingStep+CurrentDtStep <= EndStep) {
            ConfigFile = ConfigFileStart+std::to_string(CurrentSamplingStep+CurrentDtStep)+".pdb";
            //std::cout << "CurrentSecondFile: " << ConfigFile << std::endl;
            if (!(MolSecond.initializePositions(ConfigFile))) return EXIT_FAILURE;  
            COMSecond = MolSecond.centerOfMassPosition(); 
            for (unsigned i = 0; i < MolFirst.NumberOfMonomers; i++) {
                msd[CurrentDt] += (MolFirst.Monomers[i].Position - COMFirst - (MolSecond.Monomers[i].Position - COMSecond)).squaredNorm();
            }
            msd_count[CurrentDt] += MolFirst.NumberOfMonomers; 
            CurrentDtStep += DeltaMSDStep;
            CurrentDt +=  DeltaMSDTime;
        }
        CurrentSamplingStep += SamplingStep;   
    }  
    
    bool printing {true}; 
    msd_iter = msd.begin(); 
    msd_count_iter = msd_count.begin(); 
    do {
        MSDoutStream.precision(8);
	    MSDoutStream << std::scientific;
        MSDoutStream << msd_iter-> first << " " << (msd_iter -> second) / msd_count_iter -> second << " \n"; 
        MSDoutStream << std::flush; 
        auto buf = msd_iter; 
        if (++buf == msd.end()) printing = false; 
        ++msd_iter; 
        ++msd_count_iter;   
    } while(printing);
        
    



    return EXIT_SUCCESS; 
}

