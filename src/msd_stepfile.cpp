#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "Molecule.h"
#include "HelperFunctions.h"

int main(int argc, char* argv[]) {

    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, StepFile{};  
    double DeltaTSim{0.01}, CurrentDt{}; 
    unsigned long long StartStep{}, CurrentDtStep{}, Monomers{}, FileStep{};   
    std::map<double,double> msd; 
    std::map<double,double> msd_x; 
    std::map<double,double> msd_y;
    std::map<double,double> msd_z;  
    std::map<double,double>::iterator msd_iter; 
    std::map<double,double> msd_count; 
    std::map<double,double>::iterator msd_count_iter; 
   
    
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    
    if (argc != 6) {
            std::cout << "usage: ./msd DIRECTORY STARTSTEP STEPFILE MONOMERS DELTAT" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    StepFile = argv[3];
    Monomers = std::stoi(argv[4]);
    DeltaTSim = std::stod(argv[5]); 
    
    std::cout << "Directory: " << Directory << " StepFile: " << StepFile << std::endl;
	std::cout << "StartStep: " << StartStep << std::endl;
    
    Molecule MolFirst(Monomers);
    Molecule MolSecond(Monomers);  
    Vector3d COMFirst; 
    Vector3d COMSecond; 
    
    initializeStepVector(StepVector, StepFile); 
    for(auto& s : StepVector){
        s -= StartStep;
        //std::cout << s << std::endl; 
    } 
    
    

    ConfigFileStart = Directory+"/configs/config";
    
    std::string MSDout = Directory+"/data/MSD_total"; 
    std::ifstream test (MSDout); 
    if (test.good()) {
        std::cout << "File " << MSDout << " already exists! Aborting ..." << std::endl; 
        return EXIT_FAILURE; 
    }
    std::ofstream MSDoutStream(MSDout, std::ios::out | std::ios::trunc); 
    
   
    Vector3d Relative{}; 
    
    ConfigFile = ConfigFileStart+std::to_string(StartStep)+".pdb";
    if (!(MolFirst.initializePositions(ConfigFile))) {
        std::cout << "File " << ConfigFile << " does not exist! " << std::endl; 
        return EXIT_FAILURE; 
    }
    COMFirst = MolFirst.centerOfMassPosition(); 
    
    StepVectorIterator = StepVector.begin(); 
    
    while (StepVectorIterator != StepVector.end()) {
        CurrentDtStep = *StepVectorIterator; 
        //std::cout << CurrentDtStep << std::endl; 
	    CurrentDt = CurrentDtStep*DeltaTSim;
	    FileStep = StartStep+CurrentDtStep;
	    ConfigFile = ConfigFileStart+std::to_string(FileStep)+".pdb";
	    if (!(MolSecond.initializePositions(ConfigFile))) {
	        std::cout << FileStep << " , reached last step" << std::endl;
	        break; 
	    }  
        COMSecond = MolSecond.centerOfMassPosition(); 
        for (unsigned i = 0; i < MolFirst.NumberOfMonomers; i++) {
            Relative = MolFirst.Monomers[i].Position - COMFirst - (MolSecond.Monomers[i].Position - COMSecond);
            msd[CurrentDt] += Relative.squaredNorm();
            msd_x[CurrentDt] += pow(Relative(0),2); 
            msd_y[CurrentDt] += pow(Relative(1),2);
            msd_z[CurrentDt] += pow(Relative(2),2);                 
        }
        msd_count[CurrentDt] += MolFirst.NumberOfMonomers; 
        StepVectorIterator++;
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
        
    MSDoutStream.close();
    MSDout = Directory+"/data/MSD_x";
    MSDoutStream.open(MSDout, std::ios::out | std::ios::trunc); 
    
    msd_iter = msd_x.begin(); 
    msd_count_iter = msd_count.begin(); 
    printing = true; 
    do {
        MSDoutStream.precision(8);
	    MSDoutStream << std::scientific;
        MSDoutStream << msd_iter-> first << " " << (msd_iter -> second) / msd_count_iter -> second << " \n"; 
        MSDoutStream << std::flush; 
        auto buf = msd_iter; 
        if (++buf == msd_x.end()) printing = false; 
        ++msd_iter; 
        ++msd_count_iter;   
    } while(printing);    
    
    MSDoutStream.close();
    MSDout = Directory+"/data/MSD_y";
    MSDoutStream.open(MSDout, std::ios::out | std::ios::trunc); 
    
    msd_iter = msd_y.begin(); 
    msd_count_iter = msd_count.begin(); 
    printing = true;
    do {
        MSDoutStream.precision(8);
	    MSDoutStream << std::scientific;
        MSDoutStream << msd_iter-> first << " " << (msd_iter -> second) / msd_count_iter -> second << " \n"; 
        MSDoutStream << std::flush; 
        auto buf = msd_iter; 
        if (++buf == msd_y.end()) printing = false; 
        ++msd_iter; 
        ++msd_count_iter;   
    } while(printing);    

    MSDoutStream.close();
    MSDout = Directory+"/data/MSD_z";
    MSDoutStream.open(MSDout, std::ios::out | std::ios::trunc); 
    
    msd_iter = msd_z.begin(); 
    msd_count_iter = msd_count.begin(); 
    printing = true;
    do {
        MSDoutStream.precision(8);
	    MSDoutStream << std::scientific;
        MSDoutStream << msd_iter-> first << " " << (msd_iter -> second) / msd_count_iter -> second << " \n"; 
        MSDoutStream << std::flush; 
        auto buf = msd_iter; 
        if (++buf == msd_z.end()) printing = false; 
        ++msd_iter; 
        ++msd_count_iter;   
    } while(printing);    
    
    MSDoutStream.close();
    
    return EXIT_SUCCESS; 
}

