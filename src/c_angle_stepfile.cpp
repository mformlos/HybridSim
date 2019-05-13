#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include "Molecule.h"
#include "HelperFunctions.h"

int main(int argc, char* argv[]) {
    std::string DirectoryStart{}, Directory{}, EvecFileName{}, ConfigFileName{}, ConfigFileStart{}, OutputFileName{}, StepFile{};  
    std::fstream EvecFile {}, OutputFile{}; 
    double StartTime{}, MaxTime{}, DTime{}, DeltaT{0.01}, Time{}; 
    unsigned N {0}, Monomers{200};
    unsigned long long StartStep{}, MaxStep{}, DStep{}, LastStep{}, SampleStep; 
    //std::vector<unsigned, std::map<double,double>> C_angle(Monomers, std::map<double,double>);
    std::map<double,double> C_angle; 
    std::map<double, unsigned> C_angle_count;
    std::map<unsigned long long, Vector3d> Evec{}; 
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    std::vector<double> A0 (Monomers, 0.0); 
    std::vector<double> At(Monomers, 0.0);
    
     
    if (argc != 6) {
            std::cout << "usage: ./c_angle_stepfile DIRECTORY STARTSTEP DELTAT STEPFILE SAMPLESTEP" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]);
    DeltaT = std::stod(argv[3]); 
    StepFile = argv[4]; 
    SampleStep = std::stoi(argv[5]);
    initializeStepVector(StepVector, StepFile); 
    
    MaxStep = StepVector.back(); 
    LastStep = StepVector.back(); 
    
    StartTime = StartStep*DeltaT; 
    MaxTime = MaxTime*DeltaT; 
    
    std::cout << "Directory: " << Directory << std::endl; 
    std::cout << "StartStep: " << StartStep << std::endl; 
    std::cout << "DeltaT: " << DeltaT << std::endl; 
    std::cout << "StepFile: " << StepFile << std::endl; 
    std::cout << "SampleStep: " << SampleStep << std::endl; 
    
    OutputFileName = Directory+"/data/c_angle_new";
    /*std::ifstream test (OutputFileName); 
    if (test.good()) {
        std::cout << "File " << OutputFileName << " already exists! Aborting ..." << std::endl; 
        return EXIT_FAILURE; 
    }*/
   
    
    EvecFileName = Directory+"/data/eigenvalues"; 
    EvecFile.open(EvecFileName, std::ios::in); 
    
    unsigned step; 
    double ex, ey, ez, d; 
    while (EvecFile >> step >> d >> d >> d >> ex >> ey >> ez >> d >> d >> d >> d >> d >> d) {
        if (step >= StartStep) {
            Vector3d ev {ex, ey, ez}; 
            if (ev(0) < 0.0) {
                ev(0) = -ev(0); 
                ev(1) = -ev(1); 
                ev(2) = -ev(2); 
            }
            Evec[step] = ev; 
        }
    }
    EvecFile.close(); 
    
    
    /*for (auto st = StepVector.begin(); st != StepVector.end(); st++) {
        if (*st >= StartStep) {
            StepVectorIterator = st; 
            break;
        }
    }
    StepVector.erase(StepVector.begin(), StepVectorIterator); 
    */
    
    Molecule MolFirst(Monomers);
    Molecule MolSecond(Monomers);  
    Vector3d COMFirst; 
    Vector3d COMSecond; 
    
    ConfigFileStart = Directory+"/configs/config";
    
    unsigned long long T0Step {StartStep}; 
    unsigned long long T1Step {};
    unsigned long long DTStep {};
    
    
    
    
    for (auto& st : StepVector) {
        C_angle[st*DeltaT] = 0.0; 
        C_angle_count[st*DeltaT] = 0; 
    }
    
    StepVectorIterator = StepVector.begin();  
    
    //std::cout << "T0Step: " << T0Step << std::endl
    Vector3d Relative{}; 
    double dotprod {}; 
    double Amean{0.0};
    while (T0Step <= LastStep) { 
        StepVectorIterator = StepVector.begin();  
        /*for (auto st = StepVector.begin(); st != StepVector.end(); st++) {
            if (*st == T0Step) {
                StepVectorIterator = st; 
                break;
            }
            
        }*/
        //StepVector.erase(StepVector.begin(), StepVectorIterator);
        ConfigFileName = ConfigFileStart+std::to_string(T0Step)+".pdb";
        if (!(MolFirst.initializePositions(ConfigFileName))) {
            std::cout << "File " << ConfigFileName << " does not exist! " << std::endl; 
            break; 
        }
        COMFirst = MolFirst.centerOfMassPosition(); 
        
        for (unsigned m = 0; m < Monomers; m++) {
            Relative = MolFirst.Monomers[m].Position - COMFirst; 
            dotprod = Relative.dot(Evec[T0Step])/(Relative.norm()*Evec[T0Step].norm());
            //A0[m] = 2*dotprod*sqrt(1-dotprod*dotprod); 
            A0[m] = sin(2*acos(dotprod)); 
            //Amean += A0[m]; 
        }
        std::cout << "T0Step: " << T0Step << std::endl; 
        std::cout << "T1Step: " << *StepVectorIterator+T0Step << std::endl; 
        while (StepVectorIterator != StepVector.end()) {
            T1Step = *StepVectorIterator+T0Step; 
            //std::cout << T1Step << std::endl; 
            DTStep = T1Step-T0Step; 
            Time = DTStep*DeltaT;
            if ( !(C_angle.find(Time) == C_angle.end()) ) {
                
                ConfigFileName = ConfigFileStart+std::to_string(T1Step)+".pdb";
                if (!(MolSecond.initializePositions(ConfigFileName))) {
                    //std::cout << "File " << ConfigFileName << " does not exist! " << std::endl; 
                    break; 
                }
                COMSecond = MolSecond.centerOfMassPosition(); 
                for (unsigned m = 0; m < Monomers; m++) {
                    Relative = MolSecond.Monomers[m].Position - COMSecond;
                    dotprod = Relative.dot(Evec[T1Step])/(Relative.norm()*Evec[T1Step].norm());
                    //At[m] = 2*dotprod*sqrt(1-dotprod*dotprod);
                    At[m] = sin(2*acos(dotprod));
                    C_angle.at(Time) += A0[m]*At[m];  
                }
                C_angle_count.at(Time)++;
            }
            StepVectorIterator++; 
        }   
        T0Step += SampleStep; 
    }
    OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc); 
    
    for (auto& A : C_angle) { 
        //A.second /= (Monomers*Amean*Amean); 
        A.second /= Monomers; 
        OutputFile << A.first << " " << A.second << " " << C_angle_count[A.first] << std::endl; 
    }
    OutputFile.close(); 
    
    return 0; 
}
