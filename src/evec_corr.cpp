#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <../Eigen/Eigen/Dense>

using namespace Eigen; 

int main(int argc, char* argv[]) {

    std::string Directory{}, EvecFileName{}, OutputFileName{};  
    std::fstream EvecFile {}, OutputFile{}; 
    double DeltaT{}, Time{}; 
    int StartStep{}, MaxStep{}, DStep{}, LastStep{}, SampleStep; 
    std::map<double,double> Corr_evec; 
    std::map<double,double>::iterator Corr_evec_iter; 
    std::map<double, unsigned> Corr_evec_count; 
    std::map<double, unsigned>::iterator Corr_evec_count_iter; 
    std::map<unsigned, Vector3d> Evec; 
    std::map<unsigned, Vector3d>::iterator Evec_iter; 
    
    if (argc != 7) {
            std::cout << "usage: ./evec_corr DIRECTORY STARTSTEP MAXSTEP DSTEP SAMPLESTEP TIMESTEP" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    MaxStep = std::stoi(argv[3]); 
    DStep = std::stoi(argv[4]); 
    SampleStep = std::stoi(argv[5]); 
    DeltaT = std::stod(argv[6]); 
    
    
    EvecFileName = Directory+"/data/eigenvalues"; 
    EvecFile.open(EvecFileName, std::ios::in);
    
    OutputFileName = Directory + "/data/corr_evec";
    std::ifstream test (OutputFileName); 
    if (test.good()) {
        std::cout << "File " << OutputFileName << " already exists! Aborting ..." << std::endl; 
        return EXIT_FAILURE; 
    }
    OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc); 
    
    int step; 
    double ex, ey, ez, bla; 
    while (EvecFile >> step >> bla >> bla >> bla >> ex >> ey >> ez >>  bla >> bla >> bla >> bla >> bla >> bla ) {
        Evec[step](0) = ex;
        Evec[step](1) = ey; 
        Evec[step](2) = ez;  
        bla = Evec[step].squaredNorm();
        if (abs(bla-1.0) > 1e-05) {
            std::cout <<  bla << " " << Evec[step].transpose() << std::endl; 
        }
        LastStep = step;  
    }
    
    LastStep = step;  

    int T0Step {StartStep}; 
    int T1Step {};
    int DTStep {};
    for (int i = 0; i < MaxStep; i+= DStep) {
        double dt {i*DeltaT}; 
        Corr_evec[dt] = 0.0; 
        Corr_evec_count[dt] = 0;
    }
    
    while (T0Step <= LastStep) {  
        T1Step = T0Step; 
        DTStep = 0; 
        //std::cout << "T0Step = " << T0Step << std::endl; 
        while (DTStep <= MaxStep && T1Step <= LastStep) {
            Time = DTStep*DeltaT;
            //std::cout << "T1Step = " << T1Step << std::endl; 
            Corr_evec[Time] += Evec[T0Step].dot(Evec[T1Step]); 
            Corr_evec_count[Time]++; 
            T1Step += DStep;
            DTStep += DStep;  
        }
        T0Step += SampleStep; 
    }
    for (auto& corr : Corr_evec) {
        corr.second /= Corr_evec_count[corr.first]; 
    }

    for (auto& corr : Corr_evec) {
        OutputFile << corr.first << " " << corr.second << " " << Corr_evec_count[corr.first] << std::endl; 
    }
    
    OutputFile.close(); 
    

    return 0; 
}
