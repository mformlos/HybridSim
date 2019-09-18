#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include "HelperFunctions.h"


int main(int argc, char* argv[]) {

    std::string DirectoryStart{}, Directory{}, RgyrFileName{}, OutputFileName{}, StepFile{};  
    std::fstream RgyrFile {}, OutputFile{}; 
    double StartTime{}, MaxTime{}, DTime{}, DeltaT{0.01}, Time{}; 
    int StartStep{}, MaxStep{}, DStep{}, LastStep{}, SampleStep; 
    std::vector<unsigned long long> StepVector{}; 
    std::vector<unsigned long long>::iterator StepVectorIterator{}; 
    unsigned Molecules{}; 
    
    if (argc != 5) {
            std::cout << "usage: ./rgyr_corr_multiple DIRECTORY STARTSTEP MOLECULES STEPFILE" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]);
    Molecules = std::stoi(argv[3]); 
    StepFile = argv[4]; 
    initializeStepVector(StepVector, StepFile); 
    
    
    MaxStep = StepVector.back(); 
    SampleStep = MaxStep; 
    //MaxStep = std::stoi(argv[3]); 
    //DStep = std::stoi(argv[4]); 
    //SampleStep = std::stoi(argv[5]); 
    
    StartTime = StartStep*DeltaT; 
    MaxTime = MaxTime*DeltaT; 
    //DTime = DStep*DeltaT;
    
    std::cout << "Directory: " << Directory << std::endl; 
    std::cout << "StartStep: " << StartStep << std::endl; 
    std::cout << "DeltaT: " << DeltaT << std::endl; 
    std::cout << "StepFile: " << StepFile << std::endl; 
    std::cout << "SampleStep: " << SampleStep << std::endl; 
   

    for (unsigned Mol = 0; Mol < Molecules; Mol++) { 
	    std::map<double,double> Corr_rg; 
	    std::map<double,double>::iterator Corr_rg_iter; 
	    std::map<double, unsigned> Corr_rg_count; 
	    std::map<double, unsigned>::iterator Corr_rg_count_iter; 
	    std::map<unsigned, double> Rg; 
	    std::map<unsigned, double>::iterator Rg_iter; 
	    double Rg_mean{0.0}, Rg_sqmean{0.0}, Rg_quadmean{0.0};
	    std::map<double,double> Corr_rg2; 
	    std::map<double, unsigned> Corr_rg2_count; 
	    unsigned N {0}; 
	    
	    OutputFileName = Directory+"/data/corrgyr-Mol-"+std::to_string(Mol);
	    std::ifstream test (OutputFileName); 
	    if (test.good()) {
		std::cout << "File " << OutputFileName << " already exists! Aborting ..." << std::endl; 
		return EXIT_FAILURE; 
	    }
	    OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc); 
	    
	    
	    RgyrFileName = Directory+"/data/stats-Mol-"+std::to_string(Mol);
	    RgyrFile.open(RgyrFileName, std::ios::in);
	    int step; 
	    double r, bla; 
	    while (RgyrFile >> step >> r >> bla >> bla >> bla >> bla  >> bla >> bla ) {
		    if (step >= StartStep && (step-StartStep)%SampleStep == 0) {
			Rg_mean += r; 
			Rg_sqmean += r*r; 
			N++;
		    }
		    Rg[step] = r;
	    }
            RgyrFile.close(); 
	    
    
	    
	    Rg_mean /= N; 
	    Rg_sqmean /= N; 
	    Rg_quadmean /= N; 
	    
	    
	    std::cout << "Rg mean = " << Rg_mean << std::endl; 
	    std::cout << "Rg^2 mean = " << Rg_sqmean << std::endl; 
	    std::cout << "N = " << N << std::endl; 

	    Corr_rg[0.0] = 0.0; 
	    for (auto& step : StepVector) {
		double dt {step*DeltaT}; 
		Corr_rg[dt] = 0.0; 
		Corr_rg_count[dt] = 0;
	    } 
	     
	     
	LastStep = step; 
	std::cout << "REPL " << Mol << std::endl;  
	std::cout << "Laststep = " << LastStep << std::endl; 
    
	int T0Step {StartStep}; 
	int T1Step {};
	int DTStep {};

    
	while (T0Step <= LastStep) { 
	    StepVectorIterator = StepVector.begin();  
	    T1Step = T0Step+(*StepVectorIterator); 
	    DTStep = *StepVectorIterator; 
	    //std::cout << "T0Step = " << T0Step << std::endl; 
	    while (StepVectorIterator != StepVector.end() && T1Step <= LastStep) {
		
		Time = DTStep*DeltaT;
		//std::cout << "T1Step = " << T1Step << std::endl; 
		try {
		    Corr_rg[Time] += Rg.at(T0Step)*Rg.at(T1Step); //(Gxx[T0Step]-Gxx_mean)*(Gyy[T1Step]-Gyy_mean); 
		}
		catch(...) {
		    std::cout << "T0Step = " << T0Step << "T1Step = " << T1Step << std::endl;
		    return EXIT_FAILURE;
		}   
		Corr_rg2[Time] += (pow(Rg[T0Step],2)- Rg_sqmean)*(pow(Rg[T1Step],2)-Rg_sqmean); 
		Corr_rg_count[Time]++; 
		StepVectorIterator++;   
 		T1Step = T0Step+(*StepVectorIterator); 
	        DTStep = *StepVectorIterator; 

	    }
	    T0Step += SampleStep; 
	}
       
    
	    
	    for (auto& corr : Corr_rg) {
		corr.second /= Corr_rg_count[corr.first]; 
		corr.second -= (Rg_mean*Rg_mean); 
		corr.second /= (Rg_sqmean-Rg_mean*Rg_mean); 
	    }
	    
	    for (auto& corr : Corr_rg2) {
		corr.second /= Corr_rg_count[corr.first]; 
		corr.second /= (Rg_quadmean - pow(Rg_sqmean,2));
	    }
	    
	    
	    
	    for (auto& corr : Corr_rg) {
		OutputFile << corr.first << " " << corr.second << std::endl; 
	    }
	    
	    OutputFile.close(); 
	    
    }
    return 0;
    
        
} 
    
