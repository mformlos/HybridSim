#include <vector>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <math.h>
#include "Particle.h"
#include "HelperFunctions.h"
#include "BoundaryConditions.h"

int main(int argc, char* argv[]) {
    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, KabsFile{}, FileName_x {}, FileName_y {}, FileName_z {}; 
    int StartStep{}, Step{}, SamplingStep{};
    unsigned NumberOfMonomers {};
    double MDStep {0.01}, Shear {}, delrx; 
  
    std::array<unsigned, 3> BoxSize {};
    std::vector<Particle> Monomers{}; 
    
    std::vector<unsigned> KabsList {}; 
    std::map<double, double> structure_factor_x {};
    std::map<double, double> structure_factor_y {};
    std::map<double, double> structure_factor_z {};
    std::map<double, int> q_count {};

    if (argc != 10) {
        std::cout << "usage: ./structure_factor DIRECTORY NUMBEROFMONOMERS BOXSIZE_X BOXSIZE_Y BOXSIZE_Z STARTSTEP SAMPLINGSTEP KABSLIST SHEAR" << std::endl;
        return EXIT_FAILURE; 
    }

    Directory=argv[1];
	NumberOfMonomers = std::stoi(argv[2]);
	BoxSize[0] = std::stoi(argv[3]); 
	BoxSize[1] = std::stoi(argv[4]); 
	BoxSize[2] = std::stoi(argv[5]); 
	StartStep = std::stoi(argv[6]);
	SamplingStep = std::stoi(argv[7]);
	KabsFile = argv[8]; 
	Shear = std::stod(argv[9]); 
	
	std::cout << "Directory: " << Directory << std::endl;
	std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl;
	std::cout << "Number of monomers: " << NumberOfMonomers << "  BoxSize: " << BoxSize[0] << " " << BoxSize[1] << " " << BoxSize[2] <<  std::endl; 
	std::cout << "Kabs File: " << KabsFile << std::endl;
	 
	for (unsigned i = 0; i < NumberOfMonomers; i++) {
        Monomers.push_back(Particle(i)); 
    }
    
    if (!initializeVector(KabsList, KabsFile)) {
        std::cout << "problem filling KabsList, exiting!" << std::endl; 
        return EXIT_FAILURE; 
    }   
    
	FileName_x = Directory+"/data/structure_factor_x";
	std::ifstream test(FileName_x);
	if (test.good()){
		std::cout << "file " << FileName_x << " already exists! Aborting..." << std::endl;
		return EXIT_FAILURE;
	}
	std::ofstream Output_x(FileName_x, std::ios::out | std::ios::trunc);
	Output_x.precision(8);
	test.close(); 
	
	FileName_y = Directory+"/data/structure_factor_y";
	test.open(FileName_y);
	if (test.good()){
		std::cout << "file " << FileName_y << " already exists! Aborting..." << std::endl;
		return EXIT_FAILURE;
	}
	std::ofstream Output_y(FileName_y, std::ios::out | std::ios::trunc);
	Output_y.precision(8);
	test.close(); 
	
	FileName_z = Directory+"/data/structure_factor_z";
	test.open(FileName_z);
	if (test.good()){
		std::cout << "file " << FileName_z << " already exists! Aborting..." << std::endl;
		return EXIT_FAILURE;
	}
	std::ofstream Output_z(FileName_z, std::ios::out | std::ios::trunc);
	Output_z.precision(8);
	
    Vector3d KVec {}; 
    unsigned count_configs{0};
    double cossq_x {}, sinsq_x {}, cossq_y {}, sinsq_y {}, cossq_z {}, sinsq_z {};
    double sum_x {}, sum_y {}, sum_z {} ;  
    double constant_x {2.*M_PI/BoxSize[0]}, constant_y {2.*M_PI/BoxSize[1]}, constant_z {2.*M_PI/BoxSize[2]}; 
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    ////// MAIN CALCULATION //////
    
    Step = StartStep;
	ConfigFileStart = Directory+"/configs/config"; 
	while (true) {
	    ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
	    //std::cout << ConfigFile << std::endl; 
	    if (!initializePositions(Monomers, ConfigFile)) {
	        std::cout << "problem with initializing monomers" << std::endl;
	        std::cout << ConfigFile << std::endl; 
		    break; 
		}
		delrx = Step*MDStep*Shear*BoxSize[1];
		delrx -= BoxSize[0]*floor(delrx/BoxSize[0]); 
		for (auto& mono : Monomers) {
	        wrap(mono, BoxSize, Shear, delrx);
		} 
		for (auto& Kabs : KabsList) {
		    cossq_x = cossq_y = cossq_z = 0.0; 
            sinsq_x = sinsq_y = sinsq_z = 0.0;   
            for (auto& mono : Monomers) {
                cossq_x += cos(Kabs*mono.Position(0)*constant_x); 
                sinsq_x += sin(Kabs*mono.Position(0)*constant_x); 
                cossq_y += cos(Kabs*mono.Position(1)*constant_y); 
                sinsq_y += sin(Kabs*mono.Position(1)*constant_y); 
                cossq_z += cos(Kabs*mono.Position(2)*constant_z); 
                sinsq_z += sin(Kabs*mono.Position(2)*constant_z); 
                cossq_x += cos(-Kabs*mono.Position(0)*constant_x); 
                sinsq_x += sin(-Kabs*mono.Position(0)*constant_x); 
                cossq_y += cos(-Kabs*mono.Position(1)*constant_y); 
                sinsq_y += sin(-Kabs*mono.Position(1)*constant_y); 
                cossq_z += cos(-Kabs*mono.Position(2)*constant_z); 
                sinsq_z += sin(-Kabs*mono.Position(2)*constant_z); 
            }     
		    cossq_x = pow(cossq_x, 2); 
            sinsq_x = pow(sinsq_x, 2); 
            cossq_y = pow(cossq_y, 2); 
            sinsq_y = pow(sinsq_y, 2); 
            cossq_z = pow(cossq_z, 2); 
            sinsq_z = pow(sinsq_z, 2); 
            sum_x = (cossq_x + sinsq_x)/NumberOfMonomers; 
            sum_y = (cossq_y + sinsq_y)/NumberOfMonomers;              
            sum_z = (cossq_z + sinsq_z)/NumberOfMonomers;
            structure_factor_x[Kabs*constant_x] += sum_x;
            structure_factor_y[Kabs*constant_y] += sum_y;
            structure_factor_z[Kabs*constant_z] += sum_z;	 
		}
		count_configs++;
		Step += SamplingStep;
	}
	
	gettimeofday(&end, NULL); 
    double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
    std::cout << "total time: " << realTime << " , time per conformation: " << realTime/count_configs << " , total conformations: " << count_configs <<std::endl;
    
    for (auto& w : structure_factor_x) {
        Output_x << w.first << " " << w.second/count_configs << std::endl; 
    }
    Output_x.close();
    
    for (auto& w : structure_factor_y) {
        Output_y << w.first << " " << w.second/count_configs << std::endl; 
    }
    Output_y.close();
     
	for (auto& w : structure_factor_z) {
        Output_z << w.first << " " << w.second/count_configs << std::endl; 
    }
    Output_z.close();
    
    return 0; 
}
