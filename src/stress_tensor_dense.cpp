#include <vector>
#include <sys/stat.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "System.h"
#include "HelperFunctions.h"
int main(int argc, char* argv[]) {

    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, BondFile{},  MoleculeFileName{}, NeighbourCellFile{};
    double StartTime{}, EndTime{}, DeltaTSim{0.01}, Shear{}, delrx{};
    int StartStep{}, EndStep{}, Step{}, SamplingStep{}, N{}, NumberOfMolecules{}; 
    Molecule Mol(200); 
    unsigned Lx{}, Ly{}, Lz{};
    Matrix3d StressAverage{Matrix3d::Zero()};
    Matrix3d StressCurrent{Matrix3d::Zero()};
    Matrix3d StressSTD{Matrix3d::Zero()}; 
    Vector3d COMPos; 
    bool ParameterInitialized{};    
    if (argc != 7) {
            std::cout << "usage: ./stress_tensor DIRECTORY STARTSTEP ENDSTEP SAMPLINGSTEP MOLECULEFILE PARAMETERFILE" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    EndStep = std::stoi(argv[3]);
    SamplingStep = std::stoi(argv[4]); 
    MoleculeFileName = argv[5];

    std::cout << "Directory: " << Directory << std::endl; 
    std::cout << "StartStep: " << StartStep << "   EndStep: " << EndStep << "   SamplingStep: " << SamplingStep << std::endl; 
   
    std::ifstream inputfile(argv[6], ios::in);
    if (!inputfile.is_open()) {
        std::cout << "could not open file '" << argv[1] << "' , exiting" << std::endl;
        return EXIT_FAILURE;
    }
 
    Lx = extractParameter<unsigned>("BoxX", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    Ly = extractParameter<unsigned>("BoxY", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    Lz = extractParameter<unsigned>("BoxZ", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    Shear = extractParameter<double>("Shear", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    std::cout << "Box size: " << Lx << " " << Ly << " " << Lz << std::endl;  
  
    System Sys(1.2, 1.5, Lx, Ly, Lz, Shear, 0.0,false, true);
    NeighbourCellFile = Directory+"/neighbourcells";

    if (!Sys.setNeighbourDirections(NeighbourCellFile)){
        return EXIT_FAILURE;
    }

    if (!Sys.addMolecules(MoleculeFileName, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    NumberOfMolecules = Sys.Molecules.size(); 
    
    std::cout << "Number of Molecules: " << NumberOfMolecules <<std::endl; 
    
    
    N = ((EndStep-StartStep)/SamplingStep);
    
    StartTime = StartStep*DeltaTSim;
    EndTime = EndStep*DeltaTSim; 
    
    ConfigFileStart = Directory+"/configs/config";
    BondFile = Directory+"/bonds"; 
    
    if (!Sys.addLinks(BondFile)) {
        std::cout << "problem with initializing bonds" << std::endl;
        return EXIT_FAILURE; 
    }
    
    std::string outputdir=Directory+"/data/stress_tensor_list";
    
    const int dir_err = mkdir(outputdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        std::cout << "Error creating directory!" << std::endl;
        exit(1);
    }
    
    std::string name; 
    name = Directory+"/data/stress_tensor_list/tau_xx"; 
    std::ifstream test(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_xx(name, std::ios::out | std::ios::trunc);  
    tau_xx.precision(8); 
    
    name = Directory+"/data/stress_tensor_list/tau_xy"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_xy(name, std::ios::out | std::ios::trunc);  
    tau_xy.precision(8);
    
    name = Directory+"/data/stress_tensor_list/tau_xz";
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    }  
    std::ofstream tau_xz(name, std::ios::out | std::ios::trunc);  
    tau_xz.precision(8);
    
    name = Directory+"/data/stress_tensor_list/tau_yx"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_yx(name, std::ios::out | std::ios::trunc);  
    tau_yx.precision(8);
    
    name = Directory+"/data/stress_tensor_list/tau_yy"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_yy(name, std::ios::out | std::ios::trunc);  
    tau_yy.precision(8);
    
    name = Directory+"/data/stress_tensor_list/tau_yz"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_yz(name, std::ios::out | std::ios::trunc);  
    tau_yz.precision(8);
    
    name = Directory+"/data/stress_tensor_list/tau_zx"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_zx(name, std::ios::out | std::ios::trunc);  
    tau_zx.precision(8);
    
    name = Directory+"/data/stress_tensor_list/tau_zy"; 
    test.open(name); 
    if (test.good()){
        std::cout << "file " << name << " already exists! Aborting..." << std::endl; 
        return EXIT_FAILURE; 
    } 
    std::ofstream tau_zy(name, std::ios::out | std::ios::trunc);  
    tau_zy.precision(8);
    
    name = Directory+"/data/stress_tensor_list/tau_zz"; 
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
        if (!Sys.initializePositionsPDB(ConfigFile)) {
	        std::cout << Step << ", problem with initializing monomers" << std::endl;
	        Step += SamplingStep; 
	        continue; 
	    }
    	delrx = (Step+1)*DeltaTSim*Shear*Ly;
        delrx -= Lx*floor(delrx/Lx);
	Sys.delrx = delrx;
	Sys.wrapMoleculesCOM();
	Vector3d COM{Vector3d::Zero()};
	for (auto& Mol : Sys.Molecules) {
            for (auto& Mono : Mol.Monomers) {
	        COM += Mono.Position;
	    }
	} 
	COM /= Sys.NumberOfParticles();
        //Mol.calculateInternalForces(); 
        StressCurrent = Matrix3d::Zero(); 
        //Sys.calculateForcesBrute();
        Sys.updateCellLists();
    	/*Sys.calculateForcesCellList();
	for (auto& Mol : Sys.Molecules) {
	    for (auto& Mono : Mol.Monomers) {
            	for (unsigned i = 0; i < 3; i++) {
            		for (unsigned j = 0; j < 3; j++) {
                		StressCurrent(i, j) += (Mono.Position(i)-COM(i))*Mono.Force(j);
			}
		}
            }

		//Mol.calculateSpringForces(); 
            //StressCurrent += Mol.StressTensor(); 
        }*/
	StressCurrent = Sys.calculateStressTensor();
        //StressAverage += StressCurrent; 
        //StressSTD += StressCurrent.cwiseProduct(StressCurrent); 
        
        tau_xx << Step << " " << StressCurrent(0,0) << std::endl;
        tau_xy << Step << " " << StressCurrent(0,1) << std::endl;
        tau_xz << Step << " " << StressCurrent(0,2) << std::endl;
        tau_yx << Step << " " << StressCurrent(1,0) << std::endl;
        tau_yy << Step << " " << StressCurrent(1,1) << std::endl;
        tau_yz << Step << " " << StressCurrent(1,2) << std::endl;
        tau_zx << Step << " " << StressCurrent(2,0) << std::endl;
        tau_zy << Step << " " << StressCurrent(2,1) << std::endl;
        tau_zz << Step << " " << StressCurrent(2,2) << std::endl;
        Step += SamplingStep;
        
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
