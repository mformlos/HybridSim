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

    std::string Directory{}, ConfigFile{}, ConfigFileStart{}, MoleculeFileName{}, OutputName{}, AsphericityFile{}, AsphBinFile{};
    double  DeltaTSim{0.01}, Shear{}, delrx{};
    int StartStep{}, EndStep{}, Step{}, SamplingStep{}, NSteps{};
    unsigned Lx{}, Ly{}, Lz{}, NSpecies{}, NumberOfMolecules{}, NumberOfMonomers{};
    double BoxX{}, BoxY{}, BoxZ{}, Entropy_spec{}, Entropy_loc{}, EntropyCOM_spec{}, EntropyCOM_loc{}, LogCells{}, LogSpecies{};
    double normalization_spec{}, p_j{}, S_spec{};
    unsigned CellsX{}, CellsY{}, CellsZ{}, BinX{}, BinY{}, BinZ{};
    bool ParameterInitialized{};
    std::vector<std::vector<std::vector<std::vector<unsigned>>>> LocalSpecies{}, LocalSpeciesCOM{};
    std::vector<double> Asphericities{};
    std::vector<double> AsphBinLimits{};  
    std::vector<double> SpeciesPop{};
    std::vector<double> SpeciesIdentifier{}; 
    if (argc != 9) {
            std::cout << "usage: ./shannon_entropy DIRECTORY STARTSTEP ENDSTEP SAMPLINGSTEP CellsX CellsY CellsZ ASPHBINFILE" << std::endl;
            return EXIT_FAILURE;
    }
    Directory = argv[1];
    StartStep = std::stoi(argv[2]);
    EndStep = std::stoi(argv[3]);
    SamplingStep = std::stoi(argv[4]);
    CellsX = (unsigned)std::stoi(argv[5]);
    CellsY = (unsigned)std::stoi(argv[6]);
    CellsZ = (unsigned)std::stoi(argv[7]);
    AsphBinFile = argv[8]; 

    MoleculeFileName = Directory+"/chain"; 
    AsphericityFile = Directory+"/asphericities";
    ConfigFileStart = Directory+"/configs/config";
    OutputName =Directory+"/data/shannon_entropy_X"+std::to_string(CellsX)+"_Y"+std::to_string(CellsY)+"_Z"+std::to_string(CellsZ);

    /*std::ifstream test(OutputName);
    if (test.good()){
        std::cout << "file " << OutputName << " already exists! Aborting..." << std::endl;
        return EXIT_FAILURE;
    }*/
    std::ofstream Output(OutputName, std::ios::out | std::ios::trunc);
    Output.precision(8);
    
    std::cout << "Directory: " << Directory << std::endl;
    std::cout << "StartStep: " << StartStep << "   EndStep: " << EndStep << "   SamplingStep: " << SamplingStep << std::endl;

    std::ifstream inputfile(Directory+"/parameters.dat", ios::in);
    if (!inputfile.is_open()) {
        std::cout << "could not open file parameters.dat , exiting" << std::endl;
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
    if (!Sys.addMolecules(MoleculeFileName, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl;
        return EXIT_FAILURE;
    }
    NumberOfMolecules = Sys.Molecules.size();
    NumberOfMonomers = Sys.NumberOfParticles(); 
    std::cout << "Number of Molecules: " << NumberOfMolecules <<std::endl;
    std::cout << "Number of Monomers: " << NumberOfMonomers <<std::endl;
    
    //// setup local prob boxes
    BoxX = (double)Lx/CellsX; 
    BoxY = (double)Ly/CellsY; 
    BoxZ = (double)Lz/CellsZ; 
    LocalSpecies=std::vector<std::vector<std::vector<std::vector<unsigned>>>>(CellsX, std::vector<std::vector<std::vector<unsigned>>>(CellsY, std::vector<std::vector<unsigned>>(CellsZ,std::vector<unsigned>()))); 
    LocalSpeciesCOM=std::vector<std::vector<std::vector<std::vector<unsigned>>>>(CellsX, std::vector<std::vector<std::vector<unsigned>>>(CellsY, std::vector<std::vector<unsigned>>(CellsZ,std::vector<unsigned>())));
    LogCells=log(CellsX*CellsY*CellsZ); 

    //// get asphericities of molecules 
    std::ifstream file (AsphericityFile, std::ios::in);
    if (!file.is_open()) {
        std::cout << "file asphericities does not exist! Exiting..." << std::endl;
        return EXIT_FAILURE;
    }
    unsigned index{}; 
    double asph{};
    while(file >> index >> asph) {
        Asphericities.push_back(asph);
    }
    if (Asphericities.size() != NumberOfMolecules) {
	std::cout << "Number of asphericity values (" << Asphericities.size() << ") does not match number of molecules (" << NumberOfMolecules << ")! Exiting... " << std::endl; 
	return EXIT_FAILURE; 
    }

    //// get asphericity bin limits to determine species 
    initializeVector(AsphBinLimits, AsphBinFile);
    NSpecies = AsphBinLimits.size(); 
    LogSpecies=log(NSpecies); 
    std::cout << "Number of asphericity bins: " << NSpecies << std::endl;  

    for (unsigned i = 0; i < NSpecies; i++) {
        SpeciesPop.push_back(0.0);
    }
    //// determine asph bin for each molecule and speciesPop 
    unsigned species {}; 
    for (unsigned i = 0; i < NumberOfMolecules; i++) {
	species = 0; 
	while (AsphBinLimits[species] < Asphericities[i]) species++; 
	SpeciesIdentifier.push_back(species);
	SpeciesPop[species]++;  
    }
  
    Output << "#mixing entropy calculation with box sizes: X= " << BoxX << ", Y= " << BoxY << ", Z= " << BoxZ << std::endl;
    Output << "#Cells in different directions: X= " << CellsX << ", Y= " << CellsY << ", Z= " << CellsZ << std::endl;   
    Output << "#time S_spec(monomers) S_loc(monomer) S_spec(center of mass) S_loc(center of mass)" << std::endl; 
    //// setup done; starting main calculation loop  

    Step = StartStep;

    while (Step <= EndStep) {
	ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
	if (!Sys.initializePositionsPDB(ConfigFile)) {
            std::cout << Step << ", problem with initializing monomers" << std::endl;
	    Step += SamplingStep;
            continue;
        }
	std::cout << "Step: " << Step << std::endl; 
        delrx = (Step+1)*DeltaTSim*Shear*Ly;
        delrx -= Lx*floor(delrx/Lx);
        Sys.delrx = delrx;
        Sys.wrapMoleculesCOM();
        Vector3d COM{Vector3d::Zero()};
        Vector3d imPos{Vector3d::Zero()};
	for (auto& sheet : LocalSpecies) {
	    for (auto& row : sheet) {
		for (auto& list : row) list.clear(); 
	    }
	}
	for (auto& sheet : LocalSpeciesCOM) {
            for (auto& row : sheet) {
                for (auto& list : row) list.clear();
            }
        }

	//// sort species in bins
	for (unsigned mol = 0; mol < NumberOfMolecules; mol++)  {
	    COM=Sys.Molecules[mol].centerOfMassPosition();
	    BinX = (unsigned)COM(0)/BoxX;  
	    BinY = (unsigned)COM(1)/BoxY;  
	    BinZ = (unsigned)COM(2)/BoxZ;
	    //std::cout << "Bin COM: " << BinX << " " << BinY << " " << BinZ << std::endl; 
	    LocalSpeciesCOM[BinX][BinY][BinZ].push_back(SpeciesIdentifier[mol]);   
	    for (auto& mono : Sys.Molecules[mol].Monomers) {
		imPos = image(mono, Sys.BoxSize, delrx); 
		BinX = (unsigned)imPos(0)/BoxX;
                BinY = (unsigned)imPos(1)/BoxY;
                BinZ = (unsigned)imPos(2)/BoxZ;   
                //std::cout << "Bin mono: " << BinX << " " << BinY << " " << BinZ << std::endl;
		LocalSpecies[BinX][BinY][BinZ].push_back(SpeciesIdentifier[mol]);
	    }
	}

	//// calculating Entropies
	Entropy_loc = 0.0;
	Entropy_spec = 0.0;
	EntropyCOM_loc = 0.0;
	EntropyCOM_spec = 0.0;  
	double p_c {};
	for (auto& sheet : LocalSpecies) {
            for (auto& row : sheet) {
                for (auto& list : row) {
		    std::vector<unsigned> SpeciesCount (NSpecies, 0); 
		    for (auto& el : list) SpeciesCount[el]++;
		    normalization_spec = 0.0;
		    p_j = 0; 
		    for (unsigned s = 0; s < NSpecies; s++) {
			normalization_spec += SpeciesCount[s]/SpeciesPop[s];
			p_j += SpeciesCount[s]; 
		    }
		    p_j /= NumberOfMonomers; 
		    S_spec = 0.0; 
	 	    for (unsigned s = 0; s < NSpecies; s++) {
			p_c = SpeciesCount[s]/(SpeciesPop[s]*normalization_spec); 
			if (p_c != 0.0) S_spec += (p_c * log(p_c)); 
		    }
		    Entropy_spec -= p_j*S_spec; 
		    if (p_j != 0.0 ) Entropy_loc -= p_j*log(p_j); 
		}
            }
	}
	for (auto& sheet : LocalSpeciesCOM) {
            for (auto& row : sheet) {
                for (auto& list : row) {
                    std::vector<unsigned> SpeciesCount (NSpecies, 0);
                    for (auto& el : list) SpeciesCount[el]++;
		    normalization_spec = 0.0; 
                    p_j = 0;
                    for (unsigned s = 0; s < NSpecies; s++) {
			 normalization_spec += SpeciesCount[s]/SpeciesPop[s];
			 p_j += SpeciesCount[s]; 
		    }
                    p_j /= NumberOfMolecules;
                    S_spec = 0.0;	
                    for (unsigned s = 0; s < NSpecies; s++) {
                        p_c = SpeciesCount[s]/(SpeciesPop[s]*normalization_spec);
			if (p_c != 0.0) S_spec += (p_c * log(p_c));
                    }
                    EntropyCOM_spec -= p_j*S_spec;
                    if (p_j != 0.0 ) EntropyCOM_loc -= p_j*log(p_j);
                }
            }
        }
	Output << Step*DeltaTSim << " " << Entropy_spec/LogSpecies << " " << Entropy_loc/LogCells << " " << EntropyCOM_spec/LogSpecies << " " << EntropyCOM_loc/LogCells << " " << std::endl; 
	Step += SamplingStep; 
	NSteps++; 
    }
    
    return 0; 

}
