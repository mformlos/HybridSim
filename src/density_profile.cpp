#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "Molecule.h"
#include "HelperFunctions.h"

int main(int argc, char* argv[]) {
    std::string DirectoryFile{}, ConfigFile{}, ConfigFileStart{}, DensityFileName{}, DistanceFileName{}, OutputDir;
    double  BinSize{}, Lx {}, Ly{}, Lz{};
    unsigned StartStep{}, Step{}, SamplingStep{}, N{}, Monomers{}, BinsX{}, BinsY{}, BinsZ{}, BX{}, BY{}, BZ{}, Bin{}, TotalBins{}; 
    
    if (argc != 10) {
            std::cout << "usage: ./density_profile DIRECTORYFILE STARTSTEP SAMPLINGSTEP MONOMERS BINSIZE LX LY LZ OUTPUTDIR" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    DirectoryFile = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    SamplingStep = std::stoi(argv[3]); 
    Monomers = std::stoi(argv[4]); 
    BinSize = std::stod(argv[5]); 
    Lx = std::stod(argv[6]);
    Ly = std::stod(argv[7]);
    Lz = std::stod(argv[8]);
    OutputDir = argv[9]; 
    
    BinsX = (unsigned)(Lx/BinSize)+1;
    BinsY = (unsigned)(Ly/BinSize)+1;
    BinsZ = (unsigned)(Lz/BinSize)+1;
    
    TotalBins = BinsX*BinsY*BinsZ; 
    
    std::vector<std::string> Directories{}; 
    initializeStringVector(Directories, DirectoryFile); 
    
    std::cout << "DirectoryFile: " << DirectoryFile << std::endl; 
    std::cout << "StartStep: " << StartStep << "   SamplingStep: " << SamplingStep << std::endl; 
    std::cout << "Monomers: " << Monomers << " , Bin size: " << BinSize << std::endl; 
    std::cout << "Lx: " << Lx << " , Ly: " << Ly << " , Lz: " << Lz << std::endl; 
    std::cout << "BinsX: " << BinsX << " , BinsY: " << BinsY <<  " , BinsZ: " << BinsZ << std::endl; 
    
    Molecule Mol(Monomers);
    Vector3d COMPos{};  
    
    std::vector<double> DensityMap(TotalBins, 0.0); 
    std::vector<double> Distance(Monomers, 0.0); 
    
    
    
    DensityFileName = OutputDir+"/density"; 
    DistanceFileName = OutputDir+"/monomer_distance"; 
    
    std::cout << DensityFileName << std::endl; 
    std::cout << DistanceFileName << std::endl; 
    
    N = 0; 
    for (auto& Directory : Directories) {
        std::cout << Directory << std::endl; 
        ConfigFileStart = Directory+"/configs/config";    
        Step = StartStep;  
        while (true) {
            ConfigFile = ConfigFileStart+std::to_string(Step)+".pdb";
            if (!Mol.initializePositions(ConfigFile)) {
                std::cout << "reached last step at " << Step << std::endl; 
                break; 
            } 
            COMPos = Mol.centerOfMassPosition();
            Mol.translate(-COMPos); 
            for (auto& mono : Mol.Monomers) {
                BX = (unsigned)((mono.Position(0)+Lx*0.5)/BinSize); 
                BY = (unsigned)((mono.Position(1)+Ly*0.5)/BinSize); 
                BZ = (unsigned)((mono.Position(2)+Lz*0.5)/BinSize); 
                Bin = BX + BY*(BinsX)+BZ*(BinsX*BinsY); 
                DensityMap[Bin]++; 
                Distance[mono.Identifier] += mono.Position.norm(); 
            } 
            Step += SamplingStep; 
            N++;     
        }
    }
    
    for (auto& value : DensityMap) value /= (Monomers*N*pow(BinSize,3)); 
    
    for (auto& value : Distance) value /= N; 
    
    std::vector<double> Density_xy(BinsX*BinsY, 0.0); 
    std::vector<double> Density_xz(BinsX*BinsZ, 0.0); 
    unsigned Bin_xy{}; 
    unsigned Bin_xz{}; 
    for (unsigned x = 0; x < BinsX; x++) {
        for (unsigned y = 0; y < BinsY; y++) {
            Bin_xy = x + y*BinsX;
            for (unsigned z = 0; z < BinsZ; z++) {
                Bin = x + y*BinsX + z*BinsX*BinsY;
                Bin_xz = x + z*BinsX;
                Density_xy[Bin_xy] += DensityMap[Bin]; 
                Density_xz[Bin_xz] += DensityMap[Bin];
            }
        }
    }
    
    double normalization {0.0};
    for (auto& value : Density_xy) {
        //value /= BinsZ;
        normalization += value; 
    }
    std::cout << "normalization xy = " << normalization << std::endl; 
    normalization = 0.0; 
    for (auto& value : Density_xz) {
        //value /= BinsY;
        normalization += value; 
    }
    
    std::cout << "normalization xz = " << normalization << std::endl; 
    
    std::ofstream DensityFile_xy(DensityFileName+"_xy", std::ios::out | std::ios::trunc);
    DensityFile_xy.precision(5);
    
    for (unsigned y = 0; y < BinsY; y++) {
        for (unsigned x = 0; x < BinsX; x++) {
            DensityFile_xy.width(14);
            DensityFile_xy << Density_xy[x+y*BinsX];     
        }
        DensityFile_xy << std::endl; 
    }
    
    DensityFile_xy.close(); 
    
    std::ofstream DensityFile_xz(DensityFileName+"_xz", std::ios::out | std::ios::trunc);
    DensityFile_xz.precision(5);
    
    for (unsigned z = 0; z < BinsZ; z++) {
        for (unsigned x = 0; x < BinsX; x++) {
            DensityFile_xz.width(14);
            DensityFile_xz << Density_xz[x+z*BinsX];     
        }
        DensityFile_xz << std::endl; 
    }
    
    DensityFile_xz.close(); 
    
    
    
    std::ofstream DensityFile(DensityFileName, std::ios::out | std::ios::trunc); 
    DensityFile.precision(5); 
    
    for (unsigned z = 0; z < BinsZ; z++) {
        for (unsigned y = 0; y < BinsY; y++) {
            for (unsigned x = 0; x < BinsX; x++) {
                Bin = x + y*BinsX + z*BinsX*BinsY; 
                DensityFile.width(14); 
                DensityFile << x*BinSize-Lx*0.5; 
                DensityFile.width(14); 
                DensityFile << y*BinSize-Ly*0.5;
                DensityFile.width(14); 
                DensityFile << z*BinSize-Lz*0.5; 
                DensityFile.width(14);
                DensityFile << DensityMap[Bin] << std::endl; 
            }
        }
    }
    
    DensityFile.close(); 
    
    std::ofstream DistanceFile(DistanceFileName, std::ios::out | std::ios::trunc); 
    DistanceFile.precision(5); 
    for (unsigned i = 0; i < Monomers; i++) {
        DistanceFile.width(14);
        DistanceFile << i; 
        DistanceFile.width(14); 
        DistanceFile << Distance[i] << std::endl;  
    }
    
    DistanceFile.close();
    
    return 0; 
}
