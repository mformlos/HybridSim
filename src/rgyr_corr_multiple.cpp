#include <vector>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>


int main(int argc, char* argv[]) {

    std::string DirectoryStart{}, Directory{}, RgyrFileName{}, OutputFileName{};  
    std::fstream RgyrFile {}, OutputFile{}; 
    double StartTime{}, MaxTime{}, DTime{}, DeltaT{0.01}, Time{}; 
    int StartStep{}, MaxStep{}, DStep{}, LastStep{}, SampleStep; 
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
    
    if (argc != 6) {
            std::cout << "usage: ./rgyr_corr_multiple DIRECTORY STARTSTEP MAXSTEP DSTEP SAMPLESTEP" << std::endl;  
            return EXIT_FAILURE; 
    }
    
    Directory = argv[1]; 
    StartStep = std::stoi(argv[2]); 
    MaxStep = std::stoi(argv[3]); 
    DStep = std::stoi(argv[4]); 
    SampleStep = std::stoi(argv[5]); 
    
    StartTime = StartStep*DeltaT; 
    MaxTime = MaxTime*DeltaT; 
    DTime = DStep*DeltaT;
    
    OutputFileName = "/scratch-new/formanek/HYBRIDSIM/analysisprogs/combined/"+Directory+"/corr_gyr";
    std::ifstream test (OutputFileName); 
    if (test.good()) {
        std::cout << "File " << OutputFileName << " already exists! Aborting ..." << std::endl; 
        return EXIT_FAILURE; 
    }
    OutputFile.open(OutputFileName, std::ios::out | std::ios::trunc); 
    
    
    for (int i = 0; i <= 10; i++) {
        RgyrFileName = Directory+"/rgyr-"+std::to_string(i);
        RgyrFile.open(RgyrFileName, std::ios::in);
        int step; 
        double r; 
        while (RgyrFile >> step >> r) {
            if (step >= StartStep && (step-StartStep)%SampleStep == 0) {
                Rg_mean += r; 
                Rg_sqmean += r*r; 
                N++;
            }
        }
        RgyrFile.close(); 
    }
    
    
    Rg_mean /= N; 
    Rg_sqmean /= N; 
    Rg_quadmean /= N; 
    
    
    std::cout << "Rg mean = " << Rg_mean << std::endl; 
    std::cout << "Rg^2 mean = " << Rg_sqmean << std::endl; 
    std::cout << "N = " << N << std::endl; 

    for (int i = 0; i < MaxStep; i+= DStep) {
        double dt {i*DeltaT}; 
        Corr_rg[dt] = 0.0; 
        Corr_rg_count[dt] = 0;
    } 
     
     
    for (int i = 0; i <= 10; i++) {
        RgyrFileName = Directory+"/rgyr-"+std::to_string(i);
        RgyrFile.open(RgyrFileName, std::ios::in);
        int step; 
        double r; 
        while (RgyrFile >> step >> r) {
            Rg[step] = r; 
        }
        LastStep = step;  
        std::cout << "Laststep = " << LastStep << std::endl; 
    
        int T0Step {StartStep}; 
        int T1Step {};
        int DTStep {};
        
    
        while (T0Step <= LastStep) {  
            T1Step = T0Step; 
            DTStep = 0; 
            //std::cout << "T0Step = " << T0Step << std::endl; 
            while (DTStep <= MaxStep && T1Step <= LastStep) {
                Time = DTStep*DeltaT;
                //std::cout << "T1Step = " << T1Step << std::endl; 
                Corr_rg[Time] += Rg[T0Step]*Rg[T1Step]; //(Gxx[T0Step]-Gxx_mean)*(Gyy[T1Step]-Gyy_mean); 
                Corr_rg2[Time] += (pow(Rg[T0Step],2)- Rg_sqmean)*(pow(Rg[T1Step],2)-Rg_sqmean); 
                Corr_rg_count[Time]++; 
                T1Step += DStep;
                DTStep += DStep;  
            }
            T0Step += SampleStep; 
        }
        RgyrFile.close(); 
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
        OutputFile << corr.first << " " << corr.second << " " << Corr_rg2[corr.first] << " " << Corr_rg_count[corr.first] << std::endl; 
    }
    
    OutputFile.close(); 
    
    
    return 0;
    
        
} 
    
