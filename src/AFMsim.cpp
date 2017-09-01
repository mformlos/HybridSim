#include <sys/time.h>
#include <omp.h>
#include <csignal>
#include "MPC.h"
#include "Velocity_Profile.h"
#include "Velocity_Hist.h"
#include "HelperFunctions.h"
#include "System.h"

int SignalCaught {}; 
void signalHandler(int signum) 
{
    SignalCaught = signum; 
    std::cout << "signal " << signum << " caught!" << std::endl; 
}

int main(int argc, char* argv[]) {
    SignalCaught = 0;
    signal(SIGINT, signalHandler); 
    unsigned Lx{}, Ly{}, Lz{}, TotalSteps {}, Seed {}, n {}, m {};  
    double MDStep{}, TotalTime{}, SimTime{}, AnchorTime{}, Temperature{}, Time{}, SurfaceEStart{}, SurfaceEEquil {}, Gamma {};
    bool ParameterInitialized{}, AdsorptionOn{false}, Anchored {false}, AFMset {false}; 
    std::string OutputStepFile{}, MoleculeFile{}, LinkFile{}, ConfigFile{}, StatisticsFile{}, ConfigOutFile{}, ForceUpdateFile {}, ExtensionFile{}; 
    std::vector<unsigned> OutputSteps; 
    std::vector<unsigned>::iterator OutputStepsIt{};
    std::vector<ForceUpdate> ForceUpdates; 
    std::vector<ForceUpdate>::iterator ForceUpdatesIt{}; 
    
    
    if (argc != 2) {
        std::cout << "usage: ./AFMsim.cpp PARAMETER-INPUT-FILE " << std::endl;  
        return EXIT_FAILURE; 
    }
    
    //////////////////////////////
    
    ////// PARAMETER READS ///////
    
    std::ifstream inputfile(argv[1], ios::in);
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
    Temperature = extractParameter<double>("Temperature", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    Gamma = extractParameter<double>("Gamma", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    MDStep = extractParameter<double>("MDStep", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    Time = extractParameter<double>("StartTime", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    SimTime = extractParameter<double>("SimTime", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    SurfaceEStart = extractParameter<double>("SurfaceEnergyStart", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    SurfaceEEquil = extractParameter<double>("SurfaceEnergyEquil", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    AnchorTime = extractParameter<double>("AnchorTime", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    Seed = extractParameter<double>("Seed", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    OutputStepFile = extractParameter<std::string>("OutputStepFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    MoleculeFile = extractParameter<std::string>("MoleculeFile", inputfile, ParameterInitialized);      
    if (!ParameterInitialized) return EXIT_FAILURE;
    LinkFile = extractParameter<std::string>("LinkFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    ConfigFile = extractParameter<std::string>("ConfigFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    StatisticsFile = extractParameter<std::string>("StatisticsFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ConfigOutFile = extractParameter<std::string>("ConfigOutFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ForceUpdateFile = extractParameter<std::string>("ForceUpdateFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ExtensionFile = extractParameter<std::string>("ExtensionFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    
    inputfile.close(); 
    
    if (SurfaceEStart >= 0.00001) AdsorptionOn = true; 
     
    TotalSteps = (unsigned) (SimTime/MDStep); 
    std::cout << "Total number of MD steps: " << TotalSteps << std::endl;  
    
    ///////////////////////////////////////// 
     
    /////// PARAMETER CONTROL OUTPUT //////// 
    std::cout << "Lx is " << Lx << std::endl;
    std::cout << "Ly is " << Ly << std::endl;
    std::cout << "Lz is " << Lz << std::endl;
    std::cout << "Temperature is " << Temperature << std::endl;
    std::cout << "Friction Coefficient is " << Gamma << std::endl;
    std::cout << "MDStep is " << MDStep << std::endl;
    std::cout << "Surface energy constant for anchoring is " << SurfaceEStart << std::endl;
    std::cout << "Surface energy constant for simulation is " << SurfaceEEquil << std::endl;
    std::cout << "RNG seed is " << Seed << std::endl;  
    std::cout << "TotalTime is " << TotalTime << std::endl;
    std::cout << "OutputStepFile is " << OutputStepFile << std::endl;
    std::cout << "MoleculeFile is " << MoleculeFile << std::endl;
    std::cout << "LinkFile is " << LinkFile << std::endl;
    std::cout << "ConfigFile is " << ConfigFile << std::endl;
    std::cout << "ForceUpdateFile is " << ForceUpdateFile << std::endl;
    std::cout << "ExtensionFile is " << ExtensionFile << std::endl;
    /////////////////////////////////////
    
    /////// SYSTEM INITIALIZATION ///////
    System Sys(Lx, Ly, Lz, 0.0, SurfaceEStart, AdsorptionOn); 
    
    if (!Sys.addMolecules(MoleculeFile, 1.0)) {
        std::cout << "MoleculeFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    if (!Sys.addLinks(LinkFile)) {
        std::cout << "LinkFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    } 
    
    if (!Sys.initializePositions(ConfigFile)) {
        std::cout << "ConfigFile does not exist or contains too little lines!" << std::endl; 
        return EXIT_FAILURE; 
    } 
 
    Sys.initializeVelocitiesRandom(Temperature); 
    
    if (Sys.NumberOfMolecules() == 1)  {
        Sys.centerMolecule(0); //Sys.setMoleculeCOM(0,newCOM); 
        double mindist{Lz}; 
        Vector3d NewCOM {Sys.Molecules[0].centerOfMassPosition()};
        for (auto& mono : Sys.Molecules[0].Monomers) {
            if (mono.Position(2) < mindist) mindist = mono.Position(2);
        }  
        NewCOM(2) = NewCOM(2) - mindist + 1.3;
        Sys.setMoleculeCOM(0,NewCOM); 
    }
    else if (Sys.NumberOfMolecules() > 1){
        std::cout << "Code is supposed to be used with only 1 Molecule" << std::endl; 
        return EXIT_FAILURE;  
    }
    
    try {
        //Sys.updateVerletLists(); 
        //Sys.calculateForces();
        Sys.calculateForcesBrute(); 
    }
    catch (const LibraryException &ex) {
        std::cout << ex.what() << std::endl; 
        std::cout << "bad initial configuration! Terminating program." << std::endl;   
        return EXIT_FAILURE; 
    }
    
    
    /////////////////////////////////////
    
    /////// OUTPUT INITIALIZATION ///////
    
    if (!initializeStepVector(OutputSteps, OutputStepFile)) {
        std::cout << "OutputStepFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    }
    std::cout << "Output will be done " << OutputSteps.size() << " times. " << std::endl; 
    OutputStepsIt = OutputSteps.begin(); 
    
    if (!initializeForceUpdateVector(ForceUpdates, ForceUpdateFile)) {
        std::cout << "ForceUpdateFile does not exist!" << std::endl; 
        return EXIT_FAILURE;     
    }
    
    std::ofstream StatisticsStream(StatisticsFile, ios::out | ios::trunc); 
    std::ofstream ExtensionStream(ExtensionFile, ios::out | ios::trunc); 
    FILE* PDBout{}; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    ////// RANDOM ENGINE SEEDING & WARMUP //////
    Rand::seed(Seed); 
    Rand::warmup(10000); 
    
    ///////////////////////////////////
    ////// MAIN SIMULATION LOOP ///////
    m = -TotalSteps; 
    for (n = 0; n <= TotalSteps; n++, m++) {
        Time += MDStep; 
        try {
            if (m == *OutputStepsIt) Sys.propagateLangevin(MDStep, Temperature, Gamma, true);
            else Sys.propagateLangevin(MDStep, Temperature, Gamma, false);
        }
        catch (const LibraryException &ex) {
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1);
            fclose(PDBout);    
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        
        if (!Anchored && Time > AnchorTime) {
            if (Sys.Molecules[0].Monomers[0].Position(2) <= 1.05 && Sys.Molecules[0].Monomers[0].Position(2) >= 0.9 ) {
                Sys.setAnchorAuto(0); 
                Sys.SurfaceEnergy = SurfaceEEquil; 
                Anchored = true; 
                Sys.setDrive(ForceUpdatesIt -> Force, Sys.Molecules[0].NumberOfMonomers-1); 
                AFMset = true; 
                Sys.calculateForcesBrute(true);
                std::cout << "Chain anchored & starting force application " << std::endl;
                ForceUpdatesIt++; 
                m = 0;  
                Time = 0.0; 
            }
        }
        
        if (AFMset && m == ForceUpdatesIt -> Step) {
            Sys.changeDrive(ForceUpdatesIt -> Force, Sys.Molecules[0].NumberOfMonomers-1);
            Sys.calculateForcesBrute(true);
            ForceUpdatesIt++; 
        } 
        
        
        ////////// OUTPUT ////////////
        if (SignalCaught) {
            std::cout << "writing job data..." << std::endl; 
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1);
            fclose(PDBout);
            fclose(PDBout);  
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        
        if (m == *OutputStepsIt) {
            Sys.printStatistics(StatisticsStream, Time); 
            Sys.printForceExtension(ExtensionStream, Time, 2); 
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
            Sys.printPDB(PDBout, n, 1); 
            fclose(PDBout);
            std::cout << Time << " " << ForceUpdatesIt -> Force <<  std::endl; 
            OutputStepsIt++;        
        }
    }    
    
    //////////////////////////////////////
    ////// MAIN SIMULATION LOOP END///////
    
    gettimeofday(&end, NULL); 
    
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/Sys.Molecules[0].NumberOfMonomers/TotalSteps << std::endl;
    
    return EXIT_SUCCESS;

}
    

