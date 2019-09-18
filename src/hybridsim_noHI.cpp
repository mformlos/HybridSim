#include <sys/time.h>
#include <csignal>
#include <list>
#include "MPC.h"
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
    unsigned Lx{}, Ly{}, Lz{}, MPCRho{}, Seed{}; 
    unsigned long long COMWrapInterval{1000}, MPCInterval {}, EquilSteps{}, TotalSteps {}, n {}, m {}; 
    double MDStep{}, MPCStep{}, Shear{}, TotalTime{}, EquilTime{}, SimTime{}, Temperature{}, Time{};
    bool ParameterInitialized{}; 
    std::string OutputStepFile{}, MoleculeFile{}, LinkFile{}, ConfigFile{}, VelocFile{}, StatisticsFile{}, ConfigOutFile{}, NeighbourCellFile{}; 
    std::vector<unsigned long long> OutputSteps; 
    std::vector<unsigned long long>::iterator OutputStepsIt{};
    
    if (argc != 2) {
        std::cout << "usage: ./hybridsim_noHI PARAMETER-INPUT-FILE " << std::endl;  
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
    Shear = extractParameter<double>("Shear", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    Temperature = extractParameter<double>("Temperature", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    MPCRho = extractParameter<unsigned>("MPCRho", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    MDStep = extractParameter<double>("MDStep", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    MPCStep = extractParameter<double>("MPCStep", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    Time = extractParameter<double>("StartTime", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    EquilTime = extractParameter<double>("EquilTime", inputfile, ParameterInitialized); 
    if (!ParameterInitialized) return EXIT_FAILURE;
    SimTime = extractParameter<double>("SimTime", inputfile, ParameterInitialized);
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
    VelocFile = extractParameter<std::string>("VelocFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    StatisticsFile = extractParameter<std::string>("StatisticsFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    ConfigOutFile = extractParameter<std::string>("ConfigOutFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    NeighbourCellFile = extractParameter<std::string>("NeighbourCellFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    
    inputfile.close(); 
    TotalTime = EquilTime + SimTime;
    
    MPCInterval = (unsigned long long) (MPCStep/MDStep); 
    if (MPCInterval*MDStep != MPCStep) {
        std::cout << "Steps for MD and MPC MUST be multiples of each other!"; 
        return EXIT_FAILURE; 
    }
    
    EquilSteps = (unsigned long long) (EquilTime/MDStep); 
    TotalSteps = (unsigned long long) (TotalTime/MDStep); 
    std::cout << "Total number of MD steps: " << TotalSteps << std::endl;  
    
    ///////////////////////////////////////// 
     
    /////// PARAMETER CONTROL OUTPUT //////// 
    std::cout << "Lx is " << Lx << std::endl;
    std::cout << "Ly is " << Ly << std::endl;
    std::cout << "Lz is " << Lz << std::endl;
    std::cout << "Temperature is " << Temperature << std::endl;
    std::cout << "Shear is " << Shear << std::endl;
    std::cout << "MPC particles per cell is " << MPCRho << std::endl;
    std::cout << "MDStep is " << MDStep << std::endl;
    std::cout << "MPCStep is " << MPCStep << std::endl;
    std::cout << "StartTime is " << Time << std::endl;  
    std::cout << "EquilTime is " << EquilTime << std::endl;    
    std::cout << "TotalTime is " << TotalTime << std::endl;
    std::cout << "RNG seed is " << Seed << std::endl; 
    std::cout << "OutputStepFile is " << OutputStepFile << std::endl;
    std::cout << "MoleculeFile is " << MoleculeFile << std::endl;
    std::cout << "LinkFile is " << LinkFile << std::endl;
    std::cout << "ConfigFile is " << ConfigFile << std::endl;
    std::cout << "VelocFile is " << VelocFile << std::endl;
    std::cout << "NeighbourCellFile is " << NeighbourCellFile << std::endl;
    
    /////////////////////////////////////
    
    ////// RANDOM ENGINE SEEDING & WARMUP //////

    Rand::seed(Seed); 
    Rand::warmup(10000); 
    
    /////////////////////////////////////
    
    /////// SYSTEM INITIALIZATION ///////
    System Sys(1.2, 1.5, Lx, Ly, Lz, Shear, 0.0, false, true); 
    
    std::cout << "Neighbourcells in each direction: " << Sys.Cells[0] << " " <<  Sys.Cells[1] << " " << Sys.Cells[2] << std::endl; 
    std::cout << "Cell side lengths: " << Sys.CellSideLength[0] << " " <<  Sys.CellSideLength[1] << " " << Sys.CellSideLength[2] << std::endl;
    
    if (!Sys.setNeighbourDirections(NeighbourCellFile)){
        return EXIT_FAILURE;
    }
    
    if (!Sys.addMolecules(MoleculeFile, (double)MPCRho)) {
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
    
    if (VelocFile == "RANDOM") {
        Sys.initializeVelocitiesRandom(Temperature); 
        std::cout << "velocities set randomly" << std::endl; 
    }
    else {
        if (!Sys.initializeVelocities(VelocFile)) {
            std::cout << "VelocFile does not exist or contains too little lines!" << std::endl; 
            return EXIT_FAILURE; 
        } 
    }
    
    if (Sys.NumberOfMolecules() == 1)  Sys.centerMolecule(0); //Sys.setMoleculeCOM(0,newCOM); 
    else if (Sys.NumberOfMolecules() > 1){
        Sys.wrapMoleculesCOM(); 
        std::cout << "Molecules' centers of mass set inside box." << std::endl; 
    }
    
    try {
        Sys.updateCellLists();
        Sys.calculateForcesCellList();
        //Sys.updateVerletLists(); 
        //Sys.calculateForces();
        //Sys.calculateForcesBrute(); 
    }
    catch (const LibraryException &ex) {
        std::cout << ex.what() << std::endl; 
        std::cout << "bad initial configuration! Terminating program." << std::endl;   
        return EXIT_FAILURE; 
    }
    
    //////////////////////////////////
    
    /////// MPC INITIALIZATION ///////
    MPC Mpc(Lx, Ly, Lz, MPCRho, Temperature, Shear, false);
    
    Mpc.initializeSoluteVector(Sys.NumberOfParticles());
    
    /////////////////////////////////////
    
    /////// OUTPUT INITIALIZATION ///////
    if (!initializeStepVector(OutputSteps, OutputStepFile)) {
        std::cout << "OutputStepFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    }
    std::cout << "Output will be done " << OutputSteps.size() << " times. " << std::endl; 
    OutputStepsIt = OutputSteps.begin(); 
    
    std::ofstream StatisticsStream(StatisticsFile, std::ios::out | std::ios::app); 
    FILE* PDBout{}; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    ///////////////////////////////////
    ////// MAIN SIMULATION LOOP ///////
    m = -EquilSteps + (unsigned long long)(Time/MDStep); 
    for (n = 0; n <= TotalSteps; n++) {
        ////// MD STEP ///////
        Time += MDStep; 
        //std::cout << n << std::endl;
        Mpc.updateBoxShift(MDStep); 
        Sys.delrx = Mpc.delrx;     
        if (n%COMWrapInterval==0) Sys.wrapMoleculesCOM(); 
        try {
            if (m == *OutputStepsIt) Sys.propagate(MDStep, true);
            else Sys.propagate(MDStep);
        }
        catch (const LibraryException &ex) {
            std::cout << ex.what() << std::endl; 
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "r");
            if (PDBout){
                std::cout << "File " << ConfigOutFile+std::to_string(m)+".pdb" << " already exists!" << std::endl; 
            }  
            else {
                PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
                Sys.printPDB(PDBout, n, 1);
                fclose(PDBout);
            }
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }  
        //////// MPC STEP /////////
        if (n%MPCInterval == 0) {
            ////// MPC SOLUTE TRANSFER & GRID SHIFT ////
            //Mpc.shiftGrid(); 
            //std::cout << "Position before: " <<  Sys.Molecules[0].Monomers[6].Position.transpose() << std::endl;
            //std::cout << "Velocity before: " <<  Sys.Molecules[0].Monomers[6].Velocity.transpose() << std::endl;
            //std::cout << "getting solute" << std::endl; 
            Mpc.getSolute(Sys.Molecules); 
            //std::list<unsigned> SoluteCells{}; 
            /*for (auto& sol : Mpc.Solute) {
                Mpc.updateSoluteCellIndexWithList(sol, SoluteCells); 
            }*/
            //SoluteCells.sort(); 
            //SoluteCells.unique();
            /*unsigned c {0};
            for (auto& Index : SoluteCells) {
                c++; 
                Mpc.drawRotation(Index); 
                Mpc.calculateCOMVelNoHI(Index); 
            }*/
            //std::cout << c++ << endl; 
            //std::cout << "rotating solute" << std::endl; 
            for (unsigned sol = 0; sol < Mpc.Solute.size(); sol++) {
                Mpc.rotateSoluteNoHI(sol); 
            }
            //std::cout << "returning solute" << std::endl;
            Mpc.returnSolute(Sys.Molecules);
            //std::cout << "returned solute" << std::endl;
            //std::cout << "Position after: " <<  Sys.Molecules[0].Monomers[6].Position.transpose() << std::endl;
            //std::cout << "Velocity after: " <<  Sys.Molecules[0].Monomers[6].Velocity.transpose() << std::endl;
        }
        ////////// MPC STEP FINISHED //////////
            
        ////////// OUTPUT ////////////
        if (SignalCaught) {
            std::cout << "writing job data..." << std::endl; 
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "r");
            if (PDBout){
                std::cout << "File " << ConfigOutFile+std::to_string(m)+".pdb" << " already exists!" << std::endl; 
            }  
            else {
                PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
                Sys.printPDB(PDBout, n, 1);
                fclose(PDBout);
            }
            std::cout << "Terminating..." << std::endl;   
            std::terminate(); 
        }
        if (m == *OutputStepsIt) {
            Sys.printStatistics(StatisticsStream, Time);
            PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "r");
            if (PDBout){
                std::cout << "File " << ConfigOutFile+std::to_string(m)+".pdb" << " already exists!" << std::endl; 
                SignalCaught = 666;
            }  
            else {
                PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
                Sys.printPDB(PDBout, n, 1);
                fclose(PDBout);
            } 
            std::cout << Time << " " << Sys.delrx <<  std::endl; 
            OutputStepsIt++;     
        }
        m++;
    }
    
    gettimeofday(&end, NULL); 
    
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per step: " << realTime/TotalSteps << std::endl;
    
    
    return EXIT_SUCCESS; 
}
