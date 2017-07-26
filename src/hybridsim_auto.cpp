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
    unsigned Lx{}, Ly{}, Lz{}, MPCRho{}, COMWrapInterval{1000}, MPCInterval {}, EquilSteps{}, SimSteps{}, TotalSteps {}, n {}, m {}; 
    int tid{}, procs{}, maxt{}, inpar{}, dynamic{}, nested{}, nthreads{}; 
    double MDStep{}, MPCStep{}, Shear{}, TotalTime{}, EquilTime{}, SimTime{}, Temperature{}, Time{};
    bool ParameterInitialized{}; 
    std::string OutputStepFile{}, MoleculeFile{}, LinkFile{}, ConfigFile{}, StatisticsFile{}, ConfigOutFile{}, VelProfFile{}, FluidFile{}, FluidInput{}; 
    std::vector<unsigned> OutputSteps; 
    std::vector<unsigned>::iterator OutputStepsIt{};
    
    
    if (argc != 2) {
        std::cout << "usage: ./hybridsim_auto PARAMETER-INPUT-FILE " << std::endl;  
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
    FluidFile = extractParameter<std::string>("FluidFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    FluidInput = extractParameter<std::string>("FluidInput", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    VelProfFile = extractParameter<std::string>("VelProfFile", inputfile, ParameterInitialized);
    
    
    inputfile.close(); 
    TotalTime = EquilTime + SimTime;
    
    MPCInterval = (unsigned) (MPCStep/MDStep); 
    if (MPCInterval*MDStep != MPCStep) {
        std::cout << "Steps for MD and MPC MUST be multiples of each other!"; 
        return EXIT_FAILURE; 
    }
    
    EquilSteps = (unsigned) (EquilTime/MDStep); 
    SimSteps = (unsigned) (SimTime/MDStep); 
    TotalSteps = (unsigned) (TotalTime/MDStep); 
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
    std::cout << "EquilTime is " << EquilTime << std::endl;    
    std::cout << "TotalTime is " << TotalTime << std::endl;
    std::cout << "OutputStepFile is " << OutputStepFile << std::endl;
    std::cout << "MoleculeFile is " << MoleculeFile << std::endl;
    std::cout << "LinkFile is " << LinkFile << std::endl;
    std::cout << "ConfigFile is " << ConfigFile << std::endl;
    std::cout << "FluidFile is " << FluidFile << std::endl;
    std::cout << "FluidInput is " << FluidInput << std::endl;
    std::cout << "VelProfFile is " << VelProfFile << std::endl;
    
    /////////////////////////////////////
    
    /////// SYSTEM INITIALIZATION ///////
    System Sys(Lx, Ly, Lz, Shear); 
    
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
 
    Sys.initializeVelocitiesRandom(Temperature); 
    
    //Vector3d newCOM(10., 18.,10.);
    if (Sys.NumberOfMolecules() == 1)  Sys.centerMolecule(0); //Sys.setMoleculeCOM(0,newCOM); 
    else if (Sys.NumberOfMolecules() > 1){
        std::cout << "Code can not handle " << Sys.NumberOfMolecules() << " yet. Please use only 1 Molecule" << std::endl; 
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
    
    
    //////////////////////////////////
    
    /////// MPC INITIALIZATION ///////
    MPC Mpc(Lx, Ly, Lz, MPCRho, Temperature, Shear); 
    if (FluidInput == "RANDOM") {
        std::cout << "random initialization of fluid" << std::endl; 
        Mpc.initializeRandom(); 
    }
    else if (FluidInput == "PROFILE") {
        std::cout << "profile initialization of fluid" << std::endl; 
        Mpc.initializeProfile(); 
    }
    else {
        if (!Mpc.initializeFile(FluidInput)) {
            std::cout << FluidInput << " does not exist or fluid was not initialized correctly" << std::endl; 
            return EXIT_FAILURE; 
        }
        else {
            std::cout << "fluid successfully initialized from file "<< FluidInput << std::endl; 
        }
    }    
    Mpc.initializeSoluteVector(Sys.NumberOfParticles()); 
    
    
    /////////////////////////////////////
    
    /////// OUTPUT INITIALIZATION ///////
    VelocityProfile VelProf{0.2};  
    //VelocityHist VelHist{0.1}; 
    
    if (!initializeStepVector(OutputSteps, OutputStepFile)) {
        std::cout << "OutputStepFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    }
    std::cout << "Output will be done " << OutputSteps.size() << " times. " << std::endl; 
    OutputStepsIt = OutputSteps.begin(); 
    
    std::ofstream StatisticsStream(StatisticsFile, ios::out | ios::trunc); 
    FILE* PDBout{}; 
    //PDBout = fopen(ConfigOutFile.c_str(), "w"); 
    
    FILE* FluidFilePointer {}; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
  
    
    ////////////////////////////////////////
    
    /////// OPENMP BEGINNING & TESTS ///////
    
    #pragma omp parallel private(tid, n, m)
    {
        /////// OMP PARAMETERS ///////
        tid = omp_get_thread_num(); 
        if (tid == 0) {
            printf("Thread %d getting environment info...\n", tid);

            /* Get environment information */
            procs = omp_get_num_procs();
            nthreads = omp_get_num_threads();
            maxt = omp_get_max_threads();
            inpar = omp_in_parallel();
            dynamic = omp_get_dynamic();
            nested = omp_get_nested();
    
            /* Print environment information */
            printf("Number of processors = %d\n", procs);
            printf("Number of threads = %d\n", nthreads);
            printf("Max threads = %d\n", maxt);
            printf("In parallel? = %d\n", inpar);
            printf("Dynamic threads enabled? = %d\n", dynamic);
            printf("Nested parallelism enabled? = %d\n", nested);
        }
        
        ////// RANDOM ENGINE SEEDING & WARMUP //////
        Rand::seed(tid); 
        Rand::warmup(10000); 
        
        ///////////////////////////////////
        ////// MAIN SIMULATION LOOP ///////
        m = -EquilSteps; 
        for (n = 0; n <= TotalSteps; n++) {
            #pragma omp single 
            {
                //std::cout << n << std::endl; 
                Time += MDStep; 
                Mpc.updateBoxShift(MDStep); 
                Sys.delrx = Mpc.delrx; 
                if (n%COMWrapInterval==0) Sys.wrapMoleculesCOM(); 
                try {
                    if (m == *OutputStepsIt) Sys.propagate(MDStep, true);
                    else Sys.propagate(MDStep);
                }
                catch (const LibraryException &ex) {
                    //std::cout << ex.what() << std::endl; 
                    PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
                    Sys.printPDB(PDBout, n, 1);
                    FluidFilePointer = fopen(FluidFile.c_str(), "w"); 
                    Mpc.printFluid(FluidFilePointer, Time); 
                    fclose(PDBout);  
                    fclose(FluidFilePointer); 
                    if (!VelProfFile.empty()) {
                        std::ofstream VelProfFileStream{}; 
                        VelProfFileStream.open(VelProfFile, ios::out | ios::trunc); 
                        VelProf.print_result(VelProfFileStream); 
                    }     
                    std::cout << "Terminating..." << std::endl;   
                    std::terminate(); 
                }  
               
                
            }
            
            
            //////// MPC STEP /////////
            if (n%MPCInterval == 0) {
                ////// MPC SOLUTE TRANSFER & GRID SHIFT ////
                #pragma omp single 
                {
                    Mpc.shiftGrid(); 
                    Mpc.getSolute(Sys.Molecules); 
                    for (auto& sol : Mpc.Solute) {
                        Mpc.updateSoluteCellIndex(sol); 
                    } 
                }
                ////// STREAMING STEP ///////
                #pragma omp for schedule(static)
                for (unsigned part = 0; part < Mpc.NumberOfParticles; part++) {
                    Mpc.streamPlusCellAssignment(Mpc.Fluid[part], MPCStep); 
                }
                ////// SORT MPC PARTICLE VECTOR TO IMPROVE DATA LOCALITY ///// 
                #pragma omp single 
                {
                    if (n%(10*MPCInterval)==MPCInterval) Mpc.sortVector(); 
                }
                ////// SORT INTO CELLS & DRAW ROTATION AXES ////
                #pragma omp sections 
                {
                    #pragma omp section
                    Mpc.sortOnly(); 
                    
                    #pragma omp section 
                    for (unsigned Index = 0; Index < Mpc.NumberOfCells; Index++) {
                        Mpc.drawRotation(Index); 
                    }
                }
                
                ////// ROTATE FLUID & SOLUTE //////
                #pragma omp for schedule(static) 
                for (unsigned Index = 0; Index < Mpc.NumberOfCells; Index++) {
                    Mpc.calculateCOMVel(Index); 
                }
                #pragma omp for schedule(static)
                for (unsigned part = 0; part < Mpc.NumberOfParticles; part++) {
                    Mpc.rotate(part); 
                }
            
                #pragma omp for schedule(static)
                for (unsigned sol = 0; sol < Mpc.Solute.size(); sol++) {
                    Mpc.rotateSolute(sol); 
                }
                
                ////// APPLY CHANGES TO SOLUTE /////
                #pragma omp single
                {
                    Mpc.returnSolute(Sys.Molecules);
                } 
                #pragma omp barrier
            } 
            
            ////////// MPC STEP FINISHED //////////
            
            ////////// OUTPUT ////////////
            #pragma omp single 
            {
                if (SignalCaught) {
                    std::cout << "writing job data..." << std::endl; 
                    PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
                    Sys.printPDB(PDBout, n, 1);
                    fclose(PDBout);
                    FluidFilePointer = fopen(FluidFile.c_str(), "w"); 
                    Mpc.printFluid(FluidFilePointer, Time); 
                    fclose(PDBout);  
                    fclose(FluidFilePointer);  
                    if (!VelProfFile.empty()) {
                        std::ofstream VelProfFileStream{}; 
                        VelProfFileStream.open(VelProfFile, ios::out | ios::trunc); 
                        VelProf.print_result(VelProfFileStream); 
                    }  
                    std::cout << "Terminating..." << std::endl;   
                    std::terminate(); 
                }
                if (m == *OutputStepsIt) {
                    Sys.printStatistics(StatisticsStream, Time); 
                    PDBout = fopen((ConfigOutFile+std::to_string(m)+".pdb").c_str(), "w");
                    Sys.printPDB(PDBout, n, 1); 
                    fclose(PDBout);
                    if (m > 10000 && !VelProfFile.empty()) Mpc(VelProf); 
                    //Mpc(VelHist);  
                    std::cout << Time << " " << Sys.delrx <<  std::endl; 
                    OutputStepsIt++; 
                    
                }
                //if (m%1000==0) {
                    //if (!VelProfFile.empty()) Mpc(VelProf);     
                //}
            }
            m++; 
        }
        
    
    
    } 
    

    //////// PARALLEL SECTION END //////////
    
    gettimeofday(&end, NULL); 
    
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/Mpc.NumberOfParticles/TotalSteps << std::endl;
    
       
    if (!VelProfFile.empty()) {
        std::ofstream VelProfFileStream{}; 
        VelProfFileStream.open(VelProfFile, ios::out | ios::trunc); 
        VelProf.print_result(VelProfFileStream); 
    }
    
    //std::ofstream VelHistFileStream{}; 
    //VelHistFileStream.open("velhist", ios::out | ios::trunc); 
    //VelHist.print_result(VelHistFileStream); 
    
    //fclose(PDBout);   

    FluidFilePointer = fopen(FluidFile.c_str(), "w"); 
    Mpc.printFluid(FluidFilePointer, Time); 
    fclose(FluidFilePointer);  
        
    return EXIT_SUCCESS;
}
