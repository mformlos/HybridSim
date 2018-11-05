#include <sys/time.h>
#include <omp.h>
#include <csignal>
#include "MPC.h"
#include "Velocity_Profile.h"
#include "Velocity_Hist.h"
#include "HelperFunctions.h"
//#include "System.h"


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
    unsigned long long EquilSteps{}, SimSteps{}, TotalSteps {}, n {}, m {}; 
    int tid{}, procs{}, maxt{}, inpar{}, dynamic{}, nested{}, nthreads{}; 
    double MPCStep{}, Shear{}, TotalTime{}, EquilTime{}, SimTime{}, Temperature{}, Time{};
    bool ParameterInitialized{}; 
    std::string OutputStepFile{}, FluidOutputStepFile{}, StatisticsFile{}, VelProfFile{}, VelHistFile{}, FluidFile{}, FluidInput{};
    std::vector<unsigned long long> OutputSteps; 
    std::vector<unsigned long long>::iterator OutputStepsIt{}; 
    
    std::vector<unsigned long long> FluidOutputSteps; 
    std::vector<unsigned long long>::iterator FluidOutputStepsIt{}; 
 
    if (argc != 2) {
        std::cout << "usage: ./fluid_only PARAMETER-INPUT-FILE " << std::endl;  
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
    FluidOutputStepFile = extractParameter<std::string>("FluidOutputStepFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE; 
    StatisticsFile = extractParameter<std::string>("StatisticsFile", inputfile,    ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    FluidFile = extractParameter<std::string>("FluidFile", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    FluidInput = extractParameter<std::string>("FluidInput", inputfile, ParameterInitialized);
    if (!ParameterInitialized) return EXIT_FAILURE;
    VelProfFile = extractParameter<std::string>("VelProfFile", inputfile, ParameterInitialized);
    VelHistFile = extractParameter<std::string>("VelHistFile", inputfile, ParameterInitialized);


    inputfile.close(); 
    TotalTime = EquilTime + SimTime;
    
    EquilSteps = (unsigned long long) (EquilTime/MPCStep); 
    SimSteps = (unsigned long long) (SimTime/MPCStep); 
    TotalSteps = (unsigned long long) (TotalTime/MPCStep); 
    std::cout << "Total number of MD steps: " << TotalSteps << std::endl;  
    
    ///////////////////////////////////////// 
     
    /////// PARAMETER CONTROL OUTPUT //////// 
    std::cout << "Lx is " << Lx << std::endl;
    std::cout << "Ly is " << Ly << std::endl;
    std::cout << "Lz is " << Lz << std::endl;
    std::cout << "Temperature is " << Temperature << std::endl;
    std::cout << "Shear is " << Shear << std::endl;
    std::cout << "MPC particles per cell is " << MPCRho << std::endl;
    std::cout << "MPCStep is " << MPCStep << std::endl;
    std::cout << "StartTime is " << Time << std::endl;  
    std::cout << "EquilTime is " << EquilTime << std::endl;    
    std::cout << "TotalTime is " << TotalTime << std::endl;
    std::cout << "RNG seed is " << Seed << std::endl; 
    std::cout << "OutputStepFile is " << OutputStepFile << std::endl;
    std::cout << "FluidFile is " << FluidFile << std::endl;
    std::cout << "FluidInput is " << FluidInput << std::endl;
    std::cout << "VelProfFile is " << VelProfFile << std::endl;
    std::cout << "VelHistFile is " << VelHistFile << std::endl;
    
    /////////////////////////////////////
    
    ////// RANDOM ENGINE SEEDING & WARMUP //////

    Rand::seed(Seed); 
    Rand::warmup(10000); 
    
    /////////////////////////////////////
    
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
    
    /////////////////////////////////////
    
    /////// OUTPUT INITIALIZATION ///////
    VelocityProfile VelProf{0.2};  
    VelocityHist VelHist{0.1};
    
    if (!initializeStepVector(OutputSteps, OutputStepFile)) {
        std::cout << "OutputStepFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    }
    std::cout << "Output will be done " << OutputSteps.size() << " times. " << std::endl; 
    OutputStepsIt = OutputSteps.begin(); 
    
    if (!initializeStepVector(FluidOutputSteps, FluidOutputStepFile)) {
        std::cout << "FluidOutputStepFile does not exist!" << std::endl; 
        return EXIT_FAILURE; 
    }
    std::cout << "Fluid output will be done " << FluidOutputSteps.size() << " times. " << std::endl; 
    FluidOutputStepsIt = FluidOutputSteps.begin(); 
    
    std::ofstream StatisticsStream(StatisticsFile, std::ios::out | std::ios::app); 
    
    FILE* FluidFilePointer {}; 
    
    timeval start {}, end {};
    gettimeofday(&start, NULL); 
    
    ////////////////////////////////////////
    
    /////// OPENMP BEGINNING & TESTS ///////
    
    #pragma omp parallel private(tid, n, m)
    {
        /////// OMP PARAMETERS ///////
        #ifdef _OPENMP
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
        #endif
        ////// RANDOM ENGINE SEEDING & WARMUP for openMP //////
        #ifdef _OPENMP
        Rand::seed(tid); 
        Rand::warmup(10000); 
        #endif
        
        ///////////////////////////////////
        ////// MAIN SIMULATION LOOP ///////
        m = -EquilSteps + (unsigned long long)(Time/MPCStep); 
        for (n = 0; n <= TotalSteps; n++) {
        ////// MPC SOLUTE TRANSFER & GRID SHIFT ////
            #pragma omp single 
            {
                Mpc.updateBoxShift(MPCStep); 
                Mpc.shiftGrid(); 
            }
            ////// STREAMING STEP ///////
            #pragma omp for schedule(static)
            for (unsigned part = 0; part < Mpc.NumberOfParticles; part++) {
                Mpc.streamPlusCellAssignment(Mpc.Fluid[part], MPCStep); 
            }
            ////// SORT MPC PARTICLE VECTOR TO IMPROVE DATA LOCALITY ///// 
            #pragma omp single 
            {
                if (n%10==0) Mpc.sortVector(); 
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
            ////////// MPC STEP FINISHED //////////
            
            ////////// OUTPUT ////////////
            #pragma omp single 
            {
                if (SignalCaught) {
                    std::cout << "writing job data..." << std::endl; 
                    FluidFilePointer = fopen((FluidFile+std::to_string(m)).c_str(), "w"); 
                    Mpc.printFluid(FluidFilePointer, Time);
                    fclose(FluidFilePointer);  
                    if (!VelProfFile.empty()) {
                        std::ofstream VelProfFileStream{}; 
                        VelProfFileStream.open(VelProfFile+std::to_string(m), ios::out | ios::trunc); 
                        VelProf.print_result(VelProfFileStream); 
                        VelProfFileStream.close(); 
                    }  
                    if (!VelHistFile.empty()) {
                        std::ofstream VelHistFileStream{};
                        VelHistFileStream.open(VelHistFile+std::to_string(m), ios::out | ios::trunc); 
                        VelHist.print_result(VelHistFileStream); 
                        VelHistFileStream.close(); 
                    }
                    std::cout << "Terminating..." << std::endl;   
                    std::terminate(); 
                }
                if (m == *OutputStepsIt) {
                    StatisticsStream << Time << " " << Mpc.virtualTemperature() << std::endl; 
                    if (!VelProfFile.empty()) Mpc(VelProf);
                    if (!VelHistFile.empty()) Mpc(VelHist);  
                    std::cout << Time << " " << Mpc.delrx <<  std::endl; 
                    OutputStepsIt++;  
                }
                if (m == *FluidOutputStepsIt) {
                    std::cout << "printing fluid... " <<  std::endl;
                    FluidFilePointer = fopen((FluidFile+std::to_string(m)).c_str(), "w"); 
                    Mpc.printFluid(FluidFilePointer, m); 
                    fclose(FluidFilePointer); 
                    if (!VelProfFile.empty()) {
                        std::ofstream VelProfFileStream{};
                        VelProfFileStream.open(VelProfFile+std::to_string(m), ios::out | ios::trunc); 
                        VelProf.print_result(VelProfFileStream); 
                        VelProfFileStream.close(); 
                    }
                    if (!VelHistFile.empty()) {
                        std::ofstream VelHistFileStream{};
                        VelHistFileStream.open(VelHistFile+std::to_string(m), ios::out | ios::trunc); 
                        VelHist.print_result(VelHistFileStream); 
                        VelHistFileStream.close(); 
                    }
                    FluidOutputStepsIt++;   
                }
            }
            Time += MPCStep; 
            m++; 
        }    
    }
    //////// PARALLEL SECTION END //////////
    
    gettimeofday(&end, NULL); 
    
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << " , time per particle and step: " << realTime/Mpc.NumberOfParticles/TotalSteps << std::endl;
    
    if (!VelProfFile.empty()) {
        std::ofstream VelProfFileStream{};
        VelProfFileStream.open(VelProfFile+std::to_string(m), ios::out | ios::trunc); 
        VelProf.print_result(VelProfFileStream); 
        VelProfFileStream.close(); 
    }
    if (!VelHistFile.empty()) {
        std::ofstream VelHistFileStream{};
        VelHistFileStream.open(VelHistFile+std::to_string(m), ios::out | ios::trunc); 
        VelHist.print_result(VelHistFileStream); 
        VelHistFileStream.close(); 
    }
    
    FluidFilePointer = fopen((FluidFile+std::to_string(m)).c_str(), "w"); 
    Mpc.printFluid(FluidFilePointer, Time); 
    fclose(FluidFilePointer);  
    
    return EXIT_SUCCESS;     
}
