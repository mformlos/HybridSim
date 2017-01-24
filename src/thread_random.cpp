#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include "Rand.h"
#include <sys/time.h>
#include <algorithm>
#include <random>

double real_uniform(std::mt19937_64& engine) {
    uniform_real_distribution<double> dis; 
    return dis(engine); 
} 
extern std::mt19937_64 my_engine; 
#pragma omp threadprivate(my_engine)

int main() {

    timeval start, end; 

    gettimeofday(&start, NULL);
    int tid; 
    //Rand::seed(1); 
    /*std::cout << "some random numbers ... "; 
    for (int i = 0; i < 3; i++) {
        std::cout << Rand::real_uniform() << " "; 
    }
    std::cout << std::endl; 
    */
 
    int bli {-2}, bla {3}; 
    int blu {};
    blu = floor((double)bli/bla); 
    std::cout << "-2%3 equals " << blu << std::endl; 
    


    #pragma omp parallel private(tid)
    { 

    tid = omp_get_thread_num(); 
    std::mt19937_64 my_engine(1); 
    my_engine.seed(tid); 
    /*std::minstd_rand0 lc_generator(tid);
    std::uint_least32_t seed_data[std::mt19937::state_size]; 
    std::generate_n(seed_data, std::mt19937::state_size, std::ref(lc_generator)); 
    std::seed_seq q(std::begin(seed_data), std::end(seed_data)); 
    */
    int seed = ((tid+1)*57)%71; 
    std::cout << "my seed: " << seed << std::endl; 
    Rand::seed(((tid+1)*13)%71); 
    Rand::warmup(10000); 
    //std::mt19937 engine;
    //engine.seed(61737); 
    //std::uniform_real_distribution<double> dist(0.0, 1.0); 
    //Rand::seed(tid);   
    #pragma omp barrier 
 
    for (int i = 0; i < 6; i++) {
        double n {Rand::real_uniform()};
        //double n {dist(engine)}; 
        std::cout << "This is process " << tid << " and my next number is: " << n << std::endl;  
    }
    }
    gettimeofday(&end, NULL);
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << std::endl; 
    return 0;
}

