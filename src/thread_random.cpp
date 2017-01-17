#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include "NewRand.h"
#include <sys/time.h>

int main() {

    timeval start, end; 

    gettimeofday(&start, NULL);
    int tid; 
    Rand::seed(1); 
    std::cout << "some random numbers ... "; 
    for (int i = 0; i < 3; i++) {
        std::cout << Rand::real_uniform() << " "; 
    }
    std::cout << std::endl; 
    
    #pragma omp parallel private(tid)
    { 
    tid = omp_get_thread_num(); 
    
    Rand::seed(tid); 
 
    for (int i = 0; i < 3; i++) {
        double n {Rand::real_uniform()};
        std::cout << "This is process " << tid << " and my next number is: " << n << std::endl;  
    }
    }
    gettimeofday(&end, NULL);
    double realTime = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    std::cout << "total time: " << realTime << std::endl; 
    return 0;
}

