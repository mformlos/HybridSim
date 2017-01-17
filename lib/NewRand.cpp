#include "NewRand.h"

using namespace std;

//Define static members
thread_local mt19937_64 Rand::generator(1);
thread_local normal_distribution<double> Rand::dis_normal;
thread_local uniform_real_distribution<double> Rand::dis_uniform;
thread_local uniform_int_distribution<int> Rand::dis_intuniform;
thread_local chi_squared_distribution<double> Rand::dis_chisquared;
thread_local gamma_distribution<double> Rand::dis_gamma;

double Rand::real_normal(){
	return Rand::dis_normal(Rand::generator);
}
double Rand::real_normal(double mean, double std){
	return (Rand::real_normal()*std) + mean;
}

double Rand::real_uniform(){
	return Rand::dis_uniform(Rand::generator);
}

double Rand::real_uniform(double max){
	return Rand::real_uniform()*max;
}
double Rand::real_uniform(double min, double max){
	return Rand::real_uniform()*(max - min) + min;
}

double Rand::real_chisquared(unsigned n){
	dis_chisquared.param((chi_squared_distribution<double>::param_type)n);
	return dis_chisquared(generator);
}

/*double Rand::real_gamma(double shape) {
	dis_gamma.param(gamma_distribution<double>::param_type(shape, 1.0));
	return dis_gamma(generator);
}*/
double Rand::real_gamma(double a, double b) {
	dis_gamma.param(gamma_distribution<double>::param_type(a, b));
	return dis_gamma(generator);
}


void Rand::seed(int s) {
    generator.seed(s);
}


