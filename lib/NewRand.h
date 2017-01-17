#ifndef NEWRAND_H_
#define NEWRAND_H_

#include <random>
#include <chrono>


using namespace std;

class Rand{
private:
	static thread_local mt19937_64 generator;
    static thread_local normal_distribution<double> dis_normal;
	static thread_local uniform_real_distribution<double> dis_uniform;
	static thread_local uniform_int_distribution<int> dis_intuniform;
	static thread_local chi_squared_distribution<double> dis_chisquared;
	static thread_local gamma_distribution<double> dis_gamma;

public:
	static thread_local double real_normal(double, double);
	static thread_local double real_normal();
	static thread_local double real_uniform();
	static thread_local double real_uniform(double);
	static thread_local double real_uniform(double, double);
	static thread_local double real_chisquared(unsigned n);
	//static double real_gamma(double shape);
	static thread_local double real_gamma(double, double);

	static void seed(int);

};

#endif /* RAND_H_ */
