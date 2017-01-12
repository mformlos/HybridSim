/*
 * Rand.h
 *
 *  Created on: Jan 4, 2017
 *      Author: maud
 */

#ifndef RAND_H_
#define RAND_H_

#include <random>
#include <chrono>


using namespace std;

class Rand{
private:
	static mt19937_64 generator;
	static normal_distribution<double> dis_normal;
	static uniform_real_distribution<double> dis_uniform;
	static uniform_int_distribution<int> dis_intuniform;
	static chi_squared_distribution<double> dis_chisquared;
	static gamma_distribution<double> dis_gamma;
public:
	static void init();
	static double real_normal(double, double);
	static double real_normal();
	static double real_uniform();
	static double real_uniform(double);
	static double real_uniform(double, double);
	static double real_chisquared(unsigned n);
	//static double real_gamma(double shape);
	static double real_gamma(double, double);
	static void seed(int);

};

#endif /* RAND_H_ */
