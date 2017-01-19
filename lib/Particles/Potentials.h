#ifndef LIB_POTENTIALS_H_
#define LIB_POTENTIALS_H_

#include <cmath>

inline double FENE_Potential(double r2) {
    double potential {}; 
    potential = -15.0*2.25*log(1.0 - r2/2.25);
    return potential; 
}

inline double FENE_Force(double r2) {
    double force {}; 
    force = -30.0/(1.0 - (r2/2.25));
    return force; 
}

inline double RLJ_Potential(double r2) {
    double potential {}; 
    if (r2 < 1.0594631) {
		double rm24 { 1.0 / r2 };
		rm24 *= rm24*rm24;
		rm24 *= rm24;
		rm24 *= rm24;
		potential = 4.0*rm24*(rm24 - 1.0) + 1.0;
	}
	return potential;
} 

inline double RLJ_Force(double r2) {
    double force{};
	if (r2 < 1.0594631) {
		double rm2 { 1.0 / r2 };
		double rm24 { rm2*rm2*rm2 };
		rm24 *= rm24;
		rm24 *= rm24;
		force = 192.0*rm2*rm24*(rm24 - 0.5);
		if (fabs(force) > 1e4 ) std::cout << "LJ!" << std::endl;
	}
	return force;
}

#endif
