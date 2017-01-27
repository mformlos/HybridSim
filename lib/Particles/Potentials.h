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
    if (fabs(force) > 1e4 ) std::cout << "FENE! force: " << force << " radius squared " << r2 << std::endl;
    return force; 
}

inline double RLJ_Potential(double r2) {
    double potential {}; 
    if (r2 < 1.25992105) {
		double rm6 { 1.0 / r2 };
		rm6 *= rm6*rm6;
		potential = 4.0*rm6*(rm6 - 1.0) + 1.0;
	}
	return potential;
} 

inline double RLJ_Force(double r2) {
    double force{};
	if (r2 <1.25992105) {
		double rm2 { 1.0 / r2 };
		double rm6 { rm2*rm2*rm2 };
		force = 24.*rm2*rm6*(2.*rm6 - 1);
		if (fabs(force) > 1e4 ) std::cout << "LJ! force: " << force << " radius squared " << r2 << std::endl;
	}
	return force;
}

#endif
