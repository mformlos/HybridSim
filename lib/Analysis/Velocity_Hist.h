
#ifndef LIB_VELOCITY_HIST_H_
#define LIB_VELOCITY_HIST_H_

#include "Analysis.h"
#include <map>
#include <iterator>
#include <fstream>

class VelocityHist: public Analysis<Particle> {
protected: 
    std::map<double, double> hist; 
    std::map<double, double>::iterator hist_iter; 
    unsigned count; 
    
    double width; 
    
public: 
    VelocityHist(double a_width = 0.1) :
        hist {}, 
        hist_iter {}, 
        count {}, 
        width(a_width) {}
        
    void operator() (const Particle& part) {
        double vel {sqrt(pow(part.Velocity(0),2)+pow(part.Velocity(1),2)+pow(part.Velocity(2),2))};
        double bin {floor(vel/ width) * width }; 
        hist[bin] += 1.; 
        count++;
        hist_iter = hist.begin(); 
    }
    
    double value() {return 1.0;}
    std::ostream& print_result(std::ostream& os) {
		bool out {true};
		do {
			os.precision(8);
			os << std::scientific;
			os << hist_iter->first << ' ';
			os << (hist_iter->second) / (count) << '\n';


			os << std::flush;
			auto hist_iter_buf = hist_iter;
			if (++hist_iter_buf == hist.end()) out = false;
			++hist_iter;

		} while(out);
		hist_iter = hist.begin();
		return os;
	}
};

#endif
