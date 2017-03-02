
#ifndef LIB_VELOCITY_PROFILE_H_
#define LIB_VELOCITY_PROFILE_H_

#include "Analysis.h"
#include <map>
#include <iterator>
#include <fstream>
//#include "Function_Output.h"
//#include <../eigen/Eigen/Dense>

struct vel_xy {
	double vel_x;
	double vel_y;
};

class VelocityProfile: public Analysis<Particle> {
protected:
	std::map<double, vel_xy> vel_x_average;
	std::map<double, vel_xy>::iterator vel_x_average_iter;
	std::map<double, double> vel_x_average_count;
	std::map<double, double>::iterator vel_x_average_count_iter;

	std::map<double, vel_xy> vel_y_average;
	std::map<double, vel_xy>::iterator vel_y_average_iter;
	std::map<double, double> vel_y_average_count;
	std::map<double, double>::iterator vel_y_average_count_iter;



	double width;
	//Function_Output vel_x_average;

public:
	VelocityProfile(double a_width = 0.5) :
		vel_x_average { },
		vel_x_average_iter { },
		vel_x_average_count { },
		vel_x_average_count_iter { },
		width {a_width} {}

	void initialize(string filename) {
		ifstream input {filename};
		double pos{ }, vel_x { }, vel_y { };
		int count { };
	    while(input >> pos >> vel_x >> vel_y >> count) {
	    	double x {floor(pos/width)*width};
	    	vel_x_average[x].vel_x += vel_x*count;
	    	vel_x_average[x].vel_y += vel_y*count;
	    	vel_x_average_count[x] += count;
	    }
    	vel_x_average_iter = vel_x_average.begin();
		vel_x_average_count_iter = vel_x_average_count.begin();

		while(input >> pos >> vel_x >> vel_y >> count) {
			double y {floor(pos/width)*width};
			vel_y_average[y].vel_x += vel_x*count;
			vel_y_average[y].vel_y += vel_y*count;
			vel_y_average_count[y] += count;
		}
		vel_y_average_iter = vel_y_average.begin();
		vel_y_average_count_iter = vel_y_average_count.begin();

	}


	void operator() (const Particle& part) {
		//Vector3d vel {part.Velocity};
		//vel(0) -= SimBox.COM_Pos(1)*0.003;
		double x {floor(part.Position(0)/ width) * width };
		vel_x_average[x].vel_x += part.Velocity(0);
		vel_x_average[x].vel_y += part.Velocity(1);
		//double y {floor(pos(0)/width)*width };
		//vel_x_average[y] += vel(1);
		vel_x_average_iter = vel_x_average.begin();
		vel_x_average_count[x]++;
		vel_x_average_count_iter = vel_x_average_count.begin();

		double y {floor(part.Position(1)/ width) * width };
		vel_y_average[y].vel_x += part.Velocity(0);
		vel_y_average[y].vel_y += part.Velocity(1);
		//double y {floor(pos(0)/width)*width };
		//vel_x_average[y] += vel(1);
		vel_y_average_iter = vel_y_average.begin();
		vel_y_average_count[y]++;
		vel_y_average_count_iter = vel_y_average_count.begin();

	}

	double value() {return 1.0; }
	std::ostream& print_result(std::ostream& os) {
		bool out {true};
		do {
			os.precision(8);
			os << std::scientific;
			os << vel_x_average_iter->first << ' ';
			os << (vel_x_average_iter->second.vel_x) / (vel_x_average_count_iter -> second) <<  " " << (vel_x_average_iter->second.vel_y) / (vel_x_average_count_iter -> second) << " " << vel_x_average_count_iter -> second << '\n';


			os << std::flush;
			auto vel_x_average_iter_buf = vel_x_average_iter;
			if (++vel_x_average_iter_buf == vel_x_average.end()) out = false;
			++vel_x_average_iter;
			++vel_x_average_count_iter;

		} while(out);
		vel_x_average_iter = vel_x_average.begin();
		vel_x_average_count_iter = vel_x_average_count.begin();

		os << '\n';

		out = true;
		do {
			os.precision(8);
			os << std::scientific;
			os << vel_y_average_iter->first << ' ';
			os << (vel_y_average_iter->second.vel_x) / (vel_y_average_count_iter -> second) <<  " " << (vel_y_average_iter->second.vel_y) / (vel_y_average_count_iter -> second) << " " << vel_y_average_count_iter -> second << '\n';


			os << std::flush;
			auto vel_y_average_iter_buf = vel_y_average_iter;
			if (++vel_y_average_iter_buf == vel_y_average.end()) out = false;
			++vel_y_average_iter;
			++vel_y_average_count_iter;

		} while(out);

		vel_y_average_iter = vel_y_average.begin();
		vel_y_average_count_iter = vel_y_average_count.begin();

		return os;
	}



};




#endif /* LIB_VELOCITY_X_H_ */
