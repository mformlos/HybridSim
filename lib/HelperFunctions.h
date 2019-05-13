#include <../Eigen/Eigen/Dense>
#include <sstream>
#include "Particle.h"

using namespace Eigen; 

template <typename Type>
Type extractParameter(std::string key, std::ifstream& inputfile, bool& found) {
    inputfile.clear(); 
    inputfile.seekg(0, std::ios::beg);
    std::string line{}, fkey{}, equal{};
    Type param {};   
    found = false; 
    while (getline(inputfile, line)){
        if (line.find(key) != std::string::npos) {
            std::istringstream iss(line); 
            if (iss >> fkey >> equal >> param) {
                found = true; 
                return param; 
            }
        }
    }
    std::cout << "keyword not found" << std::endl; 
    return param; 
}

template <typename Type>
bool initializeVector(std::vector<Type>& vec, std::string filename) {
    std::ifstream file (filename, std::ios::in); 
    if (!file.is_open()) {
        return false; 
    }
    Type in{}; 
    while(file >> in) {
        vec.push_back(in); 
    } 
    return true;
}

bool initializeStepVector(std::vector<unsigned long long>& vec, std::string filename) {
    std::ifstream file (filename, std::ios::in); 
    if (!file.is_open()) {
        return false; 
    }
    unsigned long long Step{}; 
    while(file >> Step) {
        vec.push_back(Step); 
    } 
    return true; 
} 

bool initializeStringVector(std::vector<std::string>& vec, std::string filename) {
	vec.clear();
	std::ifstream file (filename, std::ios::in);
	if (!file.is_open()) return false;
	std::string astring {};
	while(file >> astring) vec.push_back(astring);
	file.close(); 
	return true;
}

struct ForceUpdate {
    unsigned long long Step; 
    Vector3d Force; 
    ForceUpdate(unsigned long long s, Vector3d f) : 
        Step {s}, 
        Force {f} {}
}; 

struct ConstraintUpdate {
    unsigned long long Step; 
    double Constraint; 
    ConstraintUpdate(unsigned long long s, double c) : 
        Step {s}, 
        Constraint {c} {}
}; 

bool initializeForceUpdateVector(std::vector<ForceUpdate>& vec, std::string filename) {
    std::ifstream file (filename, std::ios::in); 
    if (!file.is_open()) {
        return false; 
    }
    unsigned long long Step{};
    double fx{}, fy{}, fz{}; 
    while(file >> Step >> fx >> fy >> fz) {
        Vector3d Force(fx, fy, fz); 
        vec.push_back(ForceUpdate(Step, Force)); 
    } 
    return true; 
} 

bool initializeConstraintUpdateVector(std::vector<ConstraintUpdate>& vec, std::string filename) {
    std::ifstream file (filename, std::ios::in); 
    if (!file.is_open()) {
        return false; 
    }
    unsigned long long Step{};
    double c {};  
    while(file >> Step >> c) {
        vec.push_back(ConstraintUpdate(Step, c)); 
    } 
    return true; 
} 

bool initializePositions(std::vector<Particle>& vec, std::string filename) {
    std::ifstream file {filename};
    std::string line{};  
    if (!file.is_open()) return false;
    std::string dump, Type;
    double x, y, z;
    unsigned count {0};
    if (file.is_open()) {
        file >> dump >> dump;
        while(getline(file, line)) {
            std::istringstream iss(line); 
            if (iss >> dump >> dump >> Type >> dump >> dump >> x >> y >> z >> dump >> dump >> dump)        
            { 
			    vec[count].Position(0) = x;
			    vec[count].Position(1) = y;
			    vec[count].Position(2) = z;
			    count++;
			}
		}
        if (count != vec.size()) {
        	std::cout << "only " << count << " monomers were initialized" << std::endl;
            return false;
        }
    }
    file.close(); 
    return true;
}
