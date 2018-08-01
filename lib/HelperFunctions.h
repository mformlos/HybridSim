#include <../Eigen/Eigen/Dense>
#include <sstream>

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
