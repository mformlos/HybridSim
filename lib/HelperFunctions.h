template <typename Type>
Type extractParameter(std::string key, std::ifstream& inputfile, bool& found) {
    inputfile.clear(); 
    inputfile.seekg(0, ios::beg);
    std::string line{}, fkey{}, equal{};
    Type param {};   
    found = false; 
    while (getline(inputfile, line)){
        if (line.find(key) != string::npos) {
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

bool initializeStepVector(std::vector<unsigned>& vec, std::string filename) {
    std::ifstream file (filename, ios::in); 
    if (!file.is_open()) {
        return false; 
    }
    unsigned Step{}; 
    while(file >> Step) {
        vec.push_back(Step); 
    } 
    return true; 
} 
