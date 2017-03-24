#include <exception> 

class LibraryException : public std::exception 
{
    public: 
        ~LibraryException() throw() {}; 
        const char* what() const throw() {return "An exception regarding the library 'hybridsim' occurred"; } 
};

class ForceException : public LibraryException 
{
        public: 
        ~ForceException() throw() {}; 
        const char* what() const noexcept {return "An overflow error occurred in a force calculation. \n"; }     

}; 

class FENEException : public ForceException 
{
    public: 
        int IndexFirst; 
        int IndexSecond; 
        double Force; 
        std::string ErrorMessage; 
        ~FENEException() throw() {};  
        const char* what() const noexcept {return ErrorMessage.c_str(); }
        FENEException(int First, int Second, double aForce): 
            IndexFirst{First}, 
            IndexSecond {Second}, 
            Force(aForce) {
                ErrorMessage = "An overflow error occurred in a FENE force calculation between the following particles: \n Particle 1: "+std::to_string(IndexFirst)+" \n Particle 2: "+std::to_string(IndexSecond)+" \n Force: "+std::to_string(Force)+"\n" ; 
            }
};    

class RLJException : public ForceException 
{
    public: 
        std::string ErrorMessage; 
        ~RLJException() throw() {};  
        const char* what() const noexcept {return ErrorMessage.c_str(); }
        RLJException(int First, int Second, double aForce) 
        {
            ErrorMessage = "An overflow error occurred in a Lennard-Jones force calculation between the following particles: \n Particle 1: "+std::to_string(First)+" \n Particle 2: "+std::to_string(Second)+" \n Force: "+std::to_string(aForce)+"\n" ; 
        }
}; 

class CellAllocationException : public std::exception {
    public: 
        std::string ErrorMessage; 
        ~CellAllocationException() throw() {}; 
        const char* what() const noexcept {return ErrorMessage.c_str(); }
        CellAllocationException(const MDParticle part, const std::array<int, 3> CellIndex) {
            ErrorMessage = "Bad allocation of Particle "+std::to_string(part.Identifier)+" with Position {"+std::to_string(part.Position(0))+", "+std::to_string(part.Position(1))+", "+std::to_string(part.Position(2))+"} at CellIndex {"+std::to_string(CellIndex[0])+", "+std::to_string(CellIndex[1])+", "+std::to_string(CellIndex[2])+"}"; 
        } 
}; 

