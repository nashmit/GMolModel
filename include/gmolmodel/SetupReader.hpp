#include "Robo.hpp"

class SetupReader{
    // Constructor
    SetupReader(std::ifstream& setupFile);

    // Destructor
    ~SetupReader();

    // Read setup function
    void ReadSetup(std::ifstream& setupFile);
};
