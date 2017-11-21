#include "Robo.hpp"

// This class implements the functionality for reading arguments from a setup
// file and give back the necessary values as vectors
class SetupReader{
public:

    // Constructor
    SetupReader(const char *FN);
    SetupReader(std::string& FN);

    // Destructor
    ~SetupReader();

    // Read setup function
    void ReadSetup(const char *FN);
    void ReadSetup(std::string& FN);

    // Print all the arguments
    void dump(void);

    // Access values by key
    std::vector<std::string> getValues(const char *argKey);
    std::vector<std::string> getValues(std::string argKey);

private:
    std::map<std::string, std::vector<std::string>> Args; // arguments
    std::map<std::string, std::vector<std::string>>::iterator ArgsIt; // arguments
};
