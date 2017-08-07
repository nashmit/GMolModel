#ifndef GMOL_CONTEXT_H_
#define GMOL_CONTEXT_H_

#include <iostream>
#include <vector>
#include "bMoleculeReader.hpp"

class InitialContext{
public:

    // Constructor

    InitialContext();

    // Destructor

    ~InitialContext();

    // Insert an atom

    void addAtom(int number);

    // Print functions

    void Print(void);

private:
    std::vector<bSpecificAtom> bAtomList;

};

#endif // GMOL_CONTEXT_H_
