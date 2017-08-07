#ifndef GMOL_CONTEXTREP_H_
#define GMOL_CONTEXTREP_H_

#include <iostream>
#include <vector>
#include "bMoleculeReader.hpp"

class InitialContextRep{
public:

    // Constructor

    InitialContextRep();

    // Destructor

    ~InitialContextRep();

    // Insert an atom

    void addAtom(int number);

    // Print functions

    void Print(void);

private:

    // Atoms list and bonds list

    std::vector<bSpecificAtom> bAtomList;
    std::vector<bBond> bonds;

};

#endif // GMOL_CONTEXTREP_H_
