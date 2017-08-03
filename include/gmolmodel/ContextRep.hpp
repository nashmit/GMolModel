#ifndef GMOL_CONTEXTREP_H_
#define GMOL_CONTEXTREP_H_

#include <iostream>
#include <vector>
#include "bMoleculeReader.hpp"

class ContextRep{
public:

    // Constructor

    ContextRep();

    // Destructor

    ~ContextRep();

    // 


    // Print functions

    void Print(void);

private:

    // Atoms list and bonds list

    std::vector<bSpecificAtom> bAtomList;
    std::vector<bBond> bonds;

};

#endif // GMOL_CONTEXTREP_H_
