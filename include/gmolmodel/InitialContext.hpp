#ifndef GMOL_CONTEXT_H_
#define GMOL_CONTEXT_H_

#include <iostream>
#include <vector>
#include "IInitialContext.hpp"
#include "bMoleculeReader.hpp"

class InitialContext : public IInitialContext{
public:

    // Constructor

    InitialContext();

    // Destructor

    ~InitialContext();

    // Insert an atom

    void addAtom(int number);

    // Print functions

    void PrintAtomList(void);

private:
    std::vector<bSpecificAtom> AtomList;

};

#endif // GMOL_CONTEXT_H_
