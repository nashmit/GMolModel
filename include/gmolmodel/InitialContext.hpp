#ifndef GMOL_CONTEXT_H_
#define GMOL_CONTEXT_H_

#include <iostream>
#include "InitialContextRep.hpp"

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

    // Private implementation

    InitialContextRep *rep;

};

#endif // GMOL_CONTEXT_H_
