#ifndef GMOL_CONTEXT_H_
#define GMOL_CONTEXT_H_

#include <iostream>
#include "ContextRep.hpp"

class Context{
public:

    // Constructor

    Context();

    // Destructor

    ~Context();

    // Print functions

    void Print(void);

private:

    // Private implementation

    ContextRep *rep;

};

#endif // GMOL_CONTEXT_H_
