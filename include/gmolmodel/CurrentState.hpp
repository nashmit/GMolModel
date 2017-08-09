#ifndef __CURRENTSTATE_HPP__
#define __CURRENTSTATE_HPP__

#include "Robo.hpp"
#include "ICurrentState.hpp"
#include "bMainResidue.hpp"

class CurrentState : public ICurrentState, public SimTK::State 
{
public:

    // Constructor

    CurrentState(bMainResidue *residue, CurrentState *gmolstate);

    // Destructor

    ~CurrentState();

};

#endif // __ICURRENTSTATE_HPP__
