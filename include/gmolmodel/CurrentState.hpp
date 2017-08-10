#ifndef __CURRENTSTATE_HPP__
#define __CURRENTSTATE_HPP__

#include "Robo.hpp"
#include "ICurrentState.hpp"
#include "bMainResidue.hpp"

class CurrentState : public ICurrentState, public SimTK::State 
{
public:

    // Constructor

    CurrentState();

    // Destructor

    ~CurrentState();

};

#endif // __ICURRENTSTATE_HPP__
