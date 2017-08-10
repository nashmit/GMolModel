#ifndef __SAMPLER_HPP__
#define __SAMPLER_HPP__

#include "Robo.hpp"
#include "bMainResidue.hpp"
#include "CurrentState.hpp"

class Sampler
{
public:
    // Constructor

    Sampler(bMainResidue *residue, CurrentState *currentState);

    // Destructor

    ~Sampler();
};


#endif // __SAMPLER_HPP__

