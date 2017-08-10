#ifndef __IMONTECARLOSAMPLER_HPP__
#define __IMONTECARLOSAMPLER_HPP__

#include "Robo.hpp"
#include "bMainResidue.hpp"
#include "CurrentState.hpp"

class IMonteCarloSampler
{
public:

    // Constructor

    IMonteCarloSampler(bMainResidue *residue, CurrentState *currentState);

    // Destructor

    virtual ~IMonteCarloSampler() = 0;

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation

    virtual void update(void) = 0;

};

#endif // __IMONTECARLOSAMPLER_HPP__
