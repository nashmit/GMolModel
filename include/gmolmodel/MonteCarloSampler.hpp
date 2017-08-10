#ifndef __MONTECARLOSAMPLER_HPP__
#define __MONTECARLOSAMPLER_HPP__

#include "Robo.hpp"
#include "IMonteCarloSampler.hpp"
#include "Sampler.hpp"
#include "bMainResidue.hpp"
#include "CurrentState.hpp"

class MonteCarloSampler : public IMonteCarloSampler, public Sampler
{
public:
    // Constructor

    MonteCarloSampler(bMainResidue *residue, CurrentState *currentState);

    // Destructor

    ~MonteCarloSampler();

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
   void update(void);

};

#endif // __MONTECARLOSAMPLER_HPP__

