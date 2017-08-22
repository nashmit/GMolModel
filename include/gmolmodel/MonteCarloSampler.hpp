#ifndef __MONTECARLOSAMPLER_HPP__
#define __MONTECARLOSAMPLER_HPP__

#include "Robo.hpp"
#include "Sampler.hpp"

class Topology;
class IState;

class MonteCarloSampler : public Sampler
{
public:
    // Constructor

    MonteCarloSampler(Topology *residue, IState *currentState);

    // Destructor

    ~MonteCarloSampler();

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
   void update(void);

};

#endif // __MONTECARLOSAMPLER_HPP__

