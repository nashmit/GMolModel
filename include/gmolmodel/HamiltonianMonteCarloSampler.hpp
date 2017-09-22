#ifndef __HAMMONTECARLOSAMPLER_HPP__
#define __HAMMONTECARLOSAMPLER_HPP__

#include "MonteCarloSampler.hpp"

class Topology;
class IState;

class HamiltonianMonteCarloSampler : public MonteCarloSampler
{
public:
    // Constructor

    HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                 SimTK::SimbodyMatterSubsystem *argMatter,
                                 Topology *argResidue,
                                 SimTK::TimeStepper *argTimeStepper);

    // Destructor

    ~HamiltonianMonteCarloSampler();

    // Assign a random conformation

    void propose(SimTK::State& advanced);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(SimTK::State&);

};

#endif // __HAMMONTECARLOSAMPLER_HPP__

