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

    MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem, SimTK::SimbodyMatterSubsystem *argMatter, Topology *argResidue);

    // Destructor

    ~MonteCarloSampler();


    void assignConfFromTVector(SimTK::State& advanced);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(void);

private:
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    Topology *residue;
    SimTK::Transform *TVector;

};

#endif // __MONTECARLOSAMPLER_HPP__

