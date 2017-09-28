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
                                 //Topology *argResidue,
                                 SimTK::Compound *argResidue,
                                 SimTK::TimeStepper *argTimeStepper);

    // Destructor

    virtual ~HamiltonianMonteCarloSampler();

    // Assign a random conformation. Time measured in picoseconds

    void propose(SimTK::State& advanced, SimTK::Real timestep, int nosteps);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation

    void update(SimTK::State&);

    // Set kinetic energies    

    void setOldKE(SimTK::Real);

protected:
       SimTK::Real ke_o; 

};

#endif // __HAMMONTECARLOSAMPLER_HPP__

