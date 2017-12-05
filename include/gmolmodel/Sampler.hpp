#ifndef __SAMPLER_HPP__
#define __SAMPLER_HPP__

#include "Robo.hpp"

class Topology;

class Sampler
{
public:
    // Constructor
    Sampler(SimTK::CompoundSystem *argCompoundSystem,
            SimTK::SimbodyMatterSubsystem *argMatter,
            SimTK::Compound *argResidue,
            SimTK::DuMMForceFieldSubsystem *argDumm,
            SimTK::GeneralForceSubsystem *forces,
            SimTK::TimeStepper *argTimeStepper);


    // Destructor
    virtual ~Sampler();

    // Compute mass matrix determinant (O(n))
    SimTK::Real calcMassDeterminant(const SimTK::State& );
    SimTK::Real calcMassDeterminant(SimTK::State& );

    // Update
    virtual void update(SimTK::State&);

    void PrintSimbodyStateCache(SimTK::State& someState);

protected:
    const SimTK::System *system;
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    SimTK::Compound *residue;
    SimTK::DuMMForceFieldSubsystem *dumm;
    SimTK::GeneralForceSubsystem *forces;
    SimTK::TimeStepper *timeStepper;

};


#endif // __SAMPLER_HPP__

