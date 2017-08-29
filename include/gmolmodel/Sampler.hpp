#ifndef __SAMPLER_HPP__
#define __SAMPLER_HPP__

#include "Robo.hpp"
//#include "Topology.hpp"
//#include "IState.hpp"

class Topology;

class Sampler
{
public:
    // Constructor

    Sampler(SimTK::CompoundSystem *argCompoundSystem,
            SimTK::SimbodyMatterSubsystem *argMatter,
            Topology *argResidue,
            SimTK::TimeStepper *argTimeStepper);


    // Destructor

    ~Sampler();

protected:
    const SimTK::System *system;
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    Topology *residue;
    SimTK::TimeStepper *timeStepper;

};


#endif // __SAMPLER_HPP__

