#include "Sampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

Sampler::Sampler(SimTK::CompoundSystem *argCompoundSystem,
                 SimTK::SimbodyMatterSubsystem *argMatter,
                 Topology *argResidue,
                 SimTK::TimeStepper *argTimeStepper)
{

    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    this->residue = argResidue;
    this->timeStepper = argTimeStepper;
    this->system = &(matter->getSystem());
}

// Destructor

Sampler::~Sampler(){
    ;
}


