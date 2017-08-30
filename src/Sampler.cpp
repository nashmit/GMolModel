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

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(SimTK::State& someState)
{
    return 0;
}

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(const SimTK::State& someState)
{
    return 0;
}

// Update - to be implemented by every specific sampler

void Sampler::update(SimTK::State& somState){}

