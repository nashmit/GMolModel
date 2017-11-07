/**@file
Implementation of Sampler class. **/

#include "Sampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

Sampler::Sampler(SimTK::CompoundSystem *argCompoundSystem,
                 SimTK::SimbodyMatterSubsystem *argMatter,
                 //Topology *argResidue,
                 SimTK::Compound *argResidue,
                 SimTK::DuMMForceFieldSubsystem *argDumm,
                 SimTK::GeneralForceSubsystem *argForces,
                 SimTK::TimeStepper *argTimeStepper)
{
    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    this->residue = argResidue;
    this->dumm = argDumm;
    this->forces = argForces;
    this->timeStepper = argTimeStepper;
    this->system = &(matter->getSystem());
}

// Destructor

Sampler::~Sampler(){
    ;
}

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real *D0 = new SimTK::Real(1.0);
    matter->calcDetM(state, V, DetV, D0);
    return *D0;
}

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(const SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real *D0 = new SimTK::Real(1.0);
    matter->calcDetM(state, V, DetV, D0);
    return *D0;
}

// Update - to be implemented by every specific sampler

void Sampler::update(SimTK::State& somState){}

