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

SimTK::Real Sampler::calcMassDeterminant(SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Matrix D0(6, 6);
    Eigen::MatrixXd EiD0(6, 6);
    matter->calcDetM(state, V, DetV, D0);

    for(int i=0; i<6; i++){
        for(int j=0; j<6; j++){
            EiD0(i, j) = D0(i, j);
        }
    }
    SimTK::Real EiDetD0 = EiD0.determinant();
    for(int i=6; i<nu; i++){
        EiDetD0 *= DetV[i];
    }
    return EiDetD0;
}

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(const SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Matrix D0(6, 6);
    Eigen::MatrixXd EiD0(6, 6);
    matter->calcDetM(state, V, DetV, D0);

    for(int i=0; i<6; i++){
        for(int j=0; j<6; j++){
            EiD0(i, j) = D0(i, j);
        }
    }
    SimTK::Real EiDetD0 = EiD0.determinant();
    for(int i=6; i<nu; i++){
        EiDetD0 *= DetV[i];
    }
    return EiDetD0;
}

// Update - to be implemented by every specific sampler

void Sampler::update(SimTK::State& somState){}

