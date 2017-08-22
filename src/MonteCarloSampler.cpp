#include "MonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

MonteCarloSampler::MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem, SimTK::SimbodyMatterSubsystem *argMatter, Topology *argResidue)
{

    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    this->residue = argResidue;
    TVector = new SimTK::Transform[matter->getNumBodies()];
}

// Destructor

MonteCarloSampler::~MonteCarloSampler()
{

    delete [] TVector;

}

// Transfer coordinates from TVector to compound
 
void MonteCarloSampler::assignConfFromTVector(SimTK::State& advanced)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(advanced, TVector[i]);
    i++;
  }
}

// Assign random conformation
 
void MonteCarloSampler::assignRandomConf(SimTK::State& advanced)
{
    eng.seed(4294653137UL);

    std::cout << "MonteCarloSampler state Qs " << advanced.getQ() << std::endl;

    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SimTK::Real rand_no = urd(eng);
        mobod.setOneQ(advanced, 0, rand_no);
    }

    std::cout << "MonteCarloSampler state Qs " << advanced.getQ() << std::endl;

}



void MonteCarloSampler::update(void){

    //;

}

