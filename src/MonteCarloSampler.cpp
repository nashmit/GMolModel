/**@file
Implementation of MonteCarloSampler class. **/

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

// Stores the configuration into an internal vector of transforms TVector

void MonteCarloSampler::setTVector(SimTK::State& advanced)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(advanced);
    TVector[i] = mobod.getMobilizerTransform(advanced);
    i++;
  }
}


// Restores configuration from the internal vector of transforms TVector

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
    //randomEngine.seed(4294653137UL); // for reproductibility
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SimTK::Real rand_no = uniformRealDistribution_0_2pi(randomEngine);
        mobod.setOneQ(advanced, 0, rand_no);
    }
}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step

void MonteCarloSampler::update(SimTK::State& advanced){
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    SimTK::Real RT = getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;

    // Get energies
    SimTK::Real pe_o = getOldPE();
    SimTK::Real pe_n = getPEFromEvaluator(); // Get potential energy from OPENMM

    // Apply Metropolis criterion
    assert(!isnan(pe_n));
    if ((pe_n < pe_o) or (rand_no < exp(-(pe_n - pe_o)/RT))){
        setTVector(advanced);
        sendConfToEvaluator(); // Insert configuratin in OPENMM
        setOldPE(pe_n);
    }else{
        assignConfFromTVector(advanced);
        setOldPE(getOldPE());
    }
}

// Get stored potential energy

SimTK::Real MonteCarloSampler::getOldPE(void){return 1.0;}

// Get stored kinetic energy - only necessary for Hamiltonian Monte Carlo

SimTK::Real MonteCarloSampler::getOldKE(void){}

// Store the potential energy

void MonteCarloSampler::setOldPE(SimTK::Real argPE){}

// Store the kintetic energy - only necessary for Hamiltonian Monte Carlo

void MonteCarloSampler::setOldKE(SimTK::Real argKE){}

// Get the potential energy from an external source as far as the sampler
// is concerned - OPENMM has to be inserted here

SimTK::Real MonteCarloSampler::getPEFromEvaluator(void){return 0.0;}

// Get the simulation temperature

SimTK::Real MonteCarloSampler::getTemperature(void){}

// Set the simulatin temperature

void MonteCarloSampler::setTemperature(SimTK::Real argTemperature){temperature = argTemperature;}

// Send configuration to an external evaluator

void MonteCarloSampler::sendConfToEvaluator(void){}







