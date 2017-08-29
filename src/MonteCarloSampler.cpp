/**@file
Implementation of MonteCarloSampler class. **/

#include "MonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

MonteCarloSampler::MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     Topology *argResidue,
                                     SimTK::TimeStepper *argTimeStepper)
    : Sampler(argCompoundSystem, argMatter, argResidue, argTimeStepper)
{

    TVector = new SimTK::Transform[matter->getNumBodies()];
}

// Destructor

MonteCarloSampler::~MonteCarloSampler()
{

    delete [] TVector;

}

// Stores the configuration into an internal vector of transforms TVector

void MonteCarloSampler::setTVector(const SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
    TVector[i] = mobod.getMobilizerTransform(someState);
    i++;
  }
}


// Restores configuration from the internal vector of transforms TVector

void MonteCarloSampler::assignConfFromTVector(SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(someState, TVector[i]);
    i++;
  }
}

// Assign random conformation
 
void MonteCarloSampler::assignRandomConf(SimTK::State& someState)
{
    //randomEngine.seed(4294653137UL); // for reproductibility

    /*
    std::cout << "State info BEFORE updQ. Time = " << someState.getTime() << std::endl;
    for(int i = 0; i < someState.getNumSubsystems(); i++){
        std::cout << someState.getSubsystemName(SimTK::SubsystemIndex(i)) << " ";
        std::cout << someState.getSubsystemStage(SimTK::SubsystemIndex(i)) << " ";
        std::cout << someState.getSubsystemVersion(SimTK::SubsystemIndex(i)) << std::endl;
    }
    */
    int i = 1;
    for (SimTK::MobilizedBodyIndex mbx(i); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SimTK::Real rand_no = uniformRealDistribution_0_2pi(randomEngine);
        for(int j=0; j<mobod.getNumQ(someState); j++){
            mobod.setOneQ(someState, j, rand_no);
            //someState.updQ()[i] = rand_no;
        }
        i++;
    }

    system->realize(someState, SimTK::Stage::Acceleration); // NECESSARY

    /*
    std::cout << "State info AFTER  updQ. Time = " << someState.getTime() << std::endl;
    for(int i = 0; i < someState.getNumSubsystems(); i++){
        std::cout << someState.getSubsystemName(SimTK::SubsystemIndex(i)) << " ";
        std::cout << someState.getSubsystemStage(SimTK::SubsystemIndex(i)) << " ";
        std::cout << someState.getSubsystemVersion(SimTK::SubsystemIndex(i)) << std::endl;
    }
    */

}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step

void MonteCarloSampler::update(SimTK::State& someState){
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    SimTK::Real RT = getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;

    // Get old energy
    SimTK::Real pe_o = getOldPE();

    // Assign random configuration

    assignRandomConf(someState);
    //timeStepper->initialize(someState); // should develop timestepper class ??

    // Send configuration to evaluator  

    sendConfToEvaluator(); // OPENMM

    // Get current potential energy from evaluator

    SimTK::Real pe_n = getPEFromEvaluator(); // OPENMM

    // Apply Metropolis criterion

    assert(!isnan(pe_n));
    if ((pe_n < pe_o) or (rand_no < exp(-(pe_n - pe_o)/RT))){ // Accept
        setTVector(someState);
        setOldPE(pe_n);
    }else{ // Reject
        assignConfFromTVector(someState);
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







