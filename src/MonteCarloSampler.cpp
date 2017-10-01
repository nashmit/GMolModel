/**@file
Implementation of MonteCarloSampler class. **/

#include "MonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

MonteCarloSampler::MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     //Topology *argResidue,
                                     SimTK::Compound *argResidue,
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

// Compute Fixman potential

SimTK::Real MonteCarloSampler::calcFixman(SimTK::State& someState)
{
    int nu = someState.getNU();
    SimTK::Real detM = 1.0;
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Matrix D0(6, 6);
    Eigen::MatrixXd EiD0(6, 6);
    matter->calcDetM(someState, V, DetV, D0);
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

// Get set for stored potential energy
SimTK::Real MonteCarloSampler::getOldPE(void)
{
    return this->pe_o;
}

//
void MonteCarloSampler::setOldPE(SimTK::Real argPE)
{
    this->pe_o = argPE;
}



// Set/get Fixman potential

void MonteCarloSampler::setOldFixman(SimTK::State& someState)
{
    this->fix_o = calcFixman(someState);
}

SimTK::Real MonteCarloSampler::getOldFixman(SimTK::State& someState)
{
    return this->fix_o;
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

// Return the old configuration
SimTK::Transform * MonteCarloSampler::getTVector(void)
{
    return this->TVector;
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
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
 
void MonteCarloSampler::propose(SimTK::State& someState)
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

    propose(someState);

    // Send configuration to evaluator  

    sendConfToEvaluator(); // OPENMM

    // Get current potential energy from evaluator

    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM

    // Apply Metropolis criterion

    assert(!isnan(pe_n));
    if ((pe_n < pe_o) or (rand_no < exp(-(pe_n - pe_o)/RT))){ // Accept
        setTVector(someState);
        setOldPE(pe_n);
    }else{ // Reject
        assignConfFromTVector(someState);
    }
}

// Get stored kinetic energy - only necessary for Hamiltonian Monte Carlo

//SimTK::Real MonteCarloSampler::getOldKE(void){SimTK::Real(!assert("Not implemented"));}

// Store the kintetic energy - only necessary for Hamiltonian Monte Carlo

void MonteCarloSampler::setOldKE(SimTK::Real argKE){}

// Get the potential energy from an external source as far as the sampler
// is concerned - OPENMM has to be inserted here

SimTK::Real MonteCarloSampler::getPEFromEvaluator(SimTK::State& someState){
    return matter->calcKineticEnergy(someState);
}

// Get the simulation temperature

SimTK::Real MonteCarloSampler::getTemperature(void){}

// Set the simulatin temperature

void MonteCarloSampler::setTemperature(SimTK::Real argTemperature){temperature = argTemperature;}

// Send configuration to an external evaluator

void MonteCarloSampler::sendConfToEvaluator(void){}







