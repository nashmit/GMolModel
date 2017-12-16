/**@file
Implementation of MonteCarloSampler class. **/

#include "MonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor
MonteCarloSampler::MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     SimTK::Compound *argResidue,
                                     SimTK::DuMMForceFieldSubsystem *argDumm,
                                     SimTK::GeneralForceSubsystem *argForces,
                                     SimTK::TimeStepper *argTimeStepper)
    : Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    TVector = new SimTK::Transform[matter->getNumBodies()];
}

// Destructor
MonteCarloSampler::~MonteCarloSampler()
{
    delete [] TVector;
}

// Compute Fixman potential
SimTK::Real MonteCarloSampler::calcFixman(SimTK::State& someState){
    int nu = someState.getNU();
    SimTK::Vector V(nu);

    //for (int i=0; i < nu; ++i){
    //   V[i] = i;
    //}
    system->realize(someState, SimTK::Stage::Position);

    // Get M
    //SimTK::Matrix M(nu, nu);
    //matter->calcM(someState, M);

    // Get detM
    SimTK::Real detM = 1.0;
    SimTK::Vector DetV(nu);
    SimTK::Real* D0 = new SimTK::Real(1.0);

    // TODO: remove the request for Dynamics stage cache in SImbody files
    //std::cout << "MonteCarloSampler::calcFixman Stage: "<< matter->getStage(someState) << std::endl;
    matter->calcDetM(someState, V, DetV, D0);

    //std::cout << "MonteCarloSampler::calcFixman D0: " << *D0 << std::endl;
    //std::cout << "MonteCarloSampler::calcFixman RT: " << RT << std::endl;
    // ---- Verify with Eigen ----------
    // Eigen M determinant
    //Eigen::MatrixXd EiM(nu, nu);
    //for(int i=0; i<nu; i++){
    //    for(int j=0; j<nu; j++){
    //        EiM(i, j) = M(i, j);
    //    }
    //}
    //SimTK::Real EiDetM = EiM.determinant();
    //std::cout << "EiDetM= " << EiDetM << std::endl;
    assert(RT > SimTK::TinyReal);
    SimTK::Real result = RT * std::log(*D0);
    delete D0;
    return result;
}

// Compute Fixman potential numerically
SimTK::Real MonteCarloSampler::calcNumFixman(SimTK::State& someState){
    // Get M
    int nu = someState.getNU();
    SimTK::Matrix M(nu, nu);
    matter->calcM(someState, M);

    // Eigen M determinant
    Eigen::MatrixXd EiM(nu, nu);
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiM(i, j) = M(i, j);
        }
    }
    SimTK::Real EiDetM = EiM.determinant();
    assert(RT > SimTK::TinyReal);
    SimTK::Real result = RT * std::log(EiDetM);
    return result;
}

// Get the stored potential energy
SimTK::Real MonteCarloSampler::getOldPE(void)
{
    return this->pe_o;
}

// Set stored potential energy
void MonteCarloSampler::setOldPE(SimTK::Real argPE)
{
    this->pe_o = argPE;
}

// Set Fixman potential
void MonteCarloSampler::setOldFixman(SimTK::Real argFixman)
{
    this->fix_o = argFixman;
}

// Get Fixman potential
SimTK::Real MonteCarloSampler::getOldFixman(void)
{
    return this->fix_o;
}

// Stores the configuration into an internal vector of transforms TVector
void MonteCarloSampler::setTVector(const SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    TVector[i] = mobod.getMobilizerTransform(someState);
    i++;
  }
}

// Get the stored configuration
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

// Get the potential energy from an external source as far as the sampler
// is concerned - OPENMM has to be inserted here
SimTK::Real MonteCarloSampler::getPEFromEvaluator(SimTK::State& someState){
    return forces->getMultibodySystem().calcPotentialEnergy(someState);
}

// Get the desired simulation temperature. Not to be confused with 
// the instant temperature
SimTK::Real MonteCarloSampler::getTemperature(void){
    return this->temperature;
}

// Set the desired simulation temperature. Not to be confused with 
// the instant temperature
void MonteCarloSampler::setTemperature(SimTK::Real argTemperature)
{
    this->temperature = argTemperature;
    RT = this->temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}

// Send configuration to an external evaluator

void MonteCarloSampler::sendConfToEvaluator(void){
    assert(!"Not implemented");
}







