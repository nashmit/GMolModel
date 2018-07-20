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
    TRACE("NEW ALLOC\n");
    TVector = new SimTK::Transform[matter->getNumBodies()];
    TRACE("NEW ALLOC\n");
    SetTVector = new SimTK::Transform[matter->getNumBodies()];
    this->residualEmbeddedPotential = 0.0;
    this->alwaysAccept = false;
    acceptedSteps = 0;
}

// Destructor
MonteCarloSampler::~MonteCarloSampler()
{
    delete [] TVector;
    delete [] SetTVector;
}

// Seed the random number generator. Set simulation temperature,
// variables that store the configuration
// and variables that store the energies, both needed for the
// acception-rejection step. Also realize velocities and initialize
// the timestepper.
void MonteCarloSampler::initialize(SimTK::State& someState, SimTK::Real argTemperature, bool argUseFixman) 
{
    // Seed the random number generator
    randomEngine.seed( std::time(0) );

    // Set the simulation temperature
    setTemperature(argTemperature); // Needed for Fixman

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
    this->useFixman = argUseFixman;

    if(useFixman){
        std::cout << "Monte Carlo sampler: using Fixman potential." << std::endl;
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

}

// Same as initialize 
void MonteCarloSampler::reinitialize(SimTK::State& someState, SimTK::Real argTemperature) 
{
    // Set the simulation temperature
    setTemperature(argTemperature); // Needed for Fixman

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }
}

void MonteCarloSampler::useFixmanPotential(void)
{
    useFixman = true;
}


// Return true if use Fixman potential
bool MonteCarloSampler::isUsingFixmanPotential(void)
{
    return useFixman;
}

// Is the sampler always accepting the proposed moves
bool MonteCarloSampler::getAlwaysAccept(void)
{
    return alwaysAccept;
}

// Is the sampler always accepting the proposed moves
void MonteCarloSampler::setAlwaysAccept(bool argAlwaysAccept)
{
    alwaysAccept = argAlwaysAccept;
}


// Compute Fixman potential (should have been calcDetMInv ??)
SimTK::Real MonteCarloSampler::calcFixman(SimTK::State& someState){
    int nu = someState.getNU();
    SimTK::Vector V(nu);

    //for (int i=0; i < nu; ++i){
    //   V[i] = i;
    //}
    system->realize(someState, SimTK::Stage::Position);
    matter->realizeArticulatedBodyInertias(someState); // Move in calcDetM ?

    // Get M
    //SimTK::Matrix M(nu, nu);
    //matter->calcM(someState, M);

    // Get detM
    SimTK::Real detM = 1.0;
    SimTK::Vector DetV(nu);
    TRACE("NEW ALLOC\n");
    SimTK::Real* D0 = new SimTK::Real(1.0);

    // TODO: remove the request for Dynamics stage cache in SImbody files
    //std::cout << "MonteCarloSampler::calcFixman Stage: "<< matter->getStage(someState) << std::endl;
    matter->calcDetM(someState, V, DetV, D0);

    //std::cout << "FixmanTorque: " << "MonteCarloSampler::calcFixman logdetM: " << std::setprecision(10) << std::log(*D0) << std::setprecision(2) << std::endl;
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
    SimTK::Real result = 0.5 * RT * std::log(*D0);
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
    SimTK::Real result = 0.5 * RT * std::log(EiDetM);
    return result;
}

// Compute mass matrix determinant numerically
// Not to be confused with the Fixman potential
SimTK::Real MonteCarloSampler::calcNumDetM(SimTK::State& someState){
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
    return EiM.determinant();
}

// Get the set potential energy
SimTK::Real MonteCarloSampler::getSetPE(void)
{
    return this->pe_set;
}

// Get the stored potential energy
SimTK::Real MonteCarloSampler::getOldPE(void)
{
    return this->pe_o;
}

// Set set potential energy
void MonteCarloSampler::setSetPE(SimTK::Real argPE)
{
    this->pe_set = argPE;
}

// Set stored potential energy
void MonteCarloSampler::setOldPE(SimTK::Real argPE)
{
    this->pe_o = argPE;
}

// Set set Fixman potential
void MonteCarloSampler::setSetFixman(SimTK::Real argFixman)
{
    this->fix_set = argFixman;
}

// Set old Fixman potential
void MonteCarloSampler::setOldFixman(SimTK::Real argFixman)
{
    this->fix_o = argFixman;
}

// Get set Fixman potential
SimTK::Real MonteCarloSampler::getSetFixman(void)
{
    return this->fix_set;
}

// Get Fixman potential
SimTK::Real MonteCarloSampler::getOldFixman(void)
{
    return this->fix_o;
}

// Set/get Residual Embedded Potential
void MonteCarloSampler::setREP(SimTK::Real inp)
{
    this->residualEmbeddedPotential = inp;
}

SimTK::Real MonteCarloSampler::getREP(void)
{
    return this->residualEmbeddedPotential;
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

// Stores the configuration into an internal vector of transforms TVector
void MonteCarloSampler::setTVector(SimTK::Transform *inpTVector)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    TVector[i] = inpTVector[i];
    i++;
  }
}

// Stores the set configuration into an internal vector of transforms TVector
void MonteCarloSampler::setSetTVector(const SimTK::State& someState)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    SetTVector[i] = mobod.getMobilizerTransform(someState);
    i++;
  }
}

// Stores the configuration into an internal vector of transforms TVector
// Get the stored configuration
SimTK::Transform * MonteCarloSampler::getTVector(void)
{
    return this->TVector;
}

// Stores the configuration into an internal vector of transforms TVector
// Get the stored configuration
SimTK::Transform * MonteCarloSampler::getSetTVector(void)
{
    return this->SetTVector;
}

// Restores configuration from the internal set vector of transforms TVector
void MonteCarloSampler::assignConfFromSetTVector(SimTK::State& someState)
{
  //someState.invalidateAll(SimTK::Stage::Instance);
  //matter->invalidateArticulatedBodyInertias(someState);
  //system->realize(someState, SimTK::Stage::Position);
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    system->realize(someState, SimTK::Stage::Position);
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(someState, SetTVector[i]);

    system->realize(someState, SimTK::Stage::Position);
    //matter->realizeArticulatedBodyInertias(someState);

    //std::cout <<  "Sampler after setQ State Cache Info: " << std::endl;
    //PrintSimbodyStateCache(someState);

    i++;
  }
  system->realize(someState, SimTK::Stage::Position);
  //matter->realizeArticulatedBodyInertias(someState);
}

// Restores configuration from the internal vector of transforms TVector
void MonteCarloSampler::assignConfFromTVector(SimTK::State& someState)
{
  //someState.invalidateAll(SimTK::Stage::Instance);
  //matter->invalidateArticulatedBodyInertias(someState);
  //system->realize(someState, SimTK::Stage::Position);
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    system->realize(someState, SimTK::Stage::Position);
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(someState, TVector[i]);

    system->realize(someState, SimTK::Stage::Position);
    //matter->realizeArticulatedBodyInertias(someState);
    i++;
  }
  system->realize(someState, SimTK::Stage::Position);
  //matter->realizeArticulatedBodyInertias(someState);
}

// Assign random conformation
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

    //system->realize(someState, SimTK::Stage::Acceleration); // NECESSARY
    system->realize(someState, SimTK::Stage::Position); // NECESSARY

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

    assert(!std::isnan(pe_n));
    if ((pe_n < pe_o) or (rand_no < exp(-(pe_n - pe_o)/RT))){ // Accept
        setTVector(someState);
        setOldPE(pe_n);
        ++acceptedSteps;
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
    std::cout << " MonteCarloSampler::setTemperature " << argTemperature << std::endl;
    this->temperature = argTemperature;
    if(this->temperature < 0){
        std::cerr << "Temperature set to " << this->temperature << std::endl;
    }
    RT = this->temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}

// Set a thermostat
void MonteCarloSampler::setThermostat(ThermostatName argThermostat){
    this->thermostat = argThermostat;
}
// Set a thermostat
void MonteCarloSampler::setThermostat(std::string argThermostat){
    std::string _thermostat;
    _thermostat.resize(argThermostat.size());
    std::transform(argThermostat.begin(), argThermostat.end(),
        _thermostat.begin(), ::tolower);

    try{

        if(_thermostat == "none"){
            this->thermostat = NONE;
        }else if(_thermostat == "andersen"){
            this->thermostat = ANDERSEN;
        }else if(_thermostat == "berendsen"){
            this->thermostat = BERENDSEN;
        }else if(_thermostat == "langevin"){
            this->thermostat = LANGEVIN;
        }else if(_thermostat == "nose_hoover"){
            this->thermostat = NOSE_HOOVER;
        }else{
            throw std::invalid_argument("Thermostat");
        }

    }catch(std::invalid_argument& ia){
        std::cerr << "Invalid argument: " << ia.what() << '\n';
    }

}

// Set a thermostat
void MonteCarloSampler::setThermostat(const char *argThermostat){
    std::string _sthermostat = argThermostat;
    std::string _thermostat;
    _thermostat.resize(_sthermostat.size());
    std::transform(_sthermostat.begin(), _sthermostat.end(),
        _thermostat.begin(), ::tolower);


    try{

        if(_thermostat == "none"){
            this->thermostat = NONE;
        }else if(_thermostat == "andersen"){
            this->thermostat = ANDERSEN;
        }else if(_thermostat == "berendsen"){
            this->thermostat = BERENDSEN;
        }else if(_thermostat == "langevin"){
            this->thermostat = LANGEVIN;
        }else if(_thermostat == "nose_hoover"){
            this->thermostat = NOSE_HOOVER;
        }else{
            throw std::invalid_argument("Thermostat");
        }

    }catch(std::invalid_argument& ia){
        std::cerr << "Invalid argument: " << ia.what() << '\n';
    }

}

// Get the name of the thermostat
ThermostatName MonteCarloSampler::getThermostat(void){
    return this->thermostat;
}




// Send configuration to an external evaluator

void MonteCarloSampler::sendConfToEvaluator(void){
    assert(!"Not implemented");
}

// Get the number of accpted conformations
int MonteCarloSampler::getAcceptedSteps(void)
{
    return acceptedSteps;
}





