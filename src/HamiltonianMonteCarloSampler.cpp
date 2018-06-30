/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "HamiltonianMonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

//** Constructor **/
HamiltonianMonteCarloSampler::HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     SimTK::Compound *argResidue,
                                     SimTK::DuMMForceFieldSubsystem *argDumm,
                                     SimTK::GeneralForceSubsystem *argForces,
                                     SimTK::TimeStepper *argTimeStepper)
    : MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->useFixman = true;  
    this->fix_n = this->fix_o = 0.0;
    this->residualEmbeddedPotential = 0.0;
    sampleNumber = 0;
    this->alwaysAccept = false;

}

/** Destructor **/
HamiltonianMonteCarloSampler::~HamiltonianMonteCarloSampler()
{
}

/** Calculate sqrt(M) using Eigen. For debug purposes. **/
void HamiltonianMonteCarloSampler::calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper)
{
    assert("!Not implemented");
}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular matrix and it is computed 
multipling a set of orthonormal vectors with the sqrt(MInv). **/
void HamiltonianMonteCarloSampler::calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvU: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = SqrtMInvV[j];
        }
        V[i] = 0;
    }


}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
This is lower triangular matrix and it is computed by multipling a set of
 orthonormal vectors with the sqrt(MInv) and transpose it. **/
void HamiltonianMonteCarloSampler::calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvL: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(j, i) = SqrtMInvV[j];
        }
        V[i] = 0;
    }

}

/** Stores the accepted kinetic energy. This should be set right after a
move is accepted. It's a component of the total energy stored. **/
void HamiltonianMonteCarloSampler::setLastAcceptedKE(SimTK::Real inpKE)
{
    this->ke_lastAccepted = inpKE;
}

/** Sets the proposed kinetic energy before the proposal. This should be
set right after the velocities are initialized. **/
void HamiltonianMonteCarloSampler::setProposedKE(SimTK::Real inpKE)
{
    this->ke_proposed = inpKE;
}

/** Seed the random number generator. Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
void HamiltonianMonteCarloSampler::initialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman) 
{
    // Seed the random number generator
    randomEngine.seed( std::time(0) );

    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    timeStepper->initialize(compoundSystem->getDefaultState());

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
        std::cout << "Hamiltonian Monte Carlo sampler: using Fixman potential." << std::endl;

        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);

    // Store kinetic energies
    setProposedKE(matter->calcKineticEnergy(someState));
    setLastAcceptedKE(getProposedKE());

    // Store total energies
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman();
    this->etot_set = this->etot_proposed;

  
}

/** Same as initialize **/
void HamiltonianMonteCarloSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) 
{
    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

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

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);

    // Store kinetic energies
    setProposedKE(matter->calcKineticEnergy(someState));
    setLastAcceptedKE(getProposedKE());

    // Store total energies
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman();
    this->etot_set = this->etot_proposed;

}

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. TODO: break in two functions:
initializeVelocities and propagate/integrate **/
void HamiltonianMonteCarloSampler::propose(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    // Seed the random number generator every move
    //randomEngine.seed(4294653137UL); // for reproductibility

    // Initialize configuration - not necessary unless we modify the
    // configuration in addition to velocities
    system->realize(someState, SimTK::Stage::Position);
    int t = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        TVector[t] = SetTVector[t];
        t++;
    }

    // TODO: change the names from Old to Proposed and Set to lastAccepted
    setOldPE(getSetPE());
    setOldFixman(getSetFixman());

    // Initialize velocities according to the Maxwell-Boltzmann distribution
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);

    // TODO: Implement this in a different function
    //SimTK::Real temperatureBoost = 2.58; // sqrt(2000/300) : brings temperature from 300 to 1000
    //SqrtMInvV *= sqrtRT * temperatureBoost; // Set stddev according to temperature

    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;

    // TODEL
    // TODO: Implement this in a different function
////
////    // Unlock all mobilizers
////    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
////        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
////        mobod.unlock(someState);
////    }
////
////    // Randomly choose what mobods to lock and store the chosen ones into
////    // a binary array
////    int BinaryArray[matter->getNumBodies()];
////    for(int i = 0; i < matter->getNumBodies(); i++){
////        BinaryArray[i] = 1;
////    }
////    for(int i = 0; i < 3; i++){
////        SimTK::Real rand_no = uniformRealDistribution(randomEngine);
////        BinaryArray[int(std::floor(matter->getNumBodies() * rand_no))] = 0;
////    }
////    
////    // Lock the chosen mobilizers
////    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
////        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
////        if(BinaryArray[int(mbx) - 1] == 0){
////            SimTK::Motion::Level motionLevel = SimTK::Motion::Level::Position;
////            mobod.lock(someState, motionLevel);
////        }
////    }
    // END TODEL

    // Store the proposed kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    setProposedKE(matter->calcKineticEnergy(someState));

    // Store the proposed total energy
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman();

    // TODEL
////    std::cout << "Qs and Us before stepTo:" << std::endl;
////    PrintBigMat(someState.getQ(), someState.getNQ(), 3, "Q");
////    PrintBigMat(someState.getU(), someState.getNU(), 3, "U");
    // END TODEL

    // Write a check file
    /*
    std::string prefix;
    prefix = std::string("pdbs/HMC");
    std::string FN;
    FN = prefix + std::to_string(this->sampleNumber) + std::string("before.pdb");
    std::cout << "Writing file " << FN << std::endl;
    std::filebuf fb;
    fb.open (FN, std::ios::out);
    std::ostream os(&fb);
    ((SimTK::Compound *)residue)->writePdb(someState, os);
    fb.close();
    //
    */

    // Integrate (propagate trajectory)
//pp    std::cout << sampleNumber << ' ';
    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));

    // TODEL
////    std::cout << "Qs and Us after stepTo:" << std::endl;
////    PrintBigMat(someState.getQ(), someState.getNQ(), 3, "Q");
////    PrintBigMat(someState.getU(), someState.getNU(), 3, "U");
    // END TODEL


}

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
void HamiltonianMonteCarloSampler::update(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    // Do a trial move
    propose(someState, timestep, nosteps);

    // Get needed energies
    SimTK::Real pe_o  = getOldPE();
    if(useFixman){
        SimTK::Real fix_o = getOldFixman();
    }
    SimTK::Real ke_proposed  = getProposedKE();

    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 0.0;
    }
    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM

    system->realize(someState, SimTK::Stage::Velocity);
    SimTK::Real ke_n = matter->calcKineticEnergy(someState);

    SimTK::Real etot_proposed, etot_n;
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n;
        etot_proposed = pe_o + ke_proposed + fix_o;
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }

    etot_proposed;
    etot_n;


    std::cout<<std::setprecision(5)<<std::fixed;
//p    std::cout << "pe_o " << pe_o << " ke_o " << ke_proposed << " fix_o " << fix_o << " rep " << getREP()
//p       << " pe_n " << pe_n  << " ke_n " << ke_n << " fix_n " << fix_n
        //<< " rand_no " << rand_no << " RT " << RT << " exp(-(etot_n - etot_proposed) " << exp(-(etot_n - etot_proposed) / RT)
        //<< " etot_n " << etot_n  + getREP() << " etot_proposed " << etot_proposed + getREP()
//p        ;

//pp     std::cout << std::setprecision(10) << std::fixed << fix_n << ' ';

    // Apply Metropolis criterion
    if ( getThermostat() == ANDERSEN ){ // MD with Andersen thermostat
//p        std::cout << " acc 1 " ;
        setSetTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setSetPE(pe_n);
        setSetFixman(fix_n);
        setLastAcceptedKE(ke_n);
        this->etot_set = getSetPE() + getSetFixman() + getProposedKE(); // TODO
    }
    else if( (!isnan(pe_n)) && 
    ((etot_n < etot_proposed) || (rand_no < exp(-(etot_n - etot_proposed)/RT))) ){ // Accept
//p        std::cout << " acc 1 " ;
        setSetTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setSetPE(pe_n);
        setSetFixman(fix_n);
        setLastAcceptedKE(ke_n);
        this->etot_set = getSetPE() + getSetFixman() + getProposedKE(); // TODO
        ++acceptedSteps;
    }else{ // Reject
//p        std::cout << " acc 0 " ;
        assignConfFromSetTVector(someState);
    }

//p    std::cout << " pe_os " << getSetPE() + getREP() << " ke_os " << getLastAcceptedKE() << " fix_os " << getSetFixman()
        //<< " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n
//p        << std:: endl;

    // Keep track of how many MC trials have been done 
    ++sampleNumber;

    // Write a check file
    /*
    //std::string prefix;
    //prefix = std::string("pdbs/HMC");
    //std::string FN;
    FN = prefix + std::to_string(this->sampleNumber) + std::string("after.pdb");
    std::cout << "Writing file " << FN << std::endl;
    //std::filebuf fb;
    fb.open (FN, std::ios::out);
    std::ostream os2(&fb);
    ((SimTK::Compound *)residue)->writePdb(someState, os2);
    fb.close();
    //
    */

    // Configuration
    //std::cout << "HMC conf: ";
    //for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    //    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    //    std::cout << " mobod " << int(mbx) << " = ";
    //    std::cout << mobod.getQAsVector(someState) << std::endl ;
    //    std::cout << " P_X_F " << mobod.getInboardFrame(someState) << " ";
    //    std::cout << " F_X_M " << mobod.getMobilizerTransform(someState) << " ";
    //    //std::cout << " P_X_F * F_X_M " << mobod.getInboardFrame(someState) * mobod.getMobilizerTransform(someState) << " ";
    //    std::cout << "; ";
    //}
    //std::cout << std::endl;

    // Calculate geometric features fast
    //std::cout << someState.getQ() << std::endl ;
    //for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    //    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        //const SimTK::MobilizedBody *p_mobod = &mobod;
        //std::cout << SimTK::DuMM::Rad2Deg * ((SimTK::MobilizedBody::Pin *)(p_mobod))->getAngle(someState)  << " " ;
        //std::cout << ((SimTK::MobilizedBody::Pin *)(p_mobod))->getAngle(someState)  << " " ;
        //std::cout << mobod.getAngle(someState) << " " ;
        //std::cout << "mobod " << int(mbx) << mobod.getQAsVector(someState) << "; ";
    //}
    //std::cout << std::endl ;

    //std::cout << "Number of times the force field was evaluated: " << dumm->getForceEvaluationCount() << std::endl;
  
    


}








