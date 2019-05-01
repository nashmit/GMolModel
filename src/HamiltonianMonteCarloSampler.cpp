/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "HamiltonianMonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

//** Constructor **/
HamiltonianMonteCarloSampler::HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem
                                     ,SimTK::SimbodyMatterSubsystem *argMatter
                                     ,SimTK::Compound *argResidue
                                     ,SimTK::DuMMForceFieldSubsystem *argDumm
                                     ,SimTK::GeneralForceSubsystem *argForces
                                     ,SimTK::TimeStepper *argTimeStepper
                                     )
    : MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->useFixman = false;  
    this->fix_n = this->fix_o = 0.0;
    this->residualEmbeddedPotential = 0.0;
    nofSamples = 0;
    //printBuffIx = 0;
    this->alwaysAccept = false;
    this->timestep = 0.002; // ps
    this->temperature = 300.0;
    this->boostT = this->temperature;
    this->reproducible = false;
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

/** Get/set the TimeStepper that manages the integrator **/
const SimTK::TimeStepper * HamiltonianMonteCarloSampler::getTimeStepper(void)
{
    return timeStepper;
}

SimTK::TimeStepper * HamiltonianMonteCarloSampler::updTimeStepper(void)
{
    return timeStepper;
}

void HamiltonianMonteCarloSampler::setTimeStepper(SimTK::TimeStepper * someTimeStepper)
{
    timeStepper = someTimeStepper;
}

/** Seed the random number generator. Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
//r void HamiltonianMonteCarloSampler::initialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman) 
void HamiltonianMonteCarloSampler::initialize(SimTK::State& someState, bool randomizeConformation ) //SimTK::Real argTemperature, bool argUseFixman) 
{
    // Seed the random number generator
    if(reproducible){
        //setSeed(1);
        randomEngine.seed( getSeed() );
        std::cout << "SEED: " << getSeed() << std::endl;
    }else{
        setSeed( std::time(0) );
        randomEngine.seed( getSeed() );
        std::cout << "SEED: " << getSeed() << std::endl;
    }

    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    timeStepper->initialize(compoundSystem->getDefaultState());

    // Set the simulation temperature
//r    setTemperature(argTemperature); // Needed for Fixman

    int nu = someState.getNU();

    // Randomize configuration
    //if(randomizeConformation == true){
    //    system->realize(someState, SimTK::Stage::Position);
    //    int nq = someState.getNQ();
    //    SimTK::Vector QV(nq);
    //    for (int j=7; j < nq; ++j){
    //        QV[j] = uniformRealDistribution_mpi_pi(randomEngine);
    //    }
    //    someState.updQ() = QV;
    //}
    //

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
//r    this->useFixman = argUseFixman;
    if(useFixman){
        std::cout << "Hamiltonian Monte Carlo sampler: using Fixman potential." << std::endl;

        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
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
//r void HamiltonianMonteCarloSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) 
void HamiltonianMonteCarloSampler::reinitialize(SimTK::State& someState/*, SimTK::Real argTemperature*/) 
{
    if(reproducible){
        //randomEngine.seed( nofSamples ); // TODO change to seed + nofSamples
        randomEngine.seed( getSeed() + nofSamples ); // TODO change to seed + nofSamples
    }

    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Set the simulation temperature
//r    setTemperature(argTemperature); // Needed for Fixman

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

//std::cout << "reinitialize "
//<< dumm->CalcFullPotEnergyIncludingRigidBodies(someState) << std::endl;

    // Store potential energies
    //setOldPE(getPEFromEvaluator(someState));
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

/** Get/Set the timestep for integration **/
float HamiltonianMonteCarloSampler::getTimestep(void)
{
    return timestep;
    //return timeStepper->updIntegrator().getPredictedNextStepSize();
}

void HamiltonianMonteCarloSampler::setTimestep(float argTimestep)
{
    timeStepper->updIntegrator().setFixedStepSize(timestep);
    this->timestep = argTimestep;
}

/** Initialize the same velocities **/
bool HamiltonianMonteCarloSampler::getReproducible(void)
{
    return reproducible;
}

void HamiltonianMonteCarloSampler::setReproducible(void)
{
    this->reproducible = true;
}

/** Get/Set boost temperature **/
SimTK::Real HamiltonianMonteCarloSampler::getBoostTemperature(void)
{
    return this->boostT;
}

void HamiltonianMonteCarloSampler::setBoostTemperature(SimTK::Real argT)
{
    this->boostT = argT;
}

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. TODO: break in two functions:
initializeVelocities and propagate/integrate **/
void HamiltonianMonteCarloSampler::propose(SimTK::State& someState, int nosteps)
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


    // BEGIN Temperature boost
    // TODO: Implement this in a different function
    //SimTK::Real temperatureBoost = 1.000; // no temperature boost
    //    SimTK::Real temperatureBoost = 1.826; // sqrt(1000/300) : brings temperature from 300 to 2000
    //    //SimTK::Real temperatureBoost = 2.582; // sqrt(2000/300) : brings temperature from 300 to 2000
    //    //SimTK::Real temperatureBoost = 3.1623; // sqrt(3000/300) : brings temperature from 300 to 3000
    // Set velocities according to the boost
    SqrtMInvV *= (sqrtRT); // Set stddev according to temperature
    // END temperature boost

    // Raise the temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);

    printf("us ");
    for(unsigned int i = 0; i < nu; i++) {
        printf("%.0f ", someState.getU()[i]);
    }
    printf("\n");

    // Store the proposed energies
    setProposedKE(matter->calcKineticEnergy(someState));
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman();


    ////////////// Verify equipartition theorem using SOAa////////////
    /*
    SimTK::Vector ThetaDot(nu);
    SimTK::Vector MThetaDot(nu);
    ThetaDot = SqrtMInvV;
    matter->multiplyByM(someState, ThetaDot, MThetaDot);

    // Print velocities and momenta
    //std::cout << "ThetaDot " << ThetaDot << std::endl;
    //std::cout << "MThetaDot " << MThetaDot << std::endl;
    // Print generated random numbers
    //for (int i=0; i < nu; ++i){
    //    std::cout << V[i] << ' ';
    //}
    //std::cout << std::endl;

    // Print pv matrix
    for(int i = 0; i < nu; i++){
        for(int j = 0; j < nu; j++){
            std::cout << ThetaDot[i] * MThetaDot[j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Print kinetic energy
    //SimTK::Real KE2 = ThetaDot.transpose() * MThetaDot;
    //std::cout << "ThetaDot * MThetaDot / 2 " << 0.5 * KE2 << std::endl;
    //SimTK::Real ThetaDotMThetaDot = SimTK::dot(ThetaDot, MThetaDot);
    //SimTK::Real ThetaDotMThetaDot = ThetaDot * MThetaDot;
    //std::cout << "ThetaDotMThetaDot " << ThetaDotMThetaDot << " RT " << RT << std::endl;
    */
    ///////////////


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
    FN = prefix + std::to_string(this->nofSamples) + std::string("before.pdb");
    std::cout << "Writing file " << FN << std::endl;
    std::filebuf fb;
    fb.open (FN, std::ios::out);
    std::ostream os(&fb);
    ((SimTK::Compound *)residue)->writePdb(someState, os);
    fb.close();
    //
    */

    // START Guidance Hamiltonian with a boost temperature

/* DEBUG BEGIN
std::cout << "before "
    << dumm->CalcFullPotEnergyIncludingRigidBodies(someState) << ' '
    << matter->calcKineticEnergy(someState) << ' '
    << calcFixman(someState)
    << std::endl;
// DEBUG END */

/* BOOST BEGIN
    boostFactor = std::sqrt(boostT / this->getTemperature());
    //std::cout << "Boosting by " << boostFactor << " (" << boostT << " / " << this->getTemperature() << ")" << std::endl;
    someState.updU() *= boostFactor;
    system->realize(someState, SimTK::Stage::Velocity);
// BOOST END */


/* DEBUG BEGIN
for(int k = 0; k < nosteps; k++){
    this->timeStepper->stepTo(someState.getTime() + (timestep));

    std::cout << "guide1 "
        << dumm->CalcFullPotEnergyIncludingRigidBodies(someState) << ' '
        << matter->calcKineticEnergy(someState) << ' '
        << calcFixman(someState)
        << " ";
        //<< std::endl;

// Alanine dipeptide DIHEDRAL 4 5 7 13 5 7 13 14
    int a1, a2, a3, a4;
    a1 = 4;
    a2 = 5;
    a3 = 7;
    a4 = 13;
    SimTK::Vec3 a1pos, a2pos, a3pos, a4pos;
    a1pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
    a3pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
    a4pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));
    std::cout << bDihedral(a1pos, a2pos, a3pos, a4pos)
        << " ";
    a1 = 5;
    a2 = 7;
    a3 = 13;
    a4 = 14;
    a1pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
    a3pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
    a4pos = residue->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));
    std::cout << bDihedral(a1pos, a2pos, a3pos, a4pos)
        << std::endl;

}
// DEBUG END */

    // Integrate (propagate trajectory)
    std::cout << "time " << someState.getTime() << " + "
        << timestep << " * " << nosteps << " ";

    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));

    std::cout << someState.getTime() << std::endl;
/*    // Configuration
    std::cout << "HMC conf: ";
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        std::cout << " mobod " << int(mbx) << " = ";
        std::cout << "Q " << mobod.getQAsVector(someState) << std::endl ;
        std::cout << " P_X_F " << mobod.getInboardFrame(someState) << " ";
        std::cout << " F_X_M " << mobod.getMobilizerTransform(someState) << " ";
        std::cout << " B_X_M " << mobod.getOutboardFrame(someState) << " ";
        //std::cout << " P_X_F * F_X_M " << mobod.getInboardFrame(someState) * mobod.getMobilizerTransform(someState) << " ";
        std::cout << "; ";
    }
    std::cout << std::endl;
*/

/* DEBOOST BEGIN
    // RESET Guidance Hamiltonian with a boost temperature
    someState.updU() *= (1.0 / boostFactor);
    system->realize(someState, SimTK::Stage::Velocity);
// DEBOOST END */

/* DEBUG BEGIN
for(int k = 0; k < (nosteps * 5); k++){
    this->timeStepper->stepTo(someState.getTime() + (timestep));

    std::cout << "guide2 "
        << dumm->CalcFullPotEnergyIncludingRigidBodies(someState) << ' '
        << matter->calcKineticEnergy(someState) << ' '
        << forces->getMultibodySystem().calcEnergy(someState) << ' '
        << calcFixman(someState)
        << std::endl;
}
// DEBUG END */

/* DEBUG BEGIN
std::cout << "after "
    << dumm->CalcFullPotEnergyIncludingRigidBodies(someState) << ' '
    << matter->calcKineticEnergy(someState) << ' '
    << calcFixman(someState)
    << std::endl;

    // END Guidance */

    // TODEL
////    std::cout << "Qs and Us after stepTo:" << std::endl;
////    PrintBigMat(someState.getQ(), someState.getNQ(), 3, "Q");
////    PrintBigMat(someState.getU(), someState.getNU(), 3, "U");
    // END TODEL


}

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
void HamiltonianMonteCarloSampler::update(SimTK::State& someState, int nosteps)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    // Do a trial move
    propose(someState, nosteps);

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

    //SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM
    //std::cout << "Multibody PE " << getPEFromEvaluator(someState) << ' ' ; // OPENMM
    pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL

    system->realize(someState, SimTK::Stage::Velocity); // TODO: before computing pe_n
    ke_n = matter->calcKineticEnergy(someState);

    SimTK::Real etot_proposed, etot_n;
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n;
        etot_proposed = pe_o + ke_proposed + fix_o;
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }

    //etot_proposed;
    //etot_n;


    std::cout<<std::setprecision(5)<<std::fixed; //p
    std::cout << "pe_o " << pe_o << " pe_n " << pe_n << " pe_nB " << getPEFromEvaluator(someState) << " ke_prop " << ke_proposed << " ke_n " << ke_n
        << " fix_o " << fix_o << " fix_n " << fix_n << " "
        << " rand_no " << rand_no << " RT " << RT << " exp(-(etot_n - etot_proposed) " << exp(-(etot_n - etot_proposed) / RT)
        << " etot_n " << etot_n  << " etot_proposed " << etot_proposed
        << std::endl;

//     std::cout << std::setprecision(10) << std::fixed << fix_n << ' ';

    // Apply Metropolis criterion
    int accepted = 0;
    if ( getThermostat() == ANDERSEN ){ // MD with Andersen thermostat
        accepted = 1;
        setSetTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setSetPE(pe_n);
        setSetFixman(fix_n);
        setLastAcceptedKE(ke_n);
        this->etot_set = getSetPE() + getSetFixman() + getProposedKE(); // TODO
        ++acceptedSteps;
    }
    else if( (!std::isnan(pe_n)) && 
    ((etot_n < etot_proposed) || (rand_no < exp(-(etot_n - etot_proposed)/RT))) ){ // Correct Acceptance-Rejection 
//    ((etot_n > etot_proposed) || (rand_no < exp((etot_n - etot_proposed)/RT))) ){ // Unfold trial
    //((fix_n < fix_o) || (rand_no < exp(-(fix_n - fix_o)/RT))) ){ // Fixman Monte Carlo
        accepted = 1;
        setSetTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setSetPE(pe_n);
        setSetFixman(fix_n);
        setLastAcceptedKE(ke_n);
        this->etot_set = getSetPE() + getSetFixman() + getProposedKE(); // TODO
        ++acceptedSteps;
    }else{ // Reject
        accepted = 0;
        assignConfFromSetTVector(someState);
    }

    //std::cout << " pe_os " << getSetPE() << " ke_os " << getLastAcceptedKE() << " fix_os " << getSetFixman() //p
    //xstd::cout << " 0 ";
    //xfor (SimTK::MobilizedBodyIndex mbx(2); mbx < matter->getNumBodies(); ++mbx){
    //x    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    //x    const SimTK::MobilizedBody *p_mobod = &mobod;
    //x    std::cout<< std::setprecision(10) << std::fixed;
    //x    std::cout << ((SimTK::MobilizedBody::Pin *)(p_mobod))->getAngle(someState)  << " " ;
    //x}
    //xstd::cout << getSetFixman()  << " " << calcNumFixman(someState) << " " ;
//p    std::cout << accepted << ' ' << getPEFromEvaluator(someState) << ' ' << getLastAcceptedKE() << ' ' << getSetFixman()  << ' ' ;

//                    std::cout << bDihedral( (argResidue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(10)),
//                                            (argResidue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(0)),
//                                            (argResidue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(3)),
//                                            (argResidue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(6)) )  << " ";


        //<< " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n << " rand " << rand_no
        //<< std:: endl; //p
/* RESTORE p
    std::cout << someState.getNU() << ' ' << accepted << ' ' 
        //<< getSetPE() + getREP() << ' ' << getLastAcceptedKE() 
        << pe_o << ' '<< getSetPE() << ' ' << getLastAcceptedKE() 
        << ' ' << getSetFixman() << ' ' << fix_o << ' ' << fix_n << ' ';
*/
    // Keep track of how many MC trials have been done 
    ++nofSamples;

    // Write a check file
    /*
    //std::string prefix;
    //prefix = std::string("pdbs/HMC");
    //std::string FN;
    FN = prefix + std::to_string(this->nofSamples) + std::string("after.pdb");
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
    //for (SimTK::MobilizedBodyIndex mbx(2); mbx < matter->getNumBodies(); ++mbx){
    //    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    //    const SimTK::MobilizedBody *p_mobod = &mobod;
        //std::cout << SimTK::DuMM::Rad2Deg * ((SimTK::MobilizedBody::Pin *)(p_mobod))->getAngle(someState)  << " " ;
    //    std::cout << ((SimTK::MobilizedBody::Pin *)(p_mobod))->getAngle(someState)  << " " ;
        //std::cout << mobod.getAngle(someState) << " " ;
        //std::cout << "mobod " << int(mbx) << mobod.getQAsVector(someState) << "; ";
    //}
    //std::cout << std::endl ;

    //std::cout << "Number of times the force field was evaluated: " << dumm->getForceEvaluationCount() << std::endl;

}

/** Modifies Q randomly
 **/
void HamiltonianMonteCarloSampler::perturbQ(SimTK::State& someState)
{
    // Perturb Q
    //SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    //SimTK::Real rand_no = uniformRealDistribution_mpi_pi(randomEngine);
    int nq = someState.getNQ();
    //SimTK::Vector V(nq);
    for (int i=7; i < nq; ++i){
        //V[i] = uniformRealDistribution_mpi_pi(randomEngine);
        someState.updQ()[i] = uniformRealDistribution_mpi_pi(randomEngine);
    }
    std::cout << "perturbQ " << someState.getQ() << std::endl;
    system->realize(someState, SimTK::Stage::Position);

    // Get needed energies
    SimTK::Real pe_o  = getOldPE();
    if(useFixman){
        SimTK::Real fix_o = getOldFixman();
    }
    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 0.0;
    }

    //SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM
    //std::cout << "Multibody PE " << getPEFromEvaluator(someState) << std::endl; // OPENMM
    pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL

    int accepted = 0;

    accepted = 1;
    setSetTVector(someState);
    setSetPE(pe_n);
    setSetFixman(fix_n);
    ++acceptedSteps;
    assignConfFromSetTVector(someState);

    std::cout << someState.getNU() << ' ' << 1 << ' ' 
      //<< getSetPE() + getREP() << ' ' << getLastAcceptedKE() 
      << getSetPE() << ' ' << 0 
      << ' ' << getSetFixman() << ' ' << fix_o << ' ' << fix_n << ' ';

    // Keep track of how many MC trials have been done 
    ++nofSamples;
}






