/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "HamiltonianMonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

HamiltonianMonteCarloSampler::HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     Topology *argResidue,
                                     SimTK::TimeStepper *argTimeStepper)
    : MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argTimeStepper)
{
}

// Destructor

HamiltonianMonteCarloSampler::~HamiltonianMonteCarloSampler()
{
}

// Set old kinetic energy

void HamiltonianMonteCarloSampler::setOldKE(SimTK::Real inpKE)
{
    this->ke_o = inpKE;
}


// Assign random conformation
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
 
void HamiltonianMonteCarloSampler::propose(SimTK::State& someState, SimTK::Real timestep, int nosteps)
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
    std::cout << "HamiltonianMonteCarloSampler::propose Q: " << someState.getQ() 
        << std::endl << std::flush;
    std::cout << "HamiltonianMonteCarloSampler::propose U: " << someState.getU()
        << std::endl << std::flush;
    std::cout << "HamiltonianMonteCarloSampler::propose trying to step to: "
        << someState.getTime() + (timestep*nosteps) << std::endl;

    // Assign velocities according to Maxwell-Boltzmann distribution
    int nu = someState.getNU();
    double kTb = SimTK_BOLTZMANN_CONSTANT_MD * 300.0;
    double sqrtkTb = std::sqrt(kTb);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }
    system->realize(someState, SimTK::Stage::Position);
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtkTb; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;



    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    system->realize(someState, SimTK::Stage::Acceleration);
    (this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);
    setOldKE(matter->calcKineticEnergy(someState));

    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));
    //timeStepper->stepTo(0.0150);


}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step

void HamiltonianMonteCarloSampler::update(SimTK::State& someState){
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    SimTK::Real RT = getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;

    std::cout << "HamiltonianMonteCarloSampler::update" << std::endl;

    // Get old energy
    SimTK::Real pe_o = getOldPE();

    // Assign random configuration

    propose(someState, 0.0015, 10);

    // Send configuration to evaluator  

    sendConfToEvaluator(); // OPENMM

    // Get current potential energy from evaluator

    SimTK::Real pe_n = getPEFromEvaluator(); // OPENMM

    // Get current kinetic energy from Simbody

    SimTK::Real ke_n = matter->calcKineticEnergy(someState);

    // Apply Metropolis criterion

    assert(!isnan(pe_n));
    if ((pe_n < pe_o) or (rand_no < exp(-(pe_n - pe_o)/RT))){ // Accept
        setTVector(someState);
        setOldPE(pe_n);
    }else{ // Reject
        assignConfFromTVector(someState);
    }
}








