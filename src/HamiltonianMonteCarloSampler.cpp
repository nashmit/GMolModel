/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "HamiltonianMonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor
HamiltonianMonteCarloSampler::HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     SimTK::Compound *argResidue,
                                     SimTK::DuMMForceFieldSubsystem *argDumm,
                                     SimTK::GeneralForceSubsystem *argForces,
                                     SimTK::TimeStepper *argTimeStepper)
    : MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->fix_n = this->fix_o = 0.0;
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

// Initialize variables (identical to setTVector)
void HamiltonianMonteCarloSampler::initialize(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
      const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
      const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
      TVector[i] = mobod.getMobilizerTransform(someState);
      i++;
  }

  //timeStepper->initialize(someState);

}


// Assign random conformation
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
void HamiltonianMonteCarloSampler::propose(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    //randomEngine.seed(4294653137UL); // for reproductibility

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
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);
    setOldKE(matter->calcKineticEnergy(someState));

    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));


}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step

void HamiltonianMonteCarloSampler::update(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    SimTK::Real RT = getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;

    // Get old energy
    setOldPE(getPEFromEvaluator(someState));
    SimTK::Real pe_o = getOldPE();
    
    // Get old Fixman potential
    system->realize(someState, SimTK::Stage::Acceleration);
    fix_o = calcFixman(someState);

    // Assign random configuration

    propose(someState, timestep, nosteps);

    //system->realize(someState, SimTK::Stage::Acceleration);
    fix_n = calcFixman(someState);

    // Send configuration to evaluator  

    sendConfToEvaluator(); // OPENMM

    // Get current potential energy from evaluator

    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM

    // Get current kinetic energy from Simbody

    SimTK::Real ke_n = matter->calcKineticEnergy(someState);

    // Apply Metropolis criterion

    assert(!isnan(pe_n));
    SimTK::Real etot_n = pe_n + ke_n + fix_n;
    SimTK::Real etot_o = pe_o + ke_o + fix_o;
    std::cout << "pe_o " << pe_o << " ke_o " << ke_o << " fix_o " << fix_o << std::endl;
    std::cout << "pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n << std::endl;
    std::cout << "rand_no " << rand_no << std::endl;
    if ((etot_n < etot_o) or (rand_no < exp(-(etot_n - etot_o)/RT))){ // Accept
        std::cout << "Move accepted" << std::endl;
        setTVector(someState);
        setOldPE(pe_n);
        setOldFixman(someState);
    }else{ // Reject
        std::cout << "Move not accepted" << std::endl;
        assignConfFromTVector(someState);
    }
}








