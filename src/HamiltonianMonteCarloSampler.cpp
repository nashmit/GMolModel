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

// Assign random conformation
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
 
void HamiltonianMonteCarloSampler::propose(SimTK::State& someState)
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
    std::cout << "HamiltonianMonteCarloSampler::propose" << std::endl;
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

void HamiltonianMonteCarloSampler::update(SimTK::State& someState){
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    SimTK::Real RT = getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;

    std::cout << "HamiltonianMonteCarloSampler::update" << std::endl;

    // Get old energy
    SimTK::Real pe_o = getOldPE();

    // Assign random configuration

    propose(someState);

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








