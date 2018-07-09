#include "Molmodel.h"

#include <iostream>
#include <exception>

#include "HamiltonianMonteCarloSampler.hpp"
#include "World.hpp"

using namespace SimTK;

int main(int argc, char **argv) {
try {
    // Load the PDB file and construct the system.
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    SimTK::GeneralForceSubsystem *forces = new SimTK::GeneralForceSubsystem(system);

    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

    PDBReader pdb(argv[1]);
    pdb.createCompounds(system);
    system.modelCompounds();

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, 0.0015) );

    system.realizeTopology();
    
    // Create an initial state for the simulation.
    
    State state = system.getDefaultState();
    pdb.createState(system, state);
    //LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);
    
    // Simulate it.
    VerletIntegrator integ(system);
    integ.setAccuracy(1e-2);
    TimeStepper ts(system, integ);
    ts.initialize(state);

    // Create and initialize sampler
    HamiltonianMonteCarloSampler *p_HMCsampler = new HamiltonianMonteCarloSampler(&system, &matter, &system.updCompound( SimTK::CompoundSystem::CompoundIndex(0) ), &forceField, forces, &ts);

    const SimTK::State& constRefState = ts.getIntegrator().getState();
    SimTK::State& integAdvancedState = ts.updIntegrator().updAdvancedState();
    //p_HMCsampler->initialize(integAdvancedState, atoi(argv[4]), SimTK::Real(atof(argv[5])) );
    p_HMCsampler->initialize(integAdvancedState, SimTK::Real(atof(argv[5])) );

    // Force field scaling
    forceField.setBondStretchGlobalScaleFactor(0.0);
    forceField.setBondBendGlobalScaleFactor(0.0);
    forceField.setBondTorsionGlobalScaleFactor(0.0);
    forceField.setAmberImproperTorsionGlobalScaleFactor(0.0);

    //forceField.setVdw12ScaleFactor(0.0);
    //forceField.setVdw13ScaleFactor(0.0);
    //forceField.setVdw14ScaleFactor(0.0);
    //forceField.setVdw15ScaleFactor(0.0);
    //forceField.setVdwGlobalScaleFactor(0.0);

    forceField.setCoulomb12ScaleFactor(0.0);
    forceField.setCoulomb13ScaleFactor(0.0);
    forceField.setCoulomb14ScaleFactor(0.0);
    forceField.setCoulomb15ScaleFactor(0.0);
    forceField.setCoulombGlobalScaleFactor(0.0);

    forceField.setGbsaGlobalScaleFactor(0.0);

    forceField.dump();
    //=====================

    std::cout << "Initial const state PE: " << std::setprecision(20)
        << forces->getMultibodySystem().calcPotentialEnergy(constRefState)
        << " integ advanced state PE: "
        << forces->getMultibodySystem().calcPotentialEnergy(integAdvancedState)
        << std::endl;

    for(int i=0; i<atoi(argv[2]); i++){
        std::cout << "=========================================" << std::endl;
        std::cout << "Q before update integAdvancedState "
                  << integAdvancedState.getQ() << std::endl;
        std::cout << "U before update integAdvancedState "
                  << integAdvancedState.getU() << std::endl;
        std::cout << "Time before update: " << ts.getTime() << std::endl;

        p_HMCsampler->update(integAdvancedState, atof(argv[4]));

        std::cout << "Q after update integAdvancedState "
                  << integAdvancedState.getQ() << std::endl;
        std::cout << "Time after update: " << ts.getTime()
                  //<< "; integAdvancedState Stage after p_HMCsampler: "
                  //<< (((SimTK::Subsystem)(matter)).getStage(integAdvancedState)).getName()
                  //<< "; integAdvancedState Stage before stepping: "
                  //<< (((SimTK::Subsystem)(matter)).getStage(integAdvancedState)).getName()
                  << std::endl;
        writePdb(system.updCompound( SimTK::CompoundSystem::CompoundIndex(0) ), integAdvancedState, "pdbs", "sb_", 8, "HMCs", i);

    }

    return 0;

} 
catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
catch(...) {
    std::cerr << "ERROR: An unknown exception was raised" << std::endl;
    return 1;
}

}
