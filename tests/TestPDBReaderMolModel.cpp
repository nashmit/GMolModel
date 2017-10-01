#include "Molmodel.h"

#include <iostream>
#include <exception>

#include "HamiltonianMonteCarloSampler.hpp"
#include "World.hpp"

using namespace SimTK;

int main() {
try {
    // Load the PDB file and construct the system.
    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

    PDBReader pdb("1AKG.pdb");
    pdb.createCompounds(system);
    system.modelCompounds();

    //system.addEventHandler(new VelocityRescalingThermostat(
    //    system, 293.15, 0.1));

    // Show me a movie
    Visualizer viz(system);
    system.addEventReporter( new Visualizer::Reporter(viz, 0.0015) );

    system.realizeTopology();
    
    // Create an initial state for the simulation.
    
    State state = system.getDefaultState();
    pdb.createState(system, state);
    LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);
    
    // Simulate it.

    
    VerletIntegrator integ(system);
    integ.setAccuracy(1e-2);
    TimeStepper ts(system, integ);
    ts.initialize(state);


	// CHANGES
	HamiltonianMonteCarloSampler *p_HMCsampler = new HamiltonianMonteCarloSampler(&system, &matter, &system.updCompound( SimTK::CompoundSystem::CompoundIndex(0) ), &ts);

    SimTK::State& integAdvancedState = integ.updAdvancedState();
    for(int i=0; i<30; i++){
        std::cout << "=========================================" << std::endl;
        std::cout << "Q before update integAdvancedState "
                  << integAdvancedState.getQ() << std::endl;
        std::cout << "U before update integAdvancedState "
                  << integAdvancedState.getU() << std::endl;
        std::cout << "Time before update: " << ts.getTime() << std::endl;

        p_HMCsampler->update(integAdvancedState, 0.0015, 3);

        std::cout << "Q after update integAdvancedState "
                  << integAdvancedState.getQ() << std::endl;
        std::cout << "Time after update: " << ts.getTime()
                  //<< "; integAdvancedState Stage after p_HMCsampler: "
                  //<< (((SimTK::Subsystem)(matter)).getStage(integAdvancedState)).getName()
                  //<< "; integAdvancedState Stage before stepping: "
                  //<< (((SimTK::Subsystem)(matter)).getStage(integAdvancedState)).getName()
                  << std::endl;
        //writePdb(*((SimTK::Compound *)(world->lig1)), integAdvancedState, "pdbs", "sb_", 8, "MCs", i);
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
