/* -------------------------------------------------------------------------- *
 *                 SimTK Molmodel Example: Simple Protein                     *
 * -------------------------------------------------------------------------- *
 * This is the first example from the Molmodel User's Guide. It creates a     *
 * small protein (a five-residue peptide), simulates it and generates a live  *
 * animation while it is running.                                             *
 *                                                                            *
 * Authors: Christopher Bruns, Michael Sherman                                *
 * -------------------------------------------------------------------------- */

#include "Molmodel.h"
#include <iostream>
#include <exception>
using namespace SimTK;

int main() {
try {
    
    // MOLMODEL CONSTRUCTORS DON'T WORK IN Realease mode:
    CompoundSystem system; // Extends MolecularMechanicsSystem (Molmodel)
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

    // Build some compound here
    MagnesiumIon m1;
    MagnesiumIon::setAmberLikeParameters(forceField);
    ChlorideIon c1;
    ChlorideIon::setAmberLikeParameters(forceField);
    Vec3 v1(0.73, -4.62, 0.001);

    system.adoptCompound(c1, v1);

    system.modelCompounds(); 
    //system.addEventReporter(new Visualizer::Reporter(system, 0.020));
    //system.addEventHandler(new VelocityRescalingThermostat(	   system,  293.15, 0.1));
    State state = system.realizeTopology();
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);
    // =====================

    ts.initialize(state);
    ts.stepTo(0.06); // 0.06ps

    return 0;
} 
catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
catch(...) {
    std::cerr << "ERROR: An unknown exception was raised" 
              << std::endl;
    return 1;
}

}


