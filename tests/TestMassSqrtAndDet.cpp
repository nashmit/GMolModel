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

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

using namespace SimTK;

int main() {
try {
    
    // MOLMODEL CONSTRUCTORS DON'T WORK IN Realease mode:
    CompoundSystem system; // Extends MolecularMechanicsSystem (Molmodel)
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();
    Protein protein("G");
    protein.assignBiotypes();
    //Ethane protein;
    system.adoptCompound(protein);

    for (unsigned int r=0 ; r<protein.getNumBonds(); r++){
        protein.setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
    }

    system.modelCompounds(); 
    system.addEventReporter(new Visualizer::Reporter(system, 0.020));
    //system.addEventHandler(new VelocityRescalingThermostat(	   system,  293.15, 0.1));
    State state = system.realizeTopology();
    LocalEnergyMinimizer::minimizeEnergy(system, state, 15.0);
    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);

    // Get sqrt(MInv)
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    for (int i=0; i < nu; ++i){
       V[i] = i;
    }

    system.realize(state, SimTK::Stage::Position);
    matter.multiplyBySqrtMInv(state, V, SqrtMInvV);
    std::cout << "V= " << V << std::endl;
    std::cout << "SqrtMInvV= " << SqrtMInvV << std::endl;
    //==========


    // Get M
    SimTK::Matrix M(nu, nu);
    matter.calcM(state, M);
    //==========


    // Get detM
    SimTK::Real detM = 1.0;
    SimTK::Vector DetV(nu);
    //SimTK::Matrix D0(6, 6);
    Real* newDetM = new Real(1.0);

    matter.calcDetM(state, V, DetV, newDetM);
    std::cout << "newDetM: " << *newDetM << std::endl;
    //==========

    // ---- Verify with Eigen ----------
    // Eigen M determinant
    Eigen::MatrixXd EiM(nu, nu);
    Eigen::MatrixXd EiD0(6, 6);
    SimTK::Real JainDetM;
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiM(i, j) = M(i, j);
        }
    }
    SimTK::Real EiDetM = EiM.determinant();
    std::cout << "EiDetM= " << EiDetM << std::endl;

    // Simulate for a bit
    ts.initialize(state);
    ts.stepTo(0.12); // 0.12ps

//    // Check Fixman torque
//    unsigned int nq = state.getNQ();
//    SimTK::Vector Qs = state.getQ();
//    SimTK::Vector copyQs(nq);
//    for(int i = 0; i < nq; i++){copyQs[i] = Qs[i];} 
//    std::cout << "Qs: " << Qs << std::endl;
//    SimTK::Real prevQ = 0, prevFixmanPotential = 0, FixmanPotential = 0;
//    SimTK::Real dQ = 0, dFixman = 0, deriv = 0;
//    SimTK::Real TinyIncrement = 0.0001;
//
//    SimTK::Vector V1(nu);
//    SimTK::Vector V2(nu); // This stores the torques.
//    for(int qIx=0; qIx<nq; qIx++){
//        for(int i=0; i<10; i++){
//            matter.calcDetM(state, V, DetV, newDetM);
//            prevFixmanPotential = FixmanPotential;
//            SimTK::Real FixmanPotential = 0.5 * std::log(*newDetM);
//            dFixman = FixmanPotential - prevFixmanPotential;
//
//            prevQ = Qs[qIx];
//            copyQs[qIx] += TinyIncrement;
//            state.setQ(copyQs);
//
//            system.realize(state, SimTK::Stage::Dynamics);
//
//            matter.calcFixmanTorque(state, V1, V2, newDetM);
//            dQ = Qs[qIx] - prevQ;
//            deriv = dFixman / dQ;
//            std::cout << "calcFixmanTorque V2:" << V2 << std::endl;
//            std::cout << "numDeriv:" << deriv << std::endl;
//        }
//    }
//
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


