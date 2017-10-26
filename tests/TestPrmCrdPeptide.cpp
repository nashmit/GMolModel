/** @file
 * This file contains a native Gmolmodel test. A 2-butanol molecule
 * is simulated.
 */

#include <string>
#include <iostream>
#include <sstream>
#include "simmain.hpp"

#include "HamiltonianMonteCarloSampler.hpp"
#include "readAmberInput.hpp"

int main(int argc, char **argv)
{
    // Declare a shared memory pointer used by all the classes in the library

    //  Declare variables for simulation parameters

    int natoms;
    int nosteps;
    int ntrials;

    // Set simulation parameters

    nosteps = 10; // RESTORE DEL
    ntrials = 10; // RESTORE DEL
    std::cout<<"main ntrials: "<<ntrials<<std::endl;
    std::cout<<"main nosteps: "<<nosteps<<std::endl;

    // Set input filenames
 
    std::string mol2FN = std::string(argv[1]) + std::string("/ligand.mol2");
    std::string rbFN = std::string(argv[1]) + std::string("/ligand.rb");
    std::string gaffFN = "gaff.dat";
    std::string frcmodFN = std::string(argv[1]) + std::string("/ligand.frcmod");
    std::string flexFN = std::string(argv[1]) + std::string("/ligand.flex");

    // Simulation type:
    // IC: Internal Coordinates Dynamics
    // TD: Torsional Dynamics
    // RR: Rigid Rings Torsional Dynamics
    // RB: Rigid Bodies

    std::string ictd = "IC";

    std::cout<<"mol2FN "<<mol2FN<<std::endl<<std::flush;
    std::cout<<"rbFN "<<rbFN<<std::endl<<std::flush;
    std::cout<<"gaffFN "<<gaffFN<<std::endl<<std::flush;
    std::cout<<"frcmodFN "<<frcmodFN<<std::endl<<std::flush;
    std::cout<<"flexFN "<<flexFN<<std::endl<<std::flush;
    std::cout<<"ictd "<<ictd<<std::endl<<std::flush;

    // Read number of atoms from mol2 file

    //std::string line;
    //std::string column;

    // Read Amber prmtop and inpcrd

    readAmberInput *amberReader = new readAmberInput();
    amberReader->readAmberFiles(std::string(argv[1]) + std::string("/ligand.inpcrd"), 
        std::string(argv[1]) + std::string("/ligand.prmtop"));
    natoms = amberReader->getNumberAtoms();
    std::cout << "natoms " << natoms << std::endl << std::flush;

    // Set the shared memory size (SHMSZ)

    int natoms3 = 3*(natoms);

    // Build Gmolmodel simulation world

    World *p_world = new World(amberReader, rbFN, flexFN, ictd //, PrmToAx_po, MMTkToPrm_po, shm
        ); 

    // Seed the random number generator 

    srand (time(NULL));

    TARGET_TYPE mytimestep;
    mytimestep = 0.0015;

    p_world->InitSimulation(mytimestep, true);

    // Initialize sampler
    HamiltonianMonteCarloSampler *p_HMCsampler = new HamiltonianMonteCarloSampler(p_world->system, p_world->matter, p_world->lig1, p_world->forceField, p_world->forces, p_world->ts);
    Context *context = new Context(p_world, p_HMCsampler);
    World *world = context->getWorld();

    // Options for mass matrix, Lennard Jones

    int pyseed = 0;
    int _massMatNumOpt = 1; // EU
    int _metroFixmanOpt = 1; // EU
    double _lj14sf = 1; //--

    // Alloc memory for saving configurations WATCHOUT TOO BIG

    double **retConfsPois = new double* [ntrials];
    for(int r=0; r<ntrials; r++){
        retConfsPois[r] = new double[3 * world->mr->natoms]; // WATCHOUT
    }
    double *retPotEsPoi = new double[ntrials];
    double *accs = new double;

    *world->pyseed = pyseed;
    world->massMatNumOpt = _massMatNumOpt;
    world->metroFixmanOpt = _metroFixmanOpt;
    world->lj14sf = _lj14sf; //--
    world->sysRetConfsPois = retConfsPois;
    world->sysRetPotEsPoi = retPotEsPoi;
    world->sysAccs = accs;

   // Simulate

    std::cout << std::fixed;
    std::cout << std::setprecision(4);
    const SimTK::State& constRefState = world->integ->getState();
    SimTK::State& integAdvancedState = world->integ->updAdvancedState();
    p_HMCsampler->initialize(integAdvancedState, atof(argv[3]), atoi(argv[4]));

    //world->forceField->dump();
    //world->forceField->dumpCForceFieldParameters(std::cout);

    world->system->realize(integAdvancedState, SimTK::Stage::Dynamics);

    SimTK::Real myPE = world->forces->getMultibodySystem().calcPotentialEnergy(constRefState);
    std::cout << "Initial const state PE: " << std::setprecision(20)
        << myPE << std::endl;
    
    std::cout << "Initial const state PE: " << std::setprecision(20)
        << world->forces->getMultibodySystem().calcPotentialEnergy(constRefState)
        << " integ advanced state PE: "
        << world->forces->getMultibodySystem().calcPotentialEnergy(integAdvancedState) 
        << std::endl;

    for(int i = 0; i<atoi(argv[2]); i++){
        // -- STEPTO -- 

        //std::cout << "Time before stepping: " << world->ts->getTime()
                  //<< "; integAdvancedState Stage before stepping: " 
                  //<< (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
                  //<< "; integAdvancedState Stage before stepping: " 
                  //<< (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
                  //<< std::endl;

        //world->ts->initialize(integAdvancedState);
        //world->ts->stepTo(timeToReach);
        //world->integ->reinitialize(SimTK::Stage::Instance, true);

        //std::cout << "Time after  stepping: " << world->ts->getTime()
                  //<< "; integAdvancedState Stage after stepping: " 
                  //<< (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
                  //<< "; integAdvancedState Stage before stepping: " 
                  //<< (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
                  //<< std::endl;


        // -- UPDATE --
        std::cout << "=========================================" << std::endl;
        //std::cout << "Q before update integAdvancedState " 
        //          << integAdvancedState.getQ() << std::endl;
        //std::cout << "U before update integAdvancedState " 
        //          << integAdvancedState.getU() << std::endl;
        //std::cout << "Time before update: " << world->ts->getTime() << std::endl;

        SimTK::Real myPE = world->forces->getMultibodySystem().calcPotentialEnergy(integAdvancedState);
        std::cout << "PE: " << std::setprecision(20)
            << myPE << std::endl;
        //p_HMCsampler->update((world->ts->updIntegrator()).updAdvancedState());
        p_HMCsampler->update(integAdvancedState, atof(argv[3]), atoi(argv[4]));

        //std::cout << "Q after update integAdvancedState " 
        //          << integAdvancedState.getQ() << std::endl;
        //std::cout << "U after update integAdvancedState " 
        //          << integAdvancedState.getU() << std::endl;
        std::cout << "Time after update: " << world->ts->getTime()  
                  << "; integAdvancedState Stage after p_HMCsampler: " 
                  << (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
                  << "; integAdvancedState Stage before stepping: " 
                  << (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
                  << std::endl;
        writePdb(*((SimTK::Compound *)(world->lig1)), integAdvancedState, "pdbs", "sb_", 8, "HMCs", i);
    }


}



