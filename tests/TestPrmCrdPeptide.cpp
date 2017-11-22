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
#include "SetupReader.hpp"

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

    // Initialize setup reader
    SetupReader setupReader(argv[1]);
    setupReader.dump();

    // Get input filenames
    std::vector<std::string> argValues;
    std::vector<std::string>::iterator argValuesIt;
    argValues = setupReader.getValues("MOLECULES");
    //for(argValuesIt = argValues.begin(); argValuesIt != argValues.end(); ++argValuesIt){
        std::string prmtopFN = argValues[0] + std::string("/ligand.prmtop");
    //    std::string frcmodFN = argValuesIt + std::string("/ligand.frcmod");
        std::string inpcrdFN = argValues[0] + std::string("/ligand.inpcrd");
        std::string rbFN = argValues[0] + std::string("/ligand.rb");
        std::string flexFN = argValues[0] + std::string("/ligand.flex");
    //    std::cout << "prmtopFN " << prmtopFN << std::endl<<std::flush;
    //    std::cout << "frcmodFN " << frcmodFN << std::endl<<std::flush;
    //    std::cout << "inpcrdFN " << inpcrdFN << std::endl<<std::flush;
    //    std::cout << "rbFN "     << rbFN     << std::endl << std::flush;
    //    std::cout << "flexFN "   << flexFN   << std::endl << std::flush;
    //}

    // Get simulation type:
    // IC: Internal Coordinates Dynamics
    // TD: Torsional Dynamics
    // RR: Rigid Rings Torsional Dynamics
    // RB: Rigid Bodies
    std::string ictd = setupReader.getValues("REGIMEN")[0];
    std::cout<<"ictd "<<ictd<<std::endl<<std::flush;

    // Read Amber prmtop and inpcrd
    readAmberInput *amberReader = new readAmberInput();
    amberReader->readAmberFiles(inpcrdFN, prmtopFN);

    natoms = amberReader->getNumberAtoms();
    std::cout << "natoms " << natoms << std::endl << std::flush;

    // Build Gmolmodel simulation world

    World *p_world = new World(std::stod(setupReader.getValues("FREE_TIMESTEP")[0])); 

    // Seed the random number generator 

    srand (time(NULL));

    p_world->AddMolecule(amberReader, rbFN, flexFN, ictd);
    p_world->InitSimulation(amberReader, rbFN, flexFN, ictd);

    // Initialize sampler
    HamiltonianMonteCarloSampler *p_HMCsampler = new HamiltonianMonteCarloSampler(p_world->system, p_world->matter, p_world->lig1, p_world->forceField, p_world->forces, p_world->ts);
    Context *context = new Context(p_world, p_HMCsampler);
    World *world = context->getWorld();
    world->forceField->setTracing(true);

    

   // Simulate
    std::cout << std::fixed;
    std::cout << std::setprecision(4);
    const SimTK::State& constRefState = world->integ->getState();
    SimTK::State& integAdvancedState = world->integ->updAdvancedState();
    p_HMCsampler->initialize( integAdvancedState, 
                              std::stod(setupReader.getValues("FREE_TIMESTEP")[0]),
                              std::stoi(setupReader.getValues("STEPS")[0]),
                              SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );
    //world->forceField->dump();
    //world->forceField->dumpCForceFieldParameters(std::cout);
    //world->system->realize(integAdvancedState, SimTK::Stage::Dynamics);
    std::cout << "Initial const state PE: " << std::setprecision(20)
        << world->forces->getMultibodySystem().calcPotentialEnergy(constRefState)
        << " integ advanced state PE: "
        << world->forces->getMultibodySystem().calcPotentialEnergy(integAdvancedState) 
        << " sampler temperature: "
        << p_HMCsampler->getTemperature()
        << std::endl;

    
    for(int i = 0; i < std::stoi(setupReader.getValues("STEPS")[0]); i++){
        // -- UPDATE --
        //std::cout << "=========================================" << std::endl;
        //std::cout << "Q before update integAdvancedState " 
        //          << integAdvancedState.getQ() << std::endl;
        //std::cout << "U before update integAdvancedState " 
        //          << integAdvancedState.getU() << std::endl;
        //std::cout << "Time before update: " << world->ts->getTime() << std::endl;

        //p_HMCsampler->update((world->ts->updIntegrator()).updAdvancedState());
        p_HMCsampler->update(integAdvancedState, 
            std::stod(setupReader.getValues("FREE_TIMESTEP")[0]),
            std::stoi(setupReader.getValues("FREE_MDSTEPS")[0]) );

        //std::cout << "Q after update integAdvancedState " 
        //          << integAdvancedState.getQ() << std::endl;
        //std::cout << "U after update integAdvancedState " 
        //          << integAdvancedState.getU() << std::endl;
        //std::cout << "Time after update: " << world->ts->getTime()  
        //          << "; integAdvancedState Stage after p_HMCsampler: " 
        //          << (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
        //          << "; integAdvancedState Stage before stepping: " 
        //          << (((SimTK::Subsystem *)(world->matter))->getStage(integAdvancedState)).getName() 
        //          << std::endl;
        writePdb(*((SimTK::Compound *)(world->lig1)), integAdvancedState, "pdbs", "sb_", 8, "HMC1s", i);
        //writePdb(*((SimTK::Compound *)(world->lig2)), integAdvancedState, "pdbs", "sb_", 8, "HMC2s", i);
    }

}



