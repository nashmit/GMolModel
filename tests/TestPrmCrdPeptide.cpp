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
    // Initialize setup reader
    SetupReader setupReader(argv[1]);

    // Build Gmolmodel simulation world
    World *p_world;
    if(setupReader.getValues("VISUAL")[0] == "TRUE"){
        p_world = new World(0, true, std::stod(setupReader.getValues("FREE_TIMESTEP")[0]));
    }else{
        p_world = new World(0, false, std::stod(setupReader.getValues("FREE_TIMESTEP")[0]));
    }

    // Get input filenames
    std::vector<std::string> argValues;
    std::vector<std::string>::iterator argValuesIt;
    argValues = setupReader.getValues("MOLECULES");
    for(argValuesIt = argValues.begin(); argValuesIt != argValues.end(); ++argValuesIt){
        // Get molecule files
        std::string prmtopFN = *argValuesIt + std::string("/ligand.prmtop");
        std::string inpcrdFN = *argValuesIt + std::string("/ligand.inpcrd");
        std::string rbFN = *argValuesIt + std::string("/ligand.rb");
        std::string flexFN = *argValuesIt + std::string("/ligand.flex");

        // Get regimen
        std::string ictd = setupReader.getValues("REGIMEN")[0];

        // Read Amber prmtop and inpcrd
        readAmberInput *amberReader = new readAmberInput();
        amberReader->readAmberFiles(inpcrdFN, prmtopFN);

        // Add molecules
        p_world->AddMolecule(amberReader, rbFN, flexFN, ictd);

        delete amberReader;
    }

    p_world->Init(0.0015);


    // Initialize sampler
    srand (time(NULL));
    HamiltonianMonteCarloSampler *p_HMCsampler = new HamiltonianMonteCarloSampler(p_world->compoundSystem, p_world->matter, (SimTK::Compound *)(&(p_world->getTopology(0))), p_world->forceField, p_world->forces, p_world->ts);

    Context *context = new Context();
    context->AddWorld(p_world, true);
    World *world = context->getWorld();
    world->forceField->setTracing(true);

   // Simulate
    std::cout << std::fixed << std::setprecision(4);
    const SimTK::State& constRefState = world->integ->getState();
    SimTK::State& integAdvancedState = world->integ->updAdvancedState();
    if( setupReader.getValues("REGIMEN")[0] == "IC" ){
        p_HMCsampler->initialize( integAdvancedState, 
           std::stod(setupReader.getValues("FREE_TIMESTEP")[0]),
           std::stoi(setupReader.getValues("STEPS")[0]),
           SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ),
           false );
    }else{
        p_HMCsampler->initialize( integAdvancedState, 
           std::stod(setupReader.getValues("FREE_TIMESTEP")[0]),
           std::stoi(setupReader.getValues("STEPS")[0]),
           SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ),
           true );
    }
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

        p_HMCsampler->update(integAdvancedState, 
            std::stod(setupReader.getValues("FREE_TIMESTEP")[0]),
            std::stoi(setupReader.getValues("FREE_MDSTEPS")[0]) );

        if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
            for(int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                writePdb( ((SimTK::Compound)(world->getTopology(mol_i))),
                    integAdvancedState, 
                    "pdbs", "sb_",
                    8, 
                    (std::string("HMC") + std::to_string(mol_i)).c_str(), i);
            }
        }

    }

}



