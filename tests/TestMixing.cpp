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
    Context *context = new Context();
    srand (time(NULL));

    // Build Gmolmodel simulation worlds
    World *p_world0;

    if(setupReader.getValues("VISUAL")[0] == "TRUE"){
        p_world0 = new World(0, true, std::stod(setupReader.getValues("FREE_TIMESTEP")[0]));
    }else{
        p_world0 = new World(0, false, std::stod(setupReader.getValues("FREE_TIMESTEP")[0]));
    }

    // Set world identifiers
    p_world0->ownWorldIndex = 0;

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
        p_world0->AddMolecule(amberReader, rbFN, flexFN, ictd);

        delete amberReader;
    }

    p_world0->Init();

    // Initialize sampler
    HamiltonianMonteCarloSampler *p_HMCsampler0 = new HamiltonianMonteCarloSampler(p_world0->compoundSystem, p_world0->matter, (SimTK::Compound *)(&(p_world0->getTopology(0))), p_world0->forceField, p_world0->forces, p_world0->ts);

    // ??
    //p_world0->forceField->setTracing(true);

    // Get default and advanced states from integrators
    const SimTK::State& constRefState0 = p_world0->integ->getState();
    SimTK::State& integAdvancedState0 = p_world0->integ->updAdvancedState();

    // Initialize samplers 
    p_HMCsampler0->initialize( integAdvancedState0, 
       std::stod(setupReader.getValues("FREE_TIMESTEP")[0]),
       std::stoi(setupReader.getValues("STEPS")[0]),
       SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ),
       true );

    // Build Gmolmodel simulation worlds
    World *p_world1;

    if(setupReader.getValues("VISUAL")[0] == "TRUE"){
        p_world1 = new World(1, true, std::stod(setupReader.getValues("CONS_TIMESTEP")[0]));
    }else{
        p_world1 = new World(1, false, std::stod(setupReader.getValues("CONS_TIMESTEP")[0]));
    }

    // Set world identifiers
    p_world1->ownWorldIndex = 1;

    // Get input filenames
    //std::vector<std::string> argValues;
    //std::vector<std::string>::iterator argValuesIt;
    argValues = setupReader.getValues("MOLECULES");
    for(argValuesIt = argValues.begin(); argValuesIt != argValues.end(); ++argValuesIt){
        // Get molecule files
        std::string prmtopFN = *argValuesIt + std::string("/ligand.prmtop");
        std::string inpcrdFN = *argValuesIt + std::string("/ligand.inpcrd");
        std::string rbFN = *argValuesIt + std::string("/ligand.rb");
        std::string flexFN = *argValuesIt + std::string("/ligand.flex");

        // Get regimen
        std::string ictd = setupReader.getValues("REGIMEN")[1];

        // Read Amber prmtop and inpcrd
        readAmberInput *amberReader = new readAmberInput();
        amberReader->readAmberFiles(inpcrdFN, prmtopFN);

        // Add molecules
        p_world1->AddMolecule(amberReader, rbFN, flexFN, ictd);

        delete amberReader;
    }

    p_world1->Init();

    // Initialize sampler
    HamiltonianMonteCarloSampler *p_HMCsampler1 = new HamiltonianMonteCarloSampler(p_world1->compoundSystem, p_world1->matter, (SimTK::Compound *)(&(p_world1->getTopology(0))), p_world1->forceField, p_world1->forces, p_world1->ts);

    // ??
    //p_world1->forceField->setTracing(true);

    // Get default and advanced states from integrators
    const SimTK::State& constRefState1 = p_world1->integ->getState();
    SimTK::State& integAdvancedState1 = p_world1->integ->updAdvancedState();

    // Initialize samplers 
    p_HMCsampler1->initialize( integAdvancedState1, 
       std::stod(setupReader.getValues("CONS_TIMESTEP")[0]),
       std::stoi(setupReader.getValues("STEPS")[0]),
       SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ),
       true );

    ////////////////////////////////////////////////////////// END creating Worlds

    // Add worlds to context
    context->AddWorld(p_world0, p_HMCsampler0);
    context->AddWorld(p_world1, p_HMCsampler1);

    // Set convenient names
    World *world0 = context->getWorld(0);
    World *world1 = context->getWorld(1);

    std::cout << "World 0 Initial const state PE: " << std::setprecision(20)
        << world0->forces->getMultibodySystem().calcPotentialEnergy(constRefState0)<< std::endl;
    std::cout << "World 1 Initial const state PE: " << std::setprecision(20)
        << world1->forces->getMultibodySystem().calcPotentialEnergy(constRefState1)<< std::endl;

    // Simulate - STEPS = MIXING STEPS
    int free_mix_mcsteps = std::stoi(setupReader.getValues("FREE_MIXMCSTEPS")[0]);
    int cons_mix_mcsteps = std::stoi(setupReader.getValues("CONS_MIXMCSTEPS")[0]);
    int total_mcsteps = std::stoi(setupReader.getValues("STEPS")[0]);
    assert( (!(total_mcsteps % free_mix_mcsteps)) &&
            (!(total_mcsteps % cons_mix_mcsteps)) && 
            "Total number of steps must be divisible with free and cons MC steps within a round." );
    int round_mcsteps = free_mix_mcsteps + cons_mix_mcsteps;
    SimTK::Real free_timestep = std::stod(setupReader.getValues("FREE_TIMESTEP")[0]);
    int free_mdsteps = std::stoi(setupReader.getValues("FREE_MDSTEPS")[0]);
    SimTK::Real cons_timestep = std::stod(setupReader.getValues("CONS_TIMESTEP")[0]);
    int cons_mdsteps = std::stoi(setupReader.getValues("CONS_MDSTEPS")[0]);
    int j = -1;
    for(int i = 0; i < total_mcsteps; i += round_mcsteps){
        j += 2;
        //std::cout << "WORLD 0 +++++++++++++++++++" << std::endl;

        if(i > 0){
            integAdvancedState0 = world0->setAtomsLocationsInGround(integAdvancedState0, world1->getAtomsLocationsInGround(integAdvancedState1));
            //// Automatically accept and assign old values
            p_HMCsampler0->reinitialize( integAdvancedState0, 
                free_timestep, total_mcsteps,
                SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );
        }

        for(int k = 0; k < free_mix_mcsteps; k++){
            p_HMCsampler0->update(integAdvancedState0, free_timestep, free_mdsteps);
            if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
                world0->updateAtomLists(integAdvancedState0);
                for(int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                    (world0->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, j-1);
                    //writePdb( ((SimTK::Compound)(world0->getTopology(mol_i))), 
                    //    integAdvancedState0, "pdbs", "sb_", 10, 
                    //    (std::string("HMC") + std::to_string(mol_i)).c_str(), j-1);
                }
            }
        }

        //std::cout << "WORLD 1" << std::endl;
        //std::cout << "integAdvancedState0.getQ() before: " << integAdvancedState0.getQ() << std::endl;
        //std::cout << "integAdvancedState1.getQ() before" << integAdvancedState1.getQ() << std::endl;
        //std::cout << "MOD (WORLD1) PE before: " << world1->forces->getMultibodySystem().calcPotentialEnergy(integAdvancedState1) << std::endl;

        integAdvancedState1 = world1->setAtomsLocationsInGround(integAdvancedState1, world0->getAtomsLocationsInGround(integAdvancedState0));

        //// Automatically accept and assign old values
        p_HMCsampler1->reinitialize( integAdvancedState1, 
            cons_timestep, total_mcsteps,
            SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );
        //std::cout << "integAdvancedState0.getQ() after: " << integAdvancedState0.getQ() << std::endl;
        //std::cout << "integAdvancedState1.getQ() after" << integAdvancedState1.getQ() << std::endl;
        //std::cout << "MOD (WORLD1) PE after: " << world1->forces->getMultibodySystem().calcPotentialEnergy(integAdvancedState1) << std::endl;

        for(int k = 0; k < cons_mix_mcsteps; k++){
            p_HMCsampler1->update(integAdvancedState1, cons_timestep, cons_mdsteps);
            if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
                world1->updateAtomLists(integAdvancedState1);
                for(int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                    (world1->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, j);
                    //writePdb( ((SimTK::Compound)(world1->getTopology(mol_i))),
                    //    integAdvancedState1, "pdbs", "sb_", 10, 
                    //    (std::string("HMC") + std::to_string(mol_i)).c_str(), j);
                }
            }
        }
        //std::cout << "END ============================" << std::endl;

    } // for i in MC steps

} // END MAIN ////////////



