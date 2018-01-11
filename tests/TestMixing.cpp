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
    int nofRegimens = setupReader.getValues("REGIMENS").size();
    std::vector<int> worldIndexes;
    std::vector<World *> p_worlds;
    std::vector<HamiltonianMonteCarloSampler *> p_samplers;
    
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        if(setupReader.getValues("VISUAL")[0] == "TRUE"){
            p_worlds.push_back( new World(worldIx, true, std::stod(setupReader.getValues("TIMESTEPS")[worldIx])) );
        }else{
            p_worlds.push_back( new World(worldIx, false, std::stod(setupReader.getValues("TIMESTEPS")[worldIx])) );
        }

        // Set world identifiers
        (p_worlds[worldIx])->ownWorldIndex = worldIx;
        worldIndexes.push_back(worldIx);
    
        // Get input filenames and add molecules
        std::vector<std::string> argValues;
        std::vector<std::string>::iterator argValuesIt;
        argValues = setupReader.getValues("MOLECULES");
        for(argValuesIt = argValues.begin(); argValuesIt != argValues.end(); ++argValuesIt){
            std::string prmtopFN = *argValuesIt + std::string("/ligand.prmtop");
            std::string inpcrdFN = *argValuesIt + std::string("/ligand.inpcrd");
            std::string rbFN = *argValuesIt + std::string("/ligand.rb");
            std::string flexFN = *argValuesIt + std::string("/ligand.flex");
    
            std::string ictd = setupReader.getValues("REGIMENS")[worldIx];
    
            readAmberInput *amberReader = new readAmberInput();
            amberReader->readAmberFiles(inpcrdFN, prmtopFN);
            (p_worlds[worldIx])->AddMolecule(amberReader, rbFN, flexFN, ictd);
    
            delete amberReader;
        }
    
        // Initialize worlds
        (p_worlds[worldIx])->Init( std::stod(setupReader.getValues("TIMESTEPS")[worldIx]) );
    
        // Initialize sampler
        p_samplers.push_back( new HamiltonianMonteCarloSampler(
            (p_worlds[worldIx])->compoundSystem, (p_worlds[worldIx])->matter, 
            (SimTK::Compound *)(&((p_worlds[worldIx])->getTopology(0))),
            p_worlds[worldIx]->forceField, (p_worlds[worldIx])->forces, (p_worlds[worldIx])->ts ) );
    
        // Initialize samplers
        bool useFixman = true;
        if(setupReader.getValues("REGIMENS")[worldIx] == "IC"){
            useFixman = false;
        }
        (p_samplers[worldIx])->initialize( (p_worlds[worldIx])->integ->updAdvancedState(), 
             std::stod(setupReader.getValues("TIMESTEPS")[worldIx]),
             std::stoi(setupReader.getValues("STEPS")[0]),
             SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ),
             useFixman );
    
    }
    ////////////////////////////////////////////////////////// END creating Worlds

    // Add worlds to context
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){    
        context->AddWorld(p_worlds[worldIx], p_samplers[worldIx]);
    }

    // Set convenient names
    //World *world0 = context->getWorld(0);
    //World *world1 = context->getWorld(1);
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){    
        std::cout << "World " << worldIx << " initial const state PE: " << std::setprecision(20)
            << (p_worlds[worldIx])->forces->getMultibodySystem().calcPotentialEnergy((p_worlds[worldIx])->integ->updAdvancedState()) 
            << " useFixman = " << p_samplers[worldIx]->isUsingFixman()
            << std::endl;
    }

    // Get simulation parameters
    int total_mcsteps = std::stoi(setupReader.getValues("STEPS")[0]);
    int mix_mcsteps[nofRegimens];
    int mdsteps[nofRegimens];
    float timesteps[nofRegimens];
    int remainders[nofRegimens];
    int currentWorldIx = 0;
    int round_mcsteps = 0;

    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){    
        mix_mcsteps[worldIx] = std::stoi(setupReader.getValues("MIXMCSTEPS")[worldIx]);
        assert( (!(total_mcsteps % mix_mcsteps[worldIx])) &&
            "Total number of steps must be divisible with each regimen MC steps." );
        round_mcsteps += mix_mcsteps[worldIx];
        mdsteps[worldIx] = std::stoi(setupReader.getValues("MDSTEPS")[worldIx]);
        timesteps[worldIx] = std::stod(setupReader.getValues("TIMESTEPS")[worldIx]);
    }

    // Simulate the two worlds
    int mc_step = -1;
    // Update one round for the first regimen
    currentWorldIx = worldIndexes.front();
    SimTK::State& advancedState = (p_worlds[currentWorldIx])->integ->updAdvancedState();
    //SimTK::State& defaultState = (p_worlds[currentWorldIx])->integ->updDefaultState();
    std::cout << "Sampler " << currentWorldIx << " updating " << std::endl;
    for(int k = 0; k < mix_mcsteps[currentWorldIx]; k++){
        ++mc_step; // Increment mc_step
        p_samplers[currentWorldIx]->update(advancedState,
            timesteps[currentWorldIx], mdsteps[currentWorldIx]);
    }
    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (p_worlds[currentWorldIx])->updateAtomLists(advancedState);
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((p_worlds[currentWorldIx])->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
        }
    }

    while(mc_step < total_mcsteps){

        // Get coordinates from last world and accept them automatically
        std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());
        currentWorldIx = worldIndexes.front();

        // Transfer coordinates from last world to current
        std::cout << "main: Sending configuration from " << worldIndexes.back() << " to " << currentWorldIx 
            << " at MC step " << mc_step << std::endl;
        const SimTK::State& lastConstState = (p_worlds[worldIndexes.back()])->integ->getState();
        const SimTK::State& currentConstState = (p_worlds[currentWorldIx])->integ->getState();
        SimTK::State& lastAdvancedState = (p_worlds[worldIndexes.back()])->integ->updAdvancedState();
        SimTK::State& currentAdvancedState = (p_worlds[currentWorldIx])->integ->updAdvancedState();

        // CHECK COORD ASS Write source pdb
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            (p_worlds[worldIndexes.back()])->updateAtomLists(lastAdvancedState);
            ((p_worlds[worldIndexes.back()])->getTopology(mol_i)).writePdb("pdbs", "sb_src", ".pdb", 10, mc_step);
            std::cout << "world " << worldIndexes.back() << " ";
            (p_worlds[worldIndexes.back()])->printPoss((p_worlds[worldIndexes.back()])->getTopology(mol_i), lastAdvancedState);
        }
        currentAdvancedState = (p_worlds[currentWorldIx])->setAtomsLocationsInGround(
            currentAdvancedState, (p_worlds[worldIndexes.back()])->getAtomsLocationsInGround( lastAdvancedState ));
        // CHECK COORD ASS Write dest pdb
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((p_worlds[currentWorldIx])->getTopology(mol_i)).writePdb("pdbs", "sb_dest", ".pdb", 10, mc_step);
            std::cout << "world " << currentWorldIx << " ";
            (p_worlds[currentWorldIx])->printPoss((p_worlds[currentWorldIx])->getTopology(mol_i), currentAdvancedState);
        }

        //std::cout << "Advanced states info: "<< lastAdvancedState << " " << currentAdvancedState << std::endl;
        // Reinitialize current sampler
        //std::cout << "World " << worldIndexes.back() << " State Cache Info: " << std::endl;
        //(p_worlds[worldIndexes.back()])->PrintSimbodyStateCache(lastAdvancedState);
        //std::cout << "Test World " << currentWorldIx << " before reinitialize State Cache Info before reinitialize: " << std::endl;
        //(p_worlds[currentWorldIx])->PrintSimbodyStateCache(currentAdvancedState);
        //std::cout << "integrator addresses: "
        //    << &(p_samplers[currentWorldIx]->timeStepper->updIntegrator()) << " " << &(p_worlds[currentWorldIx]->integ) << " "
        //    << &(p_samplers[worldIndexes.back()]->timeStepper->updIntegrator()) << " " << &(p_worlds[worldIndexes.back()]->integ) << std::endl;
        p_samplers[currentWorldIx]->reinitialize( currentAdvancedState, 
            timesteps[currentWorldIx], total_mcsteps,
            SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );

        // Set residual energy for current sampler
        (p_samplers[currentWorldIx])->setREP( (p_samplers[worldIndexes.back()])->getSetPE() - (p_samplers[currentWorldIx])->getSetPE() );

        std::cout << "Sampler " << currentWorldIx << " updating " << std::endl;
        for(int k = 0; k < mix_mcsteps[currentWorldIx]; k++){
            ++mc_step; // Increment mc_step
            p_samplers[currentWorldIx]->update(currentAdvancedState, timesteps[currentWorldIx], mdsteps[currentWorldIx]);
        }
        if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
            (p_worlds[currentWorldIx])->updateAtomLists(currentAdvancedState);
            for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                ((p_worlds[currentWorldIx])->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
            }
        }

    } // for i in MC steps

    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        delete p_worlds[worldIx];
        delete p_samplers[worldIx];
    }
    delete context;

} // END MAIN ////////////



