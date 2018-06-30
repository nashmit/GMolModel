/** @file
 * This file contains a native Gmolmodel test. 
 */

#include <string>
#include <iostream>
#include <sstream>
#include "simmain.hpp"
#include "Robo.hpp"

#include "HamiltonianMonteCarloSampler.hpp"
#include "readAmberInput.hpp"
#include "SetupReader.hpp"

int main(int argc, char **argv)
{
    // Initialize setup reader
    SetupReader setupReader(argv[1]);
    TRACE("NEW ALLOC\n");
    Context *context = new Context();
    srand (time(NULL));

    // Build Gmolmodel simulation worlds
    int nofRegimens = setupReader.getValues("REGIMENS").size();
    std::vector<int> worldIndexes;
    std::vector<World *> p_worlds;
    std::vector<HamiltonianMonteCarloSampler *> p_samplers;
    
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        if(setupReader.getValues("VISUAL")[0] == "TRUE"){
            TRACE("NEW ALLOC\n");
            p_worlds.push_back( new World(worldIx, true, std::stod(setupReader.getValues("TIMESTEPS")[worldIx])) );
        }else{
            TRACE("NEW ALLOC\n");
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
    
            TRACE("NEW ALLOC\n");
            readAmberInput *amberReader = new readAmberInput();
            amberReader->readAmberFiles(inpcrdFN, prmtopFN);
            (p_worlds[worldIx])->AddMolecule(amberReader, rbFN, flexFN, ictd);
    
            delete amberReader;
        }
    
        // Initialize worlds
        if(setupReader.getValues("FIXMAN_TORQUE")[worldIx] == "TRUE"){
            (p_worlds[worldIx])->Init( std::stod(setupReader.getValues("TIMESTEPS")[worldIx]), true );
        }else{
            (p_worlds[worldIx])->Init( std::stod(setupReader.getValues("TIMESTEPS")[worldIx]), false );
        }

        (p_worlds[worldIx])->setTemperature( std::stod(setupReader.getValues("TEMPERATURE")[worldIx]) );
    
        // Initialize sampler
        TRACE("NEW ALLOC\n");
        p_samplers.push_back( new HamiltonianMonteCarloSampler(
            (p_worlds[worldIx])->compoundSystem, (p_worlds[worldIx])->matter, 
            (SimTK::Compound *)(&((p_worlds[worldIx])->getTopology(0))),
            p_worlds[worldIx]->forceField, (p_worlds[worldIx])->forces, (p_worlds[worldIx])->ts ) );
    
        // Initialize samplers
        bool useFixmanPotential = false;
        if(setupReader.getValues("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
            useFixmanPotential = true;
        }

        // We don't need Fixman potential for IC !! WATCH OUT !!
        //if(setupReader.getValues("REGIMENS")[worldIx] == "IC"){
        //    useFixmanPotential = false;
        //}

        // Set thermostats
        (p_samplers[worldIx])->setThermostat(setupReader.getValues("THERMOSTAT")[worldIx]);

        (p_samplers[worldIx])->initialize( (p_worlds[worldIx])->integ->updAdvancedState(), 
             std::stod(setupReader.getValues("TIMESTEPS")[worldIx]),
             std::stoi(setupReader.getValues("STEPS")[0]),
             SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[worldIx]) ),
             useFixmanPotential );
    
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
            << " useFixmanPotential = " << p_samplers[worldIx]->isUsingFixman()
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

    // Calculate geometric features fast
    const SimTK::Compound * p_compounds[nofRegimens];
    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
            p_compounds[worldIx] = &((p_worlds[worldIx])->getTopology(0));
        }
    }
    //

    // Simulate the two worlds
    int mc_step = -1;

    // Update one round for the first regimen
    currentWorldIx = worldIndexes.front();
    SimTK::State& advancedState = (p_worlds[currentWorldIx])->integ->updAdvancedState();

    // Update
    //std::cout << "Sampler " << currentWorldIx << " updating initially" << std::endl;
    for(int k = 0; k < mix_mcsteps[currentWorldIx]; k++){
        ++mc_step; // Increment mc_step
        p_samplers[currentWorldIx]->update(advancedState, timesteps[currentWorldIx], mdsteps[currentWorldIx]);
    }

    // Write pdb
    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (p_worlds[currentWorldIx])->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((p_worlds[currentWorldIx])->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
        }
    }

    while(mc_step < total_mcsteps){

        // Rotate worlds indeces (translate from right to left) 
        std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());
        currentWorldIx = worldIndexes.front();

        // Transfer coordinates from last world to current
        //std::cout << "main: Sending configuration from " << worldIndexes.back() << " to " << currentWorldIx 
        //    << " at MC step " << mc_step << std::endl;
        SimTK::State& lastAdvancedState = (p_worlds[worldIndexes.back()])->integ->updAdvancedState();
        SimTK::State& currentAdvancedState = (p_worlds[currentWorldIx])->integ->updAdvancedState();

        currentAdvancedState = (p_worlds[currentWorldIx])->setAtomsLocationsInGround(
            currentAdvancedState, (p_worlds[worldIndexes.back()])->getAtomsLocationsInGround( lastAdvancedState ));

        // Reinitialize current sampler
        p_samplers[currentWorldIx]->reinitialize( currentAdvancedState, 
            timesteps[currentWorldIx], total_mcsteps,
            SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );

        // Set residual energy for the rest of the samplers
        if(setupReader.getValues("REGIMENS")[worldIndexes.back()] == "IC"){
            for(int i = 0; i < nofRegimens - 1; i++){
                int restIx = worldIndexes[i];
                int backIx = worldIndexes.back();
                SimTK::Real diffPE = (p_samplers[backIx])->getSetPE() - (p_samplers[restIx])->getSetPE();
                std::cout << "Setting sampler " << restIx << " REP to " << (p_samplers[backIx])->getSetPE() << " - " << (p_samplers[restIx])->getSetPE() << " = " << diffPE << std::endl;
                (p_samplers[restIx])->setREP( diffPE );
            }
        }

        // Update
        //std::cout << "Sampler " << currentWorldIx << " updating " << std::endl;
        for(int k = 0; k < mix_mcsteps[currentWorldIx]; k++){
            ++mc_step; // Increment mc_step
            p_samplers[currentWorldIx]->update(currentAdvancedState, timesteps[currentWorldIx], mdsteps[currentWorldIx]);
        }

        // Write pdb
        if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
            //if(!((mc_step+1) % 20)){
            if(1){
                (p_worlds[currentWorldIx])->updateAtomLists(currentAdvancedState);
                std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
                for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                    ((p_worlds[currentWorldIx])->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
                }
            }
        }

        // Calculate geomtric features fast
        if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
            int dihedralIx[setupReader.getValues("DIHEDRAL").size()];
            for(unsigned int i = 0; i < setupReader.getValues("DIHEDRAL").size(); i++){
                dihedralIx[i] = atoi(setupReader.getValues("DIHEDRAL")[i].c_str());
            }

            for(int ai = 0; ai < setupReader.getValues("DIHEDRAL").size(); ai += 4){
                //std::cout << " atoms = " 
                //    << dihedralIx[ai + 0] << " " 
                //    << dihedralIx[ai + 1] << " " 
                //    << dihedralIx[ai + 2] << " " 
                //    << dihedralIx[ai + 3] << "; names = "
                //    << (p_compounds[currentWorldIx])->getAtomName(SimTK::Compound::AtomIndex(dihedralIx[ai + 0])) << " "
                //    << (p_compounds[currentWorldIx])->getAtomName(SimTK::Compound::AtomIndex(dihedralIx[ai + 1])) << " "
                //    << (p_compounds[currentWorldIx])->getAtomName(SimTK::Compound::AtomIndex(dihedralIx[ai + 2])) << " "
                //    << (p_compounds[currentWorldIx])->getAtomName(SimTK::Compound::AtomIndex(dihedralIx[ai + 3])) 
                //    << "; dihedral = " ;
                std::cout << bDihedral( (p_compounds[currentWorldIx])->calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(dihedralIx[ai + 0])), 
                                        (p_compounds[currentWorldIx])->calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(dihedralIx[ai + 1])), 
                                        (p_compounds[currentWorldIx])->calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(dihedralIx[ai + 2])), 
                                        (p_compounds[currentWorldIx])->calcAtomLocationInGroundFrame(currentAdvancedState, SimTK::Compound::AtomIndex(dihedralIx[ai + 3])) )  << " ";
            }
            std::cout << std::endl;
        }


    } // for i in MC steps

    // Free the memory
    for(int wIx = 0; wIx < nofRegimens; wIx++){
        delete p_worlds[wIx];
        delete p_samplers[wIx];
    }

    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        //delete[] p_compounds;
    }

    delete context;

} // END MAIN ////////////



