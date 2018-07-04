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

    context->LoadWorldsFromSetup(setupReader);

    context->Run(setupReader);

    /*
    // Build Gmolmodel simulation worlds
    int nofRegimens = setupReader.getValues("REGIMENS").size();
    std::vector<int> worldIndexes;

    // Create a vector of world indeces   
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        worldIndexes.push_back(worldIx);
    }

    // Set convenient names
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){    
        std::cout << "World " << worldIx << " initial const state PE: " << std::setprecision(20)
            << (context->updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((context->updWorld(worldIx))->integ->updAdvancedState()) 
            << " useFixmanPotential = " << context->updWorld(worldIx)->updSampler(0)->isUsingFixman()
            << std::endl;
    }

    // Get simulation parameters
    int total_mcsteps = std::stoi(setupReader.getValues("STEPS")[0]);
    int nofSamplesPerRound[nofRegimens];
    int nofMDSteps[nofRegimens];
    float timesteps[nofRegimens];


    int currentWorldIx = 0;
    int round_mcsteps = 0;

    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){    
        nofSamplesPerRound[worldIx] = std::stoi(setupReader.getValues("MIXMCSTEPS")[worldIx]);
        assert( (!(total_mcsteps % nofSamplesPerRound[worldIx])) &&
            "Total number of steps must be divisible with each regimen MC steps." );
        round_mcsteps += nofSamplesPerRound[worldIx];
        nofMDSteps[worldIx] = std::stoi(setupReader.getValues("MDSTEPS")[worldIx]);
        timesteps[worldIx] = std::stod(setupReader.getValues("TIMESTEPS")[worldIx]);
    }

    // Calculate geometric features fast
    const SimTK::Compound * p_compounds[nofRegimens];
    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
            p_compounds[worldIx] = &((context->updWorld(worldIx))->getTopology(0));
        }
    }
    //

    // Simulate the two worlds
    int mc_step = -1;

    // Update one round for the first regimen
    currentWorldIx = worldIndexes.front();
    SimTK::State& advancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();

    // Update
    //std::cout << "Sampler " << currentWorldIx << " updating initially" << std::endl;
    for(int k = 0; k < nofSamplesPerRound[currentWorldIx]; k++){
        ++mc_step; // Increment mc_step
        context->updWorld(currentWorldIx)->updSampler(0)->update(advancedState, timesteps[currentWorldIx], nofMDSteps[currentWorldIx]);
    }

    // Write pdb
    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (context->updWorld(currentWorldIx))->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
        }
    }

    while(mc_step < total_mcsteps){

        // Rotate worlds indeces (translate from right to left) 
        std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());
        currentWorldIx = worldIndexes.front();

        // Transfer coordinates from last world to current
        //std::cout << "main: Sending configuration from " << worldIndexes.back() << " to " << currentWorldIx 
        //    << " at MC step " << mc_step << std::endl;
        SimTK::State& lastAdvancedState = (context->updWorld(worldIndexes.back()))->integ->updAdvancedState();
        SimTK::State& currentAdvancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();

        currentAdvancedState = (context->updWorld(currentWorldIx))->setAtomsLocationsInGround(
            currentAdvancedState, (context->updWorld(worldIndexes.back()))->getAtomsLocationsInGround( lastAdvancedState ));

        // Reinitialize current sampler
        context->updWorld(currentWorldIx)->updSampler(0)->reinitialize( currentAdvancedState, 
            timesteps[currentWorldIx], total_mcsteps,
            SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );

        // Set residual energy for the rest of the samplers
        if(setupReader.getValues("REGIMENS")[worldIndexes.back()] == "IC"){
            for(int i = 0; i < nofRegimens - 1; i++){
                int restIx = worldIndexes[i];
                int backIx = worldIndexes.back();
                SimTK::Real diffPE = (context->updWorld(backIx))->updSampler(0)->getSetPE() - (context->updWorld(restIx))->updSampler(0)->getSetPE();
                std::cout << "Setting sampler " << restIx << " REP to " << (context->updWorld(backIx))->updSampler(0)->getSetPE() << " - " << (context->updWorld(restIx))->updSampler(0)->getSetPE() << " = " << diffPE << std::endl;
                (context->updWorld(restIx))->updSampler(0)->setREP( diffPE );
            }
        }

        // Update
        //std::cout << "Sampler " << currentWorldIx << " updating " << std::endl;
        for(int k = 0; k < nofSamplesPerRound[currentWorldIx]; k++){
            ++mc_step; // Increment mc_step
            context->updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState, timesteps[currentWorldIx], nofMDSteps[currentWorldIx]);
        }

        // Write pdb
        if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
            //if(!((mc_step+1) % 20)){
            if(1){
                (context->updWorld(currentWorldIx))->updateAtomLists(currentAdvancedState);
                std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
                for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                    ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
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

    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        //delete[] p_compounds;
    }
    */


    delete context;

} // END MAIN ////////////



