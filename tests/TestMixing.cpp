/** @file
 * This file tests the Context class
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
    std::cout << "SETUP" << std::endl;
    setupReader.dump();

    Context *context = new Context();


    // Add Worlds
    unsigned int nofWorlds = setupReader.getValues("WORLDS").size();
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        if(setupReader.getValues("VISUAL")[worldIx] == "TRUE"){
            context->AddWorld(true);
        }else{
            context->AddWorld(false);
        }
    }

    // Add molecules to worlds
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for(unsigned int molIx = 0; molIx < setupReader.getValues("MOLECULES").size(); molIx++){
            context->loadTopologyFile( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.prmtop") );
            context->loadCoordinatesFile( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.inpcrd") );
            context->loadRigidBodiesSpecs( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.rb") );
            context->loadFlexibleBondsSpecs( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.flex") );
            context->setRegimen( worldIx, molIx,
                setupReader.getValues("WORLDS")[worldIx] );
        }
    }
    context->loadMolecules();

    // Set worlds force field scale factors
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        // Set force field scale factors.
        if(setupReader.getValues("FFSCALE")[worldIx] == "AMBER"){
            context->setAmberForceFieldScaleFactors(worldIx);
        }else{
            context->setGlobalForceFieldScaleFactor(worldIx, std::stod(setupReader.getValues("FFSCALE")[worldIx]));
        }
        // Set world GBSA scale factor
        context->setGbsaGlobalScaleFactor(worldIx, std::stod(setupReader.getValues("GBSA")[worldIx]));
    }

    // Do we use Fixman torque
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        if(setupReader.getValues("FIXMAN_TORQUE")[worldIx] == "TRUE"){
            context->useFixmanTorque(worldIx);
        }
    }

    // Set the number of threads
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        context->setNumThreadsRequested(worldIx, std::stod(setupReader.getValues("THREADS")[worldIx]));
    }

    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        std::cout << "main context->updWorld(" << worldIx<< ")->getNumThreadsRequested() " 
            << context->updWorld(worldIx)->updForceField()->getNumThreadsRequested() << std::endl;
        std::cout << "main context->updWorld(" << worldIx<< ")->getNumThreadsInUse() " 
            << context->updWorld(worldIx)->updForceField()->getNumThreadsInUse() << std::endl;
    }

    // Model topologies
    context->modelTopologies();


    // Set Fixman torques scale factors
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        if(setupReader.getValues("FIXMAN_TORQUE")[worldIx] == "TRUE"){
            std::cout << "main call context->setFixmanTorqueScaleFactor(" << worldIx << " -1 " << std::endl;
            context->setFixmanTorqueScaleFactor(worldIx, -1.0);
        }
    }

    // Add samplers to the worlds
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        context->addSampler(worldIx, HMC);
    }

    // Set sampler parameters and initialize
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for (unsigned int samplerIx = 0; samplerIx < context->getWorld(worldIx)->getNofSamplers(); samplerIx++){
            // Set timesteps
            context->setTimestep(worldIx, samplerIx, std::stod(setupReader.getValues("TIMESTEPS")[worldIx]) );

            // Set thermostats
            context->updWorld(worldIx)->updSampler(samplerIx)->setThermostat(setupReader.getValues("THERMOSTAT")[worldIx]);

            // Activate Fixman potential if needed
            if(setupReader.getValues("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
                 context->useFixmanPotential(worldIx, samplerIx);
            }
        }
    }

    // This loop is just for check purposes (should be removed)
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for (unsigned int samplerIx = 0; samplerIx < context->getWorld(worldIx)->getNofSamplers(); samplerIx++){
            std::cout << "After setThermo world " << worldIx << " sampler " << samplerIx << "getThermostat: " ;
            std::cout << context->updWorld(worldIx)->updSampler(samplerIx)->getThermostat() ;
            std::cout << std::endl;
        }
    }


    // Make the simulation reproducible 
    if(setupReader.getValues("REPRODUCIBLE")[0] == "TRUE"){
        context->setReproducible();
        srand (0);
    }else{
        srand (time(NULL));
    }

    for(int worldIx = 0; worldIx < setupReader.getValues("WORLDS").size(); worldIx++){
        context->setTemperature(worldIx, std::stof(setupReader.getValues("TEMPERATURE")[worldIx]));
        context->setNofSamplesPerRound(worldIx, std::stoi(setupReader.getValues("SAMPLES_PER_ROUND")[worldIx]));
        context->setNofMDStepsPerSample(worldIx, std::stoi(setupReader.getValues("MDSTEPS")[worldIx]));
    }

    // Add samplers to the worlds
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for (unsigned int samplerIx = 0; samplerIx < context->getWorld(worldIx)->getNofSamplers(); samplerIx++){
            context->initializeSampler(worldIx, samplerIx);
        }
    }

    // Print thermodynamics
    for(unsigned int worldIx = 0; worldIx < setupReader.getValues("WORLDS").size(); worldIx++){
        std::cout << "MAIN World " << worldIx << " temperature = " << context->getWorld(worldIx)->getTemperature() << std::endl;
        if(setupReader.getValues("FIXMAN_TORQUE")[worldIx] == "TRUE"){
            std::cout << "MAIN World " << worldIx << " FixmanTorque temperature = " << context->updWorld(worldIx)->updFixmanTorque()->getTemperature() << std::endl;
        }
        for (unsigned int samplerIx = 0; samplerIx < context->getWorld(worldIx)->getNofSamplers(); samplerIx++){
            std::cout << "MAIN World " << worldIx << " Sampler " << samplerIx 
                << " temperature = " << context->updWorld()->updSampler(samplerIx)->getTemperature()
                << " initial const state PE: " << std::setprecision(20)
                //<< (context->updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState())
                << (context->updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy(context->updAdvancedState(worldIx, samplerIx))
                << " useFixmanPotential = " << context->updWorld(worldIx)->updSampler(0)->isUsingFixmanPotential()
                << std::endl;
        }

    }


        
    int currentWorldIx = 0;
    int round_mcsteps = 0;

    context->setNofRounds(std::stoi(setupReader.getValues("ROUNDS")[0]));

    for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){    
        //context->getNofSamplesPerRound(worldIx) = std::stoi(setupReader.getValues("MIXMCSTEPS")[worldIx]);
        round_mcsteps += context->getNofSamplesPerRound(worldIx);
        //nofMDStepsPerSample[worldIx] = std::stoi(setupReader.getValues("MDSTEPS")[worldIx]);
        //timesteps[worldIx] = std::stod(setupReader.getValues("TIMESTEPS")[worldIx]);
    }
    int total_mcsteps = round_mcsteps * context->getNofRounds();

    // Calculate geometric features
    SimTK::Real dihedrals[setupReader.getValues("DIHEDRAL").size() / 4];
    SimTK::Real dihMeans[setupReader.getValues("DIHEDRAL").size() / 4];
    SimTK::Real dihVars[setupReader.getValues("DIHEDRAL").size() / 4];
    SimTK::Real distances[setupReader.getValues("DISTANCE").size() / 2];

    const SimTK::Compound * p_compounds[context->getNofWorlds()];
    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){
            p_compounds[worldIx] = &((context->updWorld(worldIx))->getTopology(0));
        }
    }
    //

    // Simulate the two worlds
    int mc_step = -1;

    // Update one round for the first regimen
    currentWorldIx = context->worldIndexes.front();
    SimTK::State& advancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();

    // Update
    //std::cout << "Sampler " << currentWorldIx << " updating initially" << std::endl;
    for(int k = 0; k < context->getNofSamplesPerRound(currentWorldIx); k++){
        ++mc_step; // Increment mc_step
        context->updWorld(currentWorldIx)->updSampler(0)->update(advancedState, context->getNofMDStepsPerSample(currentWorldIx));
    }

    // Write pdb
    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (context->updWorld(currentWorldIx))->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
        }
    }

    for(int round = 0; round < context->getNofRounds(); round++){
        for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){
    
            // Rotate worlds indeces (translate from right to left) 
            std::rotate(context->worldIndexes.begin(), context->worldIndexes.begin() + 1, context->worldIndexes.end());

            currentWorldIx = context->worldIndexes.front();
    
            // Transfer coordinates from last world to current
            //std::cout << "main: Sending configuration from " << context->worldIndexes.back() << " to " << currentWorldIx 
            //    << " at round " << round << std::endl;
            SimTK::State& lastAdvancedState = (context->updWorld(context->worldIndexes.back()))->integ->updAdvancedState();
            SimTK::State& currentAdvancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();
    
            currentAdvancedState = (context->updWorld(currentWorldIx))->setAtomsLocationsInGround(
                currentAdvancedState, (context->updWorld(context->worldIndexes.back()))->getAtomsLocationsInGround( lastAdvancedState ));
    
            // Reinitialize current sampler
            context->updWorld(currentWorldIx)->updSampler(0)->reinitialize( currentAdvancedState
//r                , SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) 
            );
    
            // INCORRECT !!
            if(setupReader.getValues("WORLDS")[context->worldIndexes.back()] == "IC"){
                for(int i = 0; i < context->getNofWorlds() - 1; i++){
                    int restIx = context->worldIndexes[i];
                    int backIx = context->worldIndexes.back();
                    SimTK::Real diffPE = (context->updWorld(backIx))->updSampler(0)->getSetPE() - (context->updWorld(restIx))->updSampler(0)->getSetPE();
                    //std::cout << "Setting sampler " << restIx << " REP to " << (context->updWorld(backIx))->updSampler(0)->getSetPE() << " - " << (context->updWorld(restIx))->updSampler(0)->getSetPE() << " = " << diffPE << std::endl;
                    (context->updWorld(restIx))->updSampler(0)->setREP( diffPE );
                }
            }
    
            // Update
            //std::cout << "Sampler " << currentWorldIx << " updating " << std::endl;
            for(int k = 0; k < context->getNofSamplesPerRound(currentWorldIx); k++){
                context->updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState, context->getNofMDStepsPerSample(currentWorldIx));
            }
    
            // Write pdb
            if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
                //if(!((mc_step+1) % 20)){
                if(1){
                    (context->updWorld(currentWorldIx))->updateAtomLists(currentAdvancedState);
                    //std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
                    for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                        ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, round);
                    }
                }
            }
    
            // Calculate geomtric features 
            if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){

                int distanceIx[setupReader.getValues("DISTANCE").size()];
                for(unsigned int i = 0; i < setupReader.getValues("DISTANCE").size(); i++){
                    distanceIx[i] = atoi(setupReader.getValues("DISTANCE")[i].c_str());
                }
    
                for(int ai = 0; ai < (setupReader.getValues("DISTANCE").size() / 2); ai++){
                    distances[ai] = context->Distance(currentWorldIx, 0, 0, distanceIx[2*ai + 0], distanceIx[2*ai + 1]);
                    std::cout << std::setprecision(4) << distances[ai] << " ";
                }

                int dihedralIx[setupReader.getValues("DIHEDRAL").size()];
                for(unsigned int i = 0; i < setupReader.getValues("DIHEDRAL").size(); i++){
                    dihedralIx[i] = atoi(setupReader.getValues("DIHEDRAL")[i].c_str());
                }
    
                for(int ai = 0; ai < (setupReader.getValues("DIHEDRAL").size() / 4); ai++){
                    dihedrals[ai] = context->Dihedral(currentWorldIx, 0, 0, dihedralIx[4*ai + 0], dihedralIx[4*ai + 1], dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]);
                    std::cout 
                        << std::setprecision(4) << dihedrals[ai] << " ";
                }
                std::cout << std::endl;
            }else{
                std::cout << std::endl;
            }
    
        } // for i in worlds
    } // for i in rounds


    delete context;

} // END MAIN ////////////







