/** @file
 * This file tests the Context class
 */

#include <string>
#include <iostream>
#include <sstream>
#include <chrono>
#include "simmain.hpp"
#include "Robo.hpp"

#include "HamiltonianMonteCarloSampler.hpp"
#include "readAmberInput.hpp"
#include "SetupReader.hpp"

//#ifndef ROBO_DEBUG_LEVEL01
//#define ROBO_DEBUG_LEVEL01
//#endif

int main(int argc, char **argv)
{
    // Initialize setup reader
    SetupReader setupReader(argv[1]);
    std::cout << "SETUP" << std::endl;
    setupReader.dump();

    Context *context = new Context();


    // Variables
    // Output option
    unsigned long int printFreq = std::stoi(setupReader.getValues("PRINT_FREQ")[0]);

    int currentWorldIx = 0;
    int lastWorldIx = 0;
    int round_mcsteps = 0;
    int total_mcsteps = 0;

    int mc_step = -1;
    int alt_mc_step = 0;
    int restore_mc_step = 0;
    SimTK::Real convFunc = 99.0;

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
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("PRMTOP")[molIx]
            );

            context->loadCoordinatesFile( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("INPCRD")[molIx]
            );

            context->loadRigidBodiesSpecs( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("RBFILE")[molIx]
            );

            context->loadFlexibleBondsSpecs( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("FLEXFILE")[molIx]
            );

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
            context->useFixmanTorque(worldIx, std::stof(setupReader.getValues("TEMPERATURE")[worldIx]));
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

    // Use OpenMM
    if(setupReader.getValues("OPENMM")[0] == "TRUE"){
        context->setUseOpenMMAcceleration(true);
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

            //if(setupReader.getValues("WORLDS")[worldIx] == "TD"){
            //    context->initializeSampler(worldIx, samplerIx, true);
            //}else{
                context->initializeSampler(worldIx, samplerIx);
            //}

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
                << " temperature = " << context->updWorld(worldIx)->updSampler(samplerIx)->getTemperature()
                << " initial const state PE: " << std::setprecision(20)
                //<< (context->updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState())
                << (context->updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy(context->updAdvancedState(worldIx, samplerIx))
                << " useFixmanPotential = " << context->updWorld(worldIx)->updSampler(samplerIx)->isUsingFixmanPotential()
                << std::endl;
        }

    }


    currentWorldIx = 0;
    lastWorldIx = 0;
    round_mcsteps = 0;

    context->setNofRounds(std::stoi(setupReader.getValues("ROUNDS")[0]));

    for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){    
        round_mcsteps += context->getNofSamplesPerRound(worldIx);
    }

    total_mcsteps = round_mcsteps * context->getNofRounds();

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

    // Helper variables for step arithmetic
    mc_step = -1;
    alt_mc_step = 0;
    restore_mc_step = 0;
    convFunc = 99.0;

    // Update one round for the first regimen
    currentWorldIx = context->worldIndexes.front();
    SimTK::State& advancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();

    // Update
    //std::cout << "Sampler " << currentWorldIx << " updating initially" << std::endl;
    for(int k = 0; k < context->getNofSamplesPerRound(currentWorldIx); k++){
        ++mc_step; // Increment mc_step
        context->updWorld(currentWorldIx)->updSampler(0)->update(advancedState, context->getNofMDStepsPerSample(currentWorldIx));
        //context->updWorld(currentWorldIx)->updSampler(0)->perturbQ(advancedState);
    }

    // Write pdb
    std::string pdbPrefix = setupReader.getValues("MOLECULES")[0] 
        + std::to_string(context->updWorld(currentWorldIx)->updSampler(0)->getSeed());

    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (context->updWorld(currentWorldIx))->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb." + pdbPrefix + ".", ".pdb", 10, mc_step);
        }
    }

    // Output option
    //unsigned long int printFreq = std::stoi(setupReader.getValues("PRINT_FREQ")[0]);

    // Adapt variables
    bool samplesAdapted = false;
    bool restored = false;
    SimTK::Real cumU = 0.0;
    SimTK::Real avgU = 0.0;
    SimTK::Real prevAvgU = 0.0;

    if(setupReader.getValues("ADAPT_SAMPLE_RATIO")[0] == "TRUE"){
        std::cout << "Adaptive samples ratios ON." << std::endl;
    }

    // Heating
    //for(unsigned int round = 0; round < context->getNofHeatingRounds; round++){
    //}

    // Production
    for(int round = 0; round < context->getNofRounds(); round++){ // Iterate rounds

        for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){ // Iterate worlds
    
            // Rotate worlds indeces (translate from right to left) 
            std::rotate(context->worldIndexes.begin(), context->worldIndexes.begin() + 1, context->worldIndexes.end());

            // Get indeces
            currentWorldIx = context->worldIndexes.front();
            lastWorldIx = context->worldIndexes.back();
    
            // Transfer coordinates from last world to current
            SimTK::State& lastAdvancedState = (context->updWorld(context->worldIndexes.back()))->integ->updAdvancedState();
            SimTK::State& currentAdvancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();
    
            currentAdvancedState = (context->updWorld(currentWorldIx))->setAtomsLocationsInGround(
               currentAdvancedState, (context->updWorld(context->worldIndexes.back()))->getAtomsLocationsInGround( lastAdvancedState ));
    
            // Set old potential energy of the new world
            (context->updWorld(currentWorldIx))->updSampler(0)->setOldPE(
                (context->updWorld(context->worldIndexes.back()))
                ->updSampler(0)->getSetPE() );
    
            // Reinitialize current sampler
            context->updWorld(currentWorldIx)->updSampler(0)->reinitialize( currentAdvancedState);
    
            // Update
            for(int k = 0; k < context->getNofSamplesPerRound(currentWorldIx); k++){ // Iterate through samples
                context->updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState, context->getNofMDStepsPerSample(currentWorldIx));
                
                #ifdef ROBO_DEBUG_LEVEL01
                std::this_thread::sleep_for(std::chrono::milliseconds(2000)); 
                #endif

                ++mc_step;
                ++restore_mc_step;
   
            } // END for samples

            // Print energy and geometric features
            if( !(round % printFreq) ){
                // ndofs accs pe_o pe_set ke_o ke_n fix_o fix_set fix_n
                context->PrintSamplerData(currentWorldIx);
                context->PrintGeometry(setupReader, currentWorldIx);
            }
    
            // Write pdb
            if( std::stoi(setupReader.getValues("WRITEPDBS")[0]) != 0){
                if(((mc_step) % std::stoi(setupReader.getValues("WRITEPDBS")[0])) == 0){
                    (context->updWorld(currentWorldIx))->updateAtomLists(currentAdvancedState);
                    for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                        ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb." + pdbPrefix + ".", ".pdb", 10, mc_step);
                    }
                }
            } // if write pdbs

        } // for i in worlds
    } // for i in rounds

    // Write final pdbs
    for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
        ((context->updWorld( context->worldIndexes.front() ))->getTopology(mol_i)).writePdb("pdbs", "final." + pdbPrefix + ".", ".pdb", 10, context->getNofRounds());
    }
    //

    delete context;

} // END MAIN ////////////







