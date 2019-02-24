/** @file
 * This file tests the Context class
 */

#include <string>
#include <iostream>
#include <sstream>
#include <chrono>
#include <cassert>
#include <sys/stat.h>
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
    unsigned int nofWorlds = setupReader.getValues("WORLDS").size();
    std::cout << "SETUP" << std::endl;
    setupReader.dump();

    // Assert minimal requirements
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for(unsigned int molIx = 0; molIx < setupReader.getValues("MOLECULES").size(); molIx++){
            assert(SimTK::Pathname::fileExists(
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("PRMTOP")[molIx]) );
            assert(SimTK::Pathname::fileExists(
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("INPCRD")[molIx]) );
            assert(SimTK::Pathname::fileExists(
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("RBFILE")[molIx]) );
            assert(SimTK::Pathname::fileExists(
                setupReader.getValues("MOLECULES")[molIx] + std::string("/")
                + setupReader.getValues("FLEXFILE")[molIx]) );

        }
    }
    assert(SimTK::Pathname::fileExists(setupReader.getValues("OUTPUT_DIR")[0]));

    // Create pdbs directory if necessary
    if( !(SimTK::Pathname::fileExists(setupReader.getValues("OUTPUT_DIR")[0] + "/pdbs")) ){

        const int dir_err = mkdir(
            (setupReader.getValues("OUTPUT_DIR")[0] + "/pdbs").c_str(),
            S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        if (-1 == dir_err){
            std::cout << "Error creating " << setupReader.getValues("OUTPUT_DIR")[0] + "/pdbs" << std::endl;
            exit(1);
        }
    }

    // Variables
    int currentWorldIx = 0;
    int lastWorldIx = 0;
    int round_mcsteps = 0;
    int total_mcsteps = 0;

    int mc_step = -1;

    // Create context
    std::string logFilename;
    if( setupReader.find("SEED") ){
        if( !(setupReader.getValues("SEED").empty()) ){
            logFilename = std::string("log.") + setupReader.getValues("SEED")[0];
        }
    }else{
        logFilename = "x";
    }
    Context *context = new Context(logFilename);

    // Add Worlds
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
            context->useFixmanTorque(worldIx, std::stof(setupReader.getValues("TEMPERATURE_INI")[worldIx]));
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

    // Set thermodynamics
    for(int worldIx = 0; worldIx < setupReader.getValues("WORLDS").size(); worldIx++){
        context->setTemperature(worldIx, std::stof(setupReader.getValues("TEMPERATURE_INI")[worldIx]));
        context->setNofSamplesPerRound(worldIx, std::stoi(setupReader.getValues("SAMPLES_PER_ROUND")[worldIx]));
        context->setNofMDStepsPerSample(worldIx, std::stoi(setupReader.getValues("MDSTEPS")[worldIx]));
    }

    // Set the seeds for reproducibility
    if( setupReader.find("SEED") ){
        if( !(setupReader.getValues("SEED").empty()) ){
            for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){
                context->setSeed(worldIx, 0, std::stoi( setupReader.getValues("SEED")[worldIx] ));
            }
        }
    }

    // Initialize samplers
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
                << " temperature = " << context->updWorld(worldIx)->updSampler(samplerIx)->getTemperature()
                << " initial const state PE: " << std::setprecision(20)
                //<< (context->updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState())
                << (context->updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy(context->updAdvancedState(worldIx, samplerIx))
                << " useFixmanPotential = " << context->updWorld(worldIx)->updSampler(samplerIx)->isUsingFixmanPotential()
                << std::endl;
        }

    }

    // Simulation parameters
    currentWorldIx = 0;
    lastWorldIx = 0;
    round_mcsteps = 0;

    context->setNofRounds(std::stoi(setupReader.getValues("ROUNDS")[0]));

    for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){    
        round_mcsteps += context->getNofSamplesPerRound(worldIx);
    }

    total_mcsteps = round_mcsteps * context->getNofRounds();

    // Set pdb writing frequency
    context->setPdbRestartFreq( std::stoi(setupReader.getValues("WRITEPDBS")[0]) );

    // Get atom indeces for geometry calculations
    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){
            int distanceIx[setupReader.getValues("DISTANCE").size()];
            for(unsigned int i = 0; i < setupReader.getValues("DISTANCE").size(); i++){
                distanceIx[i] = atoi(setupReader.getValues("DISTANCE")[i].c_str());
            }
            // Get distances
            for(int ai = 0; ai < (setupReader.getValues("DISTANCE").size() / 2); ai++){
                context->addDistance(worldIx, 0, distanceIx[2*ai + 0], distanceIx[2*ai + 1]);
            }
    
            // Get dihedrals indeces
            int dihedralIx[setupReader.getValues("DIHEDRAL").size()];
            for(unsigned int i = 0; i < setupReader.getValues("DIHEDRAL").size(); i++){
                dihedralIx[i] = atoi(setupReader.getValues("DIHEDRAL")[i].c_str());
            }
            // Get dihedrals
            for(int ai = 0; ai < (setupReader.getValues("DIHEDRAL").size() / 4); ai++){
                context->addDihedral(worldIx, 0,
                    dihedralIx[4*ai + 0], dihedralIx[4*ai + 1],
                    dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]);
            }
        }
    }

    // wth is this?
    const SimTK::Compound * p_compounds[context->getNofWorlds()];
    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){
            p_compounds[worldIx] = &((context->updWorld(worldIx))->getTopology(0));
        }
    }
    //

    // Helper variables for step arithmetic
    mc_step = -1;

    // Update one round for the first regimen
    currentWorldIx = context->worldIndexes.front();
    SimTK::State& advancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();

    // Update
//    for(int k = 0; k < context->getNofSamplesPerRound(currentWorldIx); k++){
//        ++mc_step; // Increment mc_step
//        context->updWorld(currentWorldIx)->updSampler(0)->update(advancedState, context->getNofMDStepsPerSample(currentWorldIx));
//    }

    // Write pdb
    context->setOutputDir( setupReader.getValues("OUTPUT_DIR")[0] );
    context->setPdbPrefix( setupReader.getValues("MOLECULES")[0]
        + std::to_string(context->updWorld(currentWorldIx)->updSampler(0)->getSeed()) 
        );

    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (context->updWorld(currentWorldIx))->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb(context->getOutputDir(), "/pdbs/sb." + context->getPdbPrefix() + ".", ".pdb", 10, mc_step);
        }
    }


    if(SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT") == ""){
        std::cout << "CUDA_ROOT not set." << std::endl;
    }else{
        std::cout << "CUDA_ROOT set to " 
            << SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT") << std::endl;
    }

    // Get output printing frequency
    context->setPrintFreq( std::stoi(setupReader.getValues("PRINT_FREQ")[0]) );

    // -- Run --
    context->Run(context->getNofRounds(), 
        std::stof(setupReader.getValues("TEMPERATURE_INI")[0]),
        std::stof(setupReader.getValues("TEMPERATURE_FIN")[0]));

    // Write final pdbs
    for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
        ((context->updWorld( context->worldIndexes.front() ))->getTopology(mol_i)).writePdb(context->getOutputDir(), "/pdbs/final." + context->getPdbPrefix() + ".", ".pdb", 10, context->getNofRounds());
    }
    //

    delete context;

} // END MAIN ////////////







