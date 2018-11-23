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
                << " useFixmanPotential = " << context->updWorld(worldIx)->updSampler(0)->isUsingFixmanPotential()
                << std::endl;
        }

    }


    //for(unsigned int worldIx = 0; worldIx < setupReader.getValues("WORLDS").size(); worldIx++){
    //    std::cout << "Addresses" << " World " << << "Samplers adress "<<  << std::endl;
    //}



    int currentWorldIx = 0;
    int lastWorldIx = 0;
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
    int alt_mc_step = 0;
    int restore_mc_step = 0;
    SimTK::Real convFunc = 99.0;

    // Update one round for the first regimen
    currentWorldIx = context->worldIndexes.front();
    SimTK::State& advancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();

    // Update
    //std::cout << "Sampler " << currentWorldIx << " updating initially" << std::endl;
    for(int k = 0; k < context->getNofSamplesPerRound(currentWorldIx); k++){
        ++mc_step; // Increment mc_step
        context->updWorld(currentWorldIx)->updSampler(0)->update(advancedState, context->getNofMDStepsPerSample(currentWorldIx));
    }

    // Randomize structure
    //unsigned int nq = advancedState.getNQ();
    //SimTK::Vector V(nq);
    //advancedState.updQ() = V;
    //std::cout << "TestMixing updQ " << advancedState.getQ() << std::endl << std::flush;
    //

    // Write pdb
    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (context->updWorld(currentWorldIx))->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((context->updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
        }
    }

    // Adapt variables
    bool samplesAdapted = false;
    bool restored = false;
    SimTK::Real cumU = 0.0;
    SimTK::Real avgU = 0.0;
    SimTK::Real prevAvgU = 0.0;

    if(setupReader.getValues("ADAPT_SAMPLE_RATIO")[0] == "TRUE"){
        std::cout << "Adaptive samples ratios ON." << std::endl;
    }

    for(int round = 0; round < context->getNofRounds(); round++){ // Iterate rounds

        if(mc_step == total_mcsteps){
            break;
        }

        for(unsigned int worldIx = 0; worldIx < context->getNofWorlds(); worldIx++){ // Iterate worlds
    
            // Rotate worlds indeces (translate from right to left) 
            std::rotate(context->worldIndexes.begin(), context->worldIndexes.begin() + 1, context->worldIndexes.end());

            currentWorldIx = context->worldIndexes.front();
            lastWorldIx = context->worldIndexes.back();
    
            // Transfer coordinates from last world to current
            //std::cout << "main: Sending configuration from " << context->worldIndexes.back() << " to " << currentWorldIx 
            //    << " at round " << round << std::endl;
            SimTK::State& lastAdvancedState = (context->updWorld(context->worldIndexes.back()))->integ->updAdvancedState();
            SimTK::State& currentAdvancedState = (context->updWorld(currentWorldIx))->integ->updAdvancedState();
    
            currentAdvancedState = (context->updWorld(currentWorldIx))->setAtomsLocationsInGround(
                currentAdvancedState, (context->updWorld(context->worldIndexes.back()))->getAtomsLocationsInGround( lastAdvancedState ));
    
            // Set old potential energy of the new world
            (context->updWorld(currentWorldIx))->updSampler(0)->setOldPE( (context->updWorld(context->worldIndexes.back()))->updSampler(0)->getSetPE() );
    
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
            for(int k = 0; k < context->getNofSamplesPerRound(currentWorldIx); k++){ // Iterate through samples
                context->updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState, context->getNofMDStepsPerSample(currentWorldIx));
                ++mc_step;
                ++restore_mc_step;
    
                // Calculate geomtric features 
                if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
                    // Get distances indeces
                    int distanceIx[setupReader.getValues("DISTANCE").size()];
                    for(unsigned int i = 0; i < setupReader.getValues("DISTANCE").size(); i++){
                        distanceIx[i] = atoi(setupReader.getValues("DISTANCE")[i].c_str());
                    }
                    // Get distances
                    for(int ai = 0; ai < (setupReader.getValues("DISTANCE").size() / 2); ai++){
                        distances[ai] = context->Distance(currentWorldIx, 0, 0, distanceIx[2*ai + 0], distanceIx[2*ai + 1]);
                        std::cout << std::setprecision(4) << distances[ai] << " ";
                    }
    
                    // Get dihedrals indeces
                    int dihedralIx[setupReader.getValues("DIHEDRAL").size()];
                    for(unsigned int i = 0; i < setupReader.getValues("DIHEDRAL").size(); i++){
                        dihedralIx[i] = atoi(setupReader.getValues("DIHEDRAL")[i].c_str());
                    }
                    // Get dihedrals
                    for(int ai = 0; ai < (setupReader.getValues("DIHEDRAL").size() / 4); ai++){
                        dihedrals[ai] = context->Dihedral(currentWorldIx, 0, 0, dihedralIx[4*ai + 0], dihedralIx[4*ai + 1], dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]);
                        std::cout << std::setprecision(4) << dihedrals[ai] << " ";
                    }
                    std::cout << std::endl;
                }else{
                    std::cout << std::endl;
                }

    
            } // END for samples

            // Compute potential statistics
            SimTK::Real U, dU;
            U = context->updWorld(currentWorldIx)->updForceField()->CalcFullPotEnergyIncludingRigidBodies(context->updAdvancedState(currentWorldIx, 0)) ;
            cumU += U ;

            if((!(mc_step % 100))){
                prevAvgU = avgU;
                avgU = cumU / mc_step;
                dU = avgU - prevAvgU;

                if(mc_step != 0 ){
                    convFunc = std::abs( (avgU * std::log(mc_step)) - (prevAvgU * std::log(mc_step - 1)) );
                }

                // Adaptive mixing
                if((setupReader.getValues("ADAPT_SAMPLE_RATIO")[0] == "TRUE")){
                    std::cout << std::setprecision(4) << std::fixed;
                    std::cout << "avgU " << avgU << ' ' << dU << ' ' << convFunc << ' ';
                    if(samplesAdapted == true){ 
                        ++alt_mc_step;
                    }
                    //if( ((dihedrals[0] > -0.698) && (dihedrals[0] < 0.523)) && (samplesAdapted == false) ){
                    if(convFunc < 0.2){
                        if(samplesAdapted == false){
                            alt_mc_step = 0;
                            //std::cout << "nan adapt samples ratio to " << 1 << ' ' << 5 << std::endl;
                            std::cout << " (1.1.2) 2 " << std::endl;
    
                            context->setNofSamplesPerRound(0, 1);
                            context->setNofSamplesPerRound(1, 1);
                            samplesAdapted = true;
                            restored = false;
                            //context->setNofMDStepsPerSample(0, );
                            context->setNofMDStepsPerSample(1, 150);
                        }else{ // samples not adapted
                            std::cout << " (1.1.1) 0 " << std::endl;
                        }
                    }else if(alt_mc_step > 10){ // if dU > 1.0
                        if (restored == false){
                            std::cout << " (1.2.1.2) 1 " << std::endl;
                            restore_mc_step = 0;
                            //std::cout << "nan restore samples ratio to " ;
    
                            for(int wIx = 0; wIx < setupReader.getValues("WORLDS").size(); wIx++){
                                //std::cout << std::stoi(setupReader.getValues("SAMPLES_PER_ROUND")[wIx]) << ' '  ;
                                context->setNofSamplesPerRound(wIx, std::stoi(setupReader.getValues("SAMPLES_PER_ROUND")[wIx]));
                                context->setNofMDStepsPerSample(wIx, std::stoi(setupReader.getValues("MDSTEPS")[wIx]));
                            }
                            restored = true;
                            samplesAdapted = false;
                        }else{ // if restored already
                            std::cout << " (1.2.1.1) 0" << std::endl;
                        }
                    }else{ // if dU > 1 and alt_mc_step < 100 
                        std::cout << " (1.2.2) 0 " << std::endl;
                    }

                }else{ // if ADAPT false
                    std::cout << "avgU " << avgU << ' ' << dU << ' ';
                    std::cout << " (2) 0 " << std::endl;
                }
            } // END if mc_step 100
            // END Adaptive mixing

    
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
            } // if write pdbs
        } // for i in worlds
    } // for i in rounds

    // Write final pdb
    ((context->updWorld(0))->getTopology(0)).writePdb("pdbs", "final", ".pdb", 10, total_mcsteps);

    // Write final pdbs
    for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
        ((context->updWorld( context->worldIndexes.front() ))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, context->getNofRounds());
    }
    //

    delete context;

} // END MAIN ////////////







