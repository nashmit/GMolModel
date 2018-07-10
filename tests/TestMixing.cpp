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

/*
        // Set force field scale factors.
        if(setupReader.getValues("FFSCALE")[worldIx] == "AMBER"){
            context->setAmberForceFieldScaleFactors(worldIx);
        }else{
            context->setGlobalForceFieldScaleFactor(worldIx, std::stod(setupReader.getValues("FFSCALE")[worldIx]));
        }

*/	
	
//	context -> updForceField( worldIx ) -> setAllGlobalScaleFactors( std::stof(setupReader.getValues("setAllGlobalScaleFactors")[worldIx] ) );

	context -> updForceField( worldIx ) -> setVdwGlobalScaleFactor( std::stof(setupReader.getValues("setVdwGlobalScaleFactor")[worldIx] ) );

	context -> updForceField( worldIx ) -> setVdw12ScaleFactor( std::stof(setupReader.getValues("setVdw12ScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setVdw13ScaleFactor( std::stof(setupReader.getValues("setVdw13ScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setVdw14ScaleFactor( std::stof(setupReader.getValues("setVdw14ScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setVdw15ScaleFactor( std::stof(setupReader.getValues("setVdw15ScaleFactor")[worldIx] ) );
	

	context -> updForceField( worldIx ) -> setCoulombGlobalScaleFactor( std::stof(setupReader.getValues("setCoulombGlobalScaleFactor")[worldIx] ) );


	context -> updForceField( worldIx ) -> setCoulomb12ScaleFactor( std::stof(setupReader.getValues("setCoulomb12ScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setCoulomb13ScaleFactor( std::stof(setupReader.getValues("setCoulomb13ScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setCoulomb14ScaleFactor( std::stof(setupReader.getValues("setCoulomb14ScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setCoulomb15ScaleFactor( std::stof(setupReader.getValues("setCoulomb15ScaleFactor")[worldIx] ) );


	context -> updForceField( worldIx ) -> setGbsaGlobalScaleFactor( std::stof(setupReader.getValues("setGbsaGlobalScaleFactor")[worldIx] ) );

	context -> updForceField( worldIx ) -> setBondStretchGlobalScaleFactor( std::stof(setupReader.getValues("setBondStretchGlobalScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setBondBendGlobalScaleFactor( std::stof(setupReader.getValues("setBondBendGlobalScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setBondTorsionGlobalScaleFactor( std::stof(setupReader.getValues("setBondTorsionGlobalScaleFactor")[worldIx] ) );
	context -> updForceField( worldIx ) -> setAmberImproperTorsionGlobalScaleFactor( std::stof(setupReader.getValues("setAmberImproperTorsionGlobalScaleFactor")[worldIx] ) );

//	context -> updForceField( worldIx ) -> setCustomBondStretchGlobalScaleFactor( std::stof(setupReader.getValues("setCustomBondStretchGlobalScaleFactor")[worldIx] ) );
//	context -> updForceField( worldIx ) -> setCustomBondBendGlobalScaleFactor( std::stof(setupReader.getValues("setCustomBondBendGlobalScaleFactor")[worldIx] ) );
//	context -> updForceField( worldIx ) -> setCustomBondTorsionGlobalScaleFactor( std::stof(setupReader.getValues("setCustomBondTorsionGlobalScaleFactor")[worldIx] ) );




        // Set world GBSA scale factor
        // context->setGbsaGlobalScaleFactor(worldIx, std::stod(setupReader.getValues("GBSA")[worldIx]));
    }

<<<<<<< HEAD
    
=======
    // Model topologies
    context->modelTopologies();
>>>>>>> origin/master

    // To be removed
    context->LoadWorldsFromSetup(setupReader);

    if(setupReader.getValues("REPRODUCIBLE")[0] == "TRUE"){
        context->setReproducible();
        srand (0);
    }else{
        srand (time(NULL));
    }

    for(int worldIx = 0; worldIx < setupReader.getValues("WORLDS").size(); worldIx++){
        context->setTemperature(worldIx, std::stof(setupReader.getValues("TEMPERATURE")[worldIx]));

        context->setNofSamplesPerRound(worldIx, std::stoi(setupReader.getValues("MIXMCSTEPS")[worldIx]));
        context->setNofMDStepsPerSample(worldIx, std::stoi(setupReader.getValues("MDSTEPS")[worldIx]));
    }

    context->Run(setupReader);



    delete context;

} // END MAIN ////////////



