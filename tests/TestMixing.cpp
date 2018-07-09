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

    srand (time(NULL));

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

	
	context -> updForceField( worldIx ) -> setVdw12ScaleFactor( float(setupReader.getValues("setVdw12ScaleFactor")[worldIx] ) )
	context -> updForceField( worldIx ) -> setVdw13ScaleFactor( float(setupReader.getValues("setVdw13ScaleFactor")[worldIx] ) )
	context -> updForceField( worldIx ) -> setVdw14ScaleFactor( float(setupReader.getValues("setVdw14ScaleFactor")[worldIx] ) )
	context -> updForceField( worldIx ) -> setVdw15ScaleFactor( float(setupReader.getValues("setVdw15ScaleFactor")[worldIx] ) )
	
	context -> updForceField( worldIx ) -> setCoulomb12ScaleFactor( float(setupReader.getValues("setCoulomb12ScaleFactor")[worldIx] ) )
	context -> updForceField( worldIx ) -> setCoulomb13ScaleFactor( float(setupReader.getValues("setCoulomb13ScaleFactor")[worldIx] ) )
	context -> updForceField( worldIx ) -> setCoulomb14ScaleFactor( float(setupReader.getValues("setCoulomb14ScaleFactor")[worldIx] ) )
	context -> updForceField( worldIx ) -> setCoulomb15ScaleFactor( float(setupReader.getValues("setCoulomb15ScaleFactor")[worldIx] ) )





        // Set world GBSA scale factor
        context->setGbsaGlobalScaleFactor(worldIx, std::stod(setupReader.getValues("GBSA")[worldIx]));
    }

    

    // To be removed
    context->LoadWorldsFromSetup(setupReader);

    for(int worldIx = 0; worldIx < setupReader.getValues("WORLDS").size(); worldIx++){
        context->setTemperature(worldIx, std::stof(setupReader.getValues("TEMPERATURE")[worldIx]));

        context->setNofSamplesPerRound(worldIx, std::stoi(setupReader.getValues("MIXMCSTEPS")[worldIx]));
        context->setNofMDStepsPerSample(worldIx, std::stoi(setupReader.getValues("MDSTEPS")[worldIx]));
    }

    context->Run(setupReader);



    delete context;

} // END MAIN ////////////



