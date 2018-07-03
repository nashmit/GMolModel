#include "Context.hpp"

// Includes to get the structure of additional classes

#include "World.hpp"
#include "Sampler.hpp"

// Default constructor
Context::Context(void){
}

// Constructor
Context::Context(World *inp_p_world){
    worlds.push_back(inp_p_world);
}

// Add another world and a sampler to the context
World * Context::AddWorld(World *inp_p_world){
    worlds.push_back(inp_p_world);
    return worlds.back();
}

// Destructor
Context::~Context(){
    worlds.clear();
}

// Get world
World * Context::getWorld(void) const{
    return worlds.back();
}

// Get a specific world
World * Context::getWorld(int which) const{
    return worlds[which];
}

// Get the last mutable world
World * Context::updWorld(void){
    return worlds.back();
}

// Get a mutable specific world
World * Context::updWorld(int which){
    return worlds[which];
}

// Use a SetupReader Object to read worlds information from a file
void Context::LoadWorldsFromSetup(SetupReader& setupReader)
{
    // Build Gmolmodel simulation worlds
    int nofRegimens = setupReader.getValues("REGIMENS").size();
    std::vector<int> worldIndexes;

    // Create a vector of world indeces
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        worldIndexes.push_back(worldIx);
    }

    // Iterate thorugh regimens and add worlds to the context
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        if(setupReader.getValues("VISUAL")[0] == "TRUE"){
            TRACE("NEW ALLOC\n");
            World * pWorld = new World(worldIx, true, std::stod(setupReader.getValues("TIMESTEPS")[worldIx]));
            AddWorld(pWorld);
        }else{
            TRACE("NEW ALLOC\n");
            World * pWorld = new World(worldIx, false, std::stod(setupReader.getValues("TIMESTEPS")[worldIx]));
            AddWorld(pWorld);
        }

        // Set world identifiers
        (updWorld(worldIx))->ownWorldIndex = worldIx;
    
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
            (updWorld(worldIx))->AddMolecule(amberReader, rbFN, flexFN, ictd);
    
            delete amberReader;
        }
   
        // Set force field scale factors.
        if(setupReader.getValues("FFSCALE")[worldIx] == "AMBER"){ 
            (updWorld(worldIx))->setAmberForceFieldScaleFactors();
        }else{
            (updWorld(worldIx))->setGlobalForceFieldScaleFactor(std::stod(setupReader.getValues("FFSCALE")[worldIx]));
        }

        // Set world GBSA implicit solvent scale factor
        (updWorld(worldIx))->setGbsaGlobalScaleFactor(std::stod(setupReader.getValues("GBSA")[worldIx]));
    
        // Initialize worlds
        if(setupReader.getValues("FIXMAN_TORQUE")[worldIx] == "TRUE"){
            (updWorld(worldIx))->Init( std::stod(setupReader.getValues("TIMESTEPS")[worldIx]), true );
        }else{
            (updWorld(worldIx))->Init( std::stod(setupReader.getValues("TIMESTEPS")[worldIx]), false );
        }

        (updWorld(worldIx))->setTemperature( std::stod(setupReader.getValues("TEMPERATURE")[worldIx]) );

        // Add samplers
        TRACE("NEW ALLOC\n");

        updWorld(worldIx)->addSampler(HMC);
    
        // Do we use Fixman potential
        bool useFixmanPotential = false;
        if(setupReader.getValues("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
            useFixmanPotential = true;
        }

        // Set thermostats
        (updWorld(worldIx))->updSampler(0)->setThermostat(setupReader.getValues("THERMOSTAT")[worldIx]);

        // Initialize samplers
        (updWorld(worldIx))->updSampler(0)->initialize( (updWorld(worldIx))->integ->updAdvancedState(), 
             std::stod(setupReader.getValues("TIMESTEPS")[worldIx]),
             std::stoi(setupReader.getValues("STEPS")[0]),
             SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[worldIx]) ),
             useFixmanPotential );
    
    }
    ////////////////////////////////////////////////////////// END creating Worlds

}







