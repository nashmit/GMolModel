#include "Context.hpp"

// Includes to get the structure of additional classes

#include "World.hpp"
#include "Sampler.hpp"

// Default constructor
Context::Context(void){
    total_mcsteps = 0;
}

// Constructor
Context::Context(World *inp_p_world){
    total_mcsteps = 0;
    worlds.push_back(inp_p_world);
    worldIndexes.push_back(0);

    nofSamplesPerRound.push_back(1);
    nofMDSteps.push_back(1);
    timesteps.push_back(0.002); // ps
}

// Add another world and a sampler to the context
World * Context::AddWorld(World *inp_p_world, bool visual){
    worlds.push_back(inp_p_world);
    worldIndexes.push_back(0);

    nofSamplesPerRound.push_back(1);
    nofMDSteps.push_back(1);
    timesteps.push_back(0.002); // ps

    return worlds.back();
}

// Destructor
Context::~Context(){
    worlds.clear();
    worldIndexes.clear();
    
    nofSamplesPerRound.clear();
    nofMDSteps.clear();
    timesteps.clear(); 
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

    // Create a vector of world indeces
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        worldIndexes.push_back(worldIx);
    }

    // Iterate thorugh regimens and add worlds to the context
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){
        if(setupReader.getValues("VISUAL")[0] == "TRUE"){
            TRACE("NEW ALLOC\n");
            World * pWorld = new World(worldIx, true, std::stod(setupReader.getValues("TIMESTEPS")[worldIx]));
            AddWorld(pWorld, true);
        }else{
            TRACE("NEW ALLOC\n");
            World * pWorld = new World(worldIx, false, std::stod(setupReader.getValues("TIMESTEPS")[worldIx]));
            AddWorld(pWorld, false);
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

// --- Thermodynamics ---a
// Get/set the main temperature (acc/rej temperature for MC)
float Context::getTemperature(int whichWorld){
    return worlds[whichWorld]->temperature;
}

void  Context::setTemperature(int whichWorld, float someTemperature){
    worlds[whichWorld]->setTemperature(someTemperature);
}

    // If HMC, get/set the guidance Hamiltonian temperature
float Context::getGuidanceTemperature(void)
{
   assert(!"Not implemented"); 
}

void Context::setGuidanceTemperature(float someTemperature)
{
   assert(!"Not implemented"); 
}
//------------

// --- Simulation parameters ---
// If HMC, get/set the number of MD steps
int Context::getNofMDSteps(int whichWorld){
   assert(!"Not implemented"); 
}

void Context::setNofMDSteps(int whichWorld, int nofMDSteps){assert(!"Not implemented");}

// If HMC, get/set timestep forMD
float Context::getTimeStep(int whichWorld){assert(!"Not implemented");}
void Context::setTimeStep(int whichWorld, float timeStep){assert(!"Not implemented");}
//------------

// --- Mixing parameters ---
// Another way to do it is setting the number of rounds
int Context::getNofRounds(int nofRounds){assert(!"Not implemented");}
void Context::setNofRounds(int nofRounds){assert(!"Not implemented");}

int Context::getNofSamplesPerRound(int whichWorld){assert(!"Not implemented");}
void Context::setNofSamplesPerRound(int whichWorld, int MCStepsPerRound){assert(!"Not implemented");}

// Return the world index in position 'which'. To be used when rotationg
int Context::getWorldIndex(int which)
{
    return worldIndexes[which];
}

// --- Arrange different mixing parameters ---
void Context::initializeMixingParamters(void){assert(!"Not implemented");}
//------------

// --- Mix ---
void Context::RotateWorlds(void){assert(!"Not implemented");}
//------------

// -- Main ---
void Context::Run(SetupReader& setupReader)
{
    // Build Gmolmodel simulation worlds
    int nofRegimens = setupReader.getValues("REGIMENS").size();

    // Set convenient names
    for(int worldIx = 0; worldIx < nofRegimens; worldIx++){    
        std::cout << "World " << worldIx << " initial const state PE: " << std::setprecision(20)
            << (updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState()) 
            << " useFixmanPotential = " << updWorld(worldIx)->updSampler(0)->isUsingFixman()
            << std::endl;
    }

    // Get simulation parameters
    total_mcsteps = std::stoi(setupReader.getValues("STEPS")[0]);
    nofSamplesPerRound[nofRegimens];
    nofMDSteps[nofRegimens];
    timesteps[nofRegimens];


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
            p_compounds[worldIx] = &((updWorld(worldIx))->getTopology(0));
        }
    }
    //

    // Simulate the two worlds
    int mc_step = -1;

    // Update one round for the first regimen
    currentWorldIx = worldIndexes.front();
    SimTK::State& advancedState = (updWorld(currentWorldIx))->integ->updAdvancedState();

    // Update
    //std::cout << "Sampler " << currentWorldIx << " updating initially" << std::endl;
    for(int k = 0; k < nofSamplesPerRound[currentWorldIx]; k++){
        ++mc_step; // Increment mc_step
        updWorld(currentWorldIx)->updSampler(0)->update(advancedState, timesteps[currentWorldIx], nofMDSteps[currentWorldIx]);
    }

    // Write pdb
    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (updWorld(currentWorldIx))->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
        }
    }

    while(mc_step < total_mcsteps){

        // Rotate worlds indeces (translate from right to left) 
        std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());
        currentWorldIx = worldIndexes.front();

        // Transfer coordinates from last world to current
        //std::cout << "main: Sending configuration from " << worldIndexes.back() << " to " << currentWorldIx 
        //    << " at MC step " << mc_step << std::endl;
        SimTK::State& lastAdvancedState = (updWorld(worldIndexes.back()))->integ->updAdvancedState();
        SimTK::State& currentAdvancedState = (updWorld(currentWorldIx))->integ->updAdvancedState();

        currentAdvancedState = (updWorld(currentWorldIx))->setAtomsLocationsInGround(
            currentAdvancedState, (updWorld(worldIndexes.back()))->getAtomsLocationsInGround( lastAdvancedState ));

        // Reinitialize current sampler
        updWorld(currentWorldIx)->updSampler(0)->reinitialize( currentAdvancedState, 
            timesteps[currentWorldIx], total_mcsteps,
            SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );

        // Set residual energy for the rest of the samplers
        if(setupReader.getValues("REGIMENS")[worldIndexes.back()] == "IC"){
            for(int i = 0; i < nofRegimens - 1; i++){
                int restIx = worldIndexes[i];
                int backIx = worldIndexes.back();
                SimTK::Real diffPE = (updWorld(backIx))->updSampler(0)->getSetPE() - (updWorld(restIx))->updSampler(0)->getSetPE();
                std::cout << "Setting sampler " << restIx << " REP to " << (updWorld(backIx))->updSampler(0)->getSetPE() << " - " << (updWorld(restIx))->updSampler(0)->getSetPE() << " = " << diffPE << std::endl;
                (updWorld(restIx))->updSampler(0)->setREP( diffPE );
            }
        }

        // Update
        //std::cout << "Sampler " << currentWorldIx << " updating " << std::endl;
        for(int k = 0; k < nofSamplesPerRound[currentWorldIx]; k++){
            ++mc_step; // Increment mc_step
            updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState, timesteps[currentWorldIx], nofMDSteps[currentWorldIx]);
        }

        // Write pdb
        if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
            //if(!((mc_step+1) % 20)){
            if(1){
                (updWorld(currentWorldIx))->updateAtomLists(currentAdvancedState);
                std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
                for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                    ((updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
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


} // END of Run()
//------------

// --- Printing functions ---
void Context::WritePdb(int whichWorld){assert(!"Not implemented");}
void Context::Dihedral(int, int, int, int){assert(!"Not implemented");}
//------------




