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

    topFNs.push_back(std::vector<std::string>());
    crdFNs.push_back(std::vector<std::string>());
    rbSpecsFNs.push_back(std::vector<std::string>());
    flexSpecsFNs.push_back(std::vector<std::string>());
    regimens.push_back(std::vector<std::string>());

    nofSamplesPerRound.push_back(1);
    nofMDStepsPerSample.push_back(1);
    timesteps.push_back(0.002); // ps
}

// Add another world and a sampler to the context
World * Context::AddWorld(bool visual){

    worldIndexes.push_back(worldIndexes.size());

    World * inp_p_world = new World(worldIndexes.back(), visual);
    worlds.push_back(inp_p_world);

    topFNs.push_back(std::vector<std::string>());
    crdFNs.push_back(std::vector<std::string>());
    rbSpecsFNs.push_back(std::vector<std::string>());
    flexSpecsFNs.push_back(std::vector<std::string>());
    regimens.push_back(std::vector<std::string>());

    nofSamplesPerRound.push_back(1);
    nofMDStepsPerSample.push_back(1);
    timesteps.push_back(0.002); // ps

    return worlds.back();
}

// Add another world and a sampler to the context
World * Context::AddWorld(World *inp_p_world, bool visual){
    worlds.push_back(inp_p_world);
    worldIndexes.push_back(worldIndexes.size());

    topFNs.push_back(std::vector<std::string>());
    crdFNs.push_back(std::vector<std::string>());
    rbSpecsFNs.push_back(std::vector<std::string>());
    flexSpecsFNs.push_back(std::vector<std::string>());
    regimens.push_back(std::vector<std::string>());

    nofSamplesPerRound.push_back(1);
    nofMDStepsPerSample.push_back(1);
    timesteps.push_back(0.002); // ps

    return worlds.back();
}

// Destructor
Context::~Context(){
    worlds.clear();
    worldIndexes.clear();
   
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        topFNs[worldIx].clear();
        crdFNs[worldIx].clear();
        rbSpecsFNs[worldIx].clear();
        flexSpecsFNs[worldIx].clear();
        regimens[worldIx].clear();
    }
    topFNs.clear();
    crdFNs.clear();
    rbSpecsFNs.clear();
    flexSpecsFNs.clear();
    regimens.clear();
 
    nofSamplesPerRound.clear();
    nofMDStepsPerSample.clear();
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

// Input molecular files
bool Context::loadTopologyFile(int whichWorld, int whichMolecule, std::string topologyFilename)
{
    topFNs[whichWorld].push_back(topologyFilename);
}

bool Context::loadCoordinatesFile(int whichWorld, int whichMolecule, std::string coordinatesFilename)
{
    crdFNs[whichWorld].push_back(coordinatesFilename);
}

bool Context::loadRigidBodiesSpecs(int whichWorld, int whichMolecule, std::string RBSpecsFN)
{
    rbSpecsFNs[whichWorld].push_back(RBSpecsFN);
}

bool Context::loadFlexibleBondsSpecs(int whichWorld, int whichMolecule, std::string flexSpecsFN)
{
    flexSpecsFNs[whichWorld].push_back(flexSpecsFN);
}

void Context::setRegimen (int whichWorld, int whichMolecule, std::string regimen)
{
    regimens[whichWorld].push_back(regimen);
}

void Context::loadMolecules(void)
{

    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
        for(unsigned int molIx = 0; molIx < topFNs[worldIx].size(); molIx++){
            readAmberInput *amberReader = new readAmberInput();
            amberReader->readAmberFiles(crdFNs[worldIx][molIx], topFNs[worldIx][molIx]);
            (updWorld(worldIx))->AddMolecule(amberReader,
                rbSpecsFNs[worldIx][molIx], flexSpecsFNs[worldIx][molIx], regimens[worldIx][molIx]);
            delete amberReader; 
        }
    }
}


// Use a SetupReader Object to read worlds information from a file
void Context::LoadWorldsFromSetup(SetupReader& setupReader)
{
    // Build Gmolmodel simulation worlds
    unsigned int nofWorlds = setupReader.getValues("WORLDS").size();
    bool useFixmanPotential = false;

    // Create a vector of world indeces
    /*
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        worldIndexes.push_back(worldIx);
    }
            std::cout << "nofWorlds = " << nofWorlds << "; indeces vector = [" ;
            for(unsigned int Ix = 0; Ix < worldIndexes.size(); Ix++){
                std::cout << " " << worldIndexes[Ix] ;
            }
            std::cout << " ]" << std::endl;
    */

    // Add worlds
    /*
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        if(setupReader.getValues("VISUAL")[0] == "TRUE"){
            TRACE("NEW ALLOC\n");
//r            World * pWorld = new World(worldIx, true);
//r            AddWorld(pWorld, true);
            AddWorld(true);
        }else{
            TRACE("NEW ALLOC\n");
//r            World * pWorld = new World(worldIx, false);
//r            AddWorld(pWorld, false);
            AddWorld(false);
        }

        // Set world identifiers
//r        (updWorld(worldIx))->ownWorldIndex = worldIx;
    
    }
    */

    /*
    // Add molecules to worlds
    std::cout << " Context::loadMolecules " << std::endl << std::flush;
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        // Get input filenames and add molecules
        std::vector<std::string> argValues;
        std::vector<std::string>::iterator argValuesIt;
        argValues = setupReader.getValues("MOLECULES");
        for(unsigned int molIx = 0; molIx < setupReader.getValues("MOLECULES").size(); molIx++){
//r        for(argValuesIt = argValues.begin(); argValuesIt != argValues.end(); ++argValuesIt){
//r            std::string prmtopFN = *argValuesIt + std::string("/ligand.prmtop");
//r            std::string inpcrdFN = *argValuesIt + std::string("/ligand.inpcrd");
//r            std::string rbFN = *argValuesIt + std::string("/ligand.rb");
//r            std::string flexFN = *argValuesIt + std::string("/ligand.flex");
//r            std::string ictd = setupReader.getValues("WORLDS")[worldIx];
    
            loadTopologyFile( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.prmtop") );
            loadCoordinatesFile( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.inpcrd") );
            loadRigidBodiesSpecs( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.rb") );
            loadFlexibleBondsSpecs( worldIx, molIx,
                setupReader.getValues("MOLECULES")[molIx] + std::string("/ligand.flex") );
            setRegimen( worldIx, molIx,
                setupReader.getValues("WORLDS")[worldIx] );
    
            std::cout << "world " << worldIx << " molecule " << molIx 
                << " crdFN " << crdFNs[worldIx][molIx] << " topFN " << topFNs[worldIx][molIx]  
                << " rbSpecsFN " << rbSpecsFNs[worldIx][molIx] << " flexSpecsFN " << flexSpecsFNs[worldIx][molIx]
                << " regimen " << regimens[worldIx][molIx]
                << std::endl << std::flush;

            TRACE("NEW ALLOC\n");
            readAmberInput *amberReader = new readAmberInput();
//r            amberReader->readAmberFiles(inpcrdFN, prmtopFN);
            amberReader->readAmberFiles(crdFNs[worldIx][molIx], topFNs[worldIx][molIx]);

            (updWorld(worldIx))->AddMolecule(amberReader,
                rbSpecsFNs[worldIx][molIx], flexSpecsFNs[worldIx][molIx], regimens[worldIx][molIx]);
    
            delete amberReader;
        }
    } // END Ad molecules
    */


    //
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
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
            (updWorld(worldIx))->ModelTopologies( true );
        }else{
            (updWorld(worldIx))->ModelTopologies( false );
        }

        (updWorld(worldIx))->setTemperature( std::stod(setupReader.getValues("TEMPERATURE")[worldIx]) );

        // Add samplers
        TRACE("NEW ALLOC\n");
    
        // Do we use Fixman potential
        if(setupReader.getValues("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
            useFixmanPotential = true;
        }
     }

    // Set integrators timesteps
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        updWorld(worldIx)->addSampler(HMC);
        setTimestep(worldIx, 0, std::stod(setupReader.getValues("TIMESTEPS")[worldIx]) );
    }
        
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        // Set thermostats
        (updWorld(worldIx))->updSampler(0)->setThermostat(setupReader.getValues("THERMOSTAT")[worldIx]);

        // Initialize samplers
        (updWorld(worldIx))->updSampler(0)->initialize( (updWorld(worldIx))->integ->updAdvancedState(), 
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
int Context::getNofMDStepsPerSample(int whichWorld){
   return nofMDStepsPerSample[whichWorld]; 
}

void Context::setNofMDStepsPerSample(int whichWorld, int MDStepsPerSample)
{
   nofMDStepsPerSample[whichWorld] = MDStepsPerSample;
}

// If HMC, get/set timestep forMD
const float Context::getTimestep(int whichWorld, int whichSampler)
{
    return worlds[whichWorld]->updSampler(whichSampler)->getTimeStepper()->getIntegrator().getPredictedNextStepSize();

}

void Context::setTimestep(int whichWorld, int whichSampler, float argTimestep)
{
    worlds[whichWorld]->updSampler(whichSampler)->updTimeStepper()->updIntegrator().setFixedStepSize(argTimestep);
}
//------------

// --- Mixing parameters ---
// Another way to do it is setting the number of rounds
int Context::getNofRounds(int nofRounds){assert(!"Not implemented");}
void Context::setNofRounds(int nofRounds){assert(!"Not implemented");}

// Get the number of samples returned by the sampler in one round
int Context::getNofSamplesPerRound(int whichWorld)
{
    return nofSamplesPerRound[whichWorld];
}

// Set the number of samples returned by the sampler in one round
void Context::setNofSamplesPerRound(int whichWorld, int MCStepsPerRound)
{
    nofSamplesPerRound[whichWorld] = MCStepsPerRound;
}

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
    // 
    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){    
        std::cout << "World " << worldIx << " initial const state PE: " << std::setprecision(20)
            << (updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState()) 
            << " useFixmanPotential = " << updWorld(worldIx)->updSampler(0)->isUsingFixman()
            << std::endl;
    }

    // Get simulation parameters
    //total_mcsteps = std::stoi(setupReader.getValues("STEPS")[0]);

    //nofSamplesPerRound[worlds.size()];
    //nofMDStepsPerSample[worlds.size()];
    //timesteps[worlds.size()];

    int currentWorldIx = 0;
    int round_mcsteps = 0;

    nofRounds = std::stoi(setupReader.getValues("ROUNDS")[0]);

    for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){    
        //nofSamplesPerRound[worldIx] = std::stoi(setupReader.getValues("MIXMCSTEPS")[worldIx]);
        round_mcsteps += nofSamplesPerRound[worldIx];
        //nofMDStepsPerSample[worldIx] = std::stoi(setupReader.getValues("MDSTEPS")[worldIx]);
        //timesteps[worldIx] = std::stod(setupReader.getValues("TIMESTEPS")[worldIx]);
    }
    total_mcsteps = round_mcsteps * nofRounds;

    // Calculate geometric features fast
    const SimTK::Compound * p_compounds[worlds.size()];
    if(setupReader.getValues("GEOMETRY")[0] == "TRUE"){
        for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
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
        updWorld(currentWorldIx)->updSampler(0)->update(advancedState, nofMDStepsPerSample[currentWorldIx]);
    }

    // Write pdb
    if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
        (updWorld(currentWorldIx))->updateAtomLists(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
            ((updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, mc_step);
        }
    }

    for(int round = 0; round < nofRounds; round++){
        for(unsigned int worldIx = 0; worldIx < worlds.size(); worldIx++){
    
            // Rotate worlds indeces (translate from right to left) 
            std::rotate(worldIndexes.begin(), worldIndexes.begin() + 1, worldIndexes.end());

            currentWorldIx = worldIndexes.front();
    
            // Transfer coordinates from last world to current
            std::cout << "main: Sending configuration from " << worldIndexes.back() << " to " << currentWorldIx 
                << " at round " << round << std::endl;
            SimTK::State& lastAdvancedState = (updWorld(worldIndexes.back()))->integ->updAdvancedState();
            SimTK::State& currentAdvancedState = (updWorld(currentWorldIx))->integ->updAdvancedState();
    
            currentAdvancedState = (updWorld(currentWorldIx))->setAtomsLocationsInGround(
                currentAdvancedState, (updWorld(worldIndexes.back()))->getAtomsLocationsInGround( lastAdvancedState ));
    
            // Reinitialize current sampler
            updWorld(currentWorldIx)->updSampler(0)->reinitialize( currentAdvancedState, 
                SimTK::Real( std::stod(setupReader.getValues("TEMPERATURE")[0]) ) );
    
            // INCORRECT !!
            if(setupReader.getValues("WORLDS")[worldIndexes.back()] == "IC"){
                for(int i = 0; i < worlds.size() - 1; i++){
                    int restIx = worldIndexes[i];
                    int backIx = worldIndexes.back();
                    SimTK::Real diffPE = (updWorld(backIx))->updSampler(0)->getSetPE() - (updWorld(restIx))->updSampler(0)->getSetPE();
                    std::cout << "Setting sampler " << restIx << " REP to " << (updWorld(backIx))->updSampler(0)->getSetPE() << " - " << (updWorld(restIx))->updSampler(0)->getSetPE() << " = " << diffPE << std::endl;
                    (updWorld(restIx))->updSampler(0)->setREP( diffPE );
                }
            }
    
            // Update
            std::cout << "Sampler " << currentWorldIx << " updating " << std::endl;
            for(int k = 0; k < nofSamplesPerRound[currentWorldIx]; k++){
                updWorld(currentWorldIx)->updSampler(0)->update(currentAdvancedState, nofMDStepsPerSample[currentWorldIx]);
            }
    
            // Write pdb
            if(setupReader.getValues("WRITEPDBS")[0] == "TRUE"){
                //if(!((mc_step+1) % 20)){
                if(1){
                    (updWorld(currentWorldIx))->updateAtomLists(currentAdvancedState);
                    //std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
                    for(unsigned int mol_i = 0; mol_i < setupReader.getValues("MOLECULES").size(); mol_i++){
                        ((updWorld(currentWorldIx))->getTopology(mol_i)).writePdb("pdbs", "sb", ".pdb", 10, round);
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
    
    
        } // for i in worlds
    } // for i in rounds


} // END of Run()
//------------

// --- Printing functions ---
void Context::WritePdb(int whichWorld){assert(!"Not implemented");}
void Context::Dihedral(int, int, int, int){assert(!"Not implemented");}
//------------




