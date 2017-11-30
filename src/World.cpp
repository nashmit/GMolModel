#include "World.hpp"

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  double mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  double mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime)
{
  double mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime)
{
  double mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::PdbStructure pdb, const char *FN)
{
  std::filebuf fb;
  fb.open(FN, std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); //automatically multiplies by ten (nm to A)
  fb.close();
}

void printPossVels(const SimTK::Compound& c, SimTK::State& advanced)
{
    SimTK::Vec3 vertex, vel, acc;
    std::cout<<"MidVV: Positions:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"MidVV: Velocities:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}

void shmDump(TARGET_TYPE *shm, unsigned int natms)
{
  assert(!"Not implemented");
}

////////////////////////////
////// GRID FORCE //////////
////////////////////////////
GridForce::GridForce(SimTK::CompoundSystem *compoundSystem, SimTK::SimbodyMatterSubsystem& matter
                     , int *fassno
                     , World *Caller
                    ) : matter(matter){
  this->compoundSystem = compoundSystem;
  this->fassno = fassno;
  flag = new int;
  this->Caller = Caller;
}

void GridForce::calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
  SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const {

  #ifdef DEBUG_TIME
  boost::timer GridForce_timer;
  printf("GridForce: calcForce START\n");
  #endif

  const SimTK::Compound& c = compoundSystem->getCompound(SimTK::CompoundSystem::CompoundIndex(0));

  if((*fassno>0)){ // Get forces from MMTK

    // Apply external forces
    for(int a=0; a<c.getNumAtoms(); a++){
      SimTK::Compound::AtomIndex aIx = ( (Caller->getTopology(0, 0)).bAtomList )[a].atomIndex;
      const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(c.getAtomMobilizedBodyIndex(aIx));
      SimTK::Vec3 v_check(0.0, 0.0, 0.0);
      mobod.applyForceToBodyPoint(state, c.getAtomLocationInMobilizedBodyFrame(aIx), v_check, bodyForces);
    }

  }

  (*fassno)++;
  #ifdef DEBUG_TIME
  printf("GridForce: calcForce time %.8lf\n", GridForce_timer.elapsed());
  #endif

}


// This should be carefully analyzed. Intended to be taken from somewhere else.
SimTK::Real GridForce::calcPotentialEnergy(const SimTK::State& state) const {
  SimTK::Real energy = 0.0;
  return energy;
}

bool GridForce::dependsOnlyOnPositions() const {
  return true;
}
////////////////////////////
////// END GRID FORCE //////
////////////////////////////


////////////////////////////
////// SYMBODY SYSTEM //////
////////////////////////////
World::World(int worldIndex, bool visual, SimTK::Real visualizerFrequency)
{
    ownWorldIndex = worldIndex;
    std::cout << "World::World BEGIN: ownWorldIndex: " << this->ownWorldIndex << std::endl;
    fassno = new int; // External forces
  
    compoundSystem = new SimTK::CompoundSystem;
    matter = new SimTK::SimbodyMatterSubsystem(*compoundSystem);
    forces = new SimTK::GeneralForceSubsystem(*compoundSystem);

    if(visual){
        decorations = new SimTK::DecorationSubsystem(*compoundSystem);
        SimTK::Visualizer viz(*compoundSystem);
        vizReporter = new SimTK::Visualizer::Reporter(viz, visualizerFrequency);
        compoundSystem->addEventReporter( vizReporter );
    }
  
    forceField = new SimTK::DuMMForceFieldSubsystem(*compoundSystem);
    //forceField->loadAmber99Parameters();
    integ = new SimTK::VerletIntegrator(*compoundSystem);
    ts = new SimTK::TimeStepper(*compoundSystem, *integ);

    moleculeCount = -1;
    std::cout << "World::World END: ownWorldIndex: " << this->ownWorldIndex << std::endl;
}

void World::AddMolecule(readAmberInput *amberReader, std::string rbFN, std::string flexFN, std::string ictdF)
{
    std::cout << "World::AddMolecule BEGIN: ownWorldIndex: " << this->ownWorldIndex << std::endl;
    moleculeCount++; // Used for unique names of molecules

    this->rbFN = rbFN;
    this->frcmodF = frcmodF;
    this->flexFN = flexFN;
    this->ictdF = ictdF;
 
    std::vector<bMoleculeReader> moleculeCopies; // One copy per regimen
    std::vector<Topology> topologyCopies; // One copy per regimen
 
    moleculeCopies.push_back(bMoleculeReader(amberReader, rbFN.c_str()));
    moleculeReaders.push_back(moleculeCopies);

    std::cout << "World::AddMolecule add parameters" << std::endl;
    bAddAllParams(std::string("lig") + std::to_string(moleculeCount), amberReader, *forceField, (moleculeReaders.back().back()).bAtomList, (moleculeReaders.back().back()).bonds);
  
    std::cout << "World::AddMolecule add Compound" << std::endl;
    topologyCopies.push_back(Topology(std::string("lig") + std::to_string(moleculeCount)));
    topologies.push_back(topologyCopies);
  
    std::cout << "World::AddMolecule build Compound" << std::endl;
    (topologies.back().back()).build(*forceField, (moleculeReaders.back().back()).natoms, (moleculeReaders.back().back()).bAtomList, (moleculeReaders.back().back()).nbonds, (moleculeReaders.back().back()).bonds, flexFN, ictdF);
  
    std::cout << "World::AddMolecule adopt Compound " << &(topologies.back().back()) << std::endl;
    compoundSystem->adoptCompound(topologies.back().back());
    std::cout << "World::AddMolecule realizeTopology" << std::endl;
    compoundSystem->realizeTopology();

    std::cout << "Number of included atoms in nonbonded interactions: " << forceField->getNumNonbondAtoms() << std::endl;
    std::cout << "getVdwGlobalScaleFactor() " << forceField->getVdwGlobalScaleFactor() << std::endl;
    for(int i=0; i<(topologies.back().back()).natms; i++){
        std::cout << " DuMM VdW Radius " 
            << forceField->getVdwRadius(((topologies.back().back()).bAtomList[i]).getAtomClassIndex()) 
            << " DuMM VdW Well Depth "
            << forceField->getVdwWellDepth(((topologies.back().back()).bAtomList[i]).getAtomClassIndex())
            << std::endl;
    }
  
    std::cout << "World::AddMolecule END: ownWorldIndex: " << this->ownWorldIndex << std::endl;
}

// Initialize simulation
void World::Init(void)
{
    // Only model after loading all the compounds
    std::cout << "World::Init BEGIN: ownWorldIndex: " << this->ownWorldIndex << std::endl;
    std::cout << "World::Init model Compound" << std::endl;
    compoundSystem->modelCompounds();

    // Generate a nonbonded list
    /*
    forceField->clearIncludedNonbondAtomList();
    compoundSystem->realizeTopology();
    for(SimTK::DuMM::AtomIndex dummAIxA(0); dummAIxA < forceField->getNumAtoms(); ++dummAIxA){
        SimTK::MobilizedBodyIndex mobodIxA = forceField->getAtomBody(dummAIxA);
        SimTK::Vec3 stationA = forceField->getAtomStationOnBody(dummAIxA);
        SimTK::MobilizedBody mobodA = matter->getMobilizedBody(mobodIxA);
 
        for(SimTK::DuMM::AtomIndex dummAIxB(0); dummAIxB < forceField->getNumAtoms(); ++dummAIxB){
            SimTK::MobilizedBodyIndex mobodIxB = forceField->getAtomBody(dummAIxB);
            SimTK::Vec3 stationB = forceField->getAtomStationOnBody(dummAIxB);
            SimTK::MobilizedBody mobodB = matter->getMobilizedBody(mobodIxB);

            SimTK::Real distance = mobodA.calcStationToStationDistance(system->getDefaultState(), stationA, mobodB, stationB);

            if(distance < 8.0){
                forceField->includeNonbondAtom(dummAIxA);
            }
        }
    }
    */
   
    //SimTK::State state = compoundSystem->updDefaultState();
    //ts->initialize(compoundSystem->getDefaultState());
  
    forceField->setVdw12ScaleFactor(0.0);
    forceField->setVdw13ScaleFactor(0.0);
    forceField->setVdw14ScaleFactor(0.5);
    forceField->setVdw15ScaleFactor(1.0);
    forceField->setCoulomb12ScaleFactor(0.0);
    forceField->setCoulomb13ScaleFactor(0.0);
    forceField->setCoulomb14ScaleFactor(0.8333333333);
    forceField->setCoulomb15ScaleFactor(1.0);
    forceField->setVdwMixingRule(SimTK::DuMMForceFieldSubsystem::LorentzBerthelot);
  
    //->setSpecificDuMMScaleFactor(*forceField);
    //forceField->setBondStretchGlobalScaleFactor(0.0);
    //forceField->setBondBendGlobalScaleFactor(0.0);
    //forceField->setBondTorsionGlobalScaleFactor(0.0);
    //forceField->setAmberImproperTorsionGlobalScaleFactor(0.0);
  
    //forceField->setVdw12ScaleFactor(0.0);
    //forceField->setVdw13ScaleFactor(0.0);
    //forceField->setVdw14ScaleFactor(0.0);
    //forceField->setVdw15ScaleFactor(0.0);
    //forceField->setVdwGlobalScaleFactor(0.0);
  
    //forceField->setCoulomb12ScaleFactor(0.0);
    //forceField->setCoulomb13ScaleFactor(0.0);
    //forceField->setCoulomb14ScaleFactor(0.0);
    //forceField->setCoulomb15ScaleFactor(0.0);
    //forceField->setCoulombGlobalScaleFactor(0.0);
  
    forceField->setGbsaGlobalScaleFactor(0.0);
    //
  
    *fassno = 0;
    ExtForce = new SimTK::Force::Custom(*forces, new GridForce(compoundSystem, *matter, fassno, this));
  
    #ifdef TRY_TO_USE_OPENMM
        //forceField->setUseOpenMMAcceleration(true);
    #endif
    //forceField->setTracing(true); // log OpenMM info to console
    forceField->setNumThreadsRequested(1); // don't use this unless
  
    std::cout << "World::Init POINTERS: ownWorldIndex: " << this->ownWorldIndex << ":" << std::endl
        << " System " << &(matter->getSystem()) << std::endl
        << " System::Guts " << &((matter->getSystem()).getSystemGuts()) << std::endl
        << " CompoundSystem " << compoundSystem << std::endl
        << " SimbodyMatterSubsystem " << matter << std::endl
        << " GeneralForceSubsystem " << forces << std::endl
        << " MultibodySubsystem " << &(forces->getMultibodySystem()) << std::endl
        << " DecorationSubsystem " << decorations << std::endl
        << " Visualizer " << &viz << std::endl
        << " vizReporter " << vizReporter << std::endl
        << " DuMMForceFieldSubsystem " << forceField << std::endl
        << " VerletIntegrator " << integ << std::endl
        << " TimeStepper " << ts << std::endl
        << std::endl;

    compoundSystem->realizeTopology();

    std::cout << "World::Init END: ownWorldIndex: " << this->ownWorldIndex << std::endl;
}//end of InitSimulation


// Interface
const Topology& World::getTopology(int moleculeNumber, int moleculeCopy) const{
    return topologies[moleculeNumber][moleculeCopy];
}

Topology& World::updTopology(void){
    assert(!"Not implemented");
}


std::vector<SimTK::Vec3> World::getAtomsLocationsInGround(void)
{
    assert(!"Not implemented");
}

void World::setAtomsLocations(std::vector<SimTK::Vec3>)
{
    assert(!"Not implemented");
}












// Advance
void World::Advance(int nosteps){

  TARGET_TYPE myrealtime=0;
  myrealtime = (TARGET_TYPE)nosteps * (0.0015);
  std::cout<<"myrealtime: "<<myrealtime<<std::endl;

  SimTK::State& advanced = integ->updAdvancedState();
  std::cout<<"Advance start: integ->updAdvancedState: "<< advanced.getQ() <<std::endl;

  ts->stepTo(advanced.getTime() + nosteps*0.0015);

  std::cout<<"Advance stop: integ->updAdvancedState: "<< advanced.getQ() <<std::endl;
  
}

// Destructor
World::~World(){}



