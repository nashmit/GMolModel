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
  #endif

  const SimTK::Compound& c = compoundSystem->getCompound(SimTK::CompoundSystem::CompoundIndex(0));
  int arrays_cut = 2 + 4*3*(c.getNumAtoms());

  if((*fassno>0)){ // Get forces from MMTK

    // Apply external forces
    for(int a=0; a<c.getNumAtoms(); a++){
      SimTK::Compound::AtomIndex aIx = (Caller->lig1->bAtomList)[a].atomIndex;
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
  double energy = 0.0;
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
World::World(readAmberInput *amberReader, std::string rbF, std::string flexFN,
std::string ictdF
)
{
  passno = new int;
  vassno = new int;
  fassno = new int;
  sassno = new int;
  sysTimestep = new TARGET_TYPE;
  *sysTimestep = 0.0015;
  //
  this->mol2F = mol2F;
  this->rbF = rbF;
  this->gaffF = gaffF;
  this->frcmodF = frcmodF;
  this->flexFN = flexFN;
  this->ictdF = ictdF;
  this->pyseed = new unsigned long int;
  this->lj14sf = 1;

  system = new SimTK::CompoundSystem;
  matter = new SimTK::SimbodyMatterSubsystem(*system);
  forces = new SimTK::GeneralForceSubsystem(*system);
  decorations = new SimTK::DecorationSubsystem(*system);
  SimTK::Visualizer viz(*system);
  system->addEventReporter( new SimTK::Visualizer::Reporter(viz, 0.0015));
  forceField = new SimTK::DuMMForceFieldSubsystem(*system);
  forceField->loadAmber99Parameters();
  integ = new SimTK::VerletIntegrator(*system);
  ts = new SimTK::TimeStepper(*system, *integ);

  mr = new bMoleculeReader(amberReader, rbF.c_str());

  bAddGaffParams(
    amberReader,
    *forceField,
    mr->bAtomList,
    mr->bonds
  );

}

World::World(
  string mol2F, string rbF, string gaffF, string frcmodF,
  string ictdF 
){
  passno = new int;
  vassno = new int;
  fassno = new int;
  sassno = new int;
  sysTimestep = new TARGET_TYPE;
  *sysTimestep = 0.0015; // Default
  //
  this->mol2F = mol2F;
  this->rbF = rbF;
  this->gaffF = gaffF;
  this->frcmodF = frcmodF;
  this->ictdF = ictdF;
  this->pyseed = new unsigned long int;
  this->lj14sf = 1;

  system = new SimTK::CompoundSystem;
  matter = new SimTK::SimbodyMatterSubsystem(*system);
  forces = new SimTK::GeneralForceSubsystem(*system);
  decorations = new SimTK::DecorationSubsystem(*system);
  SimTK::Visualizer viz(*system);
  system->addEventReporter( new SimTK::Visualizer::Reporter(viz, 0.0015));
  forceField = new SimTK::DuMMForceFieldSubsystem(*system);
  forceField->loadAmber99Parameters();
  integ = new SimTK::VerletIntegrator(*system); // NEW
  ts = new SimTK::TimeStepper(*system, *integ); // NEW

  mr = new bMoleculeReader(
    *forceField,
    mol2F.c_str(),
    "mol2",
    rbF.c_str()
  );
  if (mr->bAtomList == NULL){
    std::cout<<"After bMoleculeReader: NULL bAtomList"<<std::endl;
    exit(1);
  }
  bAddGaffParams(
    *forceField,
    gaffF.c_str(),
    mr->natoms,
    mr->bAtomList,
    mr->nbonds,
    mr->bonds,
    frcmodF.c_str()
  );

}//end of constructor

// Initialize simulation
void World::InitSimulation(TARGET_TYPE extTimestep, bool first_time)
{
  *passno = 0;
  *vassno = 0;
  *fassno = 0;
  *sassno = 0;

  ExtForce = new SimTK::Force::Custom(*forces, new GridForce(system, *matter 
    , fassno
    , this
  ));

  #ifdef TRY_TO_USE_OPENMM
    forceField->setUseOpenMMAcceleration(true);
  #endif
  forceField->setTracing(true); // log OpenMM info to console
  forceField->setNumThreadsRequested(1); // default

  lig1 = new Topology();

  /*
  lig1->init(
    *forceField,
    mr->natoms,
    mr->bAtomList,
    mr->nbonds,
    mr->bonds,
    first_time,
    flexFN,
    ictdF
  );
  */
  lig1->build(
    *forceField,
    mr->natoms,
    mr->bAtomList,
    mr->nbonds,
    mr->bonds,
    flexFN,
    ictdF
  );

  lig1->setSpecificDuMMScaleFactor(*forceField);

  system->adoptCompound(*lig1);
  system->modelCompounds();
  
  TVector = new SimTK::Transform[matter->getNumBodies()];

  arrays_cut = 2 + 4*3*(lig1->natms);
  *sysTimestep = 0.0015; // CHECK
  startT = 300.0; // CHECK

  #ifdef NOSETHERMOS
  thermo = new NoseHooverThermostat(*forces, *matter, startT, (*sysTimestep)*(50)); // every (x*10000)th step
  #endif
  #ifdef VELSTHERMOS
  system->updDefaultSubsystem().addEventHandler(vthermo = new VelocityRescalingThermostat(*system, startT, (*sysTimestep)*100000)); // every (100)th step
  #endif

  system->realizeTopology();


  std::cout << "Number of included atoms in nonbonded interactions: " << forceField->getNumNonbondAtoms() << std::endl;
  std::cout << "getVdwGlobalScaleFactor() " << forceField->getVdwGlobalScaleFactor() << std::endl;
  for(int i=0; i<lig1->natms; i++){
      std::cout << " DuMM VdW Radius " 
          << forceField->getVdwRadius((lig1->bAtomList[i]).getAtomClassIndex()) 
          << " DuMM VdW Well Depth"
          << forceField->getVdwWellDepth((lig1->bAtomList[i]).getAtomClassIndex())
          << std::endl;
  }

  SimTK::State state = system->updDefaultState();

  ts->initialize(system->getDefaultState());
}//end of InitSimulation


// Interface
Topology * World::getTopology(void) const{
    return lig1;
}

Topology * World::updTopology(void){
    return lig1;
}

// Advance
void World::Advance(int nosteps){

  TARGET_TYPE myrealtime=0;
  myrealtime = (TARGET_TYPE)nosteps * (*sysTimestep);
  std::cout<<"myrealtime: "<<myrealtime<<std::endl;

  SimTK::State& advanced = integ->updAdvancedState();
  std::cout<<"Advance start: integ->updAdvancedState: "<< advanced.getQ() <<std::endl;

  ts->stepTo(advanced.getTime() + nosteps*0.0015);

  std::cout<<"Advance stop: integ->updAdvancedState: "<< advanced.getQ() <<std::endl;
  
}

// Destructor
World::~World(){}



