#include "World.hpp"
#include "MidVVIntegrator.hpp"

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
  unsigned long int shm_len = 2 + 3*4*natms + 11;
  printf("Shm dump:");
  printf("%.1f, %.1f\n", shm[0], shm[1]);
  for(int a=2; a<shm_len; a++){
    printf("shm[%3d]= %7.2f ", a, shm[a]);
    if((a-1)%3 == 0){
      printf("\n");
    }
  }
  printf("\n");
}

////////////////////////////
////// GRID FORCE //////////
////////////////////////////
GridForce::GridForce(SimTK::CompoundSystem *compoundSystem, SimTK::SimbodyMatterSubsystem& matter
                     , TARGET_TYPE **indexMap, TARGET_TYPE *PrmToAx_po, TARGET_TYPE *MMTkToPrm_po
                     , TARGET_TYPE **coords, TARGET_TYPE **vels, TARGET_TYPE **grads
                     , int *fassno
                     , TARGET_TYPE *shm
                     , World *Caller
                    ) : matter(matter){
  this->indexMap = Caller->indexMap;
  this->PrmToAx_po = PrmToAx_po;
  this->MMTkToPrm_po = MMTkToPrm_po;
  this->coords = Caller->coords;
  this->vels = Caller->vels;
  this->grads = Caller->grads;
  this->compoundSystem = compoundSystem;
  this->fassno = fassno;
  flag = new int;
  // From MMTK:
  //this->pyFFEvaluatorObject = Caller->pyFFEvaluatorObject;
  //this->p_energy_po = Caller->p_energy_po;
  //this->configuration = Caller->configuration;
  //this->universe_spec = Caller->universe_spec;
  this->shm = shm;
  this->Caller = Caller;
}

void GridForce::calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
  SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const {

  #ifdef DEBUG_TIME
  boost::timer GridForce_timer;
  #endif

  const SimTK::Compound& c = compoundSystem->getCompound(SimTK::CompoundSystem::CompoundIndex(0));
  int arrays_cut = 2 + 4*3*(c.getNumAtoms());

  float conversion_factor = -1; // Gradients to forces
  int natoms3 = 3*(c.getNumAtoms());
  int i=0;
  int tx, tshm;
  Vector3 *x, *v, *f;
  int ix = 0;
  
  if((*fassno>0)){ // Get forces from MMTK

    // Apply forces from MMTK
    for(int a=0; a<c.getNumAtoms(); a++){
      SimTK::Compound::AtomIndex aIx(Caller->_indexMap[ a ][1]);
      const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(c.getAtomMobilizedBodyIndex(aIx));
      // LAUR
      SimTK::Vec3 v_check(0, 0, -10.0);
      //====
      mobod.applyForceToBodyPoint(state, c.getAtomLocationInMobilizedBodyFrame(aIx), v_check, bodyForces);
    }

    // LAUR
    shm[arrays_cut + 5] = 0.0;
    //====

  } // * fassno % nosteps

  (*fassno)++;
  #ifdef DEBUG_TIME
  printf("GridForce: calcForce time %.8lf\n", GridForce_timer.elapsed());
  #endif

}

SimTK::Real GridForce::calcPotentialEnergy(const SimTK::State& state) const {
  const SimTK::Compound& c = compoundSystem->getCompound(SimTK::CompoundSystem::CompoundIndex(0));
  // LAUR
  double energy = 0.0;
  //====
  return energy;
}

bool GridForce::dependsOnlyOnPositions() const {
  return true;
}
////////////////////////////
////// END GRID FORCE //////
////////////////////////////


/////////////////////////////
//////// MIDVV INTEGRATOR ///
/////////////////////////////
// LAUR
//#include "MidVVIntegrator.cpp"
// ====
/////////////////////////////
//////// MIDVV INTEGRATOR ///
/////////////////////////////


////////////////////////////
////// SYMBODY SYSTEM //////
////////////////////////////

World::World(
  string mol2F, string rbF, string gaffF, string frcmodF,
  string ictdF, TARGET_TYPE *PrmToAx_po, TARGET_TYPE *MMTkToPrm_po,
  TARGET_TYPE *shm
){
  std::cout<<"hop 1"<<std::endl<<std::flush;
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
  this->PrmToAx_po = PrmToAx_po;
  this->MMTkToPrm_po = MMTkToPrm_po;
  this->shm = shm;
  this->pyseed = new unsigned long int;
  this->lj14sf = 0.5;

  system = new SimTK::CompoundSystem;
  matter = new SimTK::SimbodyMatterSubsystem(*system);
  forces = new SimTK::GeneralForceSubsystem(*system);
  decorations = new SimTK::DecorationSubsystem(*system);

  // LAUR
  SimTK::Visualizer viz(*system);
  system->addEventReporter( new SimTK::Visualizer::Reporter(viz, 0.0015));
  // ====

  forceField = new SimTK::DuMMForceFieldSubsystem(*system);
  forceField->loadAmber99Parameters();

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
    mr->natms,
    mr->bAtomList,
    mr->nbnds,
    mr->bonds,
    frcmodF.c_str()
  );

}//end of constructor

void World::InitSimulation(
  TARGET_TYPE **coords,
  TARGET_TYPE **vels,
  TARGET_TYPE **inivels,
  TARGET_TYPE **indexMap,
  TARGET_TYPE **grads,
  TARGET_TYPE extTimestep,
  bool first_time
){
  this->coords = coords;
  this->vels = vels;
  this->inivels = inivels;
  this->indexMap = indexMap;
  this->grads = grads;

  *passno = 0;
  *vassno = 0;
  *fassno = 0;
  *sassno = 0;

  ExtForce = new SimTK::Force::Custom(*forces, new GridForce(system, *matter, indexMap, PrmToAx_po, MMTkToPrm_po, coords, vels, grads, fassno
    , shm
    , this
  ));

  #ifdef TRY_TO_USE_OPENMM
    forceField->setUseOpenMMAcceleration(true);
  #endif
  forceField->setTracing(true); // log OpenMM info to console
  forceField->setNumThreadsRequested(1); // default
  forceField->setBondStretchGlobalScaleFactor(0);
  forceField->setBondBendGlobalScaleFactor(0);
  forceField->setBondTorsionGlobalScaleFactor(0);
  forceField->setAmberImproperTorsionGlobalScaleFactor(0);
  forceField->setVdw12ScaleFactor(0);
  forceField->setVdw13ScaleFactor(0);
  forceField->setVdw14ScaleFactor(0);
  forceField->setVdw15ScaleFactor(0);
  forceField->setVdwGlobalScaleFactor(0);
  forceField->setCoulomb12ScaleFactor(0);
  forceField->setCoulomb13ScaleFactor(0);
  forceField->setCoulomb14ScaleFactor(0);
  forceField->setCoulomb15ScaleFactor(0);
  forceField->setCoulombGlobalScaleFactor(0);
  forceField->setGbsaGlobalScaleFactor(0);

  lig1 = new Topology();

  lig1->init(
    *forceField,
    mr->natms,
    mr->bAtomList,
    mr->nbnds,
    mr->bonds,
    coords,
    indexMap,
    PrmToAx_po,
    MMTkToPrm_po,
    first_time,
    ictdF
  );
  system->adoptCompound(*lig1);
  system->modelCompounds();
  
  // Alloc _indexMap[int]
  _indexMap = new int*[lig1->natms];
  for(int t=0; t<lig1->natms; t++){
    _indexMap[t] = new int[3];
  }

  // Assign _indexMap
  for(int t=0; t<lig1->natms; t++){
    _indexMap[t][0] = int(indexMap[t][0]);
    _indexMap[t][1] = int(indexMap[t][1]);
    _indexMap[t][2] = int(indexMap[t][2]);
  }

  QVector = new TARGET_TYPE*[matter->getNumBodies()-1];
  for(int i=0; i<(matter->getNumBodies()-1); i++){
    QVector[i] = new TARGET_TYPE[3];
  }
  for(int i=0; i<(matter->getNumBodies()-1); i++){
    QVector[i][0] = QVector[i][1] = QVector[i][2] = 0.0;
  }

  TVector = new SimTK::Transform[matter->getNumBodies()];

  // Check versus Velocity Verlet in cart coords
  double *vv_vals = new double[27];

  arrays_cut = 2 + 4*3*(lig1->natms);
  int step = (int)(round(shm[arrays_cut + 0]));
  int nosteps = (int)(round(shm[arrays_cut + 1]));
  *sysTimestep = shm[arrays_cut + 3];
  startT = shm[arrays_cut + 2];

  #ifdef NOSETHERMOS
  thermo = new NoseHooverThermostat(*forces, *matter, startT, (*sysTimestep)*(50)); // every (x*10000)th step
  #endif
  #ifdef VELSTHERMOS
  system->updDefaultSubsystem().addEventHandler(vthermo = new VelocityRescalingThermostat(*system, startT, (*sysTimestep)*100000)); // every (100)th step
  #endif

  system->realizeTopology();

  SimTK::State state = system->updDefaultState();

  // * Build a symetric matrix that describes the mobods indeces tree * //
  // Allocate
  int nmbx = matter->getNumBodies();
  int nmbx_1 = nmbx + 1;
  branchMassVec = new SimTK::Real[nmbx];
  mbxTreeMat = new int*[nmbx];
  for(int i=0; i<(nmbx); i++){
    mbxTreeMat[i] = new int[nmbx];
  }
  // Fill first row & col
  for(int i=0; i<nmbx; i++){for(int j=0; j<nmbx; j++){mbxTreeMat[i][j] = 0;}}
  for(int j=0; j<nmbx; j++){
    mbxTreeMat[j][0] = j;
  }
  // Fill the rest
  for(int i=1; i<nmbx; i++){ // rows
    for(int j=0; j<nmbx-1; j++){ // cols
      if(mbxTreeMat[i][j] != 0){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(mbxTreeMat[i][j]));
        mbxTreeMat[i][j+1] = (mobod.getParentMobilizedBody()).getMobilizedBodyIndex();
      }
    }
  }

  // Branches
  std::set<int> branches[nmbx];
  std::set<int>::iterator branchesIt;
  std::pair<std::set<int>::iterator,bool> branchesRet;
  int flag;
  // Fill branch mbxes
  for(int b=0; b<(nmbx); b++){ // branches
    flag = 0;
    for(int i=0; i<(nmbx); i++){ // rows ------------
      for(int j=0; j<nmbx-1; j++){ // cols flag
        if(mbxTreeMat[i][j] == b){
          flag = 1;
          break;
        }
      }
      if(flag == 1){ // mobod found on this row (i)
        for(int j=0; j<nmbx-1; j++){ // cols
          if(mbxTreeMat[i][j] == 0){
            flag = 0;
            break;
          }
          branches[b].insert(mbxTreeMat[i][j]);
          if(mbxTreeMat[i][j] == b){
            flag = 0;
            break;
          }
        }
      }
    } // rows ----------------------------------------
  }
  // Fill branches masses
  for(int i=0; i<nmbx; i++){ // rows
    branchMassVec[i] = 0;
    for (branchesIt=(branches[i]).begin(); branchesIt!=(branches[i]).end(); ++branchesIt){
      if(*branchesIt){
        branchMassVec[i] += (matter->getMobilizedBody(SimTK::MobilizedBodyIndex(*branchesIt))).getBodyMass(state);
      }
    }
  }

  //integ = new MidVVIntegrator(*system, (*sysTimestep), PrmToAx_po, MMTkToPrm_po, system, this);
  integ = new SimTK::VerletIntegrator(*system);
  ts = new SimTK::TimeStepper(*system, *integ);
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
  //myrealtime = (TARGET_TYPE)nosteps * (*sysTimestep);
  myrealtime = (TARGET_TYPE)nosteps * (shm[arrays_cut + 3]);
  std::cout<<"myrealtime: "<<myrealtime<<std::endl;

  SimTK::State& advanced = integ->updAdvancedState();
  std::cout<<"Advance start: integ->updAdvancedState: "<< advanced.getQ() <<std::endl;

/*
  integ->dropBeginFlag();
  integ->dropStep0Flag();
  integ->dropMetroFlag();
  integ->resetStep();
  integ->resetTotStepsInCall();

  integ->setFixedStepSize(shm[arrays_cut + 3]); // EU
  std::cout<<"Advance: step size: "<<shm[arrays_cut + 3]<<std::endl;

  integ->setTimeToReach(advanced.getTime() + myrealtime);
  std::cout<<"Advance: integ->getTimeToReach: "<<integ->getTimeToReach()<<std::endl;

  integ->setNoSteps((int)shm[arrays_cut + 1]);
  std::cout<<"Advance: integ->getNoSteps: "<<integ->getNoSteps()<<std::endl;

  integ->setStepsPerTrial(int(round(shm[arrays_cut + 9])));
  std::cout<<"Advance: integ->getStepsPerTrial: "<<integ->getStepsPerTrial()<<std::endl;

  std::cout<<"int(shm[arrays_cut + 1]): "<<int(shm[arrays_cut + 1])
    <<" int(shm[arrays_cut + 9]): "<<int(shm[arrays_cut + 9])<<std::endl;
  integ->setNtrials(int(shm[arrays_cut + 1]) / int(shm[arrays_cut + 9]));
  std::cout<<" int(shm[arrays_cut + 1]) / int(shm[arrays_cut + 9]): "
    <<int(shm[arrays_cut + 1]) / int(shm[arrays_cut + 9])<<std::endl;
  std::cout<<"Advance: integ->getNtrials: "<<integ->getNtrials()<<std::endl;

  integ->setTrial(int(shm[arrays_cut + 10]));
  std::cout<<"Advance: integ->getTrial: "<<integ->getTrial()<<std::endl;

  integ->setMassMatNumOpt(this->massMatNumOpt); // EU
  std::cout<<"Advance: integ->getMassMatNumOpt: "<<integ->getMassMatNumOpt()<<std::endl;

  integ->setMetroFixmanOpt(this->metroFixmanOpt); // EU
  std::cout<<"Advance: integ->getMetroFixmanOpt: "<<integ->getMetroFixmanOpt()<<std::endl;

  #ifdef DEBUG_TIME
  printf("Advance: integ time %.8lf\n", Advance_timer.elapsed());
  #endif
  
  #ifdef DEBUG_SPECIFIC
  //printf("ts time advanced.getTime() %lf myrealtime %lf stepTo %lf\n", advanced.getTime(), myrealtime, advanced.getTime() + myrealtime);
  #endif
*/
  ts->stepTo(advanced.getTime() + myrealtime);

  std::cout<<"Advance stop: integ->updAdvancedState: "<< advanced.getQ() <<std::endl;
  
}

// Destructor
World::~World(){}
////////////////////////////
////// END SYMBODY SYSTEM //
////////////////////////////



