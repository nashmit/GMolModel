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

void World::printPoss(const SimTK::Compound& c, SimTK::State& advanced)
{
    SimTK::Vec3 vertex;
    std::cout<<"Positions:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
    }
}

void printVels(const SimTK::Compound& c, SimTK::State& advanced)
{
    SimTK::Vec3 vel;
    std::cout<<"Velocities:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
}

void printPossVels(const SimTK::Compound& c, SimTK::State& advanced)
{
    SimTK::Vec3 vertex, vel;
    std::cout<<"Positions:"<<std::endl;
    for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
        vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
        std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
          <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"Velocities:"<<std::endl;
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
  this->Caller = Caller;
}

//void GridForce::calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
//  SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const {
//
//  #ifdef DEBUG_TIME
//  boost::timer GridForce_timer;
//  printf("GridForce: calcForce START\n");
//  #endif
//
//  const SimTK::Compound& c = compoundSystem->getCompound(SimTK::CompoundSystem::CompoundIndex(0));
//
//  if((*fassno>0)){ // Get forces from MMTK
//
//    // Apply external forces
//    for(int a=0; a<c.getNumAtoms(); a++){
//      SimTK::Compound::AtomIndex aIx = ( (Caller->getTopology(0)).bAtomList )[a].atomIndex;
//      const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(c.getAtomMobilizedBodyIndex(aIx));
//      SimTK::Vec3 v_check(0.0, 0.0, 0.0);
//      mobod.applyForceToBodyPoint(state, c.getAtomLocationInMobilizedBodyFrame(aIx), v_check, bodyForces);
//    }
//
//  }
//
//  (*fassno)++;
//  #ifdef DEBUG_TIME
//  printf("GridForce: calcForce time %.8lf\n", GridForce_timer.elapsed());
//  #endif
//
//}


void GridForce::calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                           SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const
{
    // Compute Fixman torque
    int nu = state.getNU();
    SimTK::Vector V3(nu);
    SimTK::Vector V4(nu);
    SimTK::Real* D0 = new SimTK::Real(1.0);
    matter.calcFixmanTorque(state, V3, V4, D0);
    //std::cout << "GridForce Fixman torque " ;
    //for(int i = 0; i < nu; i++){
    //    std::cout << std::setprecision(10) << V4[i] << " ";
    //}
    //std::cout << std::endl;
    delete D0;
    // end - Compute Fixman torque

    //std::cout << "Applied " ;
    int uslot = -1;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter.getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

        for(int k = 0; k < mobod.getNumU(state); k++){
            uslot++;
            mobod.applyOneMobilityForce(state, k, V4[uslot], mobilityForces);
            //std::cout << " " << std::setprecision(10) << std::fixed << V4[uslot] << " to " << int(mbx) ;
        }

    }
    //std::cout << std::endl;

    //const SimTK::Real q = knee.getOneQ(state, 0);
    //const SimTK::Real x = q < low ? q-low : (q > high ? q-high : 0);
    //knee.applyOneMobilityForce(state, 0, -k*x, mobilityForces);
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
World::World(int worldIndex, bool isVisual, SimTK::Real visualizerFrequency)
{
    ownWorldIndex = worldIndex;
    std::cout << "World::World BEGIN: ownWorldIndex: " << this->ownWorldIndex << std::endl << std::flush;
    fassno = new int; // External forcesa
    _useFixmanTorque = false;
  
    compoundSystem = new SimTK::CompoundSystem;
    matter = new SimTK::SimbodyMatterSubsystem(*compoundSystem);
    forces = new SimTK::GeneralForceSubsystem(*compoundSystem);

    this->visual = isVisual;
    if(visual){

        decorations = new SimTK::DecorationSubsystem(*compoundSystem);
        SimTK::Visualizer visualizer(*compoundSystem);

        // Also generate our decorations
        paraMolecularDecorator = new ParaMolecularDecorator;
        visualizer.addDecorationGenerator(paraMolecularDecorator);
        //ourDecorationGenerator = new SimTK::DecorationGenerator();
        //decorations->addDecorationGenerator(SimTK::Stage::Topology, ourDecorationGenerator);
        //



        visualizerReporter = new SimTK::Visualizer::Reporter(visualizer, visualizerFrequency);
        compoundSystem->addEventReporter( visualizerReporter );
    }
  
    forceField = new SimTK::DuMMForceFieldSubsystem(*compoundSystem);
    //forceField->loadAmber99Parameters();
    integ = new SimTK::VerletIntegrator(*compoundSystem);
    ts = new SimTK::TimeStepper(*compoundSystem, *integ);

    moleculeCount = -1;

    trackStep = 0;

    std::cout << "World::World END: ownWorldIndex: " << this->ownWorldIndex << std::endl << std::flush;
}

void World::AddMolecule(readAmberInput *amberReader, std::string rbFN, std::string flexFN, std::string ictdF)
{
    std::cout << "World::AddMolecule BEGIN: ownWorldIndex: " << this->ownWorldIndex << std::endl << std::flush;
    moleculeCount++; // Used for unique names of molecules

    this->rbFN = rbFN;
    this->frcmodF = frcmodF;
    this->flexFN = flexFN;
    this->ictdF = ictdF;
 
    std::cout << "World::AddMolecule moleculeReaders.push_back" << std::endl << std::flush;
    bMoleculeReader * molRead = new bMoleculeReader(amberReader, rbFN.c_str());
    moleculeReaders.push_back(molRead);

    std::cout << "World::AddMolecule add parameters" << std::endl << std::flush;
    bAddAllParams(std::string("lig") + std::to_string(moleculeCount), amberReader, *forceField, (moleculeReaders.back())->bAtomList, (moleculeReaders.back())->bonds);
  
    std::cout << "World::AddMolecule add Compound" << std::endl << std::flush;
    Topology * top = new Topology(std::string("lig") + std::to_string(moleculeCount));
    topologies.push_back(top);
  
    std::cout << "World::AddMolecule build Compound" << std::endl << std::flush;
    (topologies.back())->build(*forceField, (moleculeReaders.back())->natoms, (moleculeReaders.back())->bAtomList, (moleculeReaders.back())->nbonds, (moleculeReaders.back())->bonds, flexFN, ictdF);
  
    std::cout << "World::AddMolecule adopt Compound " << topologies.back() << std::endl << std::flush;
    compoundSystem->adoptCompound( *(topologies.back()) );
    std::cout << "World::AddMolecule realizeTopology" << std::endl << std::flush;
    compoundSystem->realizeTopology();

    std::cout << "Number of included atoms in nonbonded interactions: " << forceField->getNumNonbondAtoms() << std::endl << std::flush;
    std::cout << "getVdwGlobalScaleFactor() " << forceField->getVdwGlobalScaleFactor() << std::endl << std::flush;
    for(int i=0; i<(topologies.back())->natms; i++){
        std::cout << " DuMM VdW Radius " 
            << forceField->getVdwRadius(((topologies.back())->bAtomList[i]).getAtomClassIndex()) 
            << " DuMM VdW Well Depth "
            << forceField->getVdwWellDepth(((topologies.back())->bAtomList[i]).getAtomClassIndex())
            << std::endl << std::flush;
    }
  
    // Our decorations
    //if(visual){
    //    for(int i = 0; i < (topologies.back())->natms; i++){
    //        ourDecorativeGeometries.push_back(SimTK::DecorativeSphere(0.1));
    //    }
    //    for(int i = 0; i < (topologies.back())->nbnds; i++){
    //        ourDecorativeGeometries.push_back(SimTK::DecorativeLine());
    //    }
    //}

    std::cout << "World::AddMolecule END: ownWorldIndex: " << this->ownWorldIndex << std::endl << std::flush;
}

// Initialize simulation
void World::Init(SimTK::Real timestep, bool useFixmanTorqueOpt)
{
    // Only model after loading all the compounds
    integ->setFixedStepSize(timestep);

    std::cout << "World::Init BEGIN: ownWorldIndex: " << this->ownWorldIndex << std::endl << std::flush;
    std::cout << "World::Init model Compound" << std::endl << std::flush;
    compoundSystem->modelCompounds();

    // Load MobilizedBodyIndex vs Compound::AtomIndex maps 
    for ( unsigned int i = 0; i < this->topologies.size(); i++){
        ((this->topologies)[i])->loadMaps();
        std::cout << "Print maps topology " << i << std::endl;
        ((this->topologies)[i])->printMaps();
    }

    // Do we use Fixman torque
    _useFixmanTorque = useFixmanTorqueOpt;

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
 
    // Amber like scale factors. These should be used during simulations.
    forceField->setVdw12ScaleFactor(0.0);
    forceField->setVdw13ScaleFactor(0.0);
    forceField->setVdw14ScaleFactor(0.5);
    forceField->setVdw15ScaleFactor(1.0);
    forceField->setCoulomb12ScaleFactor(0.0);
    forceField->setCoulomb13ScaleFactor(0.0);
    forceField->setCoulomb14ScaleFactor(0.8333333333);
    forceField->setCoulomb15ScaleFactor(1.0);
    forceField->setVdwMixingRule(SimTK::DuMMForceFieldSubsystem::LorentzBerthelot);
  
    // My specific scale factors. These should be used during debugging.
    //->setSpecificDuMMScaleFactor(*forceField);
    // Solvent treatment
    forceField->setGbsaGlobalScaleFactor(0.0);
    std::cout << "GBSA solvent dielectric " << forceField->getSolventDielectric() << std::endl;
    std::cout << "GBSA solute dielectric " << forceField->getSoluteDielectric() << std::endl;
  
//    forceField->setBondStretchGlobalScaleFactor(0.0);
//    forceField->setBondBendGlobalScaleFactor(0.0);
//    forceField->setBondTorsionGlobalScaleFactor(0.0);
//    forceField->setAmberImproperTorsionGlobalScaleFactor(0.0);
//  
//    forceField->setVdw12ScaleFactor(0.0);
//    forceField->setVdw13ScaleFactor(0.0);
//    forceField->setVdw14ScaleFactor(0.0);
//    forceField->setVdw15ScaleFactor(0.0);
//    forceField->setVdwGlobalScaleFactor(0.0);
//  
//    forceField->setCoulomb12ScaleFactor(0.0);
//    forceField->setCoulomb13ScaleFactor(0.0);
//    forceField->setCoulomb14ScaleFactor(0.0);
//    forceField->setCoulomb15ScaleFactor(0.0);
//    forceField->setCoulombGlobalScaleFactor(0.0);

  
    *fassno = 0;
    if(_useFixmanTorque){
        ExtForce = new SimTK::Force::Custom(*forces, new GridForce(compoundSystem, *matter, fassno, this));
    }
    
    #ifdef TRY_TO_USE_OPENMM
        //forceField->setUseOpenMMAcceleration(true);
    #endif
    forceField->setTracing(true); // log OpenMM info to console
    //forceField->setNumThreadsRequested(1); // don't use this unless
  
    compoundSystem->realizeTopology();



    std::cout << "World::Init END: ownWorldIndex: " << this->ownWorldIndex << std::endl;
}//end of InitSimulation


// Interface
const Topology& World::getTopology(int moleculeNumber) const{
    return *(topologies[moleculeNumber]);
}

Topology& World::updTopology(void){
    return *(topologies.back());
}

// Returns: the inner vector represents the locations as a pair of bSpecificAtom
// toget various indeces from and a location. The outer vector represents the compounds.
// The two worlds have to be identical
// This iterates through all the molecules
std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > > World::getAtomsLocationsInGround(const SimTK::State & state)
{
    std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > > returnVector;

    // Iterate through topologies
    for ( unsigned int i = 0; i < this->topologies.size(); i++){
        std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
        // Iterate through atoms
        for(int j = 0; j < (topologies[i])->getNumAtoms(); j++){
            SimTK::Compound::AtomIndex compoundAtomIndex = topologies[i]->bAtomList[j].getCompoundAtomIndex();
            SimTK::Vec3 location = (topologies[i])->calcAtomLocationInGroundFrame(state, compoundAtomIndex);
            currentTopologyInfo.push_back(std::pair<bSpecificAtom *, SimTK::Vec3>(&((topologies[i])->bAtomList[j]), location));
        }
        returnVector.push_back(currentTopologyInfo);
    }

    return returnVector;
}

// Put coordinates into bAtomLists
void World::updateAtomLists(const SimTK::State & state)
{
    // Iterate through topologies
    for ( unsigned int i = 0; i < this->topologies.size(); i++){
        // Iterate through atoms
        for(int j = 0; j < topologies[i]->getNumAtoms(); j++){
            SimTK::Compound::AtomIndex compoundAtomIndex = topologies[i]->bAtomList[j].getCompoundAtomIndex();
            SimTK::Vec3 location = topologies[i]->calcAtomLocationInGroundFrame(state, compoundAtomIndex);
            topologies[i]->bAtomList[j].setX(location[0]);
            topologies[i]->bAtomList[j].setY(location[1]);
            topologies[i]->bAtomList[j].setZ(location[2]);
            //std::cout << "updateAtomList atom " << j << " " << topologies[i]->bAtomList[j].getX() << std::endl;
        }
    }
}

// 
SimTK::State& World::setAtomsLocationsInGround(SimTK::State& someState, std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > > otherWorldsAtomsLocations)
{
    //PrintSimbodyStateCache(someState);
    //someState.invalidateAll(SimTK::Stage::Topology);
    
    SimTK::Transform G_X_T;
    SimTK::Transform T_X_atom[matter->getNumBodies()]; // related to CompoundAtom.frameInMobilizedBodyFrame s
    SimTK::Transform T_X_PAt[matter->getNumBodies()];

    // Iterate through molecules/topologies
    for(unsigned int i = 0; i < otherWorldsAtomsLocations.size(); i++){

        SimTK::Transform G_X_T = topologies[i]->getTopLevelTransform();

        if(topologies[i]->getRegimen() == "IC"){
            // Create atomTargets
            std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
            std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > currentTopology = otherWorldsAtomsLocations[i];
            for(unsigned int j = 0; j < currentTopology.size(); j++){
                SimTK::Compound::AtomIndex atomIndex = ((currentTopology[j]).first)->atomIndex;
                SimTK::Vec3 location = ((currentTopology[j]).second);
                atomTargets.insert(std::pair<SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
                //const MobilizedBody& body = matter.getMobilizedBody(getAtomMobilizedBodyIndex(atomIndex));
                //SimTK::Transform X_FM;
                //body.setQToFitTransform(someState, X_FM);
            }

            // Match
            topologies[i]->matchDefaultAtomChirality(atomTargets, 0.01, false);
            topologies[i]->matchDefaultBondLengths(atomTargets);
            topologies[i]->matchDefaultBondAngles(atomTargets);
            topologies[i]->matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
            topologies[i]->matchDefaultTopLevelTransform(atomTargets);
            topologies[i]->matchDefaultConfiguration(atomTargets, SimTK::Compound::Match_Exact, true, 150.0);

            // Checked if match is done correctly
            this->trackStep += 1;
            std::string FN(std::string("pdbs/ICmatch") + std::to_string(this->trackStep) + std::string("before.pdb"));
            std::cout << "Writing file " << FN << std::endl;
            topologies[i]->writeDefaultPdb(FN.c_str(), SimTK::Transform());

            // Get transforms and locations: P_X_M, BAt_X_atom.p()
            G_X_T = topologies[i]->getTopLevelTransform();
            SimTK::Transform P_X_M[matter->getNumBodies()]; // related to X_PFs
            P_X_M[1] = G_X_T * topologies[i]->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(0));
            T_X_atom[1] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(0));
            //std::cout << "DEBUG: " << " aIx " << aIx << " mbx " << mbx << " parentMbx " << parentMbx << std::endl;
        
            // Iterate through atoms - get P_X_M for all the bodies
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin
                    // Get body, parentBody, parentAtom
                    SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        
                    // Get P_X_M
                    T_X_atom[int(mbx)] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(aIx);
                    //std::cout << "T_X_atom: "<< std::endl << T_X_atom[int(mbx)];
                    //std::cout << "Location in mobod: " << topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx);

                    P_X_M[int(mbx)] = G_X_T * T_X_atom[int(mbx)]; // MODIFIED
                    //std::cout << "P_X_M:" << std::endl << P_X_M[int(mbx)];
                }
            }
            // Set X_PF and Q - Bottleneck!
            for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
                SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
                ((SimTK::MobilizedBody::Free&)mobod).setDefaultInboardFrame(P_X_M[int(mbx)]);
            }

        // END IC regimen
        }else if(topologies[i]->getRegimen() == "RB"){ // TD and RB
            // Create atomTargets
            std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
            std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > currentTopology = otherWorldsAtomsLocations[i];
            for(unsigned int j = 0; j < currentTopology.size(); j++){
                SimTK::Compound::AtomIndex atomIndex = ((currentTopology[j]).first)->atomIndex;
                SimTK::Vec3 location = ((currentTopology[j]).second);
                atomTargets.insert(std::pair<SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
            }
    
            // Match Default Configuration ( => default state)
            topologies[i]->matchDefaultAtomChirality(atomTargets, 0.01, false);
            topologies[i]->matchDefaultBondLengths(atomTargets);
            topologies[i]->matchDefaultBondAngles(atomTargets);
            topologies[i]->matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
            topologies[i]->matchDefaultTopLevelTransform(atomTargets);
            topologies[i]->matchDefaultConfiguration(atomTargets, SimTK::Compound::Match_Exact, true, 150.0);

            // Checked if match is done correctly
            this->trackStep += 1;
            std::string FN(std::string("pdbs/RBmatch") + std::to_string(this->trackStep) + std::string("before.pdb"));
            std::cout << "Writing file " << FN << std::endl;
            topologies[i]->writeDefaultPdb(FN.c_str(), SimTK::Transform());

            // Get transforms and locations: P_X_M, BAt_X_atom.p()
            G_X_T = topologies[i]->getTopLevelTransform();
            SimTK::Transform P_X_M[matter->getNumBodies()]; // related to X_PFs
            SimTK::Angle inboardBondDihedralAngles[matter->getNumBodies()]; // related to X_FMs
            SimTK::Vec3 locs[topologies[i]->getNumAtoms()];
            P_X_M[1] = G_X_T * topologies[i]->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(0));
            T_X_atom[1] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(0));
        
            // Iterate through atoms - get P_X_M for all the bodies
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                // Get body, parentBody, parentAtom
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin
                    SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
                    const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
                    SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
                    std::cout << "RB origin atom " << aIx << " mobod (" << mbx  
                        << " parent mobod  ("<< parentMbx << ")" << std::endl;
                    T_X_atom[int(mbx)] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(aIx);
                    P_X_M[int(mbx)] = G_X_T * T_X_atom[int(mbx)];
                }
            }
        
            // Set transforms inside the bodies = BAt_X_atom.p; Set locations for everyone
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) != 0){ // atom is not at body's origin
                    SimTK::Transform atom_X_T = ~(T_X_atom[int(mbx)]);
                    SimTK::Transform T_X_B =  topologies[i]->calcDefaultAtomFrameInCompoundFrame(aIx);
                    SimTK::Transform atom_X_B = atom_X_T * T_X_B;

                    topologies[i]->bsetFrameInMobilizedBodyFrame(aIx, atom_X_B);

                    locs[int(aIx)] = atom_X_B.p();
                }
                else{
                    locs[int(aIx)] = SimTK::Vec3(0);
                }
            }
        
            // Set stations and AtomPLacements for atoms in DuMM
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    SimTK::DuMM::AtomIndex dAIx = topologies[i]->getDuMMAtomIndex(aIx);
                    forceField->bsetAtomStationOnBody( dAIx, locs[int(aIx)] );
                    forceField->updIncludedAtomStation(dAIx) = (locs[int(aIx)]);
                    forceField->bsetAtomPlacementStation(dAIx, mbx, locs[int(aIx)] );
            }
        
            // Set X_PF and Q - Bottleneck!
            for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
                SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
                ((SimTK::MobilizedBody::Free&)mobod).setDefaultInboardFrame(P_X_M[int(mbx)]);
                //((SimTK::MobilizedBody::Free&)mobod).setDefaultQ(SimTK::Vec7(0));
            }
        // END RB regimen
        }else if(topologies[i]->getRegimen() == "TD"){
            // Create atomTargets
            std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
            std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > currentTopology = otherWorldsAtomsLocations[i];
            for(unsigned int j = 0; j < currentTopology.size(); j++){
                SimTK::Compound::AtomIndex atomIndex = ((currentTopology[j]).first)->atomIndex;
                SimTK::Vec3 location = ((currentTopology[j]).second);
                atomTargets.insert(std::pair<SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
            }
    
            // Match Default Configuration ( => default state)
            topologies[i]->matchDefaultAtomChirality(atomTargets, 0.01, false);
            topologies[i]->matchDefaultBondLengths(atomTargets);
            topologies[i]->matchDefaultBondAngles(atomTargets);
            topologies[i]->matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
            topologies[i]->matchDefaultTopLevelTransform(atomTargets);
            topologies[i]->matchDefaultConfiguration(atomTargets, SimTK::Compound::Match_Exact, true, 150.0);

            // Checked if match is done correctly
            this->trackStep += 1;
            std::string FN(std::string("pdbs/TDmatch") + std::to_string(this->trackStep) + std::string("before.pdb"));
            std::cout << "Writing file " << FN << std::endl;
            topologies[i]->writeDefaultPdb(FN.c_str(), SimTK::Transform());

            // Get transforms and locations: P_X_M, BAt_X_atom.p()
            G_X_T = topologies[i]->getTopLevelTransform();
            SimTK::Transform M_X_pin = SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis); // Moves rotation from X to Z
            SimTK::Transform P_X_M[matter->getNumBodies()]; // related to X_PFs
            SimTK::Transform T_X_atom[matter->getNumBodies()]; // related to CompoundAtom.frameInMobilizedBodyFrame s
            SimTK::Transform T_X_PAt[matter->getNumBodies()];
            SimTK::Angle inboardBondDihedralAngles[matter->getNumBodies()]; // related to X_FMs
            SimTK::Vec3 locs[topologies[i]->getNumAtoms()];
            P_X_M[1] = G_X_T * topologies[i]->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(0));
            T_X_atom[1] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(0));
        
            // Iterate through atoms - get P_X_M for all the bodies
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin
                    // Get body, parentBody, parentAtom
                    SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
                    const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
                    SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

                    //std::cout << "atom " << aIx << " mobod (" << mbx  
                    //    << " parent mobod  ("<< parentMbx << ")" << std::endl;

                    if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground
                        SimTK::Compound::AtomIndex parentAIx = (topologies[i]->getMbx2aIx()).at(parentMbx);
                        T_X_PAt[int(mbx)] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(parentAIx);
                        // Get inboard dihedral angle and put in BAt_X_M0 !!!!!!!
                        inboardBondDihedralAngles[int(mbx)] = topologies[i]->bgetDefaultInboardDihedralAngle(aIx);
                        SimTK::Transform BAt_X_M0 = SimTK::Rotation(inboardBondDihedralAngles[int(mbx)], SimTK::XAxis);

                        // Get P_X_M
                       T_X_atom[int(mbx)] = topologies[i]->calcDefaultAtomFrameInCompoundFrame(aIx);
                       SimTK::Transform T_X_M0 = T_X_atom[int(mbx)] * BAt_X_M0;
                       const SimTK::Transform& T_X_PAt = topologies[i]->calcDefaultAtomFrameInCompoundFrame(parentAIx);
                        SimTK::Transform PAt_X_T = ~T_X_PAt;
                        SimTK::Transform PAt_X_M0 = PAt_X_T * T_X_M0;
                        P_X_M[int(mbx)] = PAt_X_M0;

                    } //END if parent not Ground
                }
            }
        
            // Set transforms inside the bodies = BAt_X_atom.p; Set locations for everyone
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) != 0){ // atom is not at body's origin
                    SimTK::Transform atom_X_T = ~(T_X_atom[int(mbx)]);
                    SimTK::Transform T_X_B =  topologies[i]->calcDefaultAtomFrameInCompoundFrame(aIx);
                    SimTK::Transform atom_X_B = atom_X_T * T_X_B;
                    topologies[i]->bsetFrameInMobilizedBodyFrame(aIx, atom_X_B);
                    locs[int(aIx)] = atom_X_B.p();
                }
                else{
                    locs[int(aIx)] = SimTK::Vec3(0);
                }
            }
        
            // Set stations and AtomPLacements for atoms in DuMM
            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                    SimTK::DuMM::AtomIndex dAIx = topologies[i]->getDuMMAtomIndex(aIx);
                    forceField->bsetAtomStationOnBody( dAIx, locs[int(aIx)] );
                    forceField->updIncludedAtomStation(dAIx) = (locs[int(aIx)]);
                    forceField->bsetAtomPlacementStation(dAIx, mbx, locs[int(aIx)] );
            }
        
            // Set X_PF and Q - Bottleneck!
            for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
                SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
                if(int(mbx) == 1){ // This is dangerous TODO
                    ((SimTK::MobilizedBody::Free&)mobod).setDefaultInboardFrame(P_X_M[int(mbx)]);
                }else{
                    ((SimTK::MobilizedBody::Pin&)mobod).setDefaultInboardFrame(P_X_M[int(mbx)] * M_X_pin);
                    ((SimTK::MobilizedBody::Pin&)mobod).setDefaultQ(inboardBondDihedralAngles[int(mbx)]);
                }
            }
 
        } // END TD regimen and all regimens

    } // END iterating through molecules/topologies

    this->compoundSystem->realizeTopology();

    someState = compoundSystem->updDefaultState();

    this->compoundSystem->realize(someState, SimTK::Stage::Position);

    // Write a check file after updating someState and realizing Position
    std::string prefix;
    prefix = std::string("pdbs/") + std::string(topologies[0]->getRegimen()) + std::string("match") ;
    std::string FN2(prefix + std::to_string(this->trackStep) + std::string("after.pdb"));
    std::cout << "Writing file " << FN2 << std::endl;
    std::filebuf fb;
    fb.open (FN2, std::ios::out);
    std::ostream os(&fb);
    ((SimTK::Compound *)topologies[0])->writePdb(someState, os);
    fb.close();
     //

    updateAtomLists(someState);


    // Configuration
    //std::cout << "assigned conf: ";
    //for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    //    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    //    std::cout << " mobod " << int(mbx) << " = ";
    //    std::cout << mobod.getQAsVector(someState) << std::endl ;
    //    std::cout << " P_X_F " << mobod.getInboardFrame(someState) << " ";
    //    std::cout << " F_X_M " << mobod.getMobilizerTransform(someState) << " ";
    //    std::cout << "; ";
    //}
    //std::cout << std::endl;


    return someState;
}

void World::PrintSimbodyStateCache(SimTK::State& someState){
    std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
    for(int i = 0; i < someState.getNumSubsystems(); i++){
        std::cout << " Subsystem " << i << " Name: " << someState.getSubsystemName(SimTK::SubsystemIndex(i))
            << " Stage: " << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
            << " Version: " << someState.getSubsystemVersion(SimTK::SubsystemIndex(i)) << std::endl;
    }
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
World::~World(){
    if(this->visual == true){
        //delete paraMolecularDecorator;

        delete decorations;
        //delete vizReporter;
    }
    delete fassno;
    delete ts;
    delete integ;
    delete forceField;
    if(_useFixmanTorque){
        delete ExtForce;
    }
    delete matter;
    delete forces;
    delete compoundSystem;
    for(unsigned int i = 0; i < moleculeReaders.size(); i++){
        delete moleculeReaders[i];
    }
    for(unsigned int i = 0; i < topologies.size(); i++){
        delete topologies[i];
    }
    //forceField->loadAmber99Parameters();
}



