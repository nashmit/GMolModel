#ifndef BSYSTEM_H_
#define BSYSTEM_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <time.h>
#include <thread>
#include <array>
//#include <random>
#include <math.h>

#include <Eigen/Dense>
//#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
//using Eigen::MatrixXd;

#include "Simbody.h"
#include "Molmodel.h"

#include "readAmberInput.hpp"

#include <boost/timer.hpp>

//#include "/home/lspirido/Installers/armadillo-6.400.3/include/armadillo.hpp"

#ifndef TRY_SOFT_LJ
#define TRY_SOFT_LJ
#endif

/*
*/
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif
#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif
#ifndef DEBUG_WRITEDATA
#define DEBUG_WRITEDATA
#endif

/*
#ifndef DEBUG_WRITEALLPDBS
#define DEBUG_WRITEALLPDBS
#endif
*/

/*
#ifndef DEBUG_WRITEPDBS
#define DEBUG_WRITEPDBS
#endif
#ifndef DEBUG_WRITEPFPDB
#define DEBUG_WRITEPFPDB
#endif
*/

#ifndef DEBUG_CONF
#define DEBUG_CONF
#endif

/*
#ifndef DEBUG_SPECIFIC
#define DEBUG_SPECIFIC
#endif
#ifndef DEBUG_ENERGY
#define DEBUG_ENERGY
#endif
*/
/*
#ifndef DEBUG_TIME
#define DEBUG_TIME
#endif
*/

#ifndef TRY_TO_USE_OPENMM
#define TRY_TO_USE_OPENMM
#endif


#include "bMoleculeReader.hpp"
#include "bAddParams.hpp"
#include "server.hpp"
#include "Topology.hpp"
#include "bArgParser.hpp"

typedef double Vector3[3];

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix);

void writePdb(      SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix);

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

void writePdb(      SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

void writePdb(SimTK::PdbStructure pdb, const char *FN);

class World;

//==============================================================================
//                           CLASS GridForce
//==============================================================================
/**
 * External Custom Forces Class. It sends the cartesian coordinates calculated
 * by Compound and applies the forces passed by MMTK through PyArrayObjects 
 **/
class GridForce : public SimTK::Force::Custom::Implementation {
 public:
  SimTK::CompoundSystem *compoundSystem;
  int *fassno;
  int *flag;
  World *Caller;

  GridForce(SimTK::CompoundSystem *compoundSystem, SimTK::SimbodyMatterSubsystem& matter
            , int *fassno
            , World *Caller
            );

  void calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
    SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const;

  SimTK::Real calcPotentialEnergy(const SimTK::State& state) const;

  bool dependsOnlyOnPositions() const;

 private:
  SimTK::SimbodyMatterSubsystem& matter;
};

//==============================================================================
//                           CLASS World
//==============================================================================
/**
 *  Symbols System Class. Manages Molmodel-gMolmodel-MMTK communication.
 **/
class World{
 public:
  SimTK::CompoundSystem *compoundSystem;
  SimTK::SimbodyMatterSubsystem *matter;
  SimTK::GeneralForceSubsystem *forces;
  SimTK::Force::Custom *ExtForce;
  SimTK::DecorationSubsystem *decorations;
  SimTK::Visualizer::Reporter *vizReporter;
  bool visual;
  SimTK::DuMMForceFieldSubsystem *forceField;

  std::vector<bMoleculeReader *> moleculeReaders;
  std::vector<Topology *> topologies;

  SimTK::Visualizer *viz;
  #ifdef NOSETHERMOS
  SimTK::NoseHooverThermostat *thermo;
  #endif
  #ifdef VELSTHERMOS
  SimTK::VelocityRescalingThermostat *vthermo;
  #endif
  SimTK::VerletIntegrator *integ;
  SimTK::TimeStepper *ts;
  std::string mol2F, rbFN, frcmodF, flexFN, ictdF;

  int moleculeCount;
  int ownWorldIndex;

  int *fassno;
  SimTK::Transform *TVector;
  int **mbxTreeMat;    // tree representing the bonding
  SimTK::Real *branchMassVec; // branch masses self body included
  bool _useFixmanTorque;

  World(int worldIndex, bool isVisual=true, SimTK::Real visualizerFrequency = 0.0015);

  void AddMolecule(readAmberInput *amberReader, std::string rbFN, std::string flexFN, std::string ictdF);

  void Init(SimTK::Real integTimestep, bool useFixmanTorque = false);

  std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > >  getAtomsLocationsInGround(const SimTK::State&);

  SimTK::State& setAtomsLocationsInGround(SimTK::State&, std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > > otherWorldsAtomsLocations);

  void updateAtomLists(const SimTK::State&);
 
  void PrintSimbodyStateCache(SimTK::State&);

  void printPoss(const SimTK::Compound& c, SimTK::State& someState);
 
  // Interface
  const Topology& getTopology(int moleculeNumber) const;
  Topology& updTopology(void);

  // Manages the TimeStepper actions
  void Advance(int nosteps);

  ~World(); // destructor
};

#endif /*BSYSTEM_H_*/
