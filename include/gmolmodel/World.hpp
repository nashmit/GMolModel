#ifndef WORLD_H_
#define WORLD_H_

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

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

#include "ParaMolecularDecorator.hpp"
#include "FixmanTorque.hpp"

#include <boost/timer.hpp>

//#include "/home/lspirido/Installers/armadillo-6.400.3/include/armadillo.hpp"

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

  // Our decorations
  ParaMolecularDecorator *paraMolecularDecorator;
  //

  SimTK::DecorationSubsystem *decorations;
  SimTK::Visualizer *visualizer;
  SimTK::Visualizer::Reporter *visualizerReporter;
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

  int sampleNumber;
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

#endif /*WORLD_H_*/
