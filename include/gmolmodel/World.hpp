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

//#ifndef TRY_TO_USE_OPENMM
//#define TRY_TO_USE_OPENMM
//#endif

#include "bMoleculeReader.hpp"
#include "bAddParams.hpp"
#include "server.hpp"
#include "Topology.hpp"
#include "bArgParser.hpp"
#include "HamiltonianMonteCarloSampler.hpp"

//typedef double Vector3[3];

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
 *  Contains a Symbody system and additional data that define a regimen
 **/
class World{
public: 
    // --- Structural functions ---
    // Constructor 
    World(int worldIndex, bool isVisual=true, SimTK::Real visualizerFrequency = 0.0015);
 
    // Creates a topology object and based on amberReader forcefield parameters
    //  - defines Biotypes; - adds BAT parameters to DuMM
    void AddMolecule(readAmberInput *amberReader, std::string rbFN, std::string flexFN, std::string ictdF);
 
    // Intializes simulation parameters.
    // Beside passing the integrator step and Fixman torque option,
    // it also initializes default scaling factors for Amber
    void ModelTopologies(bool useFixmanTorque = false);
    //...............

    // --- Mixing functions: Pass configurations among Worlds
    // Get the current Compound Cartesian coords
    std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > >  getAtomsLocationsInGround(const SimTK::State&);

    // Set Compound, MultibodySystem and DuMM configurations according to
    // some other World's atoms
    SimTK::State& setAtomsLocationsInGround(SimTK::State&, std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > > otherWorldsAtomsLocations);
  
    // Update Gmolmodel bSpecificAtom Cartesian coordinates according to
    // Molmodel Compound
    void updateAtomLists(const SimTK::State&);

    // Access to molecule (Topology) objects 
    // Get a readble reference of one of the molecules 
    const Topology& getTopology(int moleculeNumber) const;

    // Get a writeble reference of one of the molecules 
    Topology& updTopology(void);  
    //...............
   
    // --- Thermodynamics ---
    // Get macro temperature
    SimTK::Real getTemperature(void);

    // Set macro temperature
    void setTemperature(SimTK::Real);
    //...............

    // --- Simulation ---
    // Amber like scale factors.
    void setAmberForceFieldScaleFactors(void);

    // Set a global scaling factor for the forcefield
    void setGlobalForceFieldScaleFactor(SimTK::Real);

    // Set GBSA implicit solvent scale factor
    void setGbsaGlobalScaleFactor(SimTK::Real);

    SimTK::DuMMForceFieldSubsystem * updForceField(void);

    // Use the Fixman torque as an additional force subsystem
    void useFixmanTorque(void);
    bool isUsingFixmanTorque(void);
    //...............

    // --- Statistics ---
    // How many samples did we have so far
    int getSampleNumber(void);

    // Sampler manipulation functions
    int getNofSamplers(void);

    int addSampler(std::string);
    int addSampler(SamplerName);

    // Get a sampler based on its position in the samplers vector
    const HamiltonianMonteCarloSampler * getSampler(int which);

    // Get a writable sampler based on its position in the samplers vector
    HamiltonianMonteCarloSampler * updSampler(int which);

    SimTK::Force::Custom * updFixmanTorque(void);

    //...............
  
    // -- Debugging / helper functions ---
    // Print information about Simbody systems
    void PrintSimbodyStateCache(SimTK::State&);
  
    // Print Compound Cartesian coordinates
    void printPoss(const SimTK::Compound& c, SimTK::State& someState);
    //...............
 
    // --- To be removed ---
    void Advance(int nosteps);
    //...............
 
    // Destructor
    ~World(); // destructor

public:

    // --- Structural data ---
    // System->MultibodySystem->MolecularMechanicsSystems->CompoundSystem
    SimTK::CompoundSystem *compoundSystem;

    // Subsystem->SimbodyMatterSubsystem
    SimTK::SimbodyMatterSubsystem *matter;

    // Subsystem->ForceSubsystem->GeneralForceSubsystem
    SimTK::GeneralForceSubsystem *forces;

    // Force::Custom
    SimTK::Force::Custom *ExtForce;

    // Subsystem->ForceSubsystem->DuMMForceFieldSubsystem
    SimTK::DuMMForceFieldSubsystem *forceField;

    // Nof molecules
    int moleculeCount;

    // Filenames for files needed to build molecules (topologies)
    std::string rbFN; // rigid bodies specifications
    std::string frcmodF; // to be removed
    std::string flexFN; // flexible bonds specifications
    std::string ictdF; // regimen specification
  
    // Molecules (topologies<-Compounds) objects
    std::vector<bMoleculeReader *> moleculeReaders;
    std::vector<Topology *> topologies;

    // This vector stores a configuration if is needed for later use
    SimTK::Transform *TVector;
  
    // Topologies graphs as tables - to be removed
    int **mbxTreeMat;    // tree representing the bonding
    SimTK::Real *branchMassVec; // branch masses self body included
    //...............

    // --- Thermodynamics ---
    SimTK::Real temperature;
    //...............

    // --- Simulation ---
    SimTK::VerletIntegrator *integ;
    SimTK::TimeStepper *ts;
    std::vector<HamiltonianMonteCarloSampler *> samplers;
    bool _useFixmanTorque;
    //...............

    // --- Statistics ---
    int sampleNumber;
    //...............
  
    // --- Graphics ---
    bool visual;
    // Our decorations
    ParaMolecularDecorator *paraMolecularDecorator;

    // Decoration subsystem
    SimTK::DecorationSubsystem *decorations;

    // Visualizer
    SimTK::Visualizer *visualizer;

    // Visualizer reporter
    SimTK::Visualizer::Reporter *visualizerReporter;
    //...............

    // --- Mixing data ---
    int ownWorldIndex;
    //...............

};

#endif /*WORLD_H_*/
