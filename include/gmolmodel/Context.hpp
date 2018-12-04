#ifndef __CONTEXT_HPP__
#define __CONTEXT_HPP__

#include "Robo.hpp"

class Sampler;
class World;

class Context{

public:
    Context(World *);
    Context();
    ~Context();

    World * AddWorld(bool visual);
    //World * AddWorld(World *, bool visual);

    World * getWorld(void) const;
    World * getWorld(int which) const;

    World * updWorld(void);
    World * updWorld(int which);

    unsigned int getNofWorlds(void);

    SimTK::DuMMForceFieldSubsystem * updForceField(int whichWorld);

    // Writeble reference to a samplers advanced state
    SimTK::State& updAdvancedState(int whichWorld, int whichSampler);

    // --- Use a SetupReader Object to read worlds information from a file ---
    bool loadTopologyFile(int whichWorld, int whichMolecule, std::string topologyFilename);
    bool loadCoordinatesFile(int whichWorld, int whichMolecule, std::string coordinatesFilename);
    bool loadRigidBodiesSpecs(int whichWorld, int whichMolecule, std::string RBSpecsFN);
    bool loadFlexibleBondsSpecs(int whichWorld, int whichMolecule, std::string FlexSpecsFN);
    void setRegimen (int whichWorld, int whichMolecule, std::string regimen);
    void loadMolecules();
    void modelTopologies(void);

    void LoadWorldsFromSetup(SetupReader&);
    //------------

    // --- Thermodynamics ---a
    // Get/set the main temperature (acc/rej temperature for MC)
    float getTemperature(int whichWorld);
    void  setTemperature(int whichWorld, float someTemperature);

    // If HMC, get/set the guidance Hamiltonian temperature
    float getGuidanceTemperature(int whichWorld, int whichSampler);
    void  setGuidanceTemperature(int whichWorld, int whichSampler, float someTemperature);
    //------------

    // --- Simulation parameters ---
    int addSampler(int whichWorld, std::string whichSampler);
    int addSampler(int whichWorld, SamplerName whichSampler);
    void initializeSampler(int whichWorld, int whichSampler, bool randomizeConformation = false);

    // Amber like scale factors.
    void setAmberForceFieldScaleFactors(int whichWorld);

    // Set a global scaling factor for the forcefield
    void setGlobalForceFieldScaleFactor(int whichWorld, SimTK::Real);

    // Set GBSA implicit solvent scale factor
    void setGbsaGlobalScaleFactor(int whichWorld, SimTK::Real);

    // If HMC, get/set the number of MD steps
    int getNofMDStepsPerSample(int whichWorld);
    void setNofMDStepsPerSample(int whichWorld, int MDStepsPerSample);

    // If HMC, get/set timestep forMD
    const float getTimestep(int whichWorld, int whichSampler);
    void setTimestep(int whichWorld, int whichSampler, float timeStep);

    // Use Fixman torque as an additional force subsystem
    void useFixmanPotential(int whichWorld, int whichSampler);
    bool isUsingFixmanPotential(int whichWorld, int whichSampler);
    void useFixmanTorque(int whichWorld, SimTK::Real argTemperature);
    bool isUsingFixmanTorque(int whichWorld);
    void setFixmanTorqueScaleFactor(int whichWorld, double scaleFactor);
    void setFixmanTorqueTemperature(int whichWorld, double temperature);
    //------------

    // --- Mixing parameters ---
    // Another way to do it is setting the number of rounds
    int getNofRounds(void);
    void setNofRounds(int nofRounds);

    int getNofSamplesPerRound(int whichWorld);
    void setNofSamplesPerRound(int whichWorld, int MCStepsPerRound);

    int getWorldIndex(int which);

    // --- Arrange different mixing parameters ---
    void initializeMixingParamters(void);
    //------------

    // --- Mix ---
    void RotateWorlds(void);
    //------------

    // --- Main ---
    void Run(SetupReader&);
    void setNumThreadsRequested(int which, int howMany);

    /** Initialize the same velocities **/
    bool getReproducible(void);
    void setReproducible(void);
    //------------

    // --- Printing functions ---
    void WritePdb(int whichWorld);
    SimTK::Real Dihedral(int whichWorld, int whichCompound, int whichSampler, int a1, int a2, int a3, int a4);
    SimTK::Real Distance(int whichWorld, int whichCompound, int whichSampler, int a1, int a2);
    //------------

public:
    std::vector<int> worldIndexes;

private:
    std::vector<World *> worlds;

    // Molecules files
    std::vector<std::vector<std::string>> topFNs;
    std::vector<std::vector<std::string>> crdFNs;
    std::vector<std::vector<std::string>> rbSpecsFNs;
    std::vector<std::vector<std::string>> flexSpecsFNs;
    std::vector<std::vector<std::string>> regimens;

    // Simulation parameters
    int nofRounds;
    //int total_mcsteps;
    std::vector<int> nofSamplesPerRound;
    std::vector<int> nofMDStepsPerSample;
    std::vector<float> timesteps;

    //
    bool reproducible;
};

#endif //__CONTEXT_HPP__
