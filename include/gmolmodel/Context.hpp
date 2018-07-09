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
    World * AddWorld(World *, bool visual);

    World * getWorld(void) const;
    World * getWorld(int which) const;

    World * updWorld(void);
    World * updWorld(int which);

    SimTK::DuMMForceFieldSubsystem * updForceField(int whichWorld);

    // --- Use a SetupReader Object to read worlds information from a file ---
    bool loadTopologyFile(int whichWorld, int whichMolecule, std::string topologyFilename);
    bool loadCoordinatesFile(int whichWorld, int whichMolecule, std::string coordinatesFilename);
    bool loadRigidBodiesSpecs(int whichWorld, int whichMolecule, std::string RBSpecsFN);
    bool loadFlexibleBondsSpecs(int whichWorld, int whichMolecule, std::string FlexSpecsFN);
    void setRegimen (int whichWorld, int whichMolecule, std::string regimen);
    void loadMolecules();

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
    //------------

    // --- Mixing parameters ---
    // Another way to do it is setting the number of rounds
    int getNofRounds(int nofRounds);
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
    //------------

    // --- Printing functions ---
    void WritePdb(int whichWorld);
    void Dihedral(int, int, int, int);
    //------------

private:
    std::vector<World *> worlds;
    std::vector<int> worldIndexes;

    // Molecules files
    std::vector<std::vector<std::string>> topFNs;
    std::vector<std::vector<std::string>> crdFNs;
    std::vector<std::vector<std::string>> rbSpecsFNs;
    std::vector<std::vector<std::string>> flexSpecsFNs;
    std::vector<std::vector<std::string>> regimens;

    // Simulation parameters
    int nofRounds;
    int total_mcsteps;
    std::vector<int> nofSamplesPerRound;
    std::vector<int> nofMDStepsPerSample;
    std::vector<float> timesteps;
};

#endif //__CONTEXT_HPP__
