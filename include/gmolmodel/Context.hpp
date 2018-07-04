#ifndef __CONTEXT_HPP__
#define __CONTEXT_HPP__

#include "Robo.hpp"
//#include "SetupReader.hpp"

class Sampler;
class World;

class Context{
//private:
//    World *p_world;

public:
    Context(World *);
    Context();
    ~Context();

    World * AddWorld(World *, bool visual);

    World * getWorld(void) const;
    World * getWorld(int which) const;

    World * updWorld(void);
    World * updWorld(int which);

    // --- Use a SetupReader Object to read worlds information from a file ---
    void LoadWorldsFromSetup(SetupReader&);
    //------------

    // --- Thermodynamics ---a
    // Get/set the main temperature (acc/rej temperature for MC)
    float getTemperature(int whichWorld);
    void  setTemperature(int whichWorld, float someTemperature);

    // If HMC, get/set the guidance Hamiltonian temperature
    float getGuidanceTemperature(void);
    void  setGuidanceTemperature(float someTemperature);
    //------------

    // --- Simulation parameters ---
    // If HMC, get/set the number of MD steps
    int getNofMDSteps(int whichWorld);
    void setNofMDSteps(int whichWorld, int nofMDSteps);

    // If HMC, get/set timestep forMD
    float getTimeStep(int whichWorld);
    void setTimeStep(int whichWorld, float timeStep);
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

    // Simulation parameters
    int nofRounds;
    int total_mcsteps;
    std::vector<int> nofSamplesPerRound;
    std::vector<int> nofMDSteps;
    std::vector<float> timesteps;
};

#endif //__CONTEXT_HPP__
