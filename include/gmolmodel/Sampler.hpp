#ifndef __SAMPLER_HPP__
#define __SAMPLER_HPP__

#include "Robo.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

//#ifndef HARMONICOSCILLATOR
//#define HARMONICOSCILLATOR // for testing purposes
//#endif

#ifndef PRINT_BUFFER_SIZE
#define PRINT_BUFFER_SIZE 4096
#endif 

class Topology;

class Sampler
{
public:
    // Constructor
    Sampler(SimTK::CompoundSystem *argCompoundSystem,
            SimTK::SimbodyMatterSubsystem *argMatter,
            SimTK::Compound *argResidue,
            SimTK::DuMMForceFieldSubsystem *argDumm,
            SimTK::GeneralForceSubsystem *forces,
            SimTK::TimeStepper *argTimeStepper);


    // Destructor
    virtual ~Sampler();

    // Compute mass matrix determinant (O(n))
    SimTK::Real calcMassDeterminant(const SimTK::State& );
    SimTK::Real calcMassDeterminant(SimTK::State& );

    // Initializtion functions
    virtual void initialize(void) {};
    virtual void reinitialize(void) {}

    /*
    // Set a thermostat (even for MCMC)
    virtual void setThermostat(ThermostatName);
    virtual void setThermostat(std::string);
    virtual void setThermostat(const char *);

    // Get the name of the thermostat
    virtual ThermostatName getThermostat(void);
    */
    
    // Doc these later
    virtual void setTemperature(SimTK::Real) = 0;

    // Extract one or more samples
    virtual void update(SimTK::State&) = 0;

    /** Returns the number of samples extracted so far. **/
    int getNofSamples(void);

    // Get set the seed
    unsigned long long int getSeed(void);
    void setSeed(unsigned long long int);

    // For debugging purposes
    void PrintSimbodyStateCache(SimTK::State& someState);

    //////////////////////////////////
    // Harmonic oscillator functions
    //////////////////////////////////

    // Potential energy function - modifies HO_f
    double HarmonicOscillatorPE(double *someX);
    
    // Kinetic energy function
    double HarmonicOscillatorKE(double *someV);
   
    // Initialize velocities according to T 
    void HO_InitializeVelocity(double *someV, double T);

    // Integrate xprop using Velocity Verlet - modifies HO_xprop, HO_v, HO_f, HO_a
    void HO_VelocityVerlet(double dt, int nofSteps);
    ///////////////////////////////////

public:
    const SimTK::System *system;
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    SimTK::Compound *residue;
    SimTK::DuMMForceFieldSubsystem *dumm;
    SimTK::GeneralForceSubsystem *forces;
    SimTK::TimeStepper *timeStepper;

    // THermodynamics
    ThermostatName thermostat;

    // Sampling
    int nofSamples;
    //int printBuffIx;
    unsigned long long int seed;
  
    // Harmonic oscillator constants
    static const int HO_D = 20; // Dimensionality
    double HO_x[HO_D]; // Equilibrium position
    double HO_x0[HO_D]; // Equilibrium position
    double HO_xini[HO_D]; // Initial position
    double HO_xprop[HO_D]; // Proposed position
    double HO_v[HO_D]; // Velocity
    double HO_f[HO_D]; // Force
    double HO_m[HO_D]; // Mass
    double HO_a[HO_D]; // Acceleration
    double HO_k[HO_D]; // Force constants

    double pertAmp; // Perturbation amplitude

    double HO_PE_x; // PE at the begining of the proposal
    double HO_PE_xprop; // PE after the proposal
    double HO_PE_set; // PE set at the end of the acc-rej step

    double HO_KE_x;
    double HO_KE_xprop;
    double HO_KE_set;

    double HO_etot_x;
    double HO_etot_xprop;
    double HO_etot_set;

    boost::random::mt19937 HO_randomEngine = boost::random::mt19937();
    boost::normal_distribution<> HO_gaurand = boost::normal_distribution<>(0.0, 1.0);

};


#endif // __SAMPLER_HPP__

