#ifndef __SAMPLER_HPP__
#define __SAMPLER_HPP__

#include "Robo.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

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

    // Update
    virtual void update(SimTK::State&);

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

protected:
    const SimTK::System *system;
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    SimTK::Compound *residue;
    SimTK::DuMMForceFieldSubsystem *dumm;
    SimTK::GeneralForceSubsystem *forces;
    SimTK::TimeStepper *timeStepper;
   
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

