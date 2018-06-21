#ifndef __HAMMONTECARLOSAMPLER_HPP__
#define __HAMMONTECARLOSAMPLER_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include "MonteCarloSampler.hpp"

class Topology;
class IState;
void writePdb(      SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

/** A Generalized Coordiantes Hamiltonian Monte Carlo sampler as described in
J Chem Theory Comput. 2017 Oct 10;13(10):4649-4659. In short it consists
of the following steps:
   1. Initialize velocities from a random normal distribution with a 
      covariance of kT sqrt(M) where M is the mass matrix tensor.
   2. Propagate the trial trajectory using a symplectic integrator provided
      by Simbody
   3. An acception-rejection step which includes the Fixman potential 
      if needed.
Step 1 and 2 are implemented in the fuction propose. Step 3 is implemented
in the function update.
**/
class HamiltonianMonteCarloSampler : public MonteCarloSampler
{
public:

    /** Constructor **/
    HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                 SimTK::SimbodyMatterSubsystem *argMatter,
                                 SimTK::Compound *argResidue,
                                 SimTK::DuMMForceFieldSubsystem *argDumm,
                                 SimTK::GeneralForceSubsystem *forces,
                                 SimTK::TimeStepper *argTimeStepper);

    /** Destructor **/
    virtual ~HamiltonianMonteCarloSampler();

    /** Calculate O(n2) the square root of the mass matrix inverse
    denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
    This is lower triangular **/
    void calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv);

    /** Calculate O(n2) the square root of the mass matrix inverse
    denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular **/
    void calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv);

    /** Calculate sqrt(M) using Eigen. For debug purposes. **/
    void calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper);

    /** Seed the random number generator. Set simulation temperature, 
    velocities to desired temperature, variables that store the configuration
    and variables that store the energies, both needed for the 
    acception-rejection step. Also realize velocities and initialize 
    the timestepper. **/
    virtual void initialize(SimTK::State& advanced, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman = true) ; 

    /** Same as initialize **/
    virtual void reinitialize(SimTK::State& advanced, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) ; 

    /** It implements the proposal move in the Hamiltonian Monte Carlo
    algorithm. It essentially propagates the trajectory after it stores
    the configuration and energies. TODO: break in two functions:
    initializeVelocities and propagate/integrate **/
    void propose(SimTK::State& someState, SimTK::Real timestep, int nosteps);

    /** Main function that contains all the 3 steps of HMC.
    Implements the acception-rejection step and sets the state of the 
    compound to the appropriate conformation wether it accepted or not. **/
    void update(SimTK::State& someState, SimTK::Real timestep, int nosteps);

    /** Get the proposed kinetic energy. This is set right  after velocities
    are initialized. **/
    SimTK::Real getProposedKE(void) { return this->ke_proposed; }
    
    /** Get the stored kinetic energy. This is set rightafter a move is
    accepted. It's a component of the total energy stored. **/
    SimTK::Real getLastAcceptedKE(void) { return this->ke_lastAccepted; }
    
    /** Sets the proposed kinetic energy before the proposal. This should be
    set right after the velocities are initialized. **/
    void setProposedKE(SimTK::Real);

    /** Stores the accepted kinetic energy. This should be set right after a 
    move is accepted. It's a component of the total energy stored. **/
    void setLastAcceptedKE(SimTK::Real);

    /** Returns the number of MC trials done by this integrator. **/
    //int getSampleNumber(void);

protected:
    SimTK::Real ke_lastAccepted; // last accepted kinetic energy
    SimTK::Real ke_proposed; // proposed kinetic energy
    SimTK::Real etot_set; // stored total energy
    SimTK::Real etot_proposed; // last accepted total energ (same with stored)
    int sampleNumber; // counter for the number of MC trials 

};

#endif // __HAMMONTECARLOSAMPLER_HPP__

