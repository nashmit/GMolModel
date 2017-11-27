#ifndef __HAMMONTECARLOSAMPLER_HPP__
#define __HAMMONTECARLOSAMPLER_HPP__

#include "MonteCarloSampler.hpp"

class Topology;
class IState;
void writePdb(      SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

class HamiltonianMonteCarloSampler : public MonteCarloSampler
{
public:

    // Constructor
    HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                 SimTK::SimbodyMatterSubsystem *argMatter,
                                 SimTK::Compound *argResidue,
                                 SimTK::DuMMForceFieldSubsystem *argDumm,
                                 SimTK::GeneralForceSubsystem *forces,
                                 SimTK::TimeStepper *argTimeStepper);

    // Destructor
    virtual ~HamiltonianMonteCarloSampler();

    // Calculate O(n2) the square root of the mass matrix inverse
    void calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv);
    void calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv);
    // Calculate sqrt(M) using Eigen
    void calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper);
    // Initialize variables (like TVector)
    void initialize(SimTK::State& advanced, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman = true); 

    // Assign a random conformation. Time measured in picoseconds
    void propose(SimTK::State& someState, SimTK::Real timestep, int nosteps);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(SimTK::State& someState, SimTK::Real timestep, int nosteps);

    // Get old kinetic energy
    SimTK::Real getOldKE(void) { return this->ke_o; }
    
    // Set kinetic energies    
    void setOldKE(SimTK::Real);

protected:
    SimTK::Real ke_o; // The kinetic energy retained for acc-rej step
    int trackStep;

};

#endif // __HAMMONTECARLOSAMPLER_HPP__

