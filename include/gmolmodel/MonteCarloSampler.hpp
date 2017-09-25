#ifndef __MONTECARLOSAMPLER_HPP__
#define __MONTECARLOSAMPLER_HPP__

#include "Robo.hpp"
#include "Sampler.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

class Topology;
class IState;

class MonteCarloSampler : public Sampler
{
public:
    // Constructor

    MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                      SimTK::SimbodyMatterSubsystem *argMatter,
                      Topology *argResidue,
                      SimTK::TimeStepper *argTimeStepper);

    // Destructor

    ~MonteCarloSampler();

    // Simulation temperature related

    SimTK::Real getTemperature(void);
    void setTemperature(SimTK::Real);

    // Store/restore the configuration from the internal TVector variable

    void setTVector(const SimTK::State& advanced);
    void assignConfFromTVector(SimTK::State& advanced);

    // Assign a random conformation

    void propose(SimTK::State& advanced);

    // Store/restore potential energy

    SimTK::Real getOldPE(void);
    SimTK::Real getOldKE(void);
    void setOldPE(SimTK::Real argPE);
    void setOldKE(SimTK::Real argKE);

    // Set/get Fixman potential

    void setOldFixman(SimTK::State& someState);
    SimTK::Real getOldFixman(SimTK::State& someState);

    // Evaluate the potential energy at current state

    SimTK::Real getPEFromEvaluator(void); 

    // Compute Fixman potential

    SimTK::Real calcFixman(SimTK::State& someState);

    // Send configuration to an external evaluator

    void sendConfToEvaluator(void);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation

    void update(SimTK::State&);

protected:
    SimTK::Transform *TVector;
    SimTK::Real pe_o;
    SimTK::Real temperature;
    
    SimTK::Real fix_o, fix_n;
 
    // Random number generators - not sure if I need two
    // Needs testing

    boost::random::mt19937 randomEngine = boost::random::mt19937();

    boost::random::uniform_real_distribution<double> uniformRealDistribution_0_2pi =
        boost::random::uniform_real_distribution<double>(SimTK::Zero, 2*SimTK::Pi);

    boost::random::uniform_real_distribution<double> uniformRealDistribution_mpi_pi =
        boost::random::uniform_real_distribution<double>((-1)*SimTK::Pi, SimTK::Pi);

    boost::random::uniform_real_distribution<double> uniformRealDistribution =
        boost::random::uniform_real_distribution<double>(SimTK::Zero, SimTK::One);

    boost::random::uniform_real_distribution<double> uniformRealDistribution_m1_1 =
        boost::random::uniform_real_distribution<double>((-1)*SimTK::One, SimTK::One);

    boost::normal_distribution<> gaurand = boost::normal_distribution<>(0.0, 1.0);

};

#endif // __MONTECARLOSAMPLER_HPP__

