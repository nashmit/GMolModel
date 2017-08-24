#ifndef __MONTECARLOSAMPLER_HPP__
#define __MONTECARLOSAMPLER_HPP__

#include "Robo.hpp"
#include "Sampler.hpp"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

class Topology;
class IState;

class MonteCarloSampler : public Sampler
{
public:
    // Constructor

    MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem, SimTK::SimbodyMatterSubsystem *argMatter, Topology *argResidue);

    // Destructor

    ~MonteCarloSampler();

    // Helper functions

    void setTVector(SimTK::State& advanced);
    void assignConfFromTVector(SimTK::State& advanced);
    void assignRandomConf(SimTK::State& advanced);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(void);

private:
    SimTK::CompoundSystem *compoundSystem;
    SimTK::SimbodyMatterSubsystem *matter;
    Topology *residue;
    SimTK::Transform *TVector;

    boost::random::mt19937 eng = boost::random::mt19937();
    boost::random::uniform_real_distribution<double> urd =
        boost::random::uniform_real_distribution<double>(0.0, 3.14);


};

#endif // __MONTECARLOSAMPLER_HPP__

