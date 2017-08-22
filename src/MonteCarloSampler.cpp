#include "MonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

MonteCarloSampler::MonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem, SimTK::SimbodyMatterSubsystem *argMatter, Topology *argResidue)
{

    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    this->residue = argResidue;
    TVector = new SimTK::Transform[matter->getNumBodies()];

}

// Destructor

MonteCarloSampler::~MonteCarloSampler()
{

    delete [] TVector;

}

// Transfer coordinates from TVector to compound
 
void MonteCarloSampler::assignConfFromTVector(SimTK::State& advanced)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(advanced, TVector[i]);
    i++;
  }
}

// Assign random conformation
 
void MonteCarloSampler::assignRandomConf(SimTK::State& advanced)
{
boost::random::mt19937 eng = boost::random::mt19937();
boost::random::uniform_real_distribution<double> urd =
   boost::random::uniform_real_distribution<double>(0.0, 3.14);

    eng.seed(4294653137UL);

    std::cout << "MonteCarloSampler state Qs " << advanced.updQ() << std::endl;

    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        SimTK::Real rand_no = urd(eng);

        std::cout << "MonteCarloSampler boostQRealRand " << rand_no << std::endl;

        //mobod.updQ(advanced, rand_no);
        i++;
    }
}



void MonteCarloSampler::update(void){

    //;

}

