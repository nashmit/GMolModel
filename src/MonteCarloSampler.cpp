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

// * Set TVector of transforms from mobods * //
void MonteCarloSampler::setTVector(SimTK::State& advanced)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(advanced);
    TVector[i] = mobod.getMobilizerTransform(advanced);
    i++;
  }
}


// Transfer coordinates from TVector to compound
 
void MonteCarloSampler::assignConfFromTVector(SimTK::State& advanced)
{
  int i = 0;
  for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
    mobod.setQToFitTransform(advanced, TVector[i]);
    i++;
  }
}

// Assign random conformation
 
void MonteCarloSampler::assignRandomConf(SimTK::State& advanced)
{
    //eng.seed(4294653137UL);

    std::cout << "MonteCarloSampler state Qs before " << advanced.getQ() << std::endl;

    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SimTK::Real rand_no = urd(eng);
        mobod.setOneQ(advanced, 0, rand_no);
    }

    std::cout << "MonteCarloSampler state Qs after " << advanced.getQ() << std::endl;
}


// Acception - rejection step
void MonteCarloSampler::update(SimTK::State& advanced){
    SimTK::Real rand_no = urd01(eng);
    SimTK::Real RT = getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;

    // Get energies
    SimTK::Real pe_o = getOldPE();
    SimTK::Real pe_n = getPEFromEvaluator(); // Get potential energy from OPENMM

    // Apply Metropolis criterion
    assert(!isnan(pe_n));
    if ((pe_n < pe_o) or (rand_no < exp(-(pe_n - pe_o)/RT))){
        setTVector(advanced);
        writeConfToEvaluator(); // Insert configuratin in OPENMM
        setOldPE(pe_n);
    }else{
        assignConfFromTVector(advanced);
        setOldPE(getOldPE());
    }
}

SimTK::Real MonteCarloSampler::getOldPE(void){return 1.0;}
SimTK::Real MonteCarloSampler::getOldKE(void){}
void MonteCarloSampler::setOldPE(SimTK::Real argPE){}
void MonteCarloSampler::setOldKE(SimTK::Real argKE){}
SimTK::Real MonteCarloSampler::getPEFromEvaluator(void){return 0.0;}
SimTK::Real MonteCarloSampler::getTemperature(void){}
void MonteCarloSampler::writeConfToEvaluator(void){}

