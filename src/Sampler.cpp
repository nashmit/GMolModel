/**@file
Implementation of Sampler class. **/

#include "Sampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

Sampler::Sampler(SimTK::CompoundSystem *argCompoundSystem,
                 SimTK::SimbodyMatterSubsystem *argMatter,
                 //Topology *argResidue,
                 SimTK::Compound *argResidue,
                 SimTK::DuMMForceFieldSubsystem *argDumm,
                 SimTK::GeneralForceSubsystem *argForces,
                 SimTK::TimeStepper *argTimeStepper)
{
    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    this->residue = argResidue;
    this->dumm = argDumm;
    this->forces = argForces;
    this->timeStepper = argTimeStepper;
    this->system = &(matter->getSystem());

    // Thermodynamics
    thermostat = NONE;

    // Harmonic oscillator constants
    for(int i = 0; i < HO_D; i++){HO_x[i] = 1;}
    for(int i = 0; i < HO_D; i++){HO_x0[i] = 0;}
    for(int i = 0; i < HO_D; i++){HO_xini[i] = 1;}
    for(int i = 0; i < HO_D; i++){HO_xprop[i] = 1;}
    for(int i = 0; i < HO_D; i++){HO_v[i] = 0;}
    for(int i = 0; i < HO_D; i++){HO_f[i] = 0;}
    for(int i = 0; i < HO_D; i++){HO_m[i] = 1;}
    for(int i = 0; i < HO_D; i++){HO_a[i] = 0;}
    for(int i = 0; i < HO_D; i++){HO_k[i] = 1;}

    pertAmp = 15.0; // Perturbation amplitude

    #ifdef HARMONICOSCILLATOR
    HO_PE_x = 0.0; // PE at the begining of the proposal
    HO_PE_xprop = 0.0; // PE after the proposal
    HO_PE_set = 0.0; // PE set at the end of the acc-rej step

    HO_KE_x = 0.0;
    HO_KE_xprop = 0.0;
    HO_KE_set = 0.0;

    HO_etot_x = 0.0;
    HO_etot_xprop = 0.0;
    HO_etot_set = 0.0;
    #endif

}

// Destructor
Sampler::~Sampler(){
    ;
}


// Set a thermostat 
void Sampler::setThermostat(Thermostat argThermostat){
    this->thermostat = argThermostat;
}

// Set a thermostat 
void Sampler::setThermostat(std::string argThermostat){
    std::string _thermostat;
    _thermostat.resize(argThermostat.size());
    std::transform(argThermostat.begin(), argThermostat.end(),
        _thermostat.begin(), ::tolower);
   
    this->thermostat = NONE; 

    try{

        if(_thermostat == "andersen"){
            this->thermostat = ANDERSEN;
        }else if(_thermostat == "berendsen"){
            this->thermostat = BERENDSEN;
        }else if(_thermostat == "langevin"){
            this->thermostat = LANGEVIN;
        }else if(_thermostat == "nose_hoover"){
            this->thermostat = NOSE_HOOVER;
        }else{
            throw std::invalid_argument("Thermostat");
        }

    }catch(std::invalid_argument& ia){
        std::cerr << "Invalid argument: " << ia.what() << '\n';
    }

}

// Set a thermostat 
void Sampler::setThermostat(const char *argThermostat){
    std::string _sthermostat = argThermostat;
    std::string _thermostat;
    _thermostat.resize(_sthermostat.size());
    std::transform(_sthermostat.begin(), _sthermostat.end(),
        _thermostat.begin(), ::tolower);
   
    this->thermostat = NONE; 

    try{

        if(_thermostat == "andersen"){
            this->thermostat = ANDERSEN;
        }else if(_thermostat == "berendsen"){
            this->thermostat = BERENDSEN;
        }else if(_thermostat == "langevin"){
            this->thermostat = LANGEVIN;
        }else if(_thermostat == "nose_hoover"){
            this->thermostat = NOSE_HOOVER;
        }else{
            throw std::invalid_argument("Thermostat");
        }

    }catch(std::invalid_argument& ia){
        std::cerr << "Invalid argument: " << ia.what() << '\n';
    }

}

// Get the name of the thermostat
Thermostat Sampler::getThermostat(void){
    return this->thermostat;
}



// Compute mass matrix determinant
SimTK::Real Sampler::calcMassDeterminant(SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;
    matter->calcDetM(state, V, DetV, &D0);
    return D0;
}

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(const SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;
    matter->calcDetM(state, V, DetV, &D0);
    return D0;
}

/** Returns the number of MC trials done by this integrator. **/
int Sampler::getSampleNumber(void){
    return sampleNumber;
}


// Update - to be implemented by every specific sampler

//void Sampler::update(SimTK::State& somState){}

void Sampler::PrintSimbodyStateCache(SimTK::State& someState){
    std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
    for(int i = 0; i < someState.getNumSubsystems(); i++){
        std::cout << " Subsystem Name: " << someState.getSubsystemName(SimTK::SubsystemIndex(i))
            << " Stage: " << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
            << " Version: " << someState.getSubsystemVersion(SimTK::SubsystemIndex(i)) << std::endl;
    }
}



//////////////////////////////////// 
// Harmonic oscillator functions
////////////////////////////////////

// Potential energy function - modifies HO_f
double Sampler::HarmonicOscillatorPE(double *someX)
{
    double pe = 0.0;
    for (int i = 0; i < HO_D; i++){
        HO_f[i] =  -1.0 * HO_k[i] *   (someX[i] - HO_x0[i]); // f = -kx
        pe  += (0.5  * HO_k[i] * ( (someX[i] - HO_x0[i]) * (someX[i] - HO_x0[i]) ) ); // PE = 1/2 k x^2
    }
    return pe;
}

// Kinetic energy function
double Sampler::HarmonicOscillatorKE(double *someV)
{
    double ke = 0.0;
    for (int i = 0; i < HO_D; i++){
        ke += ( 0.5 * HO_m[i] * someV[i] * someV[i] );
    }
    return ke;
}

// Initialize velocities according to T
void Sampler::HO_InitializeVelocity(double *someV, double T)
{
    double HO_RT = T * SimTK_BOLTZMANN_CONSTANT_MD;
    // Diagonal mass matrix case
    for (int i = 0; i < HO_D; i++){
        someV[i]  = HO_gaurand(HO_randomEngine);
        someV[i] *= std::sqrt(HO_RT/HO_m[i]);
    }
}

// Integrate x using Velocity Verlet - modifies xprop
void Sampler::HO_VelocityVerlet(double dt, int nofSteps)
{
    for(int step = 0; step < nofSteps; step++){
        // Compute forces fo = fo(xo) and put it in f
        HarmonicOscillatorPE(HO_xprop);
        // Compute accelerations ao = ao(fo) and put  it in a
        for (int i = 0; i < HO_D; i++){
            HO_a[i] = HO_f[i] / HO_m[i];
        }

        // Compute vhalf  = vhalf(vo, ao) and put in v
        for (int i = 0; i < HO_D; i++){
            HO_v[i] = HO_v[i] + (0.5 * HO_a[i] * dt);
        }

        // Compute xn = xn(xo, vhalf) and put it in x
        for (int i = 0; i < HO_D; i++){
            HO_xprop[i] = HO_xprop[i] + (HO_v[i] * dt);
        }
        // Compute new forces fn = fn(xn) and replace in f
        HarmonicOscillatorPE(HO_xprop);
        // Compute new accelerations an = an(fn) and replace in a
        for (int i = 0; i < HO_D; i++){
            HO_a[i] = HO_f[i] / HO_m[i];
        }
        // Compute new velocities vn = vn(vhalf, an) and replace v
        for (int i = 0; i < HO_D; i++){
            HO_v[i] = HO_v[i] + (0.5 * HO_a[i] * dt);
        }

        //std::cout << "HO_xprop ";
        //for (int i = 0; i < HO_D; i++){
        //    std::cout << HO_xprop[i] << " ";
        //}
        //std::cout << std::endl;

    }

}

