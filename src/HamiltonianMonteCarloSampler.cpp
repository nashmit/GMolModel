/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "HamiltonianMonteCarloSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor
HamiltonianMonteCarloSampler::HamiltonianMonteCarloSampler(SimTK::CompoundSystem *argCompoundSystem,
                                     SimTK::SimbodyMatterSubsystem *argMatter,
                                     SimTK::Compound *argResidue,
                                     SimTK::DuMMForceFieldSubsystem *argDumm,
                                     SimTK::GeneralForceSubsystem *argForces,
                                     SimTK::TimeStepper *argTimeStepper)
    : MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->useFixman = true;  
    this->fix_n = this->fix_o = 0.0;
    this->residualEmbeddedPotential = 0.0;
    trackStep = 0;
}

// Destructor
HamiltonianMonteCarloSampler::~HamiltonianMonteCarloSampler()
{
}

// Compute O(n2) the square root of the mass matrix using using Eige - doesn't work !!n
void HamiltonianMonteCarloSampler::calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper)
{
    int nu = someState.getNU();
    assert((SqrtMUpper.nrow() == nu) && (SqrtMUpper.ncol() == nu) && "calcSqrtMUpper: passed matrix doesn't have nu x nu size.");

    // Calc sqrt(MInv) and put it in SqrtM
    SimTK::Matrix SqrtMInvU(nu, nu);
    this->calcSqrtMInvU(someState, SqrtMInvU);

    // Put sqrt(MInv) in Eigen
    Eigen::MatrixXd EiSqrtMInvU(nu, nu);
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiSqrtMInvU(i, j) = SqrtMInvU(i, j);
        }
    }
    std::cout << "SqrtMInvU: " << std::endl << SqrtMInvU << std::endl;
    std::cout << "EiSqrtMInvU: " << std::endl << EiSqrtMInvU << std::endl;

    // Compute the inverse of sqrt(M) = inv(sqrt(MInv)) with Eigen
    Eigen::MatrixXd EiSqrtMUpper(nu, nu);
    EiSqrtMUpper = EiSqrtMInvU.inverse();

    // Put sqrt(M) back in Simbody style matrix SqrtM
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            SqrtMUpper(i, j) = EiSqrtMUpper(i, j);
        }
    }
    std::cout << "EiSqrtMUpper: " << std::endl << EiSqrtMUpper << std::endl;
    std::cout << "SqrtMUpper: " << std::endl << SqrtMUpper << std::endl;

    // ---- Verify with Eigen ----------
    // Get M
    SimTK::Matrix M(nu, nu);
    matter->calcM(someState, M);
    //std::cout << "M: " << M << std::endl;
    SimTK::Matrix MInv(nu, nu);
    matter->calcMInv(someState, MInv);
    std::cout << "MInv: " << MInv << std::endl;

    // Transfer into Eigen
    Eigen::MatrixXd EiVer(nu, nu);
    for(int i=0; i<nu; i++){
        for(int j=0; j<nu; j++){
            EiVer(i, j) = SqrtMUpper(i, j);
        }
    }
    std::cout << "EiMInv: " << EiSqrtMInvU.transpose() * EiSqrtMInvU << std::endl;
    //std::cout << "EiM: " << EiVer.transpose() * EiVer << std::endl;
    
}

// Calculate O(n2) the square root of the mass matrix inverse
void HamiltonianMonteCarloSampler::calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvU: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = SqrtMInvV[j];
        }
        V[i] = 0;
    }

    // ---- Verify with Eigen ----------
    // Get M
    //SimTK::Matrix MInv(nu, nu);
    //matter->calcMInv(someState, MInv);
    //std::cout << "MInv: " << MInv << std::endl;

    // Transfer into Eigen
    //Eigen::MatrixXd EiMInv(nu, nu);
    //Eigen::MatrixXd EiSqrtMInv(nu, nu);
    //for(int i=0; i<nu; i++){
    //    for(int j=0; j<nu; j++){
    //        EiSqrtMInv(i, j) = SqrtMInv(i, j);
    //    }
    //}
    //std::cout << "Eigen MInv calc: " << EiSqrtMInv.transpose() * EiSqrtMInv << std::endl;

    // Diagonalization
    ////Eigen::EigenSolver<Eigen::MatrixXd> EiSoM(EiM);
    ////Eigen::MatrixXcd EiD = EiSoM.eigenvalues().asDiagonal();
    ////Eigen::MatrixXcd EiV = EiSoM.eigenvectors();

}

void HamiltonianMonteCarloSampler::calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvL: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(j, i) = SqrtMInvV[j];
        }
        V[i] = 0;
    }

    // ---- Verify with Eigen ----------
    // Get M
    //SimTK::Matrix MInv(nu, nu);
    //matter->calcMInv(someState, MInv);
    //std::cout << "MInv: " << MInv << std::endl;

    // Transfer into Eigen
    //Eigen::MatrixXd EiMInv(nu, nu);
    //Eigen::MatrixXd EiSqrtMInv(nu, nu);
    //for(int i=0; i<nu; i++){
    //    for(int j=0; j<nu; j++){
    //        EiSqrtMInv(i, j) = SqrtMInv(i, j);
    //    }
    //}
    //std::cout << "Eigen MInv calc: " << EiSqrtMInv * EiSqrtMInv.transpose() << std::endl;

}

// Set set kinetic energy
void HamiltonianMonteCarloSampler::setSetKE(SimTK::Real inpKE)
{
    this->ke_set = inpKE;
}

// Set old kinetic energy
void HamiltonianMonteCarloSampler::setOldKE(SimTK::Real inpKE)
{
    this->ke_o = inpKE;
}

// Initialize variables (identical to setTVector)
void HamiltonianMonteCarloSampler::initialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman)
{
    randomEngine.seed( std::time(0) );
    //compoundSystem->realizeTopology(); // ELIMINATE REALIZE TOPOLOGY
    //SimTK::State state = compoundSystem->updDefaultState();
    timeStepper->initialize(compoundSystem->getDefaultState());
    setTemperature(argTemperature); // Needed for Fixman

    // Initialize x
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    this->useFixman = argUseFixman;  
    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);
    setOldKE(matter->calcKineticEnergy(someState));
    setSetKE(getOldKE());
    this->etot_o = getOldPE() + getOldKE() + getOldFixman();
    this->etot_set = this->etot_o;
  
    //timeStepper->initialize(someState);
    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////
    HO_randomEngine.seed( std::time(0) );
    // Initialize x
    for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i] = 1.0;}
    HO_PE_set = HO_PE_xprop = HO_PE_x = HarmonicOscillatorPE(HO_x);
    // Initialize v
    HO_InitializeVelocity(HO_v, getTemperature());
    HO_KE_set = HO_KE_xprop = HO_KE_x = HarmonicOscillatorKE(HO_v);
    HO_etot_x = HO_PE_x + HO_KE_x;
    HO_etot_set = HO_etot_xprop = HO_etot_x;
    #endif
}

// Initialize variables (identical to setTVector)
void HamiltonianMonteCarloSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature)
{
    setTemperature(argTemperature); // Needed for Fixman

    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Initialize x
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());
    }

    // Initialize velocities to temperature
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);
    setOldKE(matter->calcKineticEnergy(someState));
    setSetKE(getOldKE());
    this->etot_o = getOldPE() + getOldKE() + getOldFixman();
    this->etot_set = this->etot_o;
  
    //std::cout <<  "Sampler after reinitialize State Cache Info: " << std::endl;
    //PrintSimbodyStateCache(someState);

    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////
    // Initialize x
    for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i] = 1.0;}
    HO_PE_set = HO_PE_xprop = HO_PE_x = HarmonicOscillatorPE(HO_x);
    // Initialize v
    HO_InitializeVelocity(HO_v, getTemperature());
    HO_KE_set = HO_KE_xprop = HO_KE_x = HarmonicOscillatorKE(HO_v);
    HO_etot_x = HO_PE_x + HO_KE_x;
    HO_etot_set = HO_etot_xprop = HO_etot_x;
    #endif
}


// Assign random conformation
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
void HamiltonianMonteCarloSampler::propose(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    //randomEngine.seed(4294653137UL); // for reproductibility

    system->realize(someState, SimTK::Stage::Position);
    // Initialize x - not necessary
    int t = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        TVector[t] = SetTVector[t];
        t++;
    }
    setOldPE(getSetPE());
    setOldFixman(getSetFixman());

    // Assign velocities according to Maxwell-Boltzmann distribution
    // and set Old kinetic energy
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }
    //std::cout << "Q: " << someState.getQ() << std::endl;
    //std::cout << "Before stepTo PE: " << forces->getMultibodySystem().calcPotentialEnergy(someState) << std::endl;
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    //std::cout << "HamiltonianMonteCarloSampler::propose SqrtMInvV: " << SqrtMInvV << std::endl;

    //SimTK::Real temperatureBoost = 2.58; // sqrt(2000/300) : brings temperature from 300 to 1000
    //SqrtMInvV *= sqrtRT * temperatureBoost; // Set stddev according to temperature
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature

    someState.updU() = SqrtMInvV;
    //std::cout << "Before stepTo U: " << someState.getU() << std::endl;

    // Set old kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    setOldKE(matter->calcKineticEnergy(someState));
    this->etot_o = getOldPE() + getOldKE() + getOldFixman();

    // Check priinciple of equipartition of energy
    //SimTK::Vector U(nu);
    //U = SqrtMInvV;
    //SimTK::Vector P(nu);
    //SimTK::Vector elementWisePU(nu);

    //matter->multiplyByM(someState, U, P);
    //for(int i=0; i<nu; i++){
    //    elementWisePU(i) = U(i) * P(i);
    //}
    //std::cout << "HamiltonianMonteCarloSampler::propose p*u: ";
    //for(int i=0; i<nu; i++){
    //    std::cout << " " << elementWisePU(i);
    //}
    //std::cout << std::endl;

    // Get M
    //SimTK::Matrix M(nu, nu);
    //matter->calcM(someState, M);
    //std::cout << "Before stepTo M:" << std::setprecision(3) << std::endl;
    //bool worry_flag = false;
    //for(int i=0; i<nu; i++){
    //    std::cout << "M: ";
    //    for(int j=0; j<nu; j++){
    //        std::cout << M(i, j) << " ";
    //        if( std::isinf(M(i, j)) ){
    //            std::cout << "M(" << i << ", " << j << ") is inf" << std::endl;
    //            worry_flag = true;
    //            return;
    //        }else if( std::isnan(M(i, j)) ){
    //            std::cout << "M(" << i << ", " << j << ") is nan" << std::endl;
    //            worry_flag = true;
    //            return;
    //        }
    //    }
    //    std::cout << std::endl;
    //}

    // Propagate through phase space (integrate)
    //std::cout << "Before stepTo time: " << someState.getTime() << std::endl;
    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));
    //std::cout <<  "Sampler after stepTo State Cache Info: " << std::endl;
    //PrintSimbodyStateCache(someState);
    //std::cout << "After  stepTo time: " << someState.getTime() << std::endl;
    //writePdb(*residue, someState, "pdbs", "sb_", 8, "HMCprop", trackStep);
    ++trackStep;

    //std::cout << "After  stepTo Q: " << someState.getQ() << std::endl;
    //std::cout << "After  stepTo U: " << someState.getU() << std::endl;
    //std::cout << "After  stepTo PE: " << forces->getMultibodySystem().calcPotentialEnergy(someState) << std::endl;

    // Get M
    //matter->calcM(someState, M);
    //std::cout << "After stepTo:" << std::setprecision(3) << std::endl;
    //for(int i=0; i<nu; i++){
    //    std::cout << "M: ";
    //    for(int j=0; j<nu; j++){
    //        std::cout << M(i, j) << " ";
    //    }
    //    std::cout << std::endl;
    //}

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////

    for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i];}
    HO_PE_x = HO_PE_set;
    HO_InitializeVelocity(HO_v, getTemperature());
    HO_KE_x = HarmonicOscillatorKE(HO_v);
    HO_etot_x = HO_PE_x + HO_KE_x;

    HO_VelocityVerlet(timestep, 100);
    HO_PE_xprop = HarmonicOscillatorPE(HO_xprop);
    HO_KE_xprop = HarmonicOscillatorKE(HO_v);
    HO_etot_xprop = HO_PE_xprop + HO_KE_xprop;
    //////////////////////////
    #endif

}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step
void HamiltonianMonteCarloSampler::update(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    propose(someState, timestep, nosteps);

    SimTK::Real pe_o  = getOldPE();
    if(useFixman){
        SimTK::Real fix_o = getOldFixman();
    }
    SimTK::Real ke_o  = getOldKE();

    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 0.0;
    }
    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM
    system->realize(someState, SimTK::Stage::Velocity);
    SimTK::Real ke_n = matter->calcKineticEnergy(someState);

    // Apply Metropolis criterion
    SimTK::Real etot_o, etot_n;
    assert(!isnan(pe_n));
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n;
        etot_o = pe_o + ke_o + fix_o;
    }else{
        etot_n = pe_n + ke_n;
        etot_o = pe_o + ke_o;
    }

    etot_o;
    etot_n;

    //std::cout <<  "Sampler after energies calculations State Cache Info: " << std::endl;
    //PrintSimbodyStateCache(someState);

    std::cout<<std::setprecision(5)<<std::fixed;
    std::cout << "pe_o " << pe_o + getREP() << " ke_o " << ke_o << " fix_o " << fix_o << " rep " << getREP()
        << " pe_n " << pe_n  + getREP() << " ke_n " << ke_n << " fix_n " << fix_n
        //<< " rand_no " << rand_no << " RT " << RT << " exp(-(etot_n - etot_o) " << exp(-(etot_n - etot_o) / RT)
        //<< " etot_n " << etot_n  + getREP() << " etot_o " << etot_o + getREP()
        ;

    //if (( (pe_n + fix_n) < (pe_o + fix_o) ) || (rand_no < exp(-( (pe_n + fix_n) - (pe_o + fix_o) ) / RT))){ // Accept
    //if (( (pe_n) < (pe_o) ) || (rand_no < exp(-( (pe_n) - (pe_o) ) / RT))){ // Accept
    //if ((etot_n < etot_o) || (rand_no < exp(-(etot_n - etot_o)/RT))){ // Accept

    //if(1){ // Always accept
    if ((etot_n < etot_o) || (rand_no < exp(-(etot_n - etot_o)/RT))){ // Accept
        std::cout << " acc 1 " ;
        setSetTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setSetPE(pe_n);
        setSetFixman(fix_n);
        setSetKE(ke_n);
        this->etot_set = getSetPE() + getSetFixman() + getSetKE();
        //someState.updU() = 0.0;
        //setOldKE(0.0);

    }else{ // Reject
        std::cout << " acc 0 " ;
        assignConfFromSetTVector(someState);
        //someState.updU() = 0.0;
        //setOldKE(0.0);

    }

    std::cout << " pe_os " << getSetPE() + getREP() << " ke_os " << getSetKE() << " fix_os " << getSetFixman()
        //<< " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n
        << std:: endl;
    //std::cout << "Number of times the force field was evaluated: " << dumm->getForceEvaluationCount() << std::endl;

    #ifdef HARMONICOSCILLATOR
    /////////////////////////
    // Harmonic oscillator
    /////////////////////////
    //if (( (pe_n) < (pe_o) ) || (rand_no < exp(-( (pe_n) - (pe_o) ) / RT))){ // Accept
    if ((HO_etot_xprop < HO_etot_x) || (rand_no < exp(-(HO_etot_xprop - HO_etot_x)/RT))){ // Accept
        for(int i = 0; i < HO_D; i++){HO_xini[i] = HO_xprop[i];}
        HO_PE_set = HO_PE_xprop;
        HO_KE_set = HO_KE_xprop;
        HO_etot_set = HO_PE_xprop + HO_KE_xprop;
        HO_etot_set = HO_etot_x;
    }else{ // Reject
        for(int i = 0; i < HO_D; i++){HO_x[i] = HO_xini[i];}
    }
    std::cout << "HO 1 pe_os " << HO_PE_set << " ke_os " << HO_KE_set << " etot_os " << HO_etot_set << std::endl;
    #endif

}








