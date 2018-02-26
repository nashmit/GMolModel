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

    //TO BE DELETED
    prevM = SimTK::Matrix(matter->getNumMobilities(), matter->getNumMobilities());
    //TO BE DELETED

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

        // A bit of study. TO BE DELETED
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        //std::cout << "Mobod " << t+1 << " :" << std::endl;
        for(int k = 0; k < mobod.getNumU(someState); k++){
            SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(k));
        }
        // A bit of study end. TO BE DELETED

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

    // TO BE RESTORED 
    //someState.updU() = SqrtMInvV;
    // TO BE RESTORED 
    //std::cout << "Before stepTo U: " << someState.getU() << std::endl;

    // TO BE DELETED
    //SimTK::Vector myV(nu);
    for(int i = 0; i < nu; i++){someState.updU()[i] = 0;}
    //myV[7] = 0.001;
    //std::cout << "FixmanTorque: " << "myV = " << myV << std::endl;
    kForTheta = matter->getNumMobilities() - 1;
    someState.updU()[kForTheta - 1] = 0.5;
    std::cout << "FixmanTorque: " << "Us = " << someState.getU() << std::endl;
    // TO BE DELETED

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

    // TO BE DELETED
//    SimTK::Vector V1(nu);
//    SimTK::Vector V2(nu);
//    SimTK::Vector V3(nu);
//    // This stores the torques. Torques 3,4,5 of the 1st body are 0
//    SimTK::Vector V4(nu); // This stores the torques.
//    SimTK::Real* D0 = new SimTK::Real(1.0);
//
//    SimTK::Real tinyModification = 1.00;
//    unsigned int nq = someState.getNQ();
//
//        for(int k = 0; k < nq; k++){
//            std::cout << "k = " << k << std::endl;
//            for(int t = 0; t < 5; t++){
//                std::cout << "t = " << t << std::endl;
//                SimTK::Vector Qs = someState.updQ();
//                system->realize(someState, SimTK::Stage::Position);
//                std::cout << "Qs(t) = " << Qs << std::endl;
//                matter->realizeArticulatedBodyInertias(someState);
//                matter->calcDetM(someState, V1, V2, D0);
//                std::cout << "lndetM(t) = " << std::log(*D0) << std::endl;
//
//
//                Qs[k] += tinyModification;
//                system->realize(someState, SimTK::Stage::Position);
//                std::cout << "Qs(t+1) = " << Qs << std::endl;
//
//                system->realize(someState, SimTK::Stage::Dynamics);
//                matter->calcFixmanTorque(someState, V3, V4, D0);
//
//                matter->realizeArticulatedBodyInertias(someState);
//                matter->calcDetM(someState, V1, V2, D0);
//                std::cout << "lndetM(t+1) = " << std::log(*D0) << std::endl;
//
//            }
//
//        //for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
//        //    std::cout << "mbx = " << mbx << std::endl;
//        //    const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
//        //    for(int nqi = 0; nqi < mobod.getNumQ(someState); nqi++){
//        //        std::cout << "nqi = " << nqi << std::endl;
//        //        SimTK::Vector oneMobodQ(someState.getNQ());
//        //        mobod.updOneFromQPartition(someState, nqi, oneMobodQ);
//        //        std::cout << "oneMobodQ = " << oneMobodQ << std::endl;
//        //        system->realize(someState, SimTK::Stage::Dynamics);
//        //        std::cout << "Qs = " << someState.getQ() << std::endl;
//        //        matter->calcFixmanTorque(someState, V1, V2, D0);
//        //    }
//        //}
//
//
//    } // k
//    delete D0;
    // TO BE DELETED

    this->timeStepper->stepTo(someState.getTime() + (timestep*nosteps));

    // TO BE DELETED
    SimTK::Vector V3(nu);
    SimTK::Vector V4(nu);
    SimTK::Real* D0 = new SimTK::Real(1.0);
    //matter->calcFixmanTorque(someState, V3, V4, D0);
    delete D0;
    // TO BE DELETED
 
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

    // TO BE DELETED
    // Eq. 2.11 and 2.14 numerical derivs
    std::cout << "FixmanTorque: " << "Qs = " << someState.getQ() << std::endl;
    int NBODIESplusG = matter->getNumBodies();
    int SPATIAL_DOFSplusG = NBODIESplusG * 6;
    int SPATIAL_DOFS = SPATIAL_DOFSplusG - 6;
    int NU = someState.getNU();
    SimTK::Matrix M(NU, NU);
    SimTK::Matrix MInv(NU, NU);
    SimTK::Matrix dM(NU, NU);
    SimTK::Matrix dMdThetaK(NU, NU);
    SimTK::Matrix MInvdMdThetaK(NU, NU);
    SimTK::Real detM;
    SimTK::Real dDetM;
    SimTK::Real dDetMdThetaK;
    SimTK::Real numDetM;
    SimTK::Real dNumDetM;
    SimTK::Real dNumDetMdThetaK;
    SimTK::Real dThetaK;
    SimTK::Real Trace = 0.0;

    SimTK::Vector V1(NU);
    SimTK::Vector V2(NU);
    SimTK::Real* D0 = new SimTK::Real(1.0);
    matter->calcDetM(someState, V1, V2, D0);
    detM = std::log(*D0);
    delete D0;
    numDetM = std::log(calcNumDetM(someState));
    matter->calcM(someState, M);

    dM = M - prevM;
    dDetM = detM - prevDetM;
    dNumDetM = numDetM - prevNumDetM;
    dThetaK = (someState.getQ()[kForTheta]) - prevThetaK;
    dDetMdThetaK = dDetM / dThetaK;
    dNumDetMdThetaK = dNumDetM / dThetaK;
    dMdThetaK = dM / dThetaK;
    //std::cout << "FixmanTorque: " << "prevThetaK = " << prevThetaK << std::endl;
    //std::cout << "FixmanTorque matrix: " << "prevM = " << prevM << std::endl;
    //std::cout << "FixmanTorque matrix: " << "M = " << M << std::endl;
    //std::cout << "FixmanTorque matrix: " << "dM = " << dM << std::endl;
    std::cout << "FixmanTorque matrix: " << "dMdThetaK = " << dMdThetaK << std::endl;

    for(int i = 0; i < NU; i++){
        for(int j = 0; j < NU; j++){prevM(i, j) = M(i, j);}
    }
    prevThetaK = someState.getQ()[kForTheta];
    prevDetM = detM;
    prevNumDetM = numDetM;

    matter->calcMInv(someState, MInv);

    MInvdMdThetaK = MInv * dMdThetaK;
    for(int i = 0; i < NU; i++){
        Trace += MInvdMdThetaK(i, i);
    }
    //std::cout << "FixmanTorque matrix: " << "MInv = " << MInv << std::endl;
    //std::cout << "FixmanTorque matrix: " << "MInvdMdThetaK = " << MInvdMdThetaK << std::endl;
    //std::cout << "FixmanTorque: " << "dNumDetMdThetaK = " << dNumDetMdThetaK << std::endl;
    //std::cout << "FixmanTorque: " << "dDetMdThetaK = " << dDetMdThetaK << std::endl;
    //std::cout << "FixmanTorque: Trace:" << Trace << std::endl;

    // Eq. 4.11
    // First get Mks and Phi based on 3.4a
    //SimTK::Matrix J(6*NU, NU);
    SimTK::Matrix Jstar;
    SimTK::Matrix J;
    SimTK::Matrix MkTot(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    matter->calcSystemJacobian(someState, J);
    Jstar = ~J;
    SimTK::Matrix M3_4a(NU, NU, 0.0);
    std::cout << "FixmanTorque matrix: " << "J = " << J << std::endl;

    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::SpatialInertia Mk = mobod.getBodySpatialInertiaInGround(someState);
        //std::cout << "FixmanTorque matrix: " << "Mk = " << Mk.toSpatialMat() << std::endl;
        int ti = -1; // 6x6
        int tj = -1; // 6x6
        for(int i = 0; i < 2; i++){ // i for SpatialMatrix
                for(int k = 0; k < 3; k++){ // k for 3x3 
            for(int j = 0; j < 2; j++){ // j for SpatialMatrix
                    for(int l = 0; l < 3; l++){ // l for 3x3
                        tj++;
                        tj = tj % 6;
                        if(!tj){ti++;}
                        if(std::isinf(Mk.toSpatialMat()[i][j][k][l]) || std::isnan(Mk.toSpatialMat()[i][j][k][l])){
                            MkTot[int(mbx) * 6 + ti][int(mbx) * 6 + tj] = 0.0;
                        }else{
                            MkTot[int(mbx) * 6 + ti][int(mbx) * 6 + tj] = Mk.toSpatialMat()[i][j][k][l];
                        }
                    }
                }
            }
        }
    }
    //std::cout << "FixmanTorque matrix: " << "MkTot = " << MkTot << std::endl;
    // Get H and Phi
    SimTK::Matrix H(NU, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Hstar(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Matrix Phi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Phistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    int tj = -1;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        for(int k = 0; k < mobod.getNumU(someState); k++){
            tj++;
            SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(k));
            //std::cout << "mobod " << int(mbx) << " HCol " << k << " = " << HCol << std::endl;
            int ti = -1;
            for(int i = 0; i < 2; i++){
                for(int j = 0; j < 3; j++){
                    ti++;
                    Hstar[int(mbx)*6 + ti][tj] = HCol[i][j];
                }
            }
        }
        
    }

    H = Hstar.transpose(); 

    // Fill Phi
    for(int phi_I = 1; phi_I < NBODIESplusG; phi_I++){
        // Fill diagonal element
        for(int i = 0; i < 6; i++){
                Phistar[phi_I*6 + i][phi_I*6 + i] = 1.0;
        }
        // Fill lower triangle
        for(int phi_J = 1; (phi_J < NBODIESplusG) && (phi_J < phi_I); phi_J++){

            const SimTK::MobilizedBody& Pmobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(phi_I));
            const SimTK::MobilizedBody& Bmobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(phi_J));

            const SimTK::Transform X_GP = Pmobod.getBodyTransform(someState);
            const SimTK::Transform X_GB = Bmobod.getBodyTransform(someState);

            const SimTK::Transform X_PB = ~X_GP * X_GB;
            //const SimTK::Transform X_PB = ~X_GB * X_GP;

            SimTK::Mat33 PhiLCross = SimTK::crossMat(X_PB.p());
            SimTK::Matrix PhiElem(6, 6, 0.0);

            for(int i = 0; i < 6; i++){
                    PhiElem[i][i] = 1.0;
            }
            for(int i = 0; i < 3; i++){
                for(int j = 3; j < 6; j++){
                    PhiElem[i][j] = PhiLCross[i][j - 3];
                }
            }

            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    Phistar[phi_I*6 + i][phi_J*6 + j] = PhiElem[i][j];
                }
            }
    
            //std::cout << "phi " << phi_I << " " << phi_J << " X_PB = " << X_PB << std::endl;
            //std::cout << "phi " << phi_I << " " << phi_J << " X_PB.p = " << X_PB.p() << std::endl;
            //std::cout << "phi " << phi_I << " " << phi_J << " PhiLCross = " << PhiLCross << std::endl;
            std::cout << "phi " << phi_I << " " << phi_J << " PhiElem = " << PhiElem << std::endl;
        }
    }
    //PrintBigMat(Phistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "FixmanTorque matrix: Phistar"); 

    //M3_4a = (~J) * MkTot * J; // Eq 3.4a
    //M3_4a = (~Hstar * ~Phistar) * MkTot * (Phistar * Hstar); // Eq 3.4a
    //std::cout << "M = " << M << std::endl;
    //std::cout << "M3_4a = " << M3_4a << std::endl;
    //std::cout << "M - M3_4a = " << M - M3_4a << std::endl;

    // Get Jain's versions of H and J

    SimTK::Matrix J_Jain(SPATIAL_DOFS, NU, 0.0);
    SimTK::Matrix Phistar_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 1.0);
    SimTK::Matrix Hstar_Jain(SPATIAL_DOFS, NU, 0.0);

    SimTK::Matrix Jstar_Jain(NU, SPATIAL_DOFS, 0.0);
    SimTK::Matrix Phi_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 1.0);
    SimTK::Matrix H_Jain(NU, SPATIAL_DOFS, 0.0);

    SimTK::Matrix MkTot_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);

    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < NU; j++){
            Hstar_Jain[i][j] = Hstar[i + 6][j];
            J_Jain[i][j] = J[i + 6][j];
        }
    }
    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < SPATIAL_DOFS; j++){
            MkTot_Jain[i][j] = MkTot[i + 6][j + 6];
            Phistar_Jain[i][j] = Phistar[i + 6][j + 6];
        }
    }
    PrintBigMat(Phistar_Jain, SPATIAL_DOFS, SPATIAL_DOFS, 3, "FixmanTorque matrix: Phistar_Jain"); 
    //std::cout << "FixmanTorque matrix: " << "J_Jain computed = " << (Phistar_Jain * Hstar_Jain) - J_Jain << std::endl;
    //std::cout << "FixmanTorque matrix: " << "J computed = " << (Phistar * Hstar) - J<< std::endl;
    //std::cout << "verify" << Hstar_Jain  * Hstar_Jain.transpose() << std::endl;

    // Transfer into Eigen
    Eigen::MatrixXd EiHstar_Jain(SPATIAL_DOFS, NU);
    Eigen::MatrixXd EiHstar_JainLeftInv(NU, SPATIAL_DOFS);
    Eigen::MatrixXd EiHstar_JainRightInv(SPATIAL_DOFS, NU);

    for(int i=0; i<SPATIAL_DOFS; i++){
        for(int j=0; j<NU; j++){
            EiHstar_Jain(i, j) = Hstar_Jain[i][j];
        }
    }

    EiHstar_JainLeftInv = (EiHstar_Jain.transpose() * EiHstar_Jain).inverse() * EiHstar_Jain.transpose();
    EiHstar_JainRightInv = EiHstar_Jain.transpose() * (EiHstar_Jain * EiHstar_Jain.transpose()).inverse();

    //std::cout << "verify right inv MJ \n" << EiMJ_Jain * EiMJ_JainRightInv << std::endl;
    //H_Jain = Hstar_Jain.transpose(); // 
    //Phi_Jain = H_Jain * Jstar_Jain;
 
    //std::cout << "Phi_Jain = " << std::setprecision(2) << std::endl;
    for(int i = 0; i < SPATIAL_DOFS; i++){
    //for(int i = 0; i < NU; i++){
        for(int j = 0; j < SPATIAL_DOFS; j++){
        //for(int j = 0; j < NU; j++){
            //std::cout << Phi_Jain[i][j] << " ";
        }
        //std::cout << std::endl;
    }

    // Eq. 4.11
    // Build Hboldz
    int currKForTheta = 0;
    SimTK::Matrix Hbold(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Hboldstar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        for(int k = 0; k < mobod.getNumU(someState); k++){
            currKForTheta++;
            if(currKForTheta == (kForTheta - 1)){
                SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(k));
                std::cout << "mobod " << int(mbx) << " HCol " << k << " = " << HCol << std::endl;
                SimTK::Mat33 crossH = SimTK::crossMat(HCol[0]);
                std::cout << "crossH " << crossH << std::endl;
                for(int i = 0; i < 3; i++){
                    for(int j = 0; j < 3; j++){
                        Hboldstar[int(mbx)*6 + i][int(mbx)*6 + j] = crossH[i][j];
                    }
                }
                for(int i = 0; i < 3; i++){
                    for(int j = 0; j < 3; j++){
                        Hboldstar[int(mbx)*6 + i + 3][int(mbx)*6 + j + 3] = crossH[i][j];
                    }
                }
            }
        }
    }
    Hbold = ~Hboldstar;

    /*
    std::cout << "H = " << std::setprecision(2) << std::endl;
    for(int i = 0; i < NU; i++){
        for(int j = 0; j < SPATIAL_DOFSplusG; j++){
            std::cout << H[i][j] << " "; 
        }
        std::cout << std::endl;
    }
    //std::cout << "Hbold = " << std::setprecision(2) << std::endl;
    std::cout << "MkTot = " << std::setprecision(2) << std::endl;
    for(int i = 0; i < SPATIAL_DOFSplusG; i++){
        for(int j = 0; j < SPATIAL_DOFSplusG; j++){
            //std::cout << Hbold[i][j] << " ";
            std::cout << MkTot[i][j] << " ";
        }
        std::cout << std::endl;
    }
    */

    SimTK::Matrix MPhistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    MPhistar = MkTot * Phistar;

    SimTK::Matrix HboldPhiM(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix MPhistarHbold(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix SqParan(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix dMdThetaK4_11(NU, NU);

    HboldPhiM = Hbold * Phi * MkTot;
    MPhistarHbold = MkTot * Phistar * Hbold;
    SqParan = HboldPhiM - MPhistarHbold;
    dMdThetaK4_11 = Jstar * SqParan * J;

    std::cout << "J = " << std::setprecision(2) << std::endl;
    for(int i = 0; i < SPATIAL_DOFSplusG; i++){
        //for(int j = 0; j < SPATIAL_DOFSplusG; j++){
        for(int j = 0; j < NU; j++){
            std::cout << J[i][j] << " ";
        }
        std::cout << std::endl;
    }

    //std::cout << "FixmanTorque matrix: " << "HboldPhiM = " << HboldPhiM << std::endl;
    //std::cout << "FixmanTorque matrix: " << "MPhistarHbold = " << MPhistarHbold << std::endl;
    //std::cout << "FixmanTorque matrix: " << "SqParan = " << SqParan << std::endl;
    //std::cout << "FixmanTorque matrix: " << "dMdThetaK4_11 = " << dMdThetaK4_11 << std::endl;

    // TO BE DELETED

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








