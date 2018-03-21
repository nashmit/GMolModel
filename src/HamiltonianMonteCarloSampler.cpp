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
    randomEngine.seed(4294653137UL); // for reproductibility

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

////    for(int i = 0; i < nu; i++){someState.updU()[i] = 0;}

    //myV[7] = 0.001;
    //std::cout << "FixmanTorque: " << "myV = " << myV << std::endl;
    //kForTheta = matter->getNumMobilities() - 1; // penultimul U

////    kForTheta = someState.getNQ() - 1; // penultimul U
////    someState.updU()[kForTheta - 1] = 10.5;
    //someState.updU() = 0.0;

////    std::cout << "FixmanTorque 0: " << "Qs = " << someState.getQ() << std::endl;
////    std::cout << "FixmanTorque: " << "Us = " << someState.getU() << std::endl;

    // TO BE DELETED

    // Set old kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    //std::cout << "FixmanTorque 1: " << "Qs = " << someState.getQ() << std::endl;
    setOldKE(matter->calcKineticEnergy(someState));
    //std::cout << "FixmanTorque 2: " << "Qs = " << someState.getQ() << std::endl;
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

////    std::cout << "FixmanTorque 3: " << "Qs = " << someState.getQ() << std::endl;
    // TO BE DELETED
////    SimTK::Vector V3(nu);
////    SimTK::Vector V4(nu);
////    SimTK::Real* D0 = new SimTK::Real(1.0);
////    matter->calcFixmanTorque(someState, V3, V4, D0);
////    std::cout << "FixmanTorque 4: " << "Qs = " << someState.getQ() << std::endl;
////    delete D0;
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
////    std::cout << "FixmanTorque 5: " << "Qs = " << someState.getQ() << std::endl;
    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM
    system->realize(someState, SimTK::Stage::Velocity);
    SimTK::Real ke_n = matter->calcKineticEnergy(someState);
////    std::cout << "FixmanTorque 6: " << "Qs = " << someState.getQ() << std::endl;

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
////    std::cout << std::setprecision(2) << "\nFixmanTorque 7: " << "Qs = " << someState.getQ() << std::endl;

    std::cout << " pe_os " << getSetPE() + getREP() << " ke_os " << getSetKE() << " fix_os " << getSetFixman()
        //<< " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n
        << std:: endl;
    //std::cout << "Number of times the force field was evaluated: " << dumm->getForceEvaluationCount() << std::endl;

    // TO BE DELETED
/*
    // Eq. 2.11 and 2.14 numerical derivs
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
    SimTK::Real ThetaK;
    SimTK::Real dThetaK;
    SimTK::Real Trace = 0.0;

    SimTK::Vector V1(NU);
    SimTK::Vector V2(NU);
    SimTK::Real* D0 = new SimTK::Real(1.0);
    matter->calcDetM(someState, V1, V2, D0);
    std::cout << "FixmanTorque 8: " << "Qs = " << someState.getQ() << std::endl;
    detM = std::log(*D0);
    delete D0;
    numDetM = std::log(calcNumDetM(someState));
    matter->calcM(someState, M);

    dM = M - prevM;
    dDetM = detM - prevDetM;
    dNumDetM = numDetM - prevNumDetM;
    ThetaK = (someState.getQ()[kForTheta]);
    dThetaK = ThetaK - prevThetaK;
    dDetMdThetaK = dDetM / dThetaK;
    dNumDetMdThetaK = dNumDetM / dThetaK;
    dMdThetaK = dM / dThetaK;
    std::cout << std::setprecision(10) << "prevThetaK = " << prevThetaK << std::endl;
    std::cout << std::setprecision(10) << "ThetaK = " << ThetaK << std::endl;
    //std::cout << "FixmanTorque matrix: " << "prevM = " << prevM << std::endl;
    //std::cout << "FixmanTorque matrix: " << "M = " << M << std::endl;
    //PrintBigMat(M, NU, NU, 2, "M"); 
    PrintBigMat(dM, NU, NU, 4, "dM");
    std::cout << std::setprecision(10) << "dThetaK " << dThetaK << std::endl << std::setprecision(2);
    //std::cout << "FixmanTorque matrix: " << "dM = " << dM << std::endl;
    //std::cout << "FixmanTorque matrix: " << "dMdThetaK = " << dMdThetaK << std::endl;

    for(int i = 0; i < NU; i++){
        for(int j = 0; j < NU; j++){prevM(i, j) = M(i, j);}
    }
    prevThetaK = ThetaK;
    prevDetM = detM;
    prevNumDetM = numDetM;

    matter->calcMInv(someState, MInv);
    PrintBigMat(MInv, NU, NU, 2, "MInv"); 

    MInvdMdThetaK = MInv * dMdThetaK;
    for(int i = 0; i < NU; i++){
        Trace += MInvdMdThetaK(i, i);
    }
    //std::cout << "FixmanTorque matrix: " << "MInv = " << MInv << std::endl;
    std::cout << "HMC FixmanTorque matrix: " << "MInvdMdThetaK = " << MInvdMdThetaK << std::endl;
    std::cout << "HMC FixmanTorque: " << std::setprecision(10) << "1/2 dNumDetMdThetaK = " << 0.5 * dNumDetMdThetaK << std::endl;
    std::cout << "HMC FixmanTorque: " << "1/2 dDetMdThetaK = " << 0.5 * dDetMdThetaK << std::endl;
    std::cout << "HMC FixmanTorque: 1/2 Trace(MInvdMdThetaK):" << 0.5 * Trace << std::endl << std::setprecision(2);
    std::cout << "FixmanTorque 9: " << "Qs = " << someState.getQ() << std::endl;

    // Eq. 4.11
    // First get Mks and Phi based on 3.4a
    //SimTK::Matrix J(6*NU, NU);
    SimTK::Matrix Jstar;
    SimTK::Matrix J;
    SimTK::Matrix MkTot(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    matter->calcSystemJacobian(someState, J);
    Jstar = ~J;
    SimTK::Matrix M3_4a(NU, NU, 0.0);
    //std::cout << "FixmanTorque matrix: " << "J = " << J << std::endl;

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
    //PrintBigMat(Hstar, SPATIAL_DOFSplusG, NU, 2, "Hstar"); 

    // Fill Phi
    for(int phi_I = 0; phi_I < NBODIESplusG; phi_I++){
        // Fill diagonal element
        for(int i = 0; i < 6; i++){
                Phistar[phi_I*6 + i][phi_I*6 + i] = 1.0;
        }
        // Fill lower triangle
        for(int phi_J = 0; (phi_J < NBODIESplusG) && (phi_J < phi_I); phi_J++){

            const SimTK::MobilizedBody& Jmobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(phi_J));
            const SimTK::MobilizedBody& Imobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(phi_I));
            //const SimTK::MobilizedBody& Imobod = Jmobod.getParentMobilizedBody();

            const SimTK::Transform X_GI = Imobod.getBodyTransform(someState);
            const SimTK::Transform X_GJ = Jmobod.getBodyTransform(someState);

            const SimTK::Transform X_IBM = Imobod.getOutboardFrame(someState);
            const SimTK::Transform X_JBM = Jmobod.getOutboardFrame(someState);

            // X_PB = X_PF(getInboardFrame) * X_FM(getMobilizerTransform) * X_MB(getOutboardFrame)
            //const SimTK::Transform X_IJ = Jmobod.getDefaultInboardFrame() * Jmobod.getMobilizerTransform(someState) * ~Jmobod.getDefaultOutboardFrame();
            //const SimTK::Transform X_IJ = Jmobod.getInboardFrame(someState) * Jmobod.getMobilizerTransform(someState) * ~Jmobod.getOutboardFrame(someState);

            const SimTK::Transform X_IFM = Imobod.getMobilizerTransform(someState);
            const SimTK::Transform X_JFM = Jmobod.getMobilizerTransform(someState);

            // X_GP.R() * X_PB.p();
            //SimTK::Mat33 PhiLCross = SimTK::crossMat( X_GI.R() * X_IJ.p() );
            //SimTK::Mat33 PhiLCross = SimTK::crossMat( X_GJ.p() - X_GI.p() );
            SimTK::Mat33 PhiLCross = SimTK::crossMat( X_GJ.p() - X_GI.p() );

            SimTK::Matrix PhiElem(6, 6, 0.0);

            for(int i = 0; i < 6; i++){
                    PhiElem[i][i] = 1.0;
            }
            for(int i = 3; i < 6; i++){
                for(int j = 0; j < 3; j++){
                    PhiElem[i][j] = PhiLCross[i - 3][j];
                }
            }

            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    Phistar[phi_I*6 + i][phi_J*6 + j] = PhiElem[i][j];
                    //Phistar[Imobod.getMobilizedBodyIndex()*6 + i][Jmobod.getLevelInMultibodyTree()*6 + j] = PhiElem[i][j];
                }
            }
    
            //std::cout << "phi " << phi_I << " " << phi_J << " X_IJ = " << X_IJ << std::endl;
            //std::cout << "phi " << phi_I << " " << phi_J << " X_IJ.p = " << X_IJ.p() << std::endl;
            //std::cout << "phi tree levels " << Imobod.getLevelInMultibodyTree() << " " << Jmobod.getLevelInMultibodyTree() << std::endl;
            //std::cout << "phi mbx " << Imobod.getMobilizedBodyIndex() << " " << Jmobod.getMobilizedBodyIndex() << std::endl;
            //std::cout << "phi I J " << phi_I << " " << phi_J << " PhiElem = " ; //<< PhiLCross << std::endl;
            //PrintBigMat(PhiElem, 6, 6, 3, "");
            //std::cout << "phi " << phi_I << " " << phi_J << " PhiElem = " << PhiElem << std::endl;
        }
    }
    //PrintBigMat(Phistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "FixmanTorque matrix: Phistar"); 
    Phi = Phistar.transpose();

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

    Phi_Jain = Phistar_Jain.transpose();
    H_Jain = Hstar_Jain.transpose();
    Jstar_Jain = H_Jain.transpose();

    //PrintBigMat(MkTot_Jain * MkTot_Jain.transpose(), SPATIAL_DOFS, SPATIAL_DOFS, 2, "MkTot * MkTotstar");
    //PrintBigMat(MkTot_Jain.transpose() * MkTot_Jain, SPATIAL_DOFS, SPATIAL_DOFS, 2, "MkTotstar * MkTot");
    //PrintBigMat(Phistar_Jain * Phi_Jain, SPATIAL_DOFS, SPATIAL_DOFS, 2, "Phistar * Phi");
    //PrintBigMat(Phi_Jain * Phistar_Jain, SPATIAL_DOFS, SPATIAL_DOFS, 2, "Phi * Phistar");
    //PrintBigMat(Hstar_Jain * H_Jain, SPATIAL_DOFS, SPATIAL_DOFS, 2, "Hstar * H");
    //PrintBigMat(H_Jain * Hstar_Jain, NU, NU, 2, "H * Hstar");

    //PrintBigMat((H_Jain * Phi_Jain * MkTot_Jain * Phistar_Jain * Hstar_Jain) - M, NU, NU, 2, "Verify 3.4a");
    //PrintBigMat((H * Phi * MkTot * Phistar * Hstar) - M, NU, NU, 2, "Verify 3.4a");
    //PrintBigMat(Phistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Phistar"); 
    //PrintBigMat(Phistar_Jain, SPATIAL_DOFS, SPATIAL_DOFS, 2, "Phistar_Jain"); 
    //std::cout << "FixmanTorque matrix: " << "J_Jain computed = " << (Phistar_Jain * Hstar_Jain) - J_Jain << std::endl;
    //std::cout << "FixmanTorque matrix: " << "J computed - J = " << (Phistar_Jain * Hstar_Jain) - J_Jain<< std::endl;
    //std::cout << "verify" << Hstar_Jain  * Hstar_Jain.transpose() << std::endl;

    // Check Eq. 417
    SimTK::Matrix P(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix P_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);
    SimTK::Matrix D(NU, NU, 0.0);
    SimTK::Matrix DInv(NU, NU, 0.0);
    SimTK::Matrix D_Jain(NU, NU, 0.0);
    SimTK::Matrix G(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Matrix G_Jain(SPATIAL_DOFS, NU, 0.0);

    for (SimTK::MobilizedBodyIndex mbx(0); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::ArticulatedInertia Pk = matter->getArticulatedBodyInertia(someState, mbx);
        SimTK::Mat<6, 6> PkMat66;
        SpatialMat2Mat66(Pk.toSpatialMat(), PkMat66);
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                P[int(mbx) * 6 + i][int(mbx) * 6 + j] = PkMat66[i][j];
            }
        }
    }
    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < SPATIAL_DOFS; j++){
            P_Jain[i][j] = P[i + 6][j + 6];
        }
    }

    D = H * P * Hstar;
    D_Jain = H_Jain * P_Jain * Hstar_Jain;
    NumericalInverse(D, DInv, NU, NU);
    G = P * Hstar * DInv;
    G_Jain = P_Jain * Hstar_Jain * DInv;
    
    //PrintBigMat(P, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "P");
    //PrintBigMat(P_Jain, SPATIAL_DOFS, SPATIAL_DOFS, 2, "P_Jain");
    //PrintBigMat(D, NU, NU, 2, "D");
    //PrintBigMat(D_Jain, NU, NU, 2, "D_Jain");
    //PrintBigMat(DInv, NU, NU, 2, "DInv");
    //PrintBigMat(G, SPATIAL_DOFSplusG, NU, 2, "G");
    //PrintBigMat(G_Jain, SPATIAL_DOFS, NU, 2, "G_Jain");

    SimTK::Matrix EtaPhistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPhistar_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);
    SimTK::Matrix EtaPhi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPhi_Jain(SPATIAL_DOFS, SPATIAL_DOFS, 0.0);

    SimTK::Matrix K(SPATIAL_DOFSplusG, NU, 0.0);
    SimTK::Matrix K_Jain(SPATIAL_DOFS, NU, 0.0);

    for(int i = 1; i < NBODIESplusG; i++){
        for(int k = 0; k < 6; k++){
            for(int l = 0; l < 6; l++){
                EtaPhistar[i*6 + k][(i*6 - 6) + l] = Phistar[i*6 + k][(i*6 - 6) + l];
            }
        }
    }
    for(int i = 0; i < SPATIAL_DOFS; i++){
        for(int j = 0; j < SPATIAL_DOFS; j++){
            EtaPhistar_Jain[i][j] = EtaPhistar[i + 6][j + 6];
        }
    }

    EtaPhi = EtaPhistar.transpose();
    EtaPhi_Jain = EtaPhistar_Jain.transpose();
    K = EtaPhi * G;
    K_Jain = EtaPhi_Jain * G_Jain;

    //PrintBigMat(EtaPhistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "EtaPhistar"); 
    //PrintBigMat(K, SPATIAL_DOFSplusG, NU, 2, "K"); 

    //SimTK::Matrix I(NU, NU, 0.0);
    //for(int i = 0; i < NU; i++){
    //    I[i][i] = 1.0;
    //}
    //SimTK::Matrix MSqParan(NU, NU, 0.0);
    //MSqParan = I + (H * Phi * K);
    //PrintBigMat( (H * Phi) * MkTot * (Phistar * Hstar), NU, NU, 2, "Mcheck 3.4a" ); // Eq 3.4a
    //PrintBigMat( MSqParan * D * MSqParan.transpose(), NU, NU, 2, "Mcheck 3.5a 1");

    // Check EtaPhistar 1992 Rodriguez eq. 0.6
    //SimTK::Matrix I(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    //SimTK::Matrix I_EtaPhistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    //SimTK::Matrix I_EtaPhistarInv(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    //for(int i = 0; i < SPATIAL_DOFSplusG; i++){
    //    I[i][i] = 1.0;
    //}
    //I_EtaPhistar = I - EtaPhistar;
    //NumericalInverse(I_EtaPhistar, I_EtaPhistarInv, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG);
    //PrintBigMat(I_EtaPhistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "I_EtaPhistar"); 
    //PrintBigMat(Phistar - I_EtaPhistarInv, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Check EtaPhistar"); 

    // Build Psi
    // don't copy for lines
    SimTK::Matrix Psi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Psistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPsi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix EtaPsistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    EtaPsi = EtaPhi - (K * H);
    EtaPsistar = EtaPsi.transpose();

    // Fill Psi
    Psi = EtaPsi;
    for(int i = 0; i < NBODIESplusG; i++){
        for(int k = 0; k < 6; k++){
            Psi[i*6 + k][i*6 + k] = 1.0;
        }
        for(int j = i + 2; j < NBODIESplusG; j++){
            SimTK::Mat66 Psi_IJ, Psi_IJm, Psi_JmJ;
            for(int k = 0; k < 6; k++){
                for(int l = 0; l < 6; l++){
                    Psi_IJ[ k][l] = Psi[  (i*6)      + k ][  (j*6)      + l];
                    Psi_IJm[k][l] = Psi[  (i*6)      + k ][ ((j*6) - 6) + l];
                    Psi_JmJ[k][l] = Psi[ ((j*6) - 6) + k ][  (j*6)      + l];
                }
            }

            Psi_IJ = Psi_IJm * Psi_JmJ;

            for(int k = 0; k < 6; k++){
                for(int l = 0; l < 6; l++){
                    Psi[ (i*6) + k ][ (j*6) + l] = Psi_IJ[k][l];
                }
            }
            
        }
    }
    Psistar = Psi.transpose();

    // Omega
    SimTK::Matrix Omega(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    Omega = Psistar * Hstar * DInv * H * Psi;

    // Check eq. 4.14b
    //PrintBigMat(Phi * MkTot * Omega, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Phi M Omega");
    //PrintBigMat((Phi - Psi) + (P * Omega), SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "(Phi - Psi) + P Omega");
    //PrintBigMat(Phi, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Phi");
    //PrintBigMat(Psistar, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Psistar");

    // Build Upsilon using Jain 1997 recursion
    SimTK::Matrix Upsilon(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Mat66 UpsK(0.0);
    SimTK::Mat66 PsistarKKm(0.0);
    SimTK::Mat66 PsiKKm(0.0);
    SimTK::Mat66 DInv0(0.0);
    SimTK::Mat11 DInvK(0.0);
    SimTK::Mat66 Hstar0(0.0);
    SimTK::Mat66 H0(0.0);
    SimTK::Mat61 HstarK(0.0);
    SimTK::Mat16 HK(0.0);

    UpsK = 0;
    for(int k = 1; k < (NBODIESplusG); k++){
        //std::cout << "BODY K = " << k << std::endl;
        // Get PsistarKKm PsiKKm DInvK
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                //std::cout << "PsistarKKm i j " << i << " " << j << " PsistarKKm[][] " << ((k + 1)*6) + i << " " << (k*6) + j << std::flush << std::endl; 
                PsistarKKm[i][j] = Psistar[ (k*6) + i ][ ((k-1)*6) + j ];
                //PsiKKm[i][j] = Psi[ ((k + 1)*6) + i ][ (k*6) + j ];
            }
        }
        PsiKKm = PsistarKKm.transpose();
        //PrintBigMat(PsistarKKm, 6, 6, 2, "PsistarKKm");
        //PrintBigMat(PsiKKm, 6, 6, 2, "PsiKKm");

        // Get HstarK and HK
        if(k == 0){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " = 0 " << std::flush << std::endl; 
                    //std::cout << "DInv0 i j " << i << " " << j << " = 0 " << std::flush << std::endl; 
                    Hstar0[i][j] =  0.0;;
                    DInv0[i][j] =  0.0;;
                }
            }
            H0 = Hstar0.transpose();
        }else if(k == 1){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " Hstar[][] " << (k*6) + i << " " << j << std::flush << std::endl; 
                    //std::cout << "DInv0 i j " << i << " " << j << " DInv[][] " << ((k-1)*6) + i << " " << j << std::flush << std::endl; 
                    Hstar0[i][j] =  Hstar[ (k*6) + i ][ j ];
                    DInv0[i][j] =  DInv[ ((k-1)*6) + i ][ j ];
                }
            }
            H0 = Hstar0.transpose();
        }else if(k > 1){
            //std::cout << "DInvK i j " << 0 << " " << 0 << " DInv " << 5 + (k - 1) << " " << 5 + (k - 1) << std::flush << std::endl;
            DInvK[0][0] =  DInv[ 5 + (k - 1) ][5 + (k - 1)];
            for(int i = 0; i < 6; i++){
                //std::cout << "HstarK i j " << i << " " << 0 << " Hstar " << (6*k) + i << " " << 5 + (k - 1) << std::flush << std::endl;
                HstarK[i][0] =  Hstar[ (6*k) + i ][5 + (k - 1)];
            }
            HK = HstarK.transpose();
        }

        if(k <= 1){
            //std::cout << "Hstar0" << Hstar0 << std::endl;
            //std::cout << "H0" << H0 << std::endl;
            //std::cout << "DInv0" << DInv0 << std::endl;
            UpsK = (PsistarKKm * UpsK * PsiKKm) + (Hstar0 * DInv0 * H0);
        }else{
            //std::cout << "HstarK" << HstarK << std::endl;
            //std::cout << "HK" << HK << std::endl;
            //std::cout << "DInvK" << DInvK << std::endl;
            UpsK = (PsistarKKm * UpsK * PsiKKm) + (HstarK * DInvK * HK);
        }

        //std::cout << "Upsilon(" << k << ") = " << std::endl;
        //std::cout << UpsK << std::endl;

        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                Upsilon[(k*6) + i][(k*6) + j] = UpsK[i][j];
            }
        }

    }

    
    // Check 4.15
    //SimTK::Matrix IdPsi(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    //for(int i = 0; i < SPATIAL_DOFSplusG; i++){IdPsi[i][i] = 1.0;}
    //SimTK::Matrix PsiBar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    //SimTK::Matrix PsiBarstar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);

    //PsiBar = Psi - IdPsi;
    //PsiBarstar = PsiBar.transpose();

    //PrintBigMat(Omega, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Omega");
    //PrintBigMat(Upsilon + (PsiBarstar * Upsilon) + (Upsilon * PsiBar), SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Omega(Upsilon)");

    // Check 4.17
    PrintBigMat(P, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "P");
    PrintBigMat(Upsilon, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Upsilon");
    PrintBigMat(Hstar, SPATIAL_DOFSplusG, NU, 2, "Hstar"); 

    SimTK::Mat66 PKK(0.0);
    SimTK::Mat66 UKK(0.0);
    SimTK::Mat66 Hi(0.0);
    SimTK::Vec3 RotVec(0.0);
    SimTK::Vec3 RotVec_G(0.0);
    SimTK::Mat33 crossH(0.0);
    SimTK::Mat66 PYH(0.0);

    for(int k = 1; k < (NBODIESplusG); k++){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(k));
        std::cout << "BODY K = " << k << std::endl;
        // Get PsistarKKm PKK
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                //std::cout << "PKK i j " << i << " " << j << " P[][] " << (k*6) + i << " " << (k*6) + j << std::flush << std::endl; 
                PKK[i][j] = P[ (k*6) + i ][ (k*6) + j ];
                UKK[i][j] = Upsilon[ (k*6) + i ][ (k*6) + j ];
            }
        }

        // Get HstarK
        if(k == 0){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " = 0 " << std::flush << std::endl; 
                    Hstar0[i][j] =  0.0;
                }
            }
        }else if(k == 1){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " Hstar[][] " << (k*6) + i << " " << j << std::flush << std::endl; 
                    Hstar0[i][j] =  Hstar[ (k*6) + i ][ j ];
                }
            }
        }else if(k > 1){
            for(int i = 0; i < 6; i++){
                //std::cout << "HstarK i j " << i << " " << 0 << " Hstar " << (6*k) + i << " " << 5 + (k - 1) << std::flush << std::endl;
                HstarK[i][0] =  Hstar[ (6*k) + i ][5 + (k - 1)];
            }
        }

        if(k <= 1){
            for(int ExtDof = 0; ExtDof < 6; ExtDof++){
                for (int i = 0; i < 3; i++){
                    RotVec[i] = Hstar0[i][ExtDof];
                }
                const SimTK::Rotation R_GB = mobod.getBodyRotation(someState);
                RotVec_G = R_GB * RotVec;
                std::cout << "RotVec_G: " << RotVec_G << std::endl;
    
                crossH = SimTK::crossMat(RotVec_G);
                std::cout << "crossH: " << crossH << std::endl;
                Hi = SimTK::Mat66(0.0);
    
                for (int i = 0; i < 3; i++){
                    for (int j = 0; j < 3; j++){
                        Hi[i][j] = Hi[i +3][j + 3] = crossH[i][j];
                    }
                }
                std::cout << "Hi: " << Hi << std::endl;
    
                PYH = PKK * UKK * Hi;
                std::cout << std::setprecision(10) << "Tr{PYH} ExtDof " << ExtDof << " = " << PYH.trace() << std::endl << std::setprecision(2);
            }
        }else{
            for (int i = 0; i < 3; i++){
                RotVec[i] = HstarK[i][0];
            }
            const SimTK::Rotation R_GB = mobod.getBodyRotation(someState);
            RotVec_G = R_GB * RotVec;
            std::cout << "RotVec_G: " << RotVec_G << std::endl;

            crossH = SimTK::crossMat(RotVec_G);
            std::cout << "crossH: " << crossH << std::endl;
            Hi = SimTK::Mat66(0.0);

            for (int i = 0; i < 3; i++){
                for (int j = 0; j < 3; j++){
                    Hi[i][j] = Hi[i +3][j + 3] = crossH[i][j];
                }
            }
            std::cout << "Hi: " << Hi << std::endl;

            PYH = PKK * UKK * Hi;
            std::cout << std::setprecision(10) << "Tr{PYH} = " << PYH.trace() << std::endl << std::setprecision(2);
        }

    }
    // Eq. 4.11
*/
    /*
    // Build Hbold
    int currKForTheta = 0;
    SimTK::Matrix Hbold(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix Hboldstar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < NBODIESplusG; ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        for(int mobodMobility = 0; mobodMobility < mobod.getNumU(someState); mobodMobility++){
            currKForTheta++;
            //std::cout << "mbx mobodMobility currK kForTheta: " << int(mbx) << " " << mobodMobility << " " << currKForTheta << " " << kForTheta << std::endl;
            if(currKForTheta == (kForTheta)){
                SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(mobodMobility));
                //std::cout << "mobod " << int(mbx) << " HCol " << mobodMobility << " = " << HCol << std::endl;
                SimTK::Mat33 crossH = SimTK::crossMat(HCol[0]);

                SimTK::SpatialMat Hi(0);
                Hi[0][0] = crossH;
                Hi[1][1] = crossH;
                std::cout << "Hi for mobod " << int(mbx) << ": " << Hi << std::endl;
                //std::cout << "crossH " << crossH << std::endl;
                for(int i = 0; i < 3; i++){
                    for(int j = 0; j < 3; j++){
                        Hbold[int(mbx)*6 + i][int(mbx)*6 + j] = crossH[i][j];
                    }
                }
                for(int i = 0; i < 3; i++){
                    for(int j = 0; j < 3; j++){
                        Hbold[(int(mbx))*6 + i + 3][(int(mbx))*6 + j + 3] = crossH[i][j];
                    }
                }
                break;
            }
        }
    }

    Hboldstar = Hbold.transpose();
    //PrintBigMat(Hbold, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "Hbold:");

    // Check Eq. 411
    SimTK::Matrix MPhistar(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    MPhistar = MkTot * Phistar;

    SimTK::Matrix HboldPhiM(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix MPhistarHbold(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix SqParan(SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 0.0);
    SimTK::Matrix dMdThetaK4_11(NU, NU);

    HboldPhiM = Hbold * Phi * MkTot;
    MPhistarHbold = MkTot * Phistar * Hbold;

    SqParan = HboldPhiM - MPhistarHbold;

    dMdThetaK4_11 = H * Phi * SqParan * Phistar * Hstar;

    //PrintBigMat(HboldPhiM, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "HboldPhiM");
    //PrintBigMat(MPhistarHbold, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "MPhistarHbold");

    //PrintBigMat(SqParan, SPATIAL_DOFSplusG, SPATIAL_DOFSplusG, 2, "SqParan");
    //PrintBigMat(dMdThetaK, NU, NU, 4, "dMdThetaK");
    //PrintBigMat(dMdThetaK4_11, NU, NU, 4, "dMdThetaK4_11");

    //SimTK::Matrix I(NU, NU, 0.0);
    //for(int i = 0; i < NU; i++){
    //    I[i][i] = 1.0;
    //}
    //SimTK::Matrix MSqParan(NU, NU, 0.0);
    //MSqParan = I - (H * Psi * K);
    //PrintBigMat( MSqParan.transpose() * DInv * MSqParan, NU, NU, 2, "MInvcheck 3.5c");
    */

    /*
    SimTK::Mat66 IK;
    for(int i = 0; i < 6; i++){IK[i][i] = 1.0;}
    SimTK::Mat66 Hstar0;
    SimTK::Mat66 H0;
    SimTK::Mat61 HstarK;
    SimTK::Mat16 HK;

    for(int k = 1; k < NBODIESplusG; k++){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(k));
        SimTK::SpatialVec HCol = mobod.getHCol(someState, SimTK::MobilizerUIndex(k));

        SimTK::Mat66 PsiK, GK;

        for(int i = 0; i < 6; i++){
            for(int j = 0; j < 6; j++){
                PsiK[i][j] = Phi[       (k * 6) + i ][ ((k - 1) * 6) + j ];
                GK[i][j]   =   G[ ((k - 1) * 6) + i ][ ((k - 1) * 6) + j ];
            }
        }
 
        if( (k - 1) <= 1){
            for(int i = 0; i < 6; i++){
                for(int j = 0; j < 6; j++){
                    //std::cout << "Hstar0 i j " << i << " " << j << " Hstar[][] " << ((k - 1) * 6) + i << " " << j << std::endl; 
                    Hstar0[i][j] =  Hstar[ ((k - 1) * 6) + i ][ j ];
                }
            }
        }else{
            for(int i = 0; i < 6; i++){
                //std::cout << "HstarK i j " << i << " " << 0 << " HstarK " << (6 * (k - 1)) + i << " " << 5 + (k - 2) << std::endl;
                HstarK[i][0] =  Hstar[ (6 * (k - 1)) + i ][5 + (k - 2)];
            }
        }

        if( (k - 1) <= 1){
            //std::cout << "Hstar0 " << Hstar0 << std::endl;
            H0 = Hstar0.transpose();       
            PsiK = PsiK * (IK - (GK * H0));
        }else{
            HK = HstarK.transpose();       
            std::cout << "HstarK " << HstarK << std::endl;
            PsiK = PsiK * (IK - (GK * HK));
        }

    }
    */



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








