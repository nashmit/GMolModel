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
// Set old kinetic energy
void HamiltonianMonteCarloSampler::setOldKE(SimTK::Real inpKE)
{
    this->ke_o = inpKE;
}

// Initialize variables (identical to setTVector)
void HamiltonianMonteCarloSampler::initialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman)
{
    compoundSystem->realizeTopology();
    //SimTK::State state = compoundSystem->updDefaultState();
    timeStepper->initialize(compoundSystem->getDefaultState());

    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }
  
    system->realize(someState, SimTK::Stage::Position);
    setTemperature(argTemperature); // Needed for Fixman
    setOldPE(getPEFromEvaluator(someState));

    this->useFixman = argUseFixman;  

    if(useFixman){
        setOldFixman(calcFixman(someState));
    }else{
        setOldFixman(1);
    }
    setOldKE(0.0);
  
    randomEngine.seed( std::time(0) );
  
    //timeStepper->initialize(someState);
    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);

}

// Initialize variables (identical to setTVector)
void HamiltonianMonteCarloSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature)
{
    //compoundSystem->realizeTopology();
    //system->realize(someState, SimTK::Stage::Instance);

    //someState.advanceSystemToStage(SimTK::Stage::Instance);
    //someState.invalidateAllCacheAtOrAbove(SimTK::Stage::Instance);

    int nu = someState.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    system->realize(someState, SimTK::Stage::Position);
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    system->realize(someState, SimTK::Stage::Acceleration);
    setOldFixman(calcFixman(someState));

    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::Vec3& vertex = mobod.getBodyOriginLocation(someState);
        TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    setTemperature(argTemperature); // Needed for Fixman

    //system->realize(someState, SimTK::Stage::Dynamics);

    setOldPE(getPEFromEvaluator(someState));

    //setOldFixman(calcFixman(someState));

    setOldKE(0.0);
  
    // After an event handler has made a discontinuous change to the 
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Velocity, false);
}


// Assign random conformation
// In torsional dynamics the first body has 7 Q variables for 6 dofs - one
// quaternion (q) and 3 Cartesian coordinates (x). updQ will return: 
// [qw, qx, qy, qz, x1, x2, x3]
void HamiltonianMonteCarloSampler::propose(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    //randomEngine.seed(4294653137UL); // for reproductibility

    // Assign velocities according to Maxwell-Boltzmann distribution
    // and set Old kinetic energy
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);

    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
        //double z = 0;
        //V[i] = boost::math::pdf(gaurand, z);
    }

    system->realize(someState, SimTK::Stage::Position);
    //std::cout << "Before stepTo Q: " << someState.getQ() << std::endl;
    //std::cout << "Before stepTo PE: " << forces->getMultibodySystem().calcPotentialEnergy(someState) << std::endl;
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    //std::cout << "HamiltonianMonteCarloSampler::propose SqrtMInvV: " << SqrtMInvV << std::endl;
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    //std::cout << "Before stepTo U: " << someState.getU() << std::endl;

    // Set old kinetic energy
    system->realize(someState, SimTK::Stage::Acceleration);
    setOldKE(matter->calcKineticEnergy(someState));

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

}

// The update step in Monte Carlo methods consists in:
// Acception - rejection step
void HamiltonianMonteCarloSampler::update(SimTK::State& someState, SimTK::Real timestep, int nosteps)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    //system->realize(someState, SimTK::Stage::Dynamics);

    propose(someState, timestep, nosteps);

    SimTK::Real pe_o  = getOldPE();
    if(useFixman){
        SimTK::Real fix_o = getOldFixman();
    }
    SimTK::Real ke_o  = getOldKE();

    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 1;
    }
    SimTK::Real pe_n = getPEFromEvaluator(someState); // OPENMM
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

    //std::cout<<std::setprecision(10)<<std::fixed;
    //std::cout << "pe_o " << pe_o << " ke_o " << ke_o << " fix_o " << fix_o
    //    << " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n
    //    //<< " rand_no " << rand_no << " RT " << RT << " exp(-(etot_n - etot_o) " << exp(-(etot_n - etot_o) / RT)
    //    << " etot_n " << etot_n << " etot_o " << etot_o;

    if ((etot_n < etot_o) || (rand_no < exp(-(etot_n - etot_o)/RT))){ // Accept
        std::cout << " 1 " ;
        setTVector(someState);
        //sendConfToEvaluator(); // OPENMM
        setOldPE(pe_n);
        setOldFixman(fix_n);
        someState.updU() = 0.0;
        setOldKE(0.0);
    }else{ // Reject
        std::cout << " 0 " ;
        assignConfFromTVector(someState);
        someState.updU() = 0.0;
        setOldKE(0.0);
    }

    //std::cout << " : pe_o " << getOldPE() << " ke_o " << getOldKE() << " fix_o " << getOldFixman()
    //    << " pe_n " << pe_n << " ke_n " << ke_n << " fix_n " << fix_n << std:: endl;
    //std::cout << "Number of times the force field was evaluated: " << dumm->getForceEvaluationCount() << std::endl;

}








