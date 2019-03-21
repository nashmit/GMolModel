/**@file
Implementation of HamiltonianMonteCarloSampler class. **/

#include "Robo.hpp"
#include "ParaMolecularDecorator.hpp"

using namespace SimTK;

ParaMolecularDecorator::ParaMolecularDecorator(SimTK::CompoundSystem *argCompoundSystem,
                 SimTK::SimbodyMatterSubsystem *argMatter,
                 SimTK::Compound *argResidue,
                 SimTK::DuMMForceFieldSubsystem *argDumm,
                 SimTK::GeneralForceSubsystem *argForces)
{
    this->compoundSystem = argCompoundSystem;
    this->matter = argMatter;
    this->residue = argResidue;
    this->dumm = argDumm;
    this->forces = argForces;
}

void ParaMolecularDecorator::loadPoint(const Vec3 point)
{
    points.push_back(point);
}

void ParaMolecularDecorator::loadLine(const Vec3 p1, const Vec3 p2)
{
    lines.push_back(std::pair<Vec3, Vec3>(p1, p2));
}

void ParaMolecularDecorator::clearPoints(void)
{
    points.clear();
}

void ParaMolecularDecorator::clearLines(void)
{
    lines.clear();
}

void ParaMolecularDecorator::generateDecorations(const State& someState,
        Array_<DecorativeGeometry>& geometry) 
{
    /*
    for (auto p = points.begin(); p != points.end(); p++) {
        geometry.push_back(DecorativePoint(*p));
    }
    for (auto p = lines.begin(); p != lines.end(); p++) {
        geometry.push_back(DecorativeLine( p->first, p->second ));
        (geometry.back()).setLineThickness(5);
        (geometry.back()).setColor( Vec3(255, 0, 255) );
    }
    */

    // Draw DuMM based geometry
    /*
    for (SimTK::DuMM::BondIndex bIx(0); bIx < dumm->getNumBonds(); ++bIx) {
        const SimTK::DuMM::AtomIndex dAIx1 = dumm->getBondAtom(bIx, 0);
        const SimTK::DuMM::AtomIndex dAIx2 = dumm->getBondAtom(bIx, 1);
        const SimTK::MobilizedBodyIndex mbx1 = dumm->getAtomBody(dAIx1);
        const SimTK::MobilizedBodyIndex mbx2 = dumm->getAtomBody(dAIx2);
        const SimTK::MobilizedBody mobod1 = matter->getMobilizedBody(mbx1);
        const SimTK::MobilizedBody mobod2 = matter->getMobilizedBody(mbx2);

        SimTK::Transform X_GB1 = mobod1.getBodyTransform(someState);
        SimTK::Transform X_GB2 = mobod2.getBodyTransform(someState);

        SimTK::Vec3 p_BS1 = dumm->getAtomStationOnBody(dAIx1);
        SimTK::Vec3 p_GS1 = X_GB1 * p_BS1;

        SimTK::Vec3 p_BS2 = dumm->getAtomStationOnBody(dAIx2);
        SimTK::Vec3 p_GS2 = X_GB2 * p_BS2;

        geometry.push_back(DecorativeLine( p_GS1, p_GS2 ));
        (geometry.back()).setLineThickness(4);
        if( mbx1 == mbx2 ){
            (geometry.back()).setColor( Vec3(0.5, 0.5, 0.5) );
        }else{
            (geometry.back()).setColor( Vec3(255, 0, 255) );
        }
    }
    */


    for (DuMM::AtomIndex daIx(0); daIx < dumm->getNumAtoms(); ++daIx) {
        const SimTK::MobilizedBodyIndex mbx = dumm->getAtomBody(daIx);
        const SimTK::MobilizedBody mobod = matter->getMobilizedBody(mbx);

        SimTK::Transform X_GB = mobod.getBodyTransform(someState);
        SimTK::Vec3 p_BS = dumm->getAtomStationOnBody(daIx);
        SimTK::Vec3 p_GS = X_GB * p_BS;
        SimTK::Transform X_BD(Rotation(), p_GS);

        Real shrink = 0.25;
        //Real shrink = 0.8;
        //Real shrink = 0.5;
        //Real opacity = dumm->getAtomElement(daIx)==1?0.5:1;
        Real opacity = 0.5;
        Real r = dumm->getAtomRadius(daIx);
        if (r<.001) r=0.1; //nm

        geometry.push_back( DecorativeSphere(shrink * r) );
        (geometry.back()).setColor(dumm->getAtomDefaultColor(daIx));
        (geometry.back()).setOpacity(opacity);
        (geometry.back()).setResolution(3);
        (geometry.back()).setTransform(X_BD);

        ///*
        // Text
        std::ostringstream streamObj;
        streamObj << std::fixed;
        streamObj << std::setprecision(3);
        streamObj //<< X_BD.p()[0] << ' ' << X_BD.p()[1] << ' ' 
            << X_BD.p()[2];

        std::string text1 = streamObj.str();
        DecorativeText decorativeText1(text1);
        decorativeText1.setTransform(X_BD);
        decorativeText1.setScaleFactors(SimTK::Vec3(0.02, 0.02, 0.02));
        decorativeText1.setColor(SimTK::Vec3(255, 0, 255));
        geometry.push_back(decorativeText1);
        //*/

    }


    // For portal purposes
    SimTK::Transform G_X_T = residue->getTopLevelTransform();
    //DecorativeBrick topDecorativeBrick;
    //topDecorativeBrick.setTransform(G_X_T);
    //topDecorativeBrick.setScaleFactors(SimTK::Vec3(0.03, 0.03, 0.03));
    //topDecorativeBrick.setColor(SimTK::Vec3(0, 0, 0));
    //geometry.push_back( topDecorativeBrick );
    for (SimTK::Compound::AtomIndex aIx(0); aIx < residue->getNumAtoms(); ++aIx){
        SimTK::Transform T_X_atom =  residue->calcDefaultAtomFrameInCompoundFrame(SimTK::Compound::AtomIndex(aIx));
        SimTK::Transform G_X_atom = G_X_T * T_X_atom;
        // Bricks
        //DecorativeBrick decorativeBrick(SimTK::Vec3(0.03, 0.03, 0.03));
        //decorativeBrick.setTransform(G_X_atom);
        //decorativeBrick.setOpacity(0.5);
        //geometry.push_back( decorativeBrick );

        ///*
        // Text
        std::ostringstream streamObj;
        streamObj << std::fixed;
        streamObj << std::setprecision(3);
        streamObj 
            //<< G_X_atom.p()[0] << ' ' << G_X_atom.p()[1] << ' ' 
            << G_X_atom.p()[2];

        std::string text1 = streamObj.str();
        DecorativeText decorativeText1(text1);
        decorativeText1.setTransform(G_X_atom);
        decorativeText1.setScaleFactors(SimTK::Vec3(0.02, 0.02, 0.02));
        decorativeText1.setColor(SimTK::Vec3(0, 0, 0));
        geometry.push_back(decorativeText1);
        //*/

    }

    /*
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();

        SimTK::Transform G_X_B = mobod.getBodyTransform(someState);
        
        // Bricks
        DecorativeBrick decorativeBrick(SimTK::Vec3(0.03, 0.03, 0.03));
        decorativeBrick.setTransform(G_X_B);
        //decorativeBrick.setColor(SimTK::Vec3(10, 0, 0));
        decorativeBrick.setOpacity(0.5);
        geometry.push_back( decorativeBrick );

        DecorativeFrame decorativeFrame;
        decorativeFrame.setTransform(G_X_B);
        decorativeFrame.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
        decorativeFrame.setLineThickness(4); 
        decorativeFrame.setColor(SimTK::Vec3(255, 0, 0));
        decorativeFrame.setRepresentation(SimTK::DecorativeGeometry::Representation::DrawPoints); 
        geometry.push_back( decorativeFrame );

        if(mbx > 0){
            SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

            // Get all possible transforms
            SimTK::Transform G_X_P = parentMobod.getBodyTransform(someState);
            SimTK::Transform P_X_F = mobod.getInboardFrame(someState);
            SimTK::Transform G_X_F = G_X_P * P_X_F;
            SimTK::Transform F_X_M = mobod.getMobilizerTransform(someState);
            SimTK::Transform G_X_M = G_X_F * F_X_M;

            // F
            DecorativeSphere decorativeSphereF(0.02);
            decorativeSphereF.setColor(SimTK::Vec3(0, 0, 255));
            decorativeSphereF.setOpacity(0.5);
            decorativeSphereF.setTransform(G_X_F);
            geometry.push_back(decorativeSphereF);

            DecorativeFrame decorativeFrameF;
            decorativeFrameF.setTransform(G_X_F);
            decorativeFrameF.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
            decorativeFrameF.setLineThickness(4);
            decorativeFrameF.setColor(SimTK::Vec3(0, 0, 255));
            geometry.push_back( decorativeFrameF );

            std::ostringstream streamObjF;
            streamObjF << std::string("F") + std::to_string(int(mbx)); 
            std::string textF = streamObjF.str();
            DecorativeText decorativeTextF(textF);
            decorativeTextF.setTransform(G_X_F);
            decorativeTextF.setScaleFactors(SimTK::Vec3(0.02, 0.02, 0.02));
            decorativeTextF.setColor(SimTK::Vec3(0, 0, 255));
            geometry.push_back(decorativeTextF);

            // Line between P and F
            DecorativeLine decorativeLinePF(G_X_P.p(), G_X_F.p());
            decorativeLinePF.setLineThickness(3);
            geometry.push_back(decorativeLinePF);

            // M
            DecorativeSphere decorativeSphereM(0.02);
            decorativeSphereM.setColor(SimTK::Vec3(255, 0, 0));
            decorativeSphereM.setOpacity(0.5);
            decorativeSphereM.setTransform(G_X_M);
            geometry.push_back(decorativeSphereM);

            DecorativeFrame decorativeFrameM;
            decorativeFrameM.setTransform(G_X_M);
            decorativeFrameM.setScaleFactors(SimTK::Vec3(0.04, 0.04, 0.04));
            decorativeFrameM.setLineThickness(4);
            decorativeFrameM.setColor(SimTK::Vec3(255, 0, 0));
            geometry.push_back( decorativeFrameM );

            std::ostringstream streamObjM;
            streamObjM << std::string("M") + std::to_string(int(mbx)); 
            std::string textM = streamObjM.str();
            DecorativeText decorativeTextM(textM);
            decorativeTextM.setTransform(G_X_M);
            decorativeTextM.setScaleFactors(SimTK::Vec3(0.02, 0.02, 0.02));
            decorativeTextM.setColor(SimTK::Vec3(255, 0, 0));
            geometry.push_back(decorativeTextM);

            // Line between P and F
            DecorativeLine decorativeLineFM(G_X_F.p(), G_X_M.p());
            decorativeLineFM.setLineThickness(3);
            decorativeLineFM.setColor(SimTK::Vec3(255, 0, 0));
            geometry.push_back(decorativeLineFM);

                //std::cout << "mbx: " << int(mbx) << " B_X_M:"
                //  << mobod.getOutboardFrame(someState);
                //std::cout << "mbx: " << int(mbx) << " M_X_pin:"
                //  << M_X_pin;
                //std::cout << "mbx: " << int(mbx) << " Q:" // << std::endl << ' '
                //  << ((SimTK::MobilizedBody::Free&)mobod).getQ(someState); // << std::endl;

        }
    }
    */

}

ParaMolecularDecorator::~ParaMolecularDecorator()
{
    points.clear();
    lines.clear();
}













