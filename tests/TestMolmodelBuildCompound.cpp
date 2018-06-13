/* -------------------------------------------------------------------------- *
 *                 SimTK Molmodel Example: Simple Protein                     *
 * -------------------------------------------------------------------------- *
 * This is the first example from the Molmodel User's Guide. It creates a     *
 * small protein (a five-residue peptide), simulates it and generates a live  *
 * animation while it is running.                                             *
 *                                                                            *
 * Authors: Christopher Bruns, Michael Sherman                                *
 * -------------------------------------------------------------------------- */

#include "Molmodel.h"
#include <iostream>
#include <exception>
using namespace SimTK;

class  QuadrivalentAtom2 : public Compound::SingleAtom {
public:
    QuadrivalentAtom2(
        const Compound::AtomName& atomName, ///< name for new atom
        const Element& element ///< chemical element for new atom
        )
        : Compound::SingleAtom(atomName, element)
    {
        static const Angle TetrahedralAngle = 109.47 * Deg2Rad;

        // BondCenter1 dihedral will be relative to BondCenter2
        // Using Rotation() constructor that takes x axis and approximate y axis
        addFirstBondCenter( "bond1", atomName );
        // bond centers 2, 3, and 4 dihedrals will be relative to bond1
        addSecondBondCenter( "bond2", atomName, TetrahedralAngle );
        addLeftHandedBondCenter( "bond3", atomName, TetrahedralAngle, TetrahedralAngle );
        addRightHandedBondCenter( "bond4", atomName, TetrahedralAngle, TetrahedralAngle );

        // Choice of inboard bond may differ from bond priority - user may change this
        setInboardBondCenter("bond1");

        // default dihedral should be left blank, so it could be
        // defined from the other side.
        // the ultimate default will be 180 degrees anyway.
        // setDefaultInboardDihedralAngle(180*Deg2Rad);

        setCompoundName("QuadrivalentAtom"); // should be overridden by derived class constructors
    }

    QuadrivalentAtom2(
        const Compound::AtomName& atomName, ///< name for new atom
        const Element& element, ///< chemical element for new atom
        Angle bond12Angle, ///< bond angle in radians between bond center 1 and bond center 2
        Angle bond13Angle, ///< bond angle in radians between bond center 1 and bond center 3
        Angle bond14Angle, ///< bond angle in radians between bond center 1 and bond center 4
        Angle dihedral3, ///< dihedral angle of bond center 3, relative to dihedral of bond center 2
        Angle dihedral4 ///< dihedral angle of bond center 4, relative to dihedral of bond center 2
        )
        : Compound::SingleAtom(atomName, element)
    {
        // BondCenter1 dihedral will be relative to BondCenter2
        // Using Rotation() constructor that takes x axis and approximate y axis
        addFirstBondCenter( "bond1", atomName );
        // bond centers 2, 3, and 4 dihedrals will be relative to bond1
        addSecondBondCenter( "bond2", atomName, bond12Angle );

        // Compute bond23Angle and bond24Angle
        SimTK::Vec3 v1(0,0,1);
        SimTK::Vec3 v2 = SimTK::Rotation(bond12Angle, ZAxis) * SimTK::Vec3(0,0,1);
        SimTK::Vec3 v3 = SimTK::Rotation(dihedral3, YAxis) * SimTK::Rotation(bond12Angle, ZAxis) * SimTK::Vec3(0,0,1);
        SimTK::Vec3 v4 = SimTK::Rotation(dihedral4, YAxis) * SimTK::Rotation(bond12Angle, ZAxis) * SimTK::Vec3(0,0,1);

        Angle bond23Angle = std::acos(SimTK::dot(v2, v3));
        Angle bond24Angle = std::acos(SimTK::dot(v2, v4));

        // restrict angle range to (-Pi,Pi)
        while (dihedral3 >   Pi) dihedral3 -= 2.0*Pi;
        while (dihedral3 <= -Pi) dihedral3 += 2.0*Pi;
        while (dihedral4 >   Pi) dihedral4 -= 2.0*Pi;
        while (dihedral4 <= -Pi) dihedral4 += 2.0*Pi;

        if (dihedral3 > 0)
                addLeftHandedBondCenter( "bond3", atomName, bond13Angle, bond23Angle );
        else
                addRightHandedBondCenter( "bond3", atomName, bond13Angle, bond23Angle );

        if (dihedral4 > 0)
                addLeftHandedBondCenter( "bond4", atomName, bond14Angle, bond24Angle );
        else
                addRightHandedBondCenter( "bond4", atomName, bond14Angle, bond24Angle );

        // Choice of inboard bond may differ from bond priority - user may change this
        setInboardBondCenter("bond1");

        // default dihedral should be left blank, so it could be
        // defined from the other side.
        // the ultimate default will be 180 degrees anyway.
        // setDefaultInboardDihedralAngle(180*Deg2Rad);

        setCompoundName("QuadrivalentAtom"); // should be overridden by derived class constructors
    }
};

class  AliphaticHydrogen2 : public UnivalentAtom {
public:
    explicit AliphaticHydrogen2(const AtomName& atomName = "H" ///< name for new atom, defaults to "H"
        )
        : UnivalentAtom(atomName, Element::Hydrogen())
    {
        setDefaultInboardBondLength(0.1112); // for bonding to aliphatic carbon

        setCompoundName("AliphaticHydrogenAtom");
    }
};


/**
 * \brief AliphaticCarbon is a tetrahedral sp3 carbon atom for bonding to four other things
 *
 * The default geometry is perfectly terahedral, which is not bad for many uses.
 */
class  AliphaticCarbon2 : public QuadrivalentAtom2 {
public:
    explicit AliphaticCarbon2(const AtomName& atomName = "C" ///< name for new atom, defaults to "C"
        )
        : QuadrivalentAtom2(
              atomName, ///< name for new atom
              Element::Carbon(), ///< chemical element for new atom
              Angle(Deg2Rad * 90.0), // Angle bond12Angle, ///< bond angle in radians between bond center 1 and bond center 2
              Angle(Deg2Rad * 90.0), ///Angle bond13Angle < bond angle in radians between bond center 1 and bond center 3
              Angle(Deg2Rad * 90.0), //Angle bond14Angle, ///< bond angle in radians between bond center 1 and bond center 4
              Angle(Deg2Rad * 90.0), // Angle dihedral3, ///< dihedral angle of bond center 3, relative to dihedral of bond center 2
              Angle(Deg2Rad * -90.0) // Angle dihedral4 ///< dihedral angle of bond center 4, relative to dihedral of bond center 2
          )
    {
        // In case this bonds to another aliphatic carbon
        setDefaultInboardBondLength(0.15620); // for bonding to another aliphatic carbon

        setCompoundName("AliphaticCarbon");
    }
};


class  MethylGroup2 : public Compound {
public:
    MethylGroup2() {
        setBaseAtom(AliphaticCarbon2("C"));
        bondAtom(AliphaticHydrogen2("H2"), "C/bond2");
        bondAtom(AliphaticHydrogen2("H3"), "C/bond3");
        bondAtom(AliphaticHydrogen2("H4"), "C/bond4");
        nameBondCenter("bond", "C/bond1");

        setDefaultInboardBondLength(0.15620); // for bonding to another aliphatic carbon
        setCompoundName("MethylGroup");
    }
};

class  Methane2 : public Molecule {
public:
    Methane2()
    {
        setBaseCompound("methyl", MethylGroup2());
        inheritAtomNames("methyl");

        // Ordinarily, methyl group bonds to aliphatic carbon,
        // and has a default bond length to match.
        // Here we turn off the methyl inboard bond, so the
        // AliphaticHydrogen will, as desired, dictate the bond length
        convertInboardBondCenterToOutboard();

        bondAtom(AliphaticHydrogen2("H1"), "methyl/bond", 0.1112);

        setBiotypeIndex( "C", Biotype::MethaneC().getIndex() );
        setBiotypeIndex( "H1", Biotype::MethaneH().getIndex() );
        setBiotypeIndex( "H2", Biotype::MethaneH().getIndex() );
        setBiotypeIndex( "H3", Biotype::MethaneH().getIndex() );
        setBiotypeIndex( "H4", Biotype::MethaneH().getIndex() );

        setCompoundName("Methane");


    }
};

int main() {
try {
    
    CompoundSystem system; // Extends MolecularMechanicsSystem (Molmodel)
    SimbodyMatterSubsystem matter(system);
    DecorationSubsystem decorations(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

    // Build some compound here
    Methane2 c1;
    forceField.loadTestMoleculeParameters();
    Vec3 v1(0.0, 0.0, 0.0);

    system.adoptCompound(c1, v1);


    system.modelCompounds(); 

    system.addEventReporter(new Visualizer::Reporter(system, 0.020));
    //system.addEventHandler(new VelocityRescalingThermostat(	   system,  293.15, 0.1));

    //State state = system.realizeTopology();

    for (SimTK::Compound::AtomIndex aIx(0); aIx < c1.getNumAtoms(); ++aIx){
        SimTK::MobilizedBodyIndex mbx = c1.getAtomMobilizedBodyIndex(aIx);
        const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

        //DecorativeText dt(std::to_string(int(aIx)));
        DecorativeText dt(c1.getAtomName(aIx));
        dt.setScaleFactors(Vec3(0.04, 0.04, 0.04));
        dt.setColor(Vec3(1, 0., 0.));

        decorations.addBodyFixedDecoration(
            mbx,
            c1.getAtomLocationInMobilizedBodyFrame(aIx),
            dt
        );

    }
  
    State state = system.realizeTopology();

    VerletIntegrator integ(system);
    TimeStepper ts(system, integ);

    ts.initialize(state);
    ts.stepTo(0.06); // 0.06ps

    return 0;
} 
catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
catch(...) {
    std::cerr << "ERROR: An unknown exception was raised" 
              << std::endl;
    return 1;
}

}


