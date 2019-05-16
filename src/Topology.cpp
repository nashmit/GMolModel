/**@file
Implementation of Topology class. **/

#include "Topology.hpp"

using namespace SimTK;

/** Default constructor.Sets the name of this molecule to 'no_name '.
The name has no particular function and is not guaranteed to be unique **/
Topology::Topology(){
    this->name = std::string("no_name");
}

/** Constructor that sets the name of the molecule. The name has no particular 
function and is not guaranteed to be unique **/
Topology::Topology(std::string nameOfThisMolecule){
    this->name = nameOfThisMolecule;
}

/** Default destructor. It deallocates bAtomType of every atom in the bAtomList
because we want to allow the valence to change during the simulation 
e.g. semi-grand canonical ensemble. **/
Topology::~Topology(){
    for(int i = 0; i < bAtomList.size(); i++){
          delete bAtomList[i].bAtomType;
    }
}

/** Set atoms properties from a reader: number, name, element, initial
 * name, force field type, charge, coordinates, mass, LJ parameters **/
void Topology::SetAtomPropertiesFromReader(readAmberInput *amberReader)
{
    // Initialize atom variables
    for(int i = 0; i < natoms; i++){
        bAtomList[i].Zero();
    }

    // Declare handy variables
    std::string str_buf;

    // Iterate through atoms and set as much as possible from amberReader
    for(int i = 0; i < natoms; i++){

        // Assign an index like in prmtop
        bAtomList[i].setNumber(i);

        // Assign element from the first letter of the name
        str_buf = amberReader->getAtomsName(i);
        unsigned int strix;
        for (strix = 0; strix < str_buf.length(); strix++){
            if(str_buf.at(strix) != ' '){
                break;
            }
        }
        bAtomList[i].setElem(str_buf.at(strix));

        bAtomList[i].setName(GetUniqueName(i));

        // Store the initial name from prmtop
        bAtomList[i].setInName(str_buf);

        // Set atom type
        str_buf = amberReader->getAtomsNameAlias(i);
        boost::trim(str_buf);
        bAtomList[i].setFftype(str_buf);

        // Set charge as it is used in Amber
        SimTK::Real chargeMultiplier = 18.2223;
        bAtomList[i].setCharge(amberReader->getAtomsCharge(i) / chargeMultiplier);

        // Set coordinates in nm
        bAtomList[i].setX(amberReader->getAtomsXcoord(i) / 10.0);
        bAtomList[i].setY(amberReader->getAtomsYcoord(i) / 10.0);
        bAtomList[i].setZ(amberReader->getAtomsZcoord(i) / 10.0);

        // Set mass
        bAtomList[i].setMass(amberReader->getAtomsMass(i));

        // Set Lennard-Jones parameters
        bAtomList[i].setVdwRadius(amberReader->getAtomsRVdW(i));
        bAtomList[i].setLJWellDepth(amberReader->getAtomsEpsilon(i));

        // Set residue name and index
        //bAtomList[i].setResidueName(amberReader->getResidueName(i));
        //bAtomList[i].setResidueIndex(amberReader->getResidueIndex(i));
        bAtomList[i].residueName = std::string("UNK");
        bAtomList[i].residueIndex = 1;

        //bAtomList[i].Print();
    } // END atom properties
}

/** Set bonds properties from reader: bond indeces, atom neighbours **/
void Topology::SetBondingPropertiesFromReader(readAmberInput *amberReader)
{

    // Iterate through bonds and get atom indeces
    for(int i=0; i<nbonds; i++){
        bonds[i].setIndex(i);
        bonds[i].i = amberReader->getBondsAtomsIndex1(i);
        bonds[i].j = amberReader->getBondsAtomsIndex2(i);
    }

    // Assign the number of bonds an atom has and set the number of freebonds
    // equal to the number of bonds for now
    for(int i = 0; i < natoms ; i++){
        bAtomList[i].nbonds = 0;
        for(int j = 0; j < nbonds; j++){
            if((bAtomList[i].number == bonds[j].i) || \
               (bAtomList[i].number == bonds[j].j)){
                ++bAtomList[i].nbonds;
                ++bAtomList[i].freebonds;
            }
        }
    }

    // Assign neighbors and bonds involved for each atom
    // which translates into pushing bSpecificAtom * and bBond *
    // into their apropriate vectors
    for(int i=0; i<nbonds; i++){
        (bAtomList[ bonds[i].i  ]).addNeighbor( &(bAtomList[ bonds[i].j  ]) );
        (bAtomList[ bonds[i].i  ]).addBond( &(bonds[i]) );

        (bAtomList[ bonds[i].j  ]).addNeighbor( &(bAtomList[ bonds[i].i  ]) );
        (bAtomList[ bonds[i].j  ]).addBond( &(bonds[i]) );
    }
}

/** Set atoms Molmodel types (Compound::SingleAtom derived) based on
 * their valence **/
void Topology::SetAtomsMolmodelTypes()
{
    // ---------------------------------------------
    // Set every atom's (SimTK::Compound::SingleAtom *) to it's
    // appropriate element and assign it's Compound::AtomName to unique name
    // Every atom is derived from SingleAtom in turn derived from
    // Compound with one atom (AtomIndex 0)
    // Also set atom forecfield type
    // TODO: Bromine and Clorine and others
    // ---------------------------------------------
    for(int i = 0; i < (natoms); i++){
        // Atoms with one bond
        if(bAtomList[i].nbonds == 1){
            if(toupper(bAtomList[i].elem) == 'H'){
                bAtomList[i].bAtomType = new UnivalentAtom(bAtomList[i].name,
                                                           SimTK::Element( 1, "Hydrogen", "H", bAtomList[i].getMass() ));
                bAtomList[i].setAtomicNumber(1);
            }
                /*else if((toupper(bAtomList[i].name[0]) == 'C') && (toupper(bAtomList[i].name[0]) == 'L')){
                    bAtomList[i].bAtomType = new
                    UnivalentAtom(bAtomList[i].name, Element(17, "Chlorine", "Cl", bAtomList[i].getMass()));
                    bAtomList[i].setAtomicNumber(17);
                }*/
            else if(toupper(bAtomList[i].elem) == 'O'){
                bAtomList[i].bAtomType = new UnivalentAtom(bAtomList[i].name,
                                                           Element(8, "Oxygen", "O", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(8);
            }
            else if(toupper(bAtomList[i].elem) == 'F'){
                bAtomList[i].bAtomType = new
                        UnivalentAtom(bAtomList[i].name, Element(9, "Fluorine", "F", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(9);
            }
                /*
                else if((toupper(bAtomList[i].name[0]) == 'B') && (toupper(bAtomList[i].name[0]) == 'R')){
                  bAtomList[i].bAtomType = new
                    UnivalentAtom(bAtomList[i].name, Element(35, "Bromine", "Br", bAtomList[i].getMass()));
                  bAtomList[i].setAtomicNumber(35);
                }
                */
            else if(toupper(bAtomList[i].elem) == 'I'){
                bAtomList[i].bAtomType = new
                        UnivalentAtom(bAtomList[i].name, Element(53, "Iodine", "I", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(53);
            }
            else if(toupper(bAtomList[i].elem) == 'N'){
                bAtomList[i].bAtomType = new
                        UnivalentAtom(bAtomList[i].name, Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(7);
            }
            bAtomList[i].bAtomType->setDefaultInboardBondLength(0.1112); // Just for initial construction
        }
            // Atoms with two bonds
        else if (bAtomList[i].nbonds == 2){
            if(toupper(bAtomList[i].elem) == 'H'){
                bAtomList[i].bAtomType = new
                        BivalentAtom(bAtomList[i].name, Element(1, "Hydrogen", "H", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(1);
            }
            else if(toupper(bAtomList[i].elem) == 'C'){
                bAtomList[i].bAtomType = new
                        BivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(6);
            }
            else if(toupper(bAtomList[i].elem) == 'O'){
                bAtomList[i].bAtomType = new
                        BivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", bAtomList[i].getMass()),
                                     109.47*Deg2Rad);
                bAtomList[i].setAtomicNumber(8);
            }
            else if(toupper(bAtomList[i].elem) == 'N'){
                bAtomList[i].bAtomType = new
                        BivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(7);
            }
            else if(toupper(bAtomList[i].elem) == 'S'){
                bAtomList[i].bAtomType = new
                        BivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", bAtomList[i].getMass()),
                                     109.47*Deg2Rad);
                bAtomList[i].setAtomicNumber(16);
            }
            bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
        }
            // Atoms with three bonds
        else if (bAtomList[i].nbonds == 3){
            if(toupper(bAtomList[i].elem) == 'C'){
                bAtomList[i].bAtomType = new
                        TrivalentAtom(bAtomList[i].name, Element(6, "Carbon", "C", bAtomList[i].getMass()),
                                      120*Deg2Rad, 120*Deg2Rad
                );
                bAtomList[i].setAtomicNumber(6);
            }
            else if(toupper(bAtomList[i].elem) == 'O'){
                bAtomList[i].bAtomType = new
                        TrivalentAtomTetra(bAtomList[i].name,  Element(8, "Oxygen", "O", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(8);
            }
            else if(toupper(bAtomList[i].elem) == 'N'){
                bAtomList[i].bAtomType = new
                        TrivalentAtomTetra(bAtomList[i].name,  Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(7);
            }
            else if(toupper(bAtomList[i].elem) == 'S'){
                bAtomList[i].bAtomType = new
                        TrivalentAtomTetra(bAtomList[i].name,  Element(16, "Sulfur", "S", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(16);
            }
            else if(toupper(bAtomList[i].elem) == 'P'){
                bAtomList[i].bAtomType = new
                        TrivalentAtomTetra(bAtomList[i].name,  Element(15, "Phosphorus", "P", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(15);
            }
            bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
        }
            // Atoms with four bonds
        else if (bAtomList[i].nbonds == 4){
            if(toupper(bAtomList[i].elem) == 'C'){
                bAtomList[i].bAtomType = new
                        QuadrivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(6);
            }
            else if(toupper(bAtomList[i].elem) == 'O'){
                bAtomList[i].bAtomType = new
                        QuadrivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(8);
            }
            else if(toupper(bAtomList[i].elem) == 'N'){
                bAtomList[i].bAtomType = new
                        QuadrivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(7);
            }
            else if(toupper(bAtomList[i].elem) == 'S'){
                bAtomList[i].bAtomType = new
                        QuadrivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(16);
            }
            else if(toupper(bAtomList[i].elem) == 'P'){
                bAtomList[i].bAtomType = new
                        QuadrivalentAtom(bAtomList[i].name,  Element(15, "Phosphorus", "P", bAtomList[i].getMass()));
                bAtomList[i].setAtomicNumber(15);
            }
            bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
        }


    } // Finish assigning Compound::SingleAtoms
}

/** Reads information from a readAmberInput object and put it in
 * bAtomList and bonds lists **/
void Topology::loadAtomAndBondInfoFromReader(readAmberInput *amberReader)
{
    // Alloc memory for atoms and bonds list
    natoms = amberReader->getNumberAtoms();
    nbonds = amberReader->getNumberBonds();
    bAtomList.resize(natoms);
    bonds.resize(nbonds);

    assert( (!bAtomList.empty()) && "Topology::loadAtomAndBondInfoFromReader: atom list empty.");

    // Set atoms properties from a reader: number, name, element, initial
    // name, force field type, charge, coordinates, mass, LJ parameters
    SetAtomPropertiesFromReader(amberReader);

    // Set bonds properties from reader: bond indeces, atom neighbours
    SetBondingPropertiesFromReader(amberReader);

    // Set atoms Molmodel types (Compound::SingleAtom derived) based on
    // their valence
    SetAtomsMolmodelTypes();
}

/** Print atom list **/
void Topology::PrintAtomList()
{
    std::cout<<"Topology::PrintAtomList\n";
    for(unsigned int i = 0; i < bAtomList.size(); i++){
        bAtomList[i].Print();
    }
    for(unsigned int i = 0; i < bAtomList.size(); i++){
        bonds[i].Print();
    }
}

/** Biotype is a Molmodel hook that is usually used to look up molecular
force field specific parameters for an atom type. Gmolmodel defines a
new Biotype for each atom. The only thing that is specified is the element
with info about name, atomic number, valence and mass. **/
void Topology::bAddBiotypes(
      std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Iterate through atoms and define Biotypes based on resname
    for(int i = 0; i < amberReader->getNumberAtoms(); i++){
        SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
              SimTK::Element(
                  bAtomList[i].getAtomicNumber(),
                  (std::to_string(bAtomList[i].getElem())).c_str(),
                  (std::to_string(bAtomList[i].getElem())).c_str(),
                  bAtomList[i].getMass()),
            bAtomList[i].getNBonds(),
            resName.c_str(),
            (bAtomList[i].getName()).c_str(),
            SimTK::Ordinality::Any
        );

        bAtomList[i].setBiotypeIndex(biotypeIndex);

        // Assign atom's forcefield type
        bZeroCharArray(bAtomList[i].biotype, 20);

        //TODO: Assign Biotype function
        sprintf(bAtomList[i].biotype, "%s_%s", \
          bAtomList[i].name, bAtomList[i].fftype);

        std::cout << " bAddBiotypes: Defined Biotype: "
            << resName << "|" << bAtomList[i].name << " "
            << SimTK::Ordinality::Any << "|"
            << " with BiotypeIndex " << bAtomList[i].getBiotypeIndex() << std::endl;
    }
}

/** It calls DuMMs defineAtomClass, defineChargedAtomTye and 
setBiotypeChargedAtomType for every atom. These Molmodel functions contain 
information regarding the force field parameters. **/
void Topology::bAddAtomClasses(
                  std::string resName
                , readAmberInput *amberReader
                , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Define AtomClasses
    SimTK::DuMM::AtomClassIndex aIx;

    // Iterate through amberReader atoms and define AtomClasses
    for(int i = 0; i < amberReader->getNumberAtoms(); i++){
        // Get an AtomClass index
        aIx = dumm.getNextUnusedAtomClassIndex();
        bAtomList[i].setAtomClassIndex(aIx);

        // Define an AtomClass name
        const char* atomClassName = (
                std::string("top")
                + resName
                + bAtomList[i].getFftype()
                + std::string("_")
                + std::to_string(bAtomList[i].getNumber()) ).c_str();

        // Define an AtomClass (has info about van der Waals)
        dumm.defineAtomClass(
            aIx,
            atomClassName,
            bAtomList[i].getAtomicNumber(), // int atomicNumber
            bAtomList[i].getNBonds(), // expected valence
            bAtomList[i].getVdwRadius() / 10.0, // nm
            bAtomList[i].getLJWellDepth() * 4.184 // kcal to kJ
        );

        std::cout << "bAddAtomClasses: Defined AtomClass: "
            << bAtomList[i].getAtomClassIndex() << std::endl;
    }

    // Define ChargedAtomTypeIndeces
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
    std::string chargedAtomTypeName;

    // Iterate through atoms and define DuMM charged atom types
    for(int k = 0; k < amberReader->getNumberAtoms(); k++){
        // Get a ChargedAtomType index
        chargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
        bAtomList[k].setChargedAtomTypeIndex(chargedAtomTypeIndex);

        // Define a chargedAtomType name
        chargedAtomTypeName =  resName;
        chargedAtomTypeName += bAtomList[k].biotype;

        // Define a ChargedAtomType (AtomClass with a charge)
        dumm.defineChargedAtomType(
          chargedAtomTypeIndex,
          chargedAtomTypeName.c_str(),
          bAtomList[k].getAtomClassIndex(),
          bAtomList[k].charge
        );

        std::cout << "bAddAtomClasses: defineChargedAtomType atomTypeIx "<< chargedAtomTypeIndex
                  << " atomTypeName " << chargedAtomTypeName.c_str() << " atomClassIx " << bAtomList[k].getAtomClassIndex()
                  << " partialChargeInE " << bAtomList[k].charge << " chargedAtomTypeIndex " << bAtomList[k].getChargedAtomTypeIndex() << std::endl;
  
        // Associate a ChargedAtomTypeIndex with a Biotype index
        dumm.setBiotypeChargedAtomType(
          bAtomList[k].getChargedAtomTypeIndex(),
          bAtomList[k].getBiotypeIndex()
        );

        std::cout << "bAddAtomClasses: setBiotypeChargedAtomType bAtomList["<< k << "].getChargedAtomTypeIndex() "
                    << bAtomList[k].getChargedAtomTypeIndex()
                    << " bAtomList[" << k << "].biotype " << bAtomList[k].biotype
                    << std::endl << std::flush;
    }

}

/** Print Molmodel specific types as introduced in Gmolmodel **/
void Topology::PrintMolmodelAndDuMMTypes()
{

    for(int i = 0; i < bAtomList.size(); i++){
        std::cout << " name " << bAtomList[i].name << " BiotypeIndex "
                  << bAtomList[i].getBiotypeIndex() << std::endl;
    }

}

/** Calls DuMM defineBondStretch to define bonds parameters. **/
void Topology::bAddBondParams(
      std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Iterate through bonds and define their parameters
    // Suppose or try to have the same order as the reader
    for(int t = 0; t < amberReader->getNumberBonds(); t++){
        dumm.defineBondStretch_KA(
            (bAtomList[bonds[t].i]).getAtomClassIndex(), 
            (bAtomList[bonds[t].j]).getAtomClassIndex(),
            amberReader->getBondsForceK(t),  //k1
            amberReader->getBondsEqval(t)   //equil1
        );

    }
}

/** Calls DuMM defineBondBend to define angle parameters. **/
void Topology::bAddAngleParams(
      std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Iterate through angles and define their parameters
    for(int t = 0; t < amberReader->getNumberAngles(); t++){
        dumm.defineBondBend_KA(
            bAtomList[amberReader->getAnglesAtomsIndex1(t)].getAtomClassIndex(),
            bAtomList[amberReader->getAnglesAtomsIndex2(t)].getAtomClassIndex(),
            bAtomList[amberReader->getAnglesAtomsIndex3(t)].getAtomClassIndex(),
            amberReader->getAnglesForceK(t),
            ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getAnglesEqval(t))
        );
    }
}

/** Calls DuMM defineBondTorsion for 1, 2 and 3 periodicities **/
void Topology::bAddTorsionParams(
      std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    std::vector<std::pair<int, int>> pairStartAndLens = amberReader->getPairStartAndLen(); 

    for(unsigned int index=0; index<pairStartAndLens.size(); index++){

        int first    = pairStartAndLens[index].first;
        int numberOf = pairStartAndLens[index].second;

        for(int t = first; t < (first + numberOf); t++){
            if(numberOf == 1){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t),   
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))
                );
            }
            else if(numberOf == 2){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t), 
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)),
                    amberReader->getDihedralsPeriod(t+1), amberReader->getDihedralsForceK(t+1), 
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1))
                );
            }
            else if(numberOf == 3){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)),
                    amberReader->getDihedralsPeriod(t+1), amberReader->getDihedralsForceK(t+1),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)),
                    amberReader->getDihedralsPeriod(t+2), amberReader->getDihedralsForceK(t+2),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2))
                );
            }
        }
    }
}

/** Adds force field parameters read by the inputReader **/
void Topology::bAddAllParams(
    readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // We don't have any residues. The whole molecule is one residue
    std::string resName = this->name;

    // Add types
    bAddBiotypes(resName, amberReader, dumm); 
    bAddAtomClasses(resName, amberReader, dumm);

    // Add parameters
    bAddBondParams(resName, amberReader, dumm); 
    bAddAngleParams(resName, amberReader, dumm); 
    bAddTorsionParams(resName, amberReader, dumm); 
}




/**
 ** Interface **
 **/

/** Get the number of atoms in the molecule **/
int Topology::getNAtoms(void) const{
    return getNumAtoms();
}

/** Get the number of bonds in the molecule **/
int Topology::getNBonds(void) const{
    assert(!"Not implemented.");
}

/** Get a pointer to an atom object in the atom list inquiring
by number **/
bSpecificAtom * Topology::getAtomByNumber(int number) const{
    assert(!"Not implemented.");
}

/** Get a pointer to an atom object in the atom list inquiring
by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
bSpecificAtom * Topology::getAtomByAtomIx(int aIx) const{
    assert(!"Not implemented.");
}

/** Get a pointer to an atom object in the atom list inquiring
by atom name. **/
bSpecificAtom * Topology::getAtomByName(std::string name) const{assert(!"Not implemented.");}

/** Get the neighbours in the graph. **/
std::vector<bSpecificAtom *> Topology::getNeighbours(int) const{assert(!"Not implemented.");}

/** **/
bBond * Topology::getBond(int, int) const{assert(!"Not implemented.");}

/** Get bond order. **/
int Topology::getBondOrder(int, int) const{assert(!"Not implemented.");}


/** The following functions are used to build the molecular graph using bonding
information from bonds list and bondsInvolved list of each atom in bAtomList.
**/


/** The actual recursive function that builds the graph **/
void Topology::process_node(bSpecificAtom *node, bSpecificAtom *previousNode)
{
    // The base atom has to be set once Molmodel
    baseSetFlag = 0;

    // Only process unvisited nodes
    if( node->visited ){
        return;
    }

    // Mark the depth of the recursivity
    ++nofProcesses;

    // Mark Gmolmodel bond and create bond in Molmodel
    for(std::vector<bBond *>::iterator bondsInvolvedIter = (node->bondsInvolved).begin();
        bondsInvolvedIter != (node->bondsInvolved).end(); ++bondsInvolvedIter)
    {
        // Check if there is a bond between prevnode and node based on bonds
        // read from amberReader
        if ((*bondsInvolvedIter)->isThisMe(node->number, previousNode->number) ){
            (*bondsInvolvedIter)->setVisited(1);

            // Skip the first step as we don't have yet two atoms
            if(nofProcesses == 1){
                ;
            }

            // The first bond is special in Molmodel and has a different way of
            else if(nofProcesses == 2){

                // Set a base atom first
                if( baseSetFlag == 0 ){
                    this->setBaseAtom( *(previousNode->bAtomType) );
                    this->setAtomBiotype(previousNode->name, (this->name), previousNode->getName());
                    this->convertInboardBondCenterToOutboard();
                    baseSetFlag = 1;
                }

                // Bond current node by the previous (Compound function)
                std::stringstream parentBondCenterPathName;
                if(previousNode->number == baseAtomNumber){
                    parentBondCenterPathName << previousNode->name << "/bond" << previousNode->freebonds;
                }else{
                    parentBondCenterPathName << previousNode->name << "/bond" << (previousNode->nbonds - previousNode->freebonds + 1);
                }

                this->bondAtom( *(node->bAtomType), (parentBondCenterPathName.str()).c_str(), 0.149, 0); // (Compound::SingleAtom&, BondCenterPathName, Length, Angle

                // Set the final Biotype
                this->setAtomBiotype(node->name, (this->name).c_str(), node->getName());

                // Set bSpecificAtom atomIndex to the last atom added to bond
                node->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ;
                // The only time we have to set atomIndex to the previous node
                previousNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 0) ;

                // Set bBond Molmodel Compound::BondIndex
                (*bondsInvolvedIter)->setBondIndex(Compound::BondIndex(getNumBonds() - 1));
                std::pair<SimTK::Compound::BondIndex, int > pairToBeInserted(Compound::BondIndex(getNumBonds() - 1), (*bondsInvolvedIter)->getIndex());
                bondIx2bond.insert(pairToBeInserted);

                // Drop the number of available bonds
                --previousNode->freebonds;
                --node->freebonds;

            }
            // The rest of the bonds are not special
            else if(nofProcesses > 2){

                // Bond current node by the previous (Compound function)
                std::stringstream parentBondCenterPathName;
                if(previousNode->number == baseAtomNumber){
                    parentBondCenterPathName << previousNode->name << "/bond" << previousNode->freebonds;
                }else{
                    parentBondCenterPathName << previousNode->name << "/bond" << (previousNode->nbonds - previousNode->freebonds + 1);
                }

                this->bondAtom( *(node->bAtomType), (parentBondCenterPathName.str()).c_str(), 0.149, 0);

                // Set the final Biotype
                this->setAtomBiotype(node->name, (this->name), node->getName());

                // Set bSpecificAtom atomIndex to the last atom added to bond
                node->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ;

                // Set bBond Molmodel Compound::BondIndex
                (*bondsInvolvedIter)->setBondIndex(Compound::BondIndex(getNumBonds() - 1));
                std::pair<SimTK::Compound::BondIndex, int > pairToBeInserted(Compound::BondIndex(getNumBonds() - 1), (*bondsInvolvedIter)->getIndex());
                bondIx2bond.insert(pairToBeInserted);

                // Drop the number of available bonds
                --previousNode->freebonds;
                --node->freebonds;

            }

            break;

        }
    }

    // Mark the node as visited
    node->visited = 1;

    // Set the previous node to this node
    previousNode = node;

    // Choose the following node from his neighbours
    for(unsigned int i = 0; i < (node->neighbors).size(); i++) {
        process_node( (node->neighbors)[i], previousNode);
    }

}

/** Construct the molecule topology **/
void Topology::buildGraph(bSpecificAtom *root)
{
    nofProcesses = 0;
    int CurrentGeneration = 0;

    CurrentGeneration += 1;
    bSpecificAtom *previousNode = root;

    baseSetFlag = 0;
    PrintStaticVars();
    process_node(root, previousNode); // NEW
    std::cout << std::endl;
}

/** Builds the molecular tree, closes the rings, matches the configuration 
on the graph using using Molmodels matchDefaultConfiguration and sets the 
general flexibility of the molecule. **/
void Topology::build(
      SimTK::DuMMForceFieldSubsystem &dumm
    , std::string flexFN
    , std::string regimenSpec
)
{
    // Set regimen
    this->regimenSpec = regimenSpec;

    // Set the name of the Compound
    this->setCompoundName((this->name));

    // Initialize atoms to unvisited
    for(int i = 0; i < natoms; i++){
        bAtomList[i].setVisited(0);
    }

    // Initialize bonds to unvisited
    for(int i = 0; i < nbonds; i++){
        bonds[i].setVisited(0);
    }

    // Find an atom to be the root. It has to have more than one bond
    std::cout << "Start building the graph" << std::endl;
    int baseAtomListIndex = 0;
    for(int i = 0; i < natoms; i++){
        if(bAtomList[i].getNBonds() > 1){
            baseAtomListIndex = i;
            break;
        }
    }

    bSpecificAtom *root = &(bAtomList[baseAtomListIndex]);
    baseAtomNumber = root->number;

    // Build the graph
    buildGraph(root);

    // Add ring closing bonds
    for(int i=0; i<nbonds; i++){
        if(bonds[i].isVisited() == 0){
            
            bSpecificAtom *leftNode  =  &(bAtomList[bonds[i].i]);
            bSpecificAtom *rightNode =  &(bAtomList[bonds[i].j]);

            std::stringstream sbuff;
            if(leftNode->number == baseAtomNumber){
                sbuff << leftNode->name << "/bond" << leftNode->freebonds;
            }else{
                sbuff << leftNode->name << "/bond" << (leftNode->nbonds - leftNode->freebonds + 1);
            }
    
            std::stringstream otsbuff;
            if(rightNode->number == baseAtomNumber){
                otsbuff << rightNode->name << "/bond" << rightNode->freebonds;
            }else{
                otsbuff << rightNode->name << "/bond" << (rightNode->nbonds - rightNode->freebonds + 1);
            }
                  
            this->addRingClosingBond(
                (sbuff.str()).c_str(),
                (otsbuff.str()).c_str(),
                0.14,
                109*Deg2Rad,
                BondMobility::Rigid // CHANGE
            );
    
            //this->setAtomBiotype(leftNode->name, (this->name), leftNode->biotype);
            //this->setAtomBiotype(rightNode->name, (this->name), rightNode->biotype);
            this->setAtomBiotype(leftNode->name, (this->name), leftNode->getName());
            this->setAtomBiotype(rightNode->name, (this->name), rightNode->getName());
    
            //leftNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 0) ; // Not necessary
            //rightNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ; // Not necessary

            --leftNode->freebonds;
            --rightNode->freebonds;

            std::cout << "Closed ring " 
                << leftNode->name << "(" << leftNode->getInName() 
                << ") " << leftNode->number << " " << (sbuff.str()).c_str() << " to " 
                << rightNode->name << "(" << rightNode->getInName() 
                << ") " << rightNode->number << " " << (otsbuff.str()).c_str() 
                << " ... " << std::flush;

            std::cout << "done." << std::endl << std::flush;
        }
    }

    std::cout << "Topology: name inName atomIndexi:" << std::endl;
    for(int ix = 0; ix < getNumAtoms(); ++ix){
        std::cout << bAtomList[ix].name << " " << bAtomList[ix].inName << " " << bAtomList[ix].atomIndex << std::endl;
    }
    std::cout << std::endl << std::flush;

    // Assign Compound coordinates by matching bAtomList coordinates
    std::cout << "atomTargets from passed coords array: " << std::endl;
    std::map<AtomIndex, Vec3> atomTargets;
    for(int ix = 0; ix < getNumAtoms(); ++ix){
        std::cout << bAtomList[ix].atomIndex << " " 
            << this->getAtomName(bAtomList[ix].atomIndex) << " " 
            << bAtomList[ix].getX() << " " << bAtomList[ix].getY() << " " << bAtomList[ix].getZ()
            << std::endl;
        Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
        atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
    }

    std::cout << "Trying matchDefaultConfiguration Match_Exact ... " << std::flush;
    matchDefaultConfiguration(atomTargets, Match_Exact, true, 150.0); //Compound::Match_Idealized
    std::cout << "done. " << std::endl << std::flush;

    PdbStructure  pdb(*this);
    std::ostringstream sstream;
    sstream << "pdbs/sb_" << this->name <<"_ini"<<".pdb";
    std::string ofilename = sstream.str();
    std::filebuf fb;
    std::cout<<"Writing pdb file: "<<ofilename<<std::endl;
    fb.open(ofilename.c_str(), std::ios::out);
    std::ostream os(&fb);
    pdb.write(os); // automatically multiplies by ten (nm to A)
    fb.close();

    /* Just checking *////////
    std::cout << "Checking after Topology" << std::endl;
    for(int i=0; i<getNumAtoms();i++){
        bAtomList[i].Print();
    }
    for(int i=0; i<getNumBonds(); i++){ // EU
        std::cout<<"bond: "<<bonds[i].i<<" "<<bonds[i].j<<std::endl;
        fflush(stdout);
    }
    ///////////////////////////

    setRegimen(regimenSpec, flexFN);

}

/** Get regimen **/
std::string Topology::getRegimen(void){
    return this->regimen;
}

/** Set regimen **/
void Topology::setRegimen(std::string argRegimen, std::string flexFN){
    
    if(argRegimen == "IC"){
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Free, Compound::BondIndex(r));
            bonds[bondIx2bond[Compound::BondIndex(r)]].setBondMobility(BondMobility::Free);
        }
        std::cout << "Changed regimen to: " << "IC" << std::endl;
    }else if(argRegimen == "TD"){
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
        }
        std::cout << "Changed regimen to: " << "TD" << std::endl;
    }else if(argRegimen == "RB"){

        // Get flexible bonds from file and put it in PrmFlexBonds
        std::string line;
        std::ifstream F(flexFN);
        std::vector<std::pair<int, int>> PrmFlexBonds;
        std::vector<std::pair<int, int>>::iterator PrmFlexBondsIt;
        std::vector<std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>> MolmodelFlexBonds;
        std::vector<std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>>::iterator MolmodelFlexBondsIt;
        int line_i = -1;
        while(F.good()){
            line_i++;
            std::getline(F, line);
            std::istringstream iss(line);
            std::string word;
            std::vector<std::string> LineWords;
          
            int word_i = -1;
            while(iss >> word){
                if(word[0] == '#'){
                    break;
                }
                word_i++;
                LineWords.push_back(std::move(word));
            }
            if(word_i > 0){
                assert((word_i >= 1) && "2 indeces needed on each line of ligand.flex.");
                PrmFlexBonds.push_back(std::pair<int, int>( std::stod(LineWords[0]), std::stod(LineWords[1]) ));
            }
        }

        // Set all bonds to rigid first
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Rigid, SimTK::Compound::BondIndex(r));
        }

        // Iterate through prmtop flexible bonds, get Molmodel Compound::AtomIndeces and put it in MolmodelFlexBonds
        for ( PrmFlexBondsIt = PrmFlexBonds.begin(); PrmFlexBondsIt != PrmFlexBonds.end(); ++PrmFlexBondsIt){
            MolmodelFlexBonds.push_back(std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>(
                SimTK::Compound::AtomIndex(bAtomList[(*PrmFlexBondsIt).first].atomIndex), 
                SimTK::Compound::AtomIndex(bAtomList[(*PrmFlexBondsIt).second].atomIndex)) );
        }

        // Print
        for ( PrmFlexBondsIt = PrmFlexBonds.begin(); PrmFlexBondsIt != PrmFlexBonds.end(); ++PrmFlexBondsIt){
            std::cout << "Topology::setRegimen: FlexBond Prm Indeces " << (*PrmFlexBondsIt).first << " " << (*PrmFlexBondsIt).second << std::endl;
        }
        for ( MolmodelFlexBondsIt = MolmodelFlexBonds.begin(); MolmodelFlexBondsIt != MolmodelFlexBonds.end(); ++MolmodelFlexBondsIt){
            std::cout << "Topology::setRegimen: FlexBond AtomIndeces " << (*MolmodelFlexBondsIt).first << " " << (*MolmodelFlexBondsIt).second << std::endl;
        }
        for (unsigned int j = 0 ; j < getNumBonds(); j++){
            std::cout << "Topology::setRegimen: Compound Bonds AtomIndeces " << getBondAtomIndex(Compound::BondIndex(j), 0) << " " << getBondAtomIndex(Compound::BondIndex(j), 1) << std::endl;
        }

        // Iterate through  Molmodel Compound bonds and match with MolmodelFlexBonds
        for (unsigned int i = 0 ; i < MolmodelFlexBonds.size(); i++){
            for (unsigned int j = 0 ; j < getNumBonds(); j++){
                if( ((getBondAtomIndex(Compound::BondIndex(j), 0) == MolmodelFlexBonds[i].first) && (getBondAtomIndex(Compound::BondIndex(j), 1) == MolmodelFlexBonds[i].second)) || 
                    ((getBondAtomIndex(Compound::BondIndex(j), 1) == MolmodelFlexBonds[i].first) && (getBondAtomIndex(Compound::BondIndex(j), 0) == MolmodelFlexBonds[i].second)) ){
                    //setBondMobility(BondMobility::Free, Compound::BondIndex(j));
                    setBondMobility(BondMobility::Torsion, Compound::BondIndex(j));
                    std::cout << "Topology::setRegimen: Bond " << j << " set to torsion" << std::endl;
                    break;
                }
            }
        }
        std::cout << "Changed regimen to: " << "RB" << std::endl;
    }
    this->regimen = argRegimen;
}

/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
void Topology::loadMaps(void){

    // Iterate through atoms and get their MobilizedBodyIndeces
    for (SimTK::Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
        // Map mbx2aIx contains only atoms at the origin of mobods
        SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);
        std::pair<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > pairToBeInserted(mbx, aIx);
        mbx2aIx.insert(pairToBeInserted);

        // Map aIx is redundant in MobilizedBodyIndeces
        aIx2mbx.insert(std::pair< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex >(aIx, mbx));
    }

}


/** Print maps **/
void Topology::printMaps(void)
{
    std::cout << "Topology printMaps" << std::endl;
    std::cout << "mbx2aIx:" << std::endl;
    for(map<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex>::const_iterator it = mbx2aIx.begin();
       it != mbx2aIx.end(); ++it)
    {
        std::cout << "mbx " << it->first << " atomIndex " << it->second << std::endl;
    }
    std::cout << "aIx2mbx:" << std::endl;
    for(map<SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex>::const_iterator it = aIx2mbx.begin();
       it != aIx2mbx.end(); ++it)
    {
        std::cout << "atomIndex " << it->first << " mbx " << it->second << std::endl;
    }

    for(map<SimTK::Compound::BondIndex, int>::const_iterator it = bondIx2bond.begin();
       it != bondIx2bond.end(); ++it)
    {
        std::cout << "Compound bondIndex " << it->first << " bBond index " << it->second << std::endl;
    }
}

/** Write a pdb with bAtomList coordinates and inNames **/
void Topology::writePdb(std::string dirname, std::string prefix, std::string sufix, int maxNofDigits, int index) const
{
    int nofDigits = (int) floor(log10(index));
    std::string zeros("");
    if(maxNofDigits > nofDigits){
        for(int i = 0; i < (maxNofDigits - nofDigits); i++){
            zeros += std::string("0");
        }
    }
    std::stringstream sstream;
    sstream << dirname << "/" << prefix << zeros << std::to_string(index) << sufix;
    string ofilename = sstream.str();
    //std::cout << "Topology writePdb to " << ofilename << std::endl;

    FILE *oF = fopen (ofilename.c_str(),"w");
    // Pdb lines
    for(int i = 0; i < getNumAtoms(); i++){
        fprintf(oF, "%-6s%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f  %4.2f%6.2f          %2s\n"
            , "ATOM"                 // record
            , i                      // index
            , bAtomList[i].inName  // name
            , "UNK"                  // residue name
            , 'A'                    // chain
            , 1                      // residue index
            , 10.0*bAtomList[i].getX()    // x in A
            , 10.0*bAtomList[i].getY()    // y in A
            , 10.0*bAtomList[i].getZ()    // z in A
            , 1.0                    // occupancy
            , 0.0                    // beta factor
            , "  ");                 // element
    }

    fclose(oF);
}

/** Not sure we need this **/
void Topology::insertAtom(bSpecificAtom *)
{
    assert(!"Not implemented.");
}

/**  **/
void Topology::insertBond(int at1, int at2, int bondOrder)
{
    assert(!"Not implemented.");
}

/** Get bAtomList coordinates coordinates **/
void Topology::getCoordinates(
        std::vector<SimTK::Real> Xs,
        std::vector<SimTK::Real> Ys,
        std::vector<SimTK::Real> Zs)
{
    assert(Xs.size() == getNumAtoms());
    assert(Ys.size() == getNumAtoms());
    assert(Zs.size() == getNumAtoms());
    for(int ix = 0; ix < getNumAtoms(); ++ix){
        Xs[ix] = bAtomList[ix].getX();
        Ys[ix] = bAtomList[ix].getY();
        Zs[ix] = bAtomList[ix].getZ();
    }
}

// Not sure they belong here
/* ==================================================
 *    Scale all DuMM force field terms by scale_factor
 * ================================================== */
void Topology::setDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm, SimTK::Real scale_factor){    
    dumm.setBondStretchGlobalScaleFactor(scale_factor);
    dumm.setBondBendGlobalScaleFactor(scale_factor);
    dumm.setBondTorsionGlobalScaleFactor(scale_factor);
    dumm.setAmberImproperTorsionGlobalScaleFactor(scale_factor);
    dumm.setVdw12ScaleFactor(scale_factor);
    dumm.setVdw13ScaleFactor(scale_factor);
    dumm.setVdw14ScaleFactor(scale_factor);
    dumm.setVdw15ScaleFactor(scale_factor);
    dumm.setVdwGlobalScaleFactor(scale_factor);
    dumm.setCoulomb12ScaleFactor(scale_factor);
    dumm.setCoulomb13ScaleFactor(scale_factor);
    dumm.setCoulomb14ScaleFactor(scale_factor);
    dumm.setCoulomb15ScaleFactor(scale_factor);
    dumm.setCoulombGlobalScaleFactor(scale_factor);
    dumm.setGbsaGlobalScaleFactor(scale_factor);
}

/* ==================================================
 *    Scale DuMM force field terms by scale_factor
 * ================================================== */
void Topology::setSpecificDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm){

    dumm.setBondStretchGlobalScaleFactor(0.0);

    dumm.setBondBendGlobalScaleFactor(0.0);

    dumm.setBondTorsionGlobalScaleFactor(0.0);

    dumm.setAmberImproperTorsionGlobalScaleFactor(0.0);

    dumm.setVdw12ScaleFactor(0.0);
    dumm.setVdw13ScaleFactor(0.0);
    dumm.setVdw14ScaleFactor(0.5);
    dumm.setVdw15ScaleFactor(1.0);
    //dumm.setVdwGlobalScaleFactor(1.0);

    dumm.setCoulomb12ScaleFactor(0.0);
    dumm.setCoulomb13ScaleFactor(0.0);
    dumm.setCoulomb14ScaleFactor(0.0);
    dumm.setCoulomb15ScaleFactor(0.0);
    dumm.setCoulombGlobalScaleFactor(0.0);

    dumm.setGbsaGlobalScaleFactor(0.0);
}


