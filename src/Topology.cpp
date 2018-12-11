/**@file
Implementation of Topology class. **/

#include "Topology.hpp"

using namespace SimTK;

/** Default constructor.Sets the name of this molecule to 'no_name '.
The name has no particular function and is not guaranteed to be unique **/
Topology::Topology(){
    setName("no_name");
}

/** Constructor that sets the name of the molecule. The name has no particular 
function and is not guaranteed to be unique **/
Topology::Topology(std::string nameOfThisMolecule){
    setName(nameOfThisMolecule);
}

/** Default destructor. It deallocates bAtomType of every atom in the bAtomList
because we want to allow the valence to change during the simulation 
e.g. semi-grand canonical ensemble. **/
Topology::~Topology(){
    for(int i = 0; i < natoms; i++){
          delete bAtomList[i].bAtomType;
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
    for(int i=0; i<amberReader->getNumberAtoms(); i++){
        SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
              SimTK::Element(
                  bAtomList[i].getAtomicNumber()
                  , (std::to_string(bAtomList[i].getElem())).c_str()
                  , (std::to_string(bAtomList[i].getElem())).c_str()
                  , bAtomList[i].getMass())
            , bAtomList[i].getNBonds()
            , resName.c_str()
            , (bAtomList[i].getName()).c_str()
            , SimTK::Ordinality::Any
        );

        bAtomList[i].setBiotypeIndex(biotypeIndex);
        std::cout << " bAddParams: Defined Biotype: " 
            << resName << bAtomList[i].name << " " << SimTK::Ordinality::Any << "|" 
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
    // Add atom classes
    SimTK::DuMM::AtomClassIndex aIx;
    for(int i=0; i<amberReader->getNumberAtoms(); i++){
        aIx = dumm.getNextUnusedAtomClassIndex();
        bAtomList[i].setAtomClassIndex(aIx);
        dumm.defineAtomClass(
            (SimTK::DuMM::AtomClassIndex)aIx,
	        ( std::string("top") + resName + bAtomList[i].getFftype() 
                + std::string("_") + std::to_string(bAtomList[i].getNumber()) ).c_str(),
            bAtomList[i].getAtomicNumber(), // int atomicNumber
            bAtomList[i].getNBonds(),
            bAtomList[i].getVdwRadius() / 10.0, // nm
            //bAtomList[i].getVdwRadius(), // A
            bAtomList[i].getLJWellDepth() * 4.184 // kcal to kJ
        );
        std::cout << "bAddParams: Defined AtomClass: " << bAtomList[i].getAtomClassIndex() << std::endl;
    }

    // Define DuMM charged atom types
    for(int k=0; k<amberReader->getNumberAtoms(); k++){
        std::string abuff =  resName;
        abuff += bAtomList[k].biotype;
  
        SimTK::DuMM::ChargedAtomTypeIndex tempChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
        bAtomList[k].setChargedAtomTypeIndex(tempChargedAtomTypeIndex);
        std::cout << "bAddParams: defineChargedAtomType atomTypeIx "<<tempChargedAtomTypeIndex
            << " atomTypeName " << abuff.c_str() << " atomClassIx " << bAtomList[k].getAtomClassIndex()
            << " partialChargeInE " << bAtomList[k].charge << " chargedAtomTypeIndex " << bAtomList[k].getChargedAtomTypeIndex() << std::endl;
        dumm.defineChargedAtomType(
          tempChargedAtomTypeIndex,
          abuff.c_str(),
          bAtomList[k].getAtomClassIndex(),
          bAtomList[k].charge
        );
  
  
        // Associate a ChargedAtomTypeIndex with a Biotype index
        std::cout << "bAddParams: setBiotypeChargedAtomType bAtomList["<< k << "].getChargedAtomTypeIndex() "
            << bAtomList[k].getChargedAtomTypeIndex()
            << " bAtomList[" << k << "].biotype " << bAtomList[k].biotype
            << std::endl << std::flush;
        dumm.setBiotypeChargedAtomType(
          bAtomList[k].getChargedAtomTypeIndex(),
          bAtomList[k].getBiotypeIndex()
        );
    }

}

/** Calls DuMM defineBondStretch. **/
void Topology::bAddBondParams(
    std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Suppose or try to have the same order as the reader
    for(int t=0; t<amberReader->getNumberBonds(); t++){
        dumm.defineBondStretch_KA(
            (bAtomList[bonds[t].i]).getAtomClassIndex(), 
            (bAtomList[bonds[t].j]).getAtomClassIndex(),
            amberReader->getBondsForceK(t),  //k1
            amberReader->getBondsEqval(t)   //equil1
        );
    }
}

/** Calls DuMM defineBondBend. **/
void Topology::bAddAngleParams(
    std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    for(int t=0; t<amberReader->getNumberAngles(); t++){
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
    std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    bAddBiotypes(resName, amberReader, dumm); 
    bAddAtomClasses(resName, amberReader, dumm); 
    bAddBondParams(resName, amberReader, dumm); 
    bAddAngleParams(resName, amberReader, dumm); 
    bAddTorsionParams(resName, amberReader, dumm); 
}



/** Reads informatin from a readAmberInput object. **/
void Topology::loadAtomAndBondInfoFromReader(readAmberInput *amberReader)
{
    std::cout << "Topology::loadAtomAndBondInfoFromReader START" << std::endl;

    // Alloc memory for atoms and bonds list
    natoms = 0;
    int noDummies = 0; 
    natoms = amberReader->getNumberAtoms();
    nbonds = amberReader->getNumberBonds();

    //RE bAtomList = new bSpecificAtom[natoms];
    bAtomList.resize(natoms); //RE

    //RE bonds = new bBond[nbonds];
    bonds.resize(nbonds); //RE

    // Declare handy variables
    std::string str_buf;
    int a=65, b=65, c=65, d=65;
    int aRest=0, bRest=0, cRest=0;

    int charpos[] = {65, 65, 65, 65};
    std::string aStr, bStr, cStr, dStr;
    int nameCounter = 0; // same as i but clearer


    // -----------------------------------------------------------------
    // Assign atom properties. Take as much as possible from amberReader
    // -----------------------------------------------------------------
    // Make sure that variables in atoms are properly initialized
    // Shouldi not be needed
    for(int i = 0; i < natoms; i++){
        bAtomList[i].Zero();
    }

    for(int i = 0; i < natoms; i++){

        // Assign an index like in prmtop
        bAtomList[i].setNumber(i);
        //bAtomList[i].setNumber(amberReader->getAtomIndex); // not implemented

        // Assign element
        str_buf = amberReader->getAtomsName(i);
        unsigned int strix;
        for (strix = 0; strix < str_buf.length(); strix++){
            if(str_buf.at(strix) != ' '){
                break;
            }
        }
        bAtomList[i].setElem(str_buf.at(strix));

        // Assign a unique name specific to Gmolmodel. There are 60 available
        // ASCII readble characters: 0-9, A-Z and a-z. This gives a 12.960.000
        // of possible 4 character combinations in a number of the form 
        // a*60^3 + b*60^2 + c*60^1 + d. However the readble characters do not
        // form a continuous interval in the ASCII table so they have to be 
        // spread.
        std::string string_name;
        nameCounter++;

        /*
        charpos[0] = int(nameCounter / 216000);
        aRest = nameCounter % int(216000);
        charpos[1] = int(aRest / 3600);
        bRest = aRest % int(3600);
        charpos[2] = int(bRest / 60);
        cRest = bRest % int(60);
        charpos[3] = int(cRest / 1);
        // Spread the numbers so that the ASCII symbol will be readble
        for (unsigned int charposIx = 0; charposIx < 4; charposIx++){
            if(charpos[charposIx] > (10 + 25)){
                charpos[charposIx] += (7 + 6);
            }else if(charpos[charposIx] > (9)){
                charpos[charposIx] += (7);
            }
            charpos[charposIx] += 48;
        }
        aStr = (char)(charpos[0]);
        bStr = (char)(charpos[1]);
        cStr = (char)(charpos[2]);
        dStr = (char)(charpos[3]);
        */

        a = int(nameCounter / std::pow(25, 3));
        aStr = (char)(a + 65);
        aRest = nameCounter % int(std::pow(25, 3));

        b = int(aRest / std::pow(25, 2));
        bStr = (char)(b + 65);
        bRest = aRest % int(std::pow(25, 2));

        c = int(bRest / std::pow(25, 1));
        cStr = (char)(c + 65);
        cRest = bRest % int(std::pow(25, 1));

        d = int(cRest / std::pow(25, 0));
        dStr = (char)(d + 65);

        string_name = aStr + bStr + cStr + dStr;
        


        /*
        assert(nameCounter < 65536);
        std::stringstream sstream;
        std::string string_name;
        if(nameCounter < 16){
            sstream << std::hex << nameCounter;
            string_name = "000" + sstream.str();
        }else if(nameCounter < 256){
            sstream << std::hex << nameCounter;
            string_name = "00" + sstream.str();
        }else if(nameCounter < 4096){
            sstream << std::hex << nameCounter;
            string_name = "0" + sstream.str();
        }else if(nameCounter < 65536){
            sstream << std::hex << nameCounter;
            string_name = sstream.str();
        }
        */
        bAtomList[i].setName(string_name);

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


    // -----------------------
    // Assign bond properties
    // -----------------------
    //  
    for(int i=0; i<nbonds; i++){
        bonds[i].setIndex(i);
        bonds[i].i = amberReader->getBondsAtomsIndex1(i);
        bonds[i].j = amberReader->getBondsAtomsIndex2(i);
    }

    // Assign the number of bonds an atom has and sets the number of freebonds
    // equal to the number of bonds for now
    for(int i=0; i<natoms ; i++){
        bAtomList[i].nbonds = 0;
        for(int j=0; j<nbonds; j++){ 
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

    // ---------------------------------------------
    // Set every atom's (SimTK::Compound::SingleAtom *) to it's 
    // appropriate element and assign it's Compound::AtomName to unique name
    // Every atom is derived from SingleAtom in turn derived from 
    // Compound with one atom (AtomIndex 0)
    // TODO: Bromine and Clorine
    // ---------------------------------------------
    for(int i = 0; i < (natoms + noDummies); i++){
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
  
        bZeroCharArray(bAtomList[i].biotype, 20);
        sprintf(bAtomList[i].biotype, "%s_%s", \
          bAtomList[i].name, bAtomList[i].fftype);

    } // Finish assigning Compound::SingleAtom

    /* Just checking *////////
    std::cout<<"Checking after bMoleculeReader\n";
    for(int i=0; i<amberReader->getNumberAtoms();i++){
        bAtomList[i].Print();
    }
    for(int i=0; i<amberReader->getNumberBonds(); i++){ // EU
        std::cout<<"bond: "<<bonds[i].i<<" "<<bonds[i].j<<std::endl;
        fflush(stdout);
    }
    ///////////////////////////

    if (bAtomList.size() == 0){
      std::cout<<"bMoleculeReader constructor exit: 0 sized bAtomList\n";fflush(stdout);
    }
    else{
      std::cout<<"bMoleculeReader constructor: bAtomList loaded\n";fflush(stdout);
    }
    
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


/**
 * Main functions *
 **/

/** Process a graph node **/
void Topology::process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode)
{
    baseSetFlag = 0;
    std::cout << " switch to " << node->number << std::endl;

    if (node->visited == CurrentGeneration) {
        return;
    }

    ++nofProcesses;

    std::cout << " update generation " << std::endl;
    std::cout << " left bond " << previousNode->number << ' ' << node->number << std::endl;    

    // Mark Gmolmodel bond and create bond in Molmodel
    for(std::vector<bBond *>::iterator it = (node->bondsInvolved).begin();
    it != (node->bondsInvolved).end(); ++it){
        if ((*it)->isThisBond(node->number, previousNode->number) ){
            (*it)->setVisited(CurrentGeneration);

            // Create bond in Molmodel
            if(nofProcesses == 1){
                std::cout << "Skipped the first step" << std::endl;
            }else if(nofProcesses == 2){
                if( baseSetFlag == 0 ){
                    std::cout << "Set base atom" << std::endl;
                    this->setBaseAtom( *(previousNode->bAtomType) );
                    //this->setAtomBiotype(previousNode->name, (this->name), previousNode->biotype);
                    this->setAtomBiotype(previousNode->name, (this->name), previousNode->getName());
                    this->convertInboardBondCenterToOutboard();
                    baseSetFlag = 1;
                }

                std::stringstream sbuff;
                if(previousNode->number == baseAtomNumber){
                    sbuff << previousNode->name << "/bond" << previousNode->freebonds;
                }else{
                    sbuff << previousNode->name << "/bond" << (previousNode->nbonds - previousNode->freebonds + 1);
                }

                std::cout << "Trying to connect " << node->name << "(" << node->getInName() << ") " 
                    << node->number << " to " << previousNode->number << "(" << previousNode->getInName() << ") "
                    << (sbuff.str()).c_str() << " ... " << std::flush;

                std::cout << "Check 1: " << std::endl << std::flush;
                std::cout << "Check 1: " << " (sbuff.str()).c_str() " << (sbuff.str()).c_str() << std::endl << std::flush;
                std::cout << "Check 1: *(node->bAtomType) = " << *(node->bAtomType) << std::endl << std::flush;

                std::cout << "Node inName: " << *(node->inName) << std::endl << std::flush;

                this->bondAtom( *(node->bAtomType), (sbuff.str()).c_str(), 0.149, 0); // (Compound::SingleAtom&, BondCenterPathName, Length, Angle
                //this->setAtomBiotype(node->name, (this->name), node->biotype);

                std::cout << "Topology::process_node: setAtomBiotype: " 
                    << " node->name " << node->name << " (this->name).c_str() " << (this->name).c_str() << " node->getName() " << node->getName() 
                    << std::endl;

                this->setAtomBiotype(node->name, (this->name).c_str(), node->getName());

                // Set bSpecificAtom atomIndex to the last atom added to bond
                node->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ; // Set bSpecificAtom atomIndex to the last atom added to bond
                previousNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 0) ; // The only time we have to set atomIndex to the previous node

                // Set bBond Molmodel Compound::BondIndex
                (*it)->setBondIndex(Compound::BondIndex(getNumBonds() - 1));
                std::pair<SimTK::Compound::BondIndex, int > pairToBeInserted(Compound::BondIndex(getNumBonds() - 1), (*it)->getIndex());
                bondIx2bond.insert(pairToBeInserted);

                --previousNode->freebonds;
                --node->freebonds;

                std::cout << "done." << std::endl << std::flush;

            }else if(nofProcesses > 2){
                std::stringstream sbuff;
                if(previousNode->number == baseAtomNumber){
                    sbuff << previousNode->name << "/bond" << previousNode->freebonds;
                }else{
                    sbuff << previousNode->name << "/bond" << (previousNode->nbonds - previousNode->freebonds + 1);
                }


                std::cout << "Trying to connect " << node->name << "(" << node->getInName() << ") " 
                    << node->number << " to " << previousNode->number << "(" << previousNode->getInName() << ") "
                    << (sbuff.str()).c_str() << " ... " << std::flush;

                std::cout << "Check 2: " << std::endl << std::flush;
                std::cout << "Check 2: " << " (sbuff.str()).c_str() " << (sbuff.str()).c_str() << std::endl << std::flush;
                std::cout << "Check 2: *(node->bAtomType) = " << *(node->bAtomType) << std::endl << std::flush;

                this->bondAtom( *(node->bAtomType), (sbuff.str()).c_str(), 0.149, 0);
                //this->setAtomBiotype(node->name, (this->name), node->biotype);
                std::cout << "Topology::process_node: setAtomBiotype: " 
                    << " node->name " << node->name << " (this->name).c_str() " << (this->name).c_str() << " node->getName() " << node->getName() 
                    << std::endl;

                this->setAtomBiotype(node->name, (this->name), node->getName());

                // Set bSpecificAtom atomIndex to the last atom added to bond
                node->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ; // Set bSpecificAtom atomIndex to the last atom added to bond

                // Set bBond Molmodel Compound::BondIndex
                (*it)->setBondIndex(Compound::BondIndex(getNumBonds() - 1));
                std::pair<SimTK::Compound::BondIndex, int > pairToBeInserted(Compound::BondIndex(getNumBonds() - 1), (*it)->getIndex());
                bondIx2bond.insert(pairToBeInserted);

                --previousNode->freebonds;
                --node->freebonds;

                std::cout << "done." << std::endl << std::flush;

            }
            break;
        }
    }
    // ========
    previousNode = node;

    node->visited = CurrentGeneration;

    std::cout << " start checking neighbors " << std::endl;
    unsigned int i;
    for(i = 0; i < (node->neighbors).size(); i++) {
        //process_node( (node->neighbors)[i], CurrentGeneration, previousNode, nofProcesses, baseAtomNumber );
        process_node( (node->neighbors)[i], CurrentGeneration, previousNode);
    }

    // At the end of the graph walking
    // Compound can no longer be a child to the geometry of another compound
    if(node->number == 0){
        //this->convertInboardBondCenterToOutboard(); 
    }
    std::cout << " end processing " << node->number << std::endl;
}

// Construct the molecule topology
//void Topology::walkGraph(bSpecificAtom *root, int baseAtomNumber)
void Topology::walkGraph(bSpecificAtom *root)
{
    nofProcesses = 0;
    int CurrentGeneration = 0;

    CurrentGeneration += 1;
    bSpecificAtom *previousNode = root;

    baseSetFlag = 0;
    PrintStaticVars();
    process_node(root, CurrentGeneration, previousNode);
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
    this->regimenSpec = regimenSpec;

    this->setCompoundName((this->name));
    for(int i=0; i<natoms; i++){
        bAtomList[i].visited = 0;
    }
    for(int i=0; i<nbonds; i++){
        bonds[i].setVisited(0);
    }

    // Walk graph
    std::cout << "Walk the graph" << std::endl;

    int baseAtomListIndex = 0;
    for(int i=0; i<natoms; i++){
        std::cout << bAtomList[i].getNBonds() << std::endl;
        if(bAtomList[i].getNBonds() > 1){
            baseAtomListIndex = i;
            break;
        }
    }

    bSpecificAtom *root = &(bAtomList[baseAtomListIndex]);
    baseAtomNumber = root->number;
    walkGraph( &(bAtomList[baseAtomListIndex]));

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

// Get regimen  
std::string Topology::getRegimen(void){
    return this->regimen;
}

// Set regimen
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
                    setBondMobility(BondMobility::Free, Compound::BondIndex(j));
                    //setBondMobility(BondMobility::Torsion, Compound::BondIndex(j));
                    std::cout << "Topology::setRegimen: Bond " << j << " set to free" << std::endl;
                    break;
                }
            }
        }
        std::cout << "Changed regimen to: " << "RB" << std::endl;
    }
    this->regimen = argRegimen;
}

// Load MobilizedBody vs Compound::AtomIndex maps
void Topology::loadMaps(void){
    // Load MobilizedBodyIndex vs CompoundAtomIndex maps. 
    // These contain only atoms at the origin of mobods !!!
    for (SimTK::Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
        SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);
        //std::cout << "check transform = " << getAtomLocationInMobilizedBodyFrame(aIx) << std::endl;
        //if(getAtomLocationInMobilizedBodyFrame(aIx) == 0){
            std::pair<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > pairToBeInserted(mbx, aIx);
            mbx2aIx.insert(pairToBeInserted);
        //}
        aIx2mbx.insert(std::pair< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex >(aIx, mbx));
        //for(int i = 0; i< getNumAtoms(); i++){
        //    if((bAtomList[i]).atomIndex == aIx){
        //        aIx2number.insert(std::pair< SimTK::Compound::AtomIndex, int >(aIx, (bAtomList[i]).number));
        //        break;
        //    }
        //}
    }
}

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
        //std::cout << "Topology writePdb atom " << i << " " << bAtomList[i].getX() << std::endl;
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

  // Not sure we need this
void Topology::insertAtom(bSpecificAtom *)
{
    assert(!"Not implemented.");
}

void Topology::insertBond(int at1, int at2, int bondOrder)
{
    assert(!"Not implemented.");
}

// Get coordinates
void Topology::getCoordinates(std::vector<SimTK::Real> Xs, std::vector<SimTK::Real> Ys, std::vector<SimTK::Real> Zs)
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


