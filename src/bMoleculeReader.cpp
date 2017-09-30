#include "bMoleculeReader.hpp"

/*
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif

#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif
*/

//using namespace std;
using namespace SimTK;

/****
 *  Trivalent Atom Class with tetrahedral geometry.
 *  Bond centers are named "bond1", "bond2", and "bond3"
 ****/
TrivalentAtomTetra::TrivalentAtomTetra(
        const Compound::AtomName& atomName,   ///< name for new atom
        const Element& element   ///< chemical element for new atom
) : Compound::SingleAtom(atomName, element)
{
  static const Angle TetrahedralAngle = 109.47 * Deg2Rad;
  // BondCenter1 dihedral will be relative to BondCenter2
  addFirstBondCenter( "bond1", atomName );
  // bond centers 2 and 3 dihedrals relative to bond center 1
  //addSecondBondCenter( "bond2", atomName,  TetrahedralAngle);
  addLeftHandedBondCenter( "bond2", atomName, TetrahedralAngle, TetrahedralAngle );
  addRightHandedBondCenter( "bond3", atomName, TetrahedralAngle, TetrahedralAngle );
  // Choice of inboard bond may differ from bond priority - user may change this
  setInboardBondCenter("bond1");
  setCompoundName("TrivalentAtomTetra"); // overridden
}


/****
 * bPDBReader
 ****/
bPDBReader::bPDBReader()
{
    ;
}
bPDBReader::~bPDBReader()
{
    ;
}


/****
 * bSpecificAtom
 ****/
bSpecificAtom::bSpecificAtom(){
    nbonds = 0;
    freebonds = 0;
    number = 0;
    bZeroCharArray(name, 5);
    bZeroCharArray(inName, 5);
    bZeroCharArray(fftype, 20);
    bZeroCharArray(biotype, 20);
    x = -999;
    y = -999;
    z = -999;
    visited = 0;
}

bSpecificAtom::~bSpecificAtom(){;}

void bSpecificAtom::Zero(void){
    nbonds = 0;
    freebonds = 0;
    number = 0;
    bZeroCharArray(name, 5);
    bZeroCharArray(inName, 5);
    bZeroCharArray(fftype, 20);
    bZeroCharArray(biotype, 20);
    x = -999;
    y = -999;
    z = -999;
    visited = 0;
    bAtomType = NULL;
}

void bSpecificAtom::Print(void)
{
  std::cout<<"bSpecificAtom Print: nbonds "<<nbonds<<" freebonds "<<freebonds<<" name "<<name<<" inName "<<inName
    <<" number "<<number<<" atomIndex  "<<atomIndex<<" elem "<<elem<<" atomicNumber "<<atomicNumber<<" x "<<x<<" y "<< y<<" z "<<z
    <<" mass "<<mass<<" vdwRadius  "<<vdwRadius<<" LJWellDepth  "<<LJWellDepth<<" fftype "<<fftype
    <<" atomClassIndex  "<<atomClassIndex<<" biotype "<<biotype<<" bAtomType "<< bAtomType 
    <<" charge "<<charge<<" mobile "<<mobile<<" visited "<<visited<<std::endl;

    // Suggested by Eliza
    //std::string residueName;
    //long int residueIndex;
    //std::string chain;
    //int moleculeIndex;


}


// bSpecificAtom Interface

//  
int bSpecificAtom::getNBonds(void){
    return this->nbonds;
}

//
int bSpecificAtom::getFreebonds(void)
{
    assert(!"Not implemented");
}

//
std::string bSpecificAtom::getName(void){assert(!"Not implemented");}
std::string bSpecificAtom::getInName(void){assert(!"Not implemented");}

//
int bSpecificAtom::getNumber(void)
{
    return this->number;
}

// Return atom element
char bSpecificAtom::getElem(void)
{
    return this->elem;
}

SimTK::Real bSpecificAtom::getX(void)
{
    return this->x;
}

SimTK::Real bSpecificAtom::getY(void)
{
    return this->y;
}

SimTK::Real bSpecificAtom::getZ(void)
{
    return this->z;
}

//
std::string bSpecificAtom::getFftype(void)
{
    return this->fftype;
}

// Get atom biotype
const char * bSpecificAtom::getBiotype(void)
{
    return this->biotype;
}
SimTK::Compound::SingleAtom * bSpecificAtom::getBAtomType(void){assert(!"Not implemented");}

// Get atom's index as held in Compound
SimTK::Compound::AtomIndex bSpecificAtom::getCompoundAtomIndex(void)
{
    return this->atomIndex;
}

SimTK::Real bSpecificAtom::getCharge(void){assert(!"Not implemented");}
int bSpecificAtom::getIsMobile(void){assert(!"Not implemented");}
int bSpecificAtom::getIsVisited(void){assert(!"Not implemented");}
void bSpecificAtom::setNbonds(int){assert(!"Not implemented");}
void bSpecificAtom::setFreebonds(int){assert(!"Not implemented");}

// Set atom unique name
void bSpecificAtom::setName(std::string inpName){
    strncpy(this->name, inpName.c_str(), 5);
}

// Set initial name

void bSpecificAtom::setInName(std::string inpInName){
    strncpy(this->inName, inpInName.c_str(), 4);
}

// Set number
void bSpecificAtom::setNumber(int inpNumber){
    this->number = inpNumber;
}

// Set element
void bSpecificAtom::setElem(char inpElem){
    this->elem = inpElem;
}

// Set the X coordinate
void bSpecificAtom::setX(SimTK::Real inpX){
    this->x = inpX;
}

// Set the Y coordinate
void bSpecificAtom::setY(SimTK::Real inpY){
    this->y = inpY;
}

// Set the Z coordinate
void bSpecificAtom::setZ(SimTK::Real inpZ){
    this->z = inpZ;
}

// Get atomic mass
SimTK::mdunits::Mass bSpecificAtom::getMass(void)
{
    return this->mass;
}

// Set atomic mass
void bSpecificAtom::setMass(SimTK::Real inpMass)
{
    this->mass = inpMass;
}

// Set force field atom type
void bSpecificAtom::setFftype(std::string inpFftype){
    //strncpy(this->fftype, inpFftype.c_str(), 20);
    sprintf(this->fftype, "%s", inpFftype.c_str());
}

// Set atom Biotype name - dangerous
void bSpecificAtom::setBiotype(std::string inpBiotype)
{
    for(int i=0; i<20; i++){
        this->biotype[i] = '\0';
    }
    std::strncpy(this->biotype, inpBiotype.c_str(), 20);
}

// Set atom Biotype name - dangerous
void bSpecificAtom::setBiotype(const char * inpBiotype)
{
    for(int i=0; i<20; i++){
        this->biotype[i] = '\0';
    }
    std::strncpy(this->biotype, inpBiotype, 20);
}

// Set 
void bSpecificAtom::setBAtomType(SimTK::Compound::SingleAtom *){assert(!"Not implemented");}

void bSpecificAtom::setCompoundAtomIndex(SimTK::Compound::AtomIndex inpAtomIndex)
{
    this->atomIndex = inpAtomIndex;
}

// Set charge
void bSpecificAtom::setCharge(SimTK::Real inpCharge){
    this->charge = inpCharge;
}

// Get the DuMM ChargedAtomTypeIndex
SimTK::DuMM::ChargedAtomTypeIndex bSpecificAtom::getChargedAtomTypeIndex(void)
{
    return this->chargedAtomTypeIndex;
}

// Set the DuMM ChargedAtomTypeIndex
void bSpecificAtom::setChargedAtomTypeIndex(SimTK::DuMM::ChargedAtomTypeIndex inpChargedAtomTypeIndex)
{
    this->chargedAtomTypeIndex = inpChargedAtomTypeIndex;
}

void bSpecificAtom::setIsMobile(int){assert(!"Not implemented");}
void bSpecificAtom::setIsVisited(int){assert(!"Not implemented");}

// Get the atom class index
DuMM::AtomClassIndex bSpecificAtom::getAtomClassIndex(void)
{
    return this->atomClassIndex;
}

// Set the atom class index
void bSpecificAtom::setAtomClassIndex(DuMM::AtomClassIndex inpAtomClassIndex)
{
    this->atomClassIndex = inpAtomClassIndex;
}

// Get the atomic number
int bSpecificAtom::getAtomicNumber(void)
{
    return this->atomicNumber;
}

// Set the atomic number
void bSpecificAtom::setAtomicNumber(int inpAtomicNumber)
{
    this->atomicNumber = inpAtomicNumber;
}

//
void bSpecificAtom::setVdwRadius(SimTK::Real inpVdwRadius)
{
    this->vdwRadius = inpVdwRadius;
}

//
SimTK::Real bSpecificAtom::getVdwRadius(void)
{
    return this->vdwRadius;
}

//
void bSpecificAtom::setLJWellDepth(SimTK::Real inpLJWellDepth)
{
    this->LJWellDepth = inpLJWellDepth;
}

//
SimTK::Real bSpecificAtom::getLJWellDepth(void)
{
    return this->LJWellDepth;
}



/********************
 *     FUNCTIONS
 * ******************/
int bAtomAssign(MolAtom *dest, const bSpecificAtom *src)
{
  dest->name = new char[5];
  strncpy(dest->name, src->name, 4);

  dest->num = src->number;

  switch(src->elem){
    case 'C': dest->type = MOL_ATOM_ELEMENT_CARBON; break;
    case 'H': dest->type = MOL_ATOM_ELEMENT_HYDROGEN; break;
    case 'N': dest->type = MOL_ATOM_ELEMENT_NITROGEN; break;
    case 'O': dest->type = MOL_ATOM_ELEMENT_OXYGEN; break;
    case 'S': dest->type = MOL_ATOM_ELEMENT_SULFUR; break;
    case 'P': dest->type = MOL_ATOM_ELEMENT_PHOSPHORUS; break;
    default: dest->type = MOL_ATOM_ELEMENT_UNKNOWN; break;
  }

  dest->pos[0] = src->x;
  dest->pos[1] = src->y;
  dest->pos[2] = src->z;

  return 0;
}

/********************************/

/****
 * intpair
 ****/
intpair::intpair(){
    i = 0; j = 0;
  }
intpair::intpair(int inI, int inJ){
    this->i = inI;
    this->j = inJ;
}
intpair::~intpair(){}

bool intpair::operator==(const intpair *other){
  return (
    ((this->i == other->i) && (this->j == other->j)) ||
    ((this->i == other->j) && (this->j == other->i))
  );
}

bool intpair::operator!=(const intpair *other){
  return (
    ((this->i != other->i) || (this->j != other->j)) &&
    ((this->i != other->j) || (this->j != other->i))
  );
}

bool intpair::isTheSameAs(const intpair *other){
  return (
    ((this->i == other->i) && (this->j == other->j)) ||
    ((this->i == other->j) && (this->j == other->i))
  );
}

void intpair::swap(void){
  int inter;
  inter = this->i;
  this->i = this->j;
  this->j = inter;
}

void intpair::dump(void){
  std::cout<<i<<' '<<j<<std::endl;
}

std::string intpair::getString(void){
  std::stringstream ret;
  ret<<i<<' '<<j;
  return ret.str();
}

/********************************/


/****
 * bBond
 ****/
bBond::bBond(void) : intpair(){
  inring = 0;
  rigid = 0;
  ring_closing = 0; // later
  ring_no = 0; // later
}

bBond::bBond(int a, int b) : intpair(a, b){
  inring = 0;
  rigid = 0;
  ring_closing = 0;
  ring_no = 0; // later
}

bBond::~bBond(void){;}

bool bBond::isInRing(void){
  if(this->inring == 1)
    return true;
  else
    return false;
}

bool bBond::isRingClosing(void){
  if(this->ring_closing == 1)
    return true;
  else
    return false;
}

bool bBond::isRigid(void){
  if(this->rigid == 1)
    return true;
  else
    return false;
}

int bBond::ringNo(void){
  return this->ring_no;
}

void bBond::setInRing(void){
  this->inring = 1;
}

void bBond::setAsRingClosing(void){
  this->ring_closing = 1;
}

void bBond::setAsRigid(void){
  this->rigid = 1;
}

void bBond::setRingNo(int rn){
  this->ring_no = rn;
}

Compound::BondIndex bBond::getBondIndex(void){
  return bondIndex;
}
  
void bBond::setBondIndex(Compound::BondIndex otherIx){
  this->bondIndex = otherIx;
}



/********************************/

/****
 * intriad
 ****/
intriad::intriad(){
  i=0; j=0; k=0;
}
intriad::intriad(int inI, int inJ, int inK){
  this->i = inI;
  this->j = inJ;
  this->k = inK;
}
intriad::~intriad(){}
//Sorry for crowding
bool intriad::operator==(const intriad *other){
  return (
    ((this->i==other->i)&&(this->j==other->j)&&(this->k==other->k))||
    ((this->i==other->j)&&(this->j==other->i)&&(this->k==other->k))||
    ((this->i==other->k)&&(this->j==other->i)&&(this->k==other->j))||
    ((this->i==other->k)&&(this->j==other->j)&&(this->k==other->i))||
    ((this->i==other->i)&&(this->j==other->k)&&(this->k==other->j))||
    ((this->i==other->j)&&(this->j==other->k)&&(this->k==other->i))
  );
}
bool intriad::isTheSameAs(const intriad *other){
  return (
    ((this->i==other->i)&&(this->j==other->j)&&(this->k==other->k))||
    ((this->i==other->k)&&(this->j==other->j)&&(this->k==other->i))
  );
}
void intriad::dump(void){
  std::cout<<i<<' '<<j<<' '<<k<<std::endl;
}
std::string intriad::getString(void){
  std::stringstream ret;
  ret<<i<<' '<<j<<' '<<k;
  return ret.str();
}

/********************************/

/****
 * bMoleculeReader PRMTOP
 ****/

bMoleculeReader::bMoleculeReader(readAmberInput *amberReader, const char *rbfilename)
{
    natoms = 0;
    bBond buffij;
    char elem = 'x';
    int noDummies = 0; 
    std::string line;
    char line_c[1000];
    natoms = amberReader->getNumberAtoms();
    nbonds = amberReader->getNumberBonds();

    bAtomList = new bSpecificAtom[natoms]; /*Alloc - 1 for dummy*/
    bonds = new bBond[nbonds];
    // Read atoms READER

    std::string str_buf;
    SimTK::Real chargeMultiplier = 18.2223;
    for(int i = 0; i < natoms; i++){
        bAtomList[i].Zero();

        bAtomList[i].setNumber(i);
        str_buf = amberReader->getAtomsName(i);
        boost::trim(str_buf);
        bAtomList[i].setElem(str_buf.at(0));
        //bAtomList[i].setName(str_buf);
        std::string elemStr = std::string(1, bAtomList[i].getElem());
        bAtomList[i].setName(elemStr + std::to_string(bAtomList[i].getNumber()));
        bAtomList[i].setInName(str_buf);


        str_buf = amberReader->getAtomsNameAlias(i);
        boost::trim(str_buf);
        bAtomList[i].setFftype(str_buf);

        bAtomList[i].setCharge(amberReader->getAtomsCharge(i) / chargeMultiplier);

        bAtomList[i].setX(amberReader->getAtomsXcoord(i) / 10.0);
        bAtomList[i].setY(amberReader->getAtomsYcoord(i) / 10.0);
        bAtomList[i].setZ(amberReader->getAtomsZcoord(i) / 10.0);

        bAtomList[i].setMass(amberReader->getAtomsMass(i));

        bAtomList[i].setVdwRadius(amberReader->getAtomsRVdW(i));
        bAtomList[i].setLJWellDepth(amberReader->getAtomsEpsilon(i));

        bAtomList[i].Print();
    }

    /*READ BONDS*/
    for(int i=0; i<nbonds; i++){
        bonds[i].i = amberReader->getBondsAtomsIndex1(i);
        bonds[i].j = amberReader->getBondsAtomsIndex2(i);
    }

    /*Assign atoms nbonds and freebonds*/
    for(int i=0; i<natoms ; i++){
      bAtomList[i].nbonds = 0;
      for(int j=0; j<nbonds; j++){ 
        if((bAtomList[i].number == bonds[j].i) ||\
           (bAtomList[i].number == bonds[j].j)){
          ++bAtomList[i].nbonds;
          ++bAtomList[i].freebonds;
        }
      }
    }


    // Every *valentAtom is derived from SingleAtom in turn derived from Compound with one atom (AtomIndex 0)
    for(int i=0; i<natoms+noDummies; i++){
      if(bAtomList[i].nbonds == 1){
        if(toupper(bAtomList[i].elem) == 'H'){
          bAtomList[i].bAtomType = new
            UnivalentAtom(bAtomList[i].name, SimTK::Element( 1, "Hydrogen", "H", bAtomList[i].getMass() )); // Prmtop mass
          bAtomList[i].setAtomicNumber(1);
        }
        else if((toupper(bAtomList[i].name[0]) == 'C') && (toupper(bAtomList[i].name[0]) == 'L')){
          bAtomList[i].bAtomType = new
            UnivalentAtom(bAtomList[i].name, Element(17, "Chlorine", "Cl", bAtomList[i].getMass()));
          bAtomList[i].setAtomicNumber(17);
        }
        else if(toupper(bAtomList[i].name[0]) == 'O'){
          bAtomList[i].bAtomType = new
            UnivalentAtom(bAtomList[i].name, Element(8, "Oxygen", "O", bAtomList[i].getMass()));
          bAtomList[i].setAtomicNumber(8);
        }
        else if(toupper(bAtomList[i].name[0]) == 'F'){
          bAtomList[i].bAtomType = new
            UnivalentAtom(bAtomList[i].name, Element(9, "Fluorine", "F", bAtomList[i].getMass()));
          bAtomList[i].setAtomicNumber(9);
        }
        else if((toupper(bAtomList[i].name[0]) == 'B') && (toupper(bAtomList[i].name[0]) == 'R')){
          bAtomList[i].bAtomType = new
            UnivalentAtom(bAtomList[i].name, Element(35, "Bromine", "Br", bAtomList[i].getMass()));
          bAtomList[i].setAtomicNumber(35);
        }
        else if(toupper(bAtomList[i].name[0]) == 'I'){
          bAtomList[i].bAtomType = new
            UnivalentAtom(bAtomList[i].name, Element(53, "Iodine", "I", bAtomList[i].getMass()));
          bAtomList[i].setAtomicNumber(53);
        }
        else if(toupper(bAtomList[i].name[0]) == 'N'){
          bAtomList[i].bAtomType = new
            UnivalentAtom(bAtomList[i].name, Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
          bAtomList[i].setAtomicNumber(7);
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.1112); // Just for initial construction
      }
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

    }

    if (bAtomList == NULL){
      std::cout<<"bMoleculeReader constructor exit: NULL bAtomList\n";fflush(stdout);
    }
    else{
      std::cout<<"bMoleculeReader constructor: bAtomList not NULL\n";fflush(stdout);
    }
    


  /*Now read rigid bodies specifications*/
  FILE *rfpo;
  rfpo = fopen(rbfilename, "r");
  if(rfpo == NULL){
    printf("Usage:\n<program> -mol2 <mol2_file> -rb <rb_file> -gaff  gaff.dat -frcmod <frcmod_file>\n");
    printf("rb_file not provided. Exiting...\n");
    exit(1);
  }

  std::string sbuff;
  std::vector<int> ring;
  char curr_char;
  int bond_found = 0;
  unsigned int boi;
  int par_open = 0;
  int ring_no = -1;
  int ring_closing = 0;

  while(fgets(line_c, MAX_LINE_LENGTH, rfpo)){
    line = line_c; //RESTORE
    if((sbuff = line.substr(0,5)) == "rings"){ // RESTORE
      sbuff = "";
      for(int i=0; i<strlen(line_c); i++){ 
        curr_char = line_c[i]; // RESTORE
        if(curr_char == ']'){
          ++ring_no;
          ring.push_back(atoi(sbuff.c_str()));
          sbuff = "";
          bond_found = 0;
          for(boi=0; boi<nbonds; boi++){ // EU
            for(unsigned int ri=0; ri<ring.size(); ri++){
              for(unsigned int rj=0; (rj<ring.size()) && (rj<ri); rj++){
                if( ((ring[ri] == bonds[boi].i) && (ring[rj] == bonds[boi].j)) ||
                  ((ring[ri] == bonds[boi].j) && (ring[rj] == bonds[boi].i))){
                  bonds[boi].setInRing();
                  bonds[boi].setRingNo(ring_no);
                  if(ring_closing == 0){
                    bonds[boi].setAsRingClosing();
                  }
                  ++ring_closing;
                  bond_found = 1;
                  break;
                }
              }
            }
          }
          ring.clear();
          ring_closing = 0;
        }
        else if(curr_char == ','){
          ring.push_back(atoi(sbuff.c_str()));
          sbuff = "";
        }
        else if((curr_char >= '0') && (curr_char <='9')){
          sbuff += curr_char;
        }
      }
    }

    else if((sbuff = line.substr(0,17)) == "non_ring_pi_bonds"){ // RESTORE
      sbuff = "";
      for(int i=0; i<line.size(); i++){
        curr_char = line_c[i]; // RESTORE
        if(curr_char == ')'){
          par_open = 0;
          ring.push_back(atoi(sbuff.c_str()));
          for(boi = 0; boi<nbonds; boi++){ // EU
            if( (bonds[boi].i == ring[0]) && (bonds[boi].j == ring[1]) ||
              (bonds[boi].i == ring[1]) && (bonds[boi].j == ring[0]) ){
              bonds[boi].setAsRigid();
              ring.clear();
              sbuff = "";
              break;
            }
          }
        }
        else if((curr_char == ',') && (par_open == 1)){
          ring.push_back(atoi(sbuff.c_str()));
          sbuff = "";
        }
        else if((curr_char >= '0') && (curr_char <='9')){
          sbuff += curr_char;
        }
        else if(curr_char == '('){
          par_open = 1;
        }
      }
    }

    else if((sbuff = line.substr(0,17)) == "rigid_bodies"){
      // TODO
    }

  }

  fclose(rfpo);

  /* Just checking *////////
  std::cout<<"Checking after bMoleculeReader\n";
  for(int i=0; i<natoms;i++){
      bAtomList[i].Print();
  }
  for(int i=0; i<nbonds; i++){ // EU
    std::cout<<"bond: "<<bonds[i].i<<" "<<bonds[i].j<<std::endl;
    fflush(stdout);
  }
  ///////////////////////////


}

bMoleculeReader::bMoleculeReader(DuMMForceFieldSubsystem& dumm,
        const char *filename,
        const char *filetype,
        const char *rbfilename){
  #ifdef DEBUG_LEVEL02
  std::cout<<"bMoleculeReader::bMoleculeReader() BEGIN"<<std::endl<<std::flush;
  #endif

  FILE *fpo;
  fpo = fopen(filename, "r");
  /*rewind(fpo);*/ /*Doesn't work for Windows files*/
  MAX_LINE_LENGTH = 1000;
  char *buff = new char[80];
  //bZeroCharArray(buff, 80);
  bZeroCharArray(buff, 80);
  natoms = 0;
  //int nbonds = 0;
  bBond buffij;
  unsigned int i = 0;
  unsigned int j = 0;
  char line_c[MAX_LINE_LENGTH];
  std::string line;
  std::string str_buf;
  char elem = 'x';
  unsigned int lno = 0;
  int noDummies = 0;

  /*+++++++++ External readers ++++++++++*/
  if(strstr(filetype, "external")){
    ;
  }

  /*+++++++++ MOL2 type ++++++++++*/
  else if(strstr(filetype, "mol2")){
    bZeroCharArray(buff, 80); // EU
    natoms = 0;
    i = 0; j = 0;
    char elem = 'x';
    unsigned int lno = 0;
 
    /*Read number of atoms*/
    bZeroCharArray(line_c, MAX_LINE_LENGTH); // EU
    while(fgets(line_c, MAX_LINE_LENGTH, fpo)){
      ++lno;
      if(lno == 3){
        bSubstr(buff, line_c, 0, 5);
        natoms = atoi(buff);
        bZeroCharArray(buff, 80);

        bSubstr(buff, line_c, 6, 6); // EU
        nbonds = atoi(buff); // EU
        bonds = new bBond[nbonds];
        bZeroCharArray(buff, 80); // EU

        break;
      }
      bZeroCharArray(line_c, MAX_LINE_LENGTH);
    }
    bZeroCharArray(line_c, MAX_LINE_LENGTH);

    bAtomList = new bSpecificAtom[natoms + noDummies]; /*Alloc - 1 for dummy*/

    /*Jump to the atom section*/
    while(fgets(line_c, MAX_LINE_LENGTH, fpo)){
      if(line_c[0] == '@'){
        break;
      }
      bZeroCharArray(line_c, MAX_LINE_LENGTH);
    }
    bZeroCharArray(line_c, MAX_LINE_LENGTH);

    /*Read position, element, fftype and charge*/
    lno = 0;
    while(fgets(line_c, MAX_LINE_LENGTH, fpo) && (lno < natoms)){ 
      ++lno;
      if(line_c != NULL){
        std::cout<<"bMoleculeReader::bMoleculeReader() line_c "<<strlen(line_c)<<" chars"<<std::endl<<std::flush;
      }
      line = line_c; // INTERFACE
      bAtomList[lno-1].Zero();

      // Create new unique name for atom
      elem = line.at(8);
      bAtomList[lno-1].setElem(elem);

      str_buf = elem;
      str_buf += std::to_string(lno);
      bAtomList[lno-1].setName(str_buf);

      bAtomList[lno-1].setNumber(lno);

      str_buf = "gaff_";
      str_buf += line.substr(47, 2);
      //str_buf.erase(std::remove_if(str_buf.begin(), str_buf.end(), ::isspace), str_buf.end());
      boost::trim(str_buf);
      bAtomList[lno-1].setFftype(str_buf);

      str_buf = line.substr(67,9);
      bAtomList[lno-1].setCharge(std::stod(str_buf));

      str_buf = line.substr(8, 4);
      bAtomList[lno-1].setInName(str_buf);

      str_buf = line.substr(17,9);
      bAtomList[lno-1].setX(std::stod(str_buf));

      str_buf = line.substr(27,9);
      bAtomList[lno-1].setY(std::stod(str_buf));

      str_buf = line.substr(37,9);
      bAtomList[lno-1].setZ(std::stod(str_buf));

      //bAtomList[lno-1].Print(); // EU
    }

    #ifdef DEBUG_LEVEL02
    //std::cout<<std::endl<<"bMoleculeReader::bMoleculeReader 2 Read Bonds down"<<std::endl<<std::flush;
    #endif
    /*READ BONDS*/
    //fgets(line_c, MAX_LINE_LENGTH, fpo); // RESTORE
    lno = 0;
    //do{ // RESTORE
    std::cout<<"bMoleculeReader hop nbonds "<<nbonds<<std::endl<<std::flush;
    while(fgets(line_c, MAX_LINE_LENGTH, fpo) && (lno < nbonds)){ // EU 
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() bond line= "<<line_c<<std::flush;
      #endif
      ++lno;
      //line = line_c; // RESTORE
      sscanf(line_c, "%d%d%d", &i, &buffij.i, &buffij.j);
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() bond line sscanf "<<std::endl<<std::flush;
      std::cout<<"bMoleculeReader::bMoleculeReader() i, buffij.i buffij.j "
        <<i<<" "<<buffij.i<<" "<<buffij.j<<std::endl<<std::flush;
      #endif
      //bonds.push_back(buffij); // RESTORE
      bonds[i-1].i = buffij.i; // EU
      bonds[i-1].j = buffij.j; // EU
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() bonds["<<i-1<<"].i "<<bonds[i-1].i<<" .j  "<<bonds[i-1].j<<std::endl<<std::flush;
      #endif
    }
    //while(fgets(line_c, MAX_LINE_LENGTH, fpo) && (line_c[0] != '@')); // RESTORE

    
    /*Assign atoms nbonds and freebonds*/
    for(i=0; i<natoms;i++){
      bAtomList[i].nbonds = 0;
      //for(j=0; j<bonds.size(); j++){ // RESTORE
      for(j=0; j<nbonds; j++){ // EU
        if((bAtomList[i].number == bonds[j].i) ||\
           (bAtomList[i].number == bonds[j].j)){
          ++bAtomList[i].nbonds;
          ++bAtomList[i].freebonds;
        }
      }
    }

    #ifdef DEBUG_LEVEL02
    std::cout<<"bMoleculeReader::bMoleculeReader() nbonds & freebonds assigned "<<line_c<<std::endl<<std::fflush;
    #endif


    /*TODO Develop TrivalentAtomTetra for adding custom angles*/
    // Every *valentAtom is derived from SingleAtom in turn derived from Compound with one atom (AtomIndex 0)
    for(i=0; i<natoms+noDummies; i++){
      #ifdef DEBUG_LEVEL02
      //std::cout<<"bMoleculeReader::bMoleculeReader 3"<<std::endl<<std::flush;
      #endif
      if(bAtomList[i].nbonds == 1){
        if(toupper(bAtomList[i].elem) == 'H'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, MMTKElement(1, "MMTKHydrogen", "MMTKH", 1.00797598)); // MMTK mass
            //UnivalentAtom(bAtomList[i].name, MMTKElement::MMTKHydrogen()); // Molmodel Mass
            UnivalentAtom(bAtomList[i].name, Element::Hydrogen()); // Molmodel Mass
        }
        else if(toupper(bAtomList[i].name[0]) == 'C'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(17, "Chlorine", "Cl", 35.45));
            UnivalentAtom(bAtomList[i].name, Element::Chlorine());
        }
        else if(toupper(bAtomList[i].name[0]) == 'O'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(8, "Oxygen", "O", 16.00));
            UnivalentAtom(bAtomList[i].name, Element::Oxygen());
        }
        else if(toupper(bAtomList[i].name[0]) == 'F'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(9, "Fluorine", "F", 19.00));
            UnivalentAtom(bAtomList[i].name, Element::Fluorine());
        }
        else if(toupper(bAtomList[i].name[0]) == 'B'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(35, "Bromine", "Br", 79.90));
            UnivalentAtom(bAtomList[i].name, Element::Bromine());
        }
        else if(toupper(bAtomList[i].name[0]) == 'I'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(53, "Iodine", "I", 126.9));
            UnivalentAtom(bAtomList[i].name, Element::Iodine());
        }
        else if(toupper(bAtomList[i].name[0]) == 'N'){
          bAtomList[i].bAtomType = new
            //UnivalentAtom(bAtomList[i].name, Element(7, "Nitrogen", "N", 14.01));
            UnivalentAtom(bAtomList[i].name, Element::Nitrogen());
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.1112);
      }
      else if (bAtomList[i].nbonds == 2){
        if(toupper(bAtomList[i].elem) == 'H'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name, Element(1, "MMTKHydrogen", "MMTKH", 1.00797598)); // MMTK mass
            //BivalentAtom(bAtomList[i].name, MMTKElement::MMTKHydrogen()); // Molmodel Mass
            BivalentAtom(bAtomList[i].name, Element::Hydrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'C'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", 12.01));
            BivalentAtom(bAtomList[i].name,  Element::Carbon());
        }
        else if(toupper(bAtomList[i].elem) == 'O'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", 16.00),
            BivalentAtom(bAtomList[i].name,  Element::Oxygen(),
            109.47*Deg2Rad);
        }
        else if(toupper(bAtomList[i].elem) == 'N'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", 14.01));
            BivalentAtom(bAtomList[i].name,  Element::Nitrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'S'){
          bAtomList[i].bAtomType = new
            //BivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", 32.06),
            BivalentAtom(bAtomList[i].name,  Element::Sulfur(),
            109.47*Deg2Rad);
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
      }
      else if (bAtomList[i].nbonds == 3){
        if(toupper(bAtomList[i].elem) == 'C'){
          bAtomList[i].bAtomType = new
            //TrivalentAtom(bAtomList[i].name, Element(6, "Carbon", "C", 12.01),
            TrivalentAtom(bAtomList[i].name, Element::Carbon(),
              120*Deg2Rad, 120*Deg2Rad
            );
        }
        else if(toupper(bAtomList[i].elem) == 'O'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(8, "Oxygen", "O", 16.00));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Oxygen());
        }
        else if(toupper(bAtomList[i].elem) == 'N'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(7, "Nitrogen", "N", 14.01));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Nitrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'S'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(16, "Sulfur", "S", 32.06));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Sulfur());
        }
        else if(toupper(bAtomList[i].elem) == 'P'){
          bAtomList[i].bAtomType = new
            //TrivalentAtomTetra(bAtomList[i].name,  Element(15, "Phosphorus", "P", 30.97));
            TrivalentAtomTetra(bAtomList[i].name,  Element::Phosphorus());
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
      }
      else if (bAtomList[i].nbonds == 4){
        if(toupper(bAtomList[i].elem) == 'C'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", 12.01103690)); // MMTK mass
            QuadrivalentAtom(bAtomList[i].name,  Element::Carbon()); // Molmodel mass
        }
        else if(toupper(bAtomList[i].elem) == 'O'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", 16.00));
            QuadrivalentAtom(bAtomList[i].name,  Element::Oxygen());
        }
        else if(toupper(bAtomList[i].elem) == 'N'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", 14.01));
            QuadrivalentAtom(bAtomList[i].name,  Element::Nitrogen());
        }
        else if(toupper(bAtomList[i].elem) == 'S'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", 32.06));
            QuadrivalentAtom(bAtomList[i].name,  Element::Sulfur());
        }
        else if(toupper(bAtomList[i].elem) == 'P'){
          bAtomList[i].bAtomType = new
            //QuadrivalentAtom(bAtomList[i].name,  Element(15, "Phosphorus", "P", 30.97));
            QuadrivalentAtom(bAtomList[i].name,  Element::Phosphorus());
        }
        bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
      }

      bZeroCharArray(bAtomList[i].biotype, 20);
      sprintf(bAtomList[i].biotype, "%s_%s", \
        bAtomList[i].name, bAtomList[i].fftype);

      #ifdef DEBUG_LEVEL02
      //std::cout<<"bMoleculeReader::bMoleculeReader 4"<<std::endl<<std::flush;
      #endif
    }

    if (bAtomList == NULL){
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader constructor exit: NULL bAtomList\n";fflush(stdout);
      #endif
    }
    else{
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader constructor: bAtomList not NULL\n";fflush(stdout);
      #endif
    }
    
  }

  /* Just checking *////////
  #ifdef DEBUG_LEVEL02
  std::cout<<"Just checking\n";
  for(i=0; i<natoms;i++){
    printf(" -- name(%s) inName(%s) number(%d) elem(%c) fftype(%s) biotype(%s) charge(%f)\n", 
      bAtomList[i].name, bAtomList[i].inName, bAtomList[i].number, 
      bAtomList[i].elem, bAtomList[i].fftype,
      bAtomList[i].biotype, bAtomList[i].charge);
    fflush(stdout);
  }
  //for(i=0; i<bonds.size(); i++){ // RESTORE
  for(i=0; i<nbonds; i++){ // EU
    std::cout<<"bond: "<<bonds[i].i<<" "<<bonds[i].j<<std::endl;
    fflush(stdout);
  }
  #endif
  ///////////////////////////


  fclose(fpo);
  #ifdef DEBUG_LEVEL02
  std::cout<<"File closed succesfully\n";fflush(stdout);
  #endif
  //////////////////////////////////////////


  /*Now read rigid bodies specifications*/
  FILE *rfpo;
  rfpo = fopen(rbfilename, "r");
  if(rfpo == NULL){
    printf("Usage:\n<program> -mol2 <mol2_file> -rb <rb_file> -gaff  gaff.dat -frcmod <frcmod_file>\n");
    printf("rb_file not provided. Exiting...\n");
    exit(1);
  }

  #ifdef DEBUG_LEVEL02
  std::cout<<"\nREAD RIGID BODIES BEGIN\n"; fflush(stdout); 
  #endif

  std::string sbuff;
  std::vector<int> ring;
  char curr_char;
  int bond_found = 0;
  unsigned int boi;
  int par_open = 0;
  int ring_no = -1;
  int ring_closing = 0;

  while(fgets(line_c, MAX_LINE_LENGTH, rfpo)){
    line = line_c; //RESTORE
    if((sbuff = line.substr(0,5)) == "rings"){ // RESTORE
      #ifdef DEBUG_LEVEL02
      std::cout<<"bMoleculeReader::bMoleculeReader() sbuff"<<sbuff<<std::endl;
      #endif
      sbuff = "";
      //for(i=0; i<line.size(); i++){ // RESTORE
      for(i=0; i<strlen(line_c); i++){ // EU
        //curr_char = line.at(i); // RESTORE
        curr_char = line_c[i]; // RESTORE
        if(curr_char == ']'){
          ++ring_no;
          ring.push_back(atoi(sbuff.c_str()));
          #ifdef DEBUG_LEVEL02
          std::cout<<"push "<<atoi(sbuff.c_str())<<std::endl;
          #endif
          #ifdef DEBUG_LEVEL02
          std::cout<<"ring "<<ring_no<<": ";
          for(int tz; tz<ring.size(); tz++){std::cout<<ring[tz]<<' ';}
          std::cout<<std::endl;
          #endif
          sbuff = "";
          bond_found = 0;
          //for(boi=0; boi<bonds.size(); boi++){ // RESTORE
          for(boi=0; boi<nbonds; boi++){ // EU
            for(unsigned int ri=0; ri<ring.size(); ri++){
              for(unsigned int rj=0; (rj<ring.size()) && (rj<ri); rj++){
                if( ((ring[ri] == bonds[boi].i) && (ring[rj] == bonds[boi].j)) ||
                  ((ring[ri] == bonds[boi].j) && (ring[rj] == bonds[boi].i))){
                  bonds[boi].setInRing();
                  bonds[boi].setRingNo(ring_no);
                  if(ring_closing == 0){
                    bonds[boi].setAsRingClosing();
                  }
                  #ifdef DEBUG_LEVEL02
                  std::cout<<"RING BOND: "<<bonds[boi].i<<' '<<bonds[boi].j<<' '
                    <<bonds[boi].ringNo()<<' '
                    <<bonds[boi].isRingClosing()<<std::endl;
                  #endif
                  ++ring_closing;
                  bond_found = 1;
                  break;
                }
              }
            }
          }
          ring.clear();
          ring_closing = 0;
        }
        else if(curr_char == ','){
          ring.push_back(atoi(sbuff.c_str()));
          #ifdef DEBUG_LEVEL02
          std::cout<<"push "<<atoi(sbuff.c_str())<<std::endl;
          #endif
          sbuff = "";
        }
        else if((curr_char >= '0') && (curr_char <='9')){
          sbuff += curr_char;
        }
      }
    }

    else if((sbuff = line.substr(0,17)) == "non_ring_pi_bonds"){ // RESTORE
      #ifdef DEBUG_LEVEL02
      std::cout<<std::endl<<sbuff<<std::endl<<std::flush;
      #endif
      sbuff = "";
      for(i=0; i<line.size(); i++){
        //curr_char = line.at(i); // RESTORE
        curr_char = line_c[i]; // RESTORE
        if(curr_char == ')'){
          par_open = 0;
          ring.push_back(atoi(sbuff.c_str()));
          // Set the bond as rigid
          //for(boi = 0; boi<bonds.size(); boi++){ // RESTORE
          for(boi = 0; boi<nbonds; boi++){ // EU
            if( (bonds[boi].i == ring[0]) && (bonds[boi].j == ring[1]) ||
              (bonds[boi].i == ring[1]) && (bonds[boi].j == ring[0]) ){
              bonds[boi].setAsRigid();
              #ifdef DEBUG_LEVEL02
              std::cout<<"RIGID BOND: "<<bonds[boi].i
                <<' '<<bonds[boi].j<<std::endl;
              fflush(stdout);
              #endif
              ring.clear();
              sbuff = "";
              break;
            }
          }
        }
        else if((curr_char == ',') && (par_open == 1)){
          ring.push_back(atoi(sbuff.c_str()));
          sbuff = "";
        }
        else if((curr_char >= '0') && (curr_char <='9')){
          sbuff += curr_char;
        }
        else if(curr_char == '('){
          par_open = 1;
        }
      }
    }

    else if((sbuff = line.substr(0,17)) == "rigid_bodies"){
      // TODO
    }

  }

  #ifdef DEBUG_LEVEL01
  std::cout<<"\nREAD RIGID BODIES END\n\n"; fflush(stdout); 
  #endif

  fclose(rfpo);
}

bMoleculeReader::~bMoleculeReader(){;}




