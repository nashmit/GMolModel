#ifndef __BMOLECULEREADER__
#define __BMOLECULEREADER__

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file
 * This defines the bMoleculeReader class and additional heloer classes
 **/

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"
#include "readAmberInput.hpp"
//#include "Molmodel.h"
//#include "mol.h"

//#include "SimTKcommon.h"
//#include "molmodel/internal/common.h"
//#include "molmodel/internal/GrinPointer.h"
//#include "molmodel/internal/units.h"

/*
*/
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif

#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif

//==============================================================================
//                           CLASS TrivalentAtomTetra
//==============================================================================
/** 
 * Trivalent Atom Class with tetrahedral geometry.
 * Bond centers are named "bond1", "bond2", and "bond3"
**/

class  TrivalentAtomTetra : public SimTK::Compound::SingleAtom {
 public:
  TrivalentAtomTetra(
    const SimTK::Compound::AtomName& atomName,   ///< name for new atom
    const SimTK::Element& element              /// element for new atom
  );
};


//==============================================================================
//                           CLASS PDBReader
//==============================================================================
/** 
 * Pdb File Reader Class. Not necessary.
**/
class bPDBReader{
 public:
  bPDBReader();
  ~bPDBReader();
};


//==============================================================================
//                           CLASS SpecificAtom
//==============================================================================
/** 
 * gMolmodel Specific Atom Type Class.
 * This incorporates additional Amber forcefield data.
**/
class bSpecificAtom{ /*Information like in sdf*/
public:
    int nbonds;
    int freebonds;
    char name[5];
    char inName[5];
    int number;
    char elem;
    int atomicNumber;
    SimTK::Real mass;
    SimTK::Real vdwRadius;
    SimTK::Real LJWellDepth; // Lennard-Jones well depth
    double charge;
    char fftype[20];
    SimTK::DuMM::AtomClassIndex atomClassIndex;
    SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;

    float x;
    float y;
    float z;
    char biotype[20];
    SimTK::BiotypeIndex biotypeIndex;
    SimTK::Compound::SingleAtom *bAtomType;
    SimTK::Compound::AtomIndex atomIndex;
    int mobile;
    int visited;
    
    // Suggested by Eliza
    std::string residueName;
    long int residueIndex;
    std::string chain;
    int moleculeIndex;
  
    bSpecificAtom();
    ~bSpecificAtom();
    void Print(void);
    void Zero(void);
  
    // Interface
  
    int getNBonds(void);
    int getFreebonds(void);
    std::string getName(void);
    std::string getInName(void);
    int getNumber(void);
    char getElem(void);
    SimTK::mdunits::Mass getMass(void);
    void setMass(SimTK::Real);
    SimTK::Real getX(void);
    SimTK::Real getY(void);
    SimTK::Real getZ(void);
    std::string getFftype(void);
    SimTK::DuMM::AtomClassIndex getAtomClassIndex(void);
    void setAtomClassIndex(SimTK::DuMM::AtomClassIndex);
    const char * getBiotype(void);
    SimTK::Compound::SingleAtom * getBAtomType(void);
    SimTK::Compound::AtomIndex getCompoundAtomIndex(void);
    SimTK::Real getCharge(void);
    int getIsMobile(void);
    int getIsVisited(void);

    int getAtomicNumber(void);
    void setAtomicNumber(int);
    SimTK::Real getVdwRadius(void);
    void setVdwRadius(SimTK::Real);
    SimTK::Real getLJWellDepth(void);
    void setLJWellDepth(SimTK::Real);
    SimTK::DuMM::ChargedAtomTypeIndex getChargedAtomTypeIndex(void);
    void setChargedAtomTypeIndex(SimTK::DuMM::ChargedAtomTypeIndex);

    void setNbonds(int);
    void setFreebonds(int);
    void setName(std::string);
    void setInName(std::string);
    void setNumber(int);
    void setElem(char);
    void setX(SimTK::Real);
    void setY(SimTK::Real);
    void setZ(SimTK::Real);
    void setFftype(std::string);
    void setBiotype(std::string);
    void setBiotype(const char *);
    void setBAtomType(SimTK::Compound::SingleAtom *);
    void setCompoundAtomIndex(SimTK::Compound::AtomIndex);
    void setCharge(SimTK::Real);
    void setIsMobile(int);
    void setIsVisited(int);

};


int bAtomAssign(MolAtom *dest, const bSpecificAtom *src);

//==============================================================================
//                           CLASS intpair
//==============================================================================
/** 
 * Intpair Class is a two int vector used for connectivity definition in MoleculeReader.
**/
class intpair{
 public:
  int i; int j; // These will correspond to bSpecificAtom.number
  intpair();
  intpair(int inI, int inJ);
  ~intpair();

  bool operator==(const intpair *other);
  bool operator!=(const intpair *other);
  bool isTheSameAs(const intpair *other);
  void swap(void);
  void dump(void);
  std::string getString(void);
};

//==============================================================================
//                           CLASS Bond
//==============================================================================
/** 
 * Bond Class used for connectivity definition in MoleculeReader.
**/
class bBond : public intpair{
 private:
  int inring;
  int ring_closing;
  int rigid;
  int ring_no;
  SimTK::Compound::BondIndex bondIndex;

 public:
  bBond(void);
  bBond(int a, int b);
  ~bBond(void);

  bool isInRing(void);
  bool isRingClosing(void);
  bool isRigid(void);
  int ringNo(void);

  void setInRing(void);
  void setAsRingClosing(void);
  void setAsRigid(void);
  void setRingNo(int rn);

  SimTK::Compound::BondIndex getBondIndex(void);
  void setBondIndex(SimTK::Compound::BondIndex otherIx);
};


//==============================================================================
//                           CLASS intrad
//==============================================================================
/** 
 * intriad Class is a three int vector used for connectivity definition in MoleculeReader.
**/
class intriad{
 public:
  int i;
  int j;
  int k;
  intriad();
  intriad(int inI, int inJ, int inK);
  ~intriad();
  //Sorry for crowding
  bool operator==(const intriad *other);
  bool isTheSameAs(const intriad *other);
  void dump(void);
  std::string getString(void);
};


//==============================================================================
//                           CLASS MoleculeReader
//==============================================================================
/** 
 * Molecule Reader Class. Creates a list of gMolmodel specific atoms bAtomList
 * from a molecular structure file (mol2) which contain all the information needed
 * to create a Compound
**/
class bMoleculeReader{
 public:
  bSpecificAtom *bAtomList;
  //std::vector<bBond> bonds; // RESTORE
  bBond *bonds;
  unsigned int natoms;
  unsigned int nbonds; // EU
  unsigned int MAX_LINE_LENGTH;

  bMoleculeReader(readAmberInput *amberReader, const char *);

  bMoleculeReader(SimTK::DuMMForceFieldSubsystem& dumm,
          const char *filename,
          const char *filetype,
          const char *rbfilename);

  ~bMoleculeReader();
};


#endif  //__BMOLECULEREADER__


