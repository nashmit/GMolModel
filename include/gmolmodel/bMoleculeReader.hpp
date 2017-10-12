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
#include "bSpecificAtom.hpp"
#include "bBond.hpp"
#include "readAmberInput.hpp"

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
//                           CLASS intrad
//==============================================================================
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
  int natoms;
  int nbonds; // EU
  unsigned int MAX_LINE_LENGTH;

  bMoleculeReader(readAmberInput *amberReader, const char *);

  bMoleculeReader(SimTK::DuMMForceFieldSubsystem& dumm,
          const char *filename,
          const char *filetype,
          const char *rbfilename);

  ~bMoleculeReader();
};


#endif  //__BMOLECULEREADER__


