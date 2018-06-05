#ifndef __BADDPARAMS__
#define __BADDPARAMS__

#include <algorithm>

#include <iostream>
#include <fstream>
#include <sstream>

#include "bgeneral.hpp"
//#include "Simbody.h"
//#include "Molmodel.h"
#include "bMoleculeReader.hpp"

#ifndef DEBUG_PARAMS_LEVEL01
#define DEBUG_PARAMS_LEVEL01
#endif

#ifndef DEBUG_PARAMS_LEVEL02
#define DEBUG_PARAMS_LEVEL02
#endif
/*
*/

// Conveninet functions (Pointers have to be allocated first).
//SimTK::Real bDihedral(SimTK::Vec3, SimTK::Vec3, SimTK::Vec3, SimTK::Vec3);

bool Type2atomType(char *Type, char *atomType, int ATOMTYPE_MAX_LEN);
bool Type2atomType(string Type, char *atomType, int ATOMTYPE_MAX_LEN);

//==============================================================================
//                           FUNCTION AddGaffParams
//==============================================================================
/**
 * Add Parameters Function. It adds the parameters read from prmtop
 * file to the DuMM force field.
 **/
void bAddAllParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
);
void bAddBiotypes(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
);
void bAddAtomClasses(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
);
void bAddBondParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
);
void bAddAngleParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
);
void bAddTorsionParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
);

/**
 * Add gaff Parameters Function. It adds the parameters read from an Amber parameters
 * file to the DuMM force field.
void bAddGaffParams(
  SimTK::DuMMForceFieldSubsystem& dumm,
  const char *filename,
  int natms,
  bSpecificAtom *bAtomList,
  //vector<bBond> bonds, // RESTORE
  int nbnds, // EU
  bBond *bonds, // EU
  string frcmodfn
);
 **/

#endif  //__BADDPARAMS__

