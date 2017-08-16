#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

#include "bMoleculeReader.hpp"
#include "bgeneral.hpp"
#include "server.hpp"

/*
*/
#ifndef MAIN_RESIDUE_DEBUG_SPECIFIC
#define MAIN_RESIDUE_DEBUG_SPECIFIC 1
#endif
#ifndef MAIN_RESIDUE_DEBUG_LEVEL01
#define MAIN_RESIDUE_DEBUG_LEVEL01
#endif

#ifndef MAIN_RESIDUE_DEBUG_LEVEL02
#define MAIN_RESIDUE_DEBUG_LEVEL02
#endif
//using namespace SimTK;

void mol_StructureChainsBuild (MolStructure *, int);

//==============================================================================
//                           CLASS MainResidue
//==============================================================================
/**
 * Main Residue Class. It represents the main compound.
 **/
class Topology : public SimTK::Compound{
public:

  MolAtom *bMolAtomList;
  MolStructure *struc;
  MolModel *model;
  bool hasBuiltSystem;
  unsigned int natms;
  bSpecificAtom *bAtomList;
  //std::vector<bBond> bonds; // RESTORE
  unsigned int nbnds; // EU
  bBond *bonds; // EU
  std::string ictdF;
  TARGET_TYPE *PrmToAx_po;
  TARGET_TYPE *MMTkToPrm_po;

  // Constructor
  Topology();
  // Destructor
  ~Topology();

  // In case we already know the graph and order
  void init(
    SimTK::DuMMForceFieldSubsystem &dumm,
    unsigned int natms,
    bSpecificAtom *bAtomList,
    unsigned int nbnds,
    //std::vector<bBond> bonds, // RESTORE
    bBond *bonds, // EU
    TARGET_TYPE **coords,
    TARGET_TYPE **indexMap,
    TARGET_TYPE *PrmToAx_po,
    TARGET_TYPE *MMTkToPrm_po,
    bool first_time=true,
    std::string ictdF="IC"
  );

  // Interface

  // Set graph

  void insertAtom(bSpecificAtom *);
  void insertBond(int, int, bondOrder);

  // Parameters

  void setDuMMAtomParams(int, vdw, well);
  void setDuMMBondParams(int, int, Real k, Real equil);
  void setDuMMAngleParams(int, int, int, Real k, Real equil);

  void setDuMMDihedralParams(int, int, int, int,
      int periodicity, Real ampInKJ, Real phaseInDegrees
  );
  void setDuMMDihedralParams(int, int, int, int, 
      int periodicity1, Real ampInKJ1, Real phaseInDegrees1,
      int periodicity2, Real ampInKJ2, Real phaseInDegrees2,
  );
  void setDuMMDihedralParams(int, int, int, int, 
      int periodicity1, Real ampInKJ1, Real phaseInDegrees1,
      int periodicity2, Real ampInKJ2, Real phaseInDegrees2,
      int periodicity3, Real ampInKJ3, Real phaseInDegrees3,
  );

  void setDuMMImproperParams(int, int, int, int,
      int periodicity, Real ampInKJ, Real phaseInDegrees
  );
  void setDuMMImproperParams(int, int, int, int, 
      int periodicity1, Real ampInKJ1, Real phaseInDegrees1,
      int periodicity2, Real ampInKJ2, Real phaseInDegrees2,
  );
  void setDuMMImproperParams(int, int, int, int, 
      int periodicity1, Real ampInKJ1, Real phaseInDegrees1,
      int periodicity2, Real ampInKJ2, Real phaseInDegrees2,
      int periodicity3, Real ampInKJ3, Real phaseInDegrees3,
  );

  // Get

  int getNAtoms(void) const;
  int getNBonds(void) const;

  bSpecificAtom * getAtomByNumber(int number) const;
  bSpecificAtom * getAtomByAtomIx(int aIx) const;
  bSpecificAtom * getAtomByName(std::string name) const;
 
  std::vector<bSpecificAtom> getNeighbours(int) const;
  bBond * getBond(int, int) const;
  int getBondOrder(int, int) const;
 
};

#endif //TOPOLOGY_H_
