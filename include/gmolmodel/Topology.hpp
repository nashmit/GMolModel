#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

#include "bMoleculeReader.hpp"
#include "server.hpp"

#ifndef MAIN_RESIDUE_DEBUG_SPECIFIC
#define MAIN_RESIDUE_DEBUG_SPECIFIC 1
#endif
#ifndef MAIN_RESIDUE_DEBUG_LEVEL01
#define MAIN_RESIDUE_DEBUG_LEVEL01
#endif

#ifndef MAIN_RESIDUE_DEBUG_LEVEL02
#define MAIN_RESIDUE_DEBUG_LEVEL02
#endif

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
  int natms;
  bSpecificAtom *bAtomList;
  int nbnds; // EU
  bBond *bonds; // EU
  std::string ictdF;

  // Constructor
  Topology();
  // Destructor
  virtual ~Topology();

  // In case we already know the graph and order
  void init(
    SimTK::DuMMForceFieldSubsystem &dumm,
    int natms,
    bSpecificAtom *bAtomList,
    int nbnds,
    bBond *bonds, // EU
    bool first_time=true,
    std::string flexFN="ligand.flex",
    std::string ictdF="IC"
  );


  // Set a MolModel and a MolStructure - to be removed

  void setMolModel(void);

  // Scale all DuMM force field terms by scale_factor

  void setDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm, SimTK::Real scale_factor);

  // Scale specific DuMM force field terms by scale_factor

  void setSpecificDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm);

  // Interface

  // Set graph

  void insertAtom(bSpecificAtom *atom);
  void insertBond(int atom_no1, int atom_no2, int bondOrder);

  // Parameters

  void setDuMMAtomParam(int, SimTK::Real vdw, SimTK::Real well);
  void setDuMMBondParam(int, int, SimTK::Real k, SimTK::Real equil);
  void setDuMMAngleParam(int, int, int, SimTK::Real k, SimTK::Real equil);

  void setDuMMDihedralParam(int, int, int, int,
      int periodicity, SimTK::Real ampInKJ, SimTK::Real phaseInDegrees
  );
  void setDuMMDihedralParam(int, int, int, int, 
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2
  );
  void setDuMMDihedralParam(int, int, int, int, 
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2,
      int periodicity3, SimTK::Real ampInKJ3, SimTK::Real phaseInDegrees3
  );

  void setDuMMImproperParam(int, int, int, int,
      int periodicity, SimTK::Real ampInKJ, SimTK::Real phaseInDegrees
  );
  void setDuMMImproperParam(int, int, int, int, 
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2
  );
  void setDuMMImproperParam(int, int, int, int, 
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2,
      int periodicity3, SimTK::Real ampInKJ3, SimTK::Real phaseInDegrees3
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
 
  // Process a graph node
  void process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode, int nofProcesses, int baseAtomNumber);

  // Construct the molecule topology
  void walkGraph(bSpecificAtom *root, int baseAtomNumber);

  // Build Molmodel Compound
  void build(
      SimTK::DuMMForceFieldSubsystem &dumm,
      int natms,
      bSpecificAtom *bAtomList,
      int nbnds,
      bBond *bonds,
      std::string flexFN="ligand.flex",
      std::string ictdF="IC"
  );

};




#endif //TOPOLOGY_H_
