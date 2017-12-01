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
  Topology(std::string argName);

  // Destructor
  virtual ~Topology();

  // Set a MolModel and a MolStructure - to be removed
  void setMolModel(void);

  // Scale all DuMM force field terms by scale_factor
  void setDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm, SimTK::Real scale_factor);

  // Scale specific DuMM force field terms by scale_factor
  void setSpecificDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm);

  // Set graph
  void insertAtom(bSpecificAtom *atom);
  void insertBond(int atom_no1, int atom_no2, int bondOrder);

  // Interface
  int getNAtoms(void) const;
  int getNBonds(void) const;

  bSpecificAtom * getAtomByNumber(int number) const;
  bSpecificAtom * getAtomByAtomIx(int aIx) const;
  bSpecificAtom * getAtomByName(std::string name) const;
 
  std::vector<bSpecificAtom> getNeighbours(int) const;
  bBond * getBond(int, int) const;
  int getBondOrder(int, int) const;
 
  // Process a graph node
  //void process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode, int nofProcesses, int baseAtomNumber);
  void process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode);

  // Construct the molecule topology
  void walkGraph(bSpecificAtom *root);

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

  void setRegimen(std::string argRegimen, std::string flexFN);

  std::string getname(void){return this->name;}
  void setName(std::string argName){this->name = argName;}

  void PrintStaticVars(void){
      std::cout << "Topology static vars:" 
          << " regimen " << regimen << " name " << name
          << " nofProcesses " << nofProcesses << " baseSetFlag " << baseSetFlag << " baseAtomNumber " << baseAtomNumber
          << std::endl;
  };

  void loadMaps(void);

  std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > getMbx2aIx(void){
      return mbx2aIx;
  }

  std::map< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex > getAIx2mbx(void){
      return aIx2mbx;
  }

private:
  std::string regimen;
  std::string name;
  int nofProcesses;
  int baseSetFlag;
  int baseAtomNumber;

  // Map mbx2aIx contains only atoms at the origin of mobods
  std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > mbx2aIx;
  std::map< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex > aIx2mbx;

};




#endif //TOPOLOGY_H_
