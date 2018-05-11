#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

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

/** Topological information (bonds graph) for one molecule.
 * It maps to one compound in Molmodel thus it is derived 
 * from Molmodel Compound class.
 * Has two main functionalities:
   - contructs the graph based on a list of bSpecificAtom objects each of 
     which already contains bonding information from the inmput files
   - defines the rigid bodies.
 * Contains a list of atoms bAtomList which consists of bSpecificAtom 
 * objects 
 **/
class Topology : public SimTK::Compound{
public:

    /** Default Constructor **/
    Topology();

    /** Constructor that sets the name of the molecule.**/
    Topology(std::string nameOfThisMolecule);

    /** Default Destructor **/
    virtual ~Topology();

    // Interface functions
    std::string getname(void){return this->name;}
    void setName(std::string nameOfThisMolecule){this->name = nameOfThisMolecule;}
    int getNAtoms(void) const;
    int getNBonds(void) const;

    bSpecificAtom * getAtomByNumber(int number) const;
    bSpecificAtom * getAtomByAtomIx(int aIx) const;
    bSpecificAtom * getAtomByName(std::string name) const;
 
    std::vector<bSpecificAtom> getNeighbours(int) const;
    bBond * getBond(int, int) const;
    int getBondOrder(int, int) const;
 
    std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > getMbx2aIx(void){
        return mbx2aIx;
    }
    std::map< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex > getAIx2mbx(void){
        return aIx2mbx;
    }

    void writePdb(std::string dirname, std::string prefix, std::string sufix, int maxNofDigits, int index) const;

    // Graph building functions
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
    std::string getRegimen(void);

    void PrintStaticVars(void){
        std::cout << "Topology static vars:" 
            << " regimen " << regimen << " name " << name
            << " nofProcesses " << nofProcesses << " baseSetFlag " << baseSetFlag << " baseAtomNumber " << baseAtomNumber
            << std::endl;
    };

    void loadMaps(void);
    void printMaps(void);

    // Not sure they belong here
    //// Scale all DuMM force field terms by scale_factor
    void setDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm, SimTK::Real scale_factor);
    
    // Scale specific DuMM force field terms by scale_factor
    void setSpecificDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm);
   
    // Not sure we need them 
    // Set graph
    void insertAtom(bSpecificAtom *atom);
    void insertBond(int atom_no1, int atom_no2, int bondOrder);

public:

    bool hasBuiltSystem;
    int natms;
    bSpecificAtom *bAtomList;
    int nbnds; // EU
    bBond *bonds; // EU
    std::string ictdF;

private:

    std::string regimen;
    std::string name;
    int nofProcesses;
    int baseSetFlag;
    int baseAtomNumber;

    // Map mbx2aIx contains only atoms at the origin of mobods
    std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex > mbx2aIx;
    std::map< SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex > aIx2mbx;
    //std::map< SimTK::Compound::AtomIndex, int > aIx2number;

};




#endif //TOPOLOGY_H_
