#ifndef TOPOLOGY_H_
#define TOPOLOGY_H_

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

//RE #include "bMoleculeReader.hpp"
#include "TrivalentAtomTetra.hpp"
#include "bSpecificAtom.hpp"
#include "bBond.hpp"
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

    /** Constructor that sets the name of the molecule. The name has no 
    particular function and is not guaranteed to be unique. **/
    Topology(std::string nameOfThisMolecule);

    /** Default Destructor **/
    virtual ~Topology();

    // Interface:

    /** Get the name of this molecule **/
    std::string getname(void){return this->name;}

    /** Set the name of this molecule **/
    void setName(std::string nameOfThisMolecule){this->name = nameOfThisMolecule;}

    /** Get the number of atoms. **/
    int getNAtoms(void) const;

    /** Reads data from a specific reader (readAmberInput for now) object **/
    void loadAtomAndBondInfoFromReader(readAmberInput *amberReader);

    /** Get the number of bonds. **/
    int getNBonds(void) const;

    /** Get a pointer to an atom object in the atom list inquiring
    by number **/
    bSpecificAtom * getAtomByNumber(int number) const;

    /** Get a pointer to an atom object in the atom list inquiring
    by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
    bSpecificAtom * getAtomByAtomIx(int aIx) const;

    /** Get a pointer to an atom object in the atom list inquiring
    by atom name **/
    bSpecificAtom * getAtomByName(std::string name) const;

    /** Get the neighbours in the graph **/
    std::vector<bSpecificAtom *> getNeighbours(int) const;

    /** **/
    bBond * getBond(int, int) const;

    /** Get bond order. **/
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
        int natoms,
        bSpecificAtom *bAtomList,
        int nbonds,
        bBond *bonds,
        std::string flexFN="ligand.flex",
        std::string regimenSpec="IC"
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

    // Get coordinates
    void getCoordinates(std::vector<SimTK::Real> Xs, std::vector<SimTK::Real> Ys, std::vector<SimTK::Real> Zs);

public:

    bool hasBuiltSystem;
    int natoms;
    bSpecificAtom *bAtomList;
    int nbonds; // EU
    bBond *bonds; // EU
    std::string regimenSpec;
    unsigned int MAX_LINE_LENGTH;

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
    std::map< SimTK::Compound::BondIndex, int > bondIx2bond;

};




#endif //TOPOLOGY_H_
