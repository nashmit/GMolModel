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

#ifndef DEBUG_PARAMS_LEVEL01
#define DEBUG_PARAMS_LEVEL01
#endif

#ifndef DEBUG_PARAMS_LEVEL02
#define DEBUG_PARAMS_LEVEL02
#endif

/** Topological information (bonds graph) for one molecule.
It maps to one compound in Molmodel thus it is derived 
from Molmodel Compound class.
It does  things:
   - loads information from input files such as Amber input prmtop / inpcrd
   - adds parameters to a DuMM force field which belongs to the World class
     because one DuMM class should be used for multiple molecules
   - contructs the graph based on a list of bSpecificAtom objects each of 
     which already contains bonding information from the input files
   - defines the rigid bodies based on imput files provided by the users.
Contains a list of atoms bAtomList which consists of bSpecificAtom 
objects **/
class Topology : public SimTK::Compound{
public:

    /** Default Constructor. Sets the name of this molecule to 'no_name '.
    The name has no particular function and is not guaranteed to be unique.**/
    Topology();

    /** Constructor that sets the name of the molecule. The name has no 
    particular function and is not guaranteed to be unique. **/
    Topology(std::string nameOfThisMolecule);

    /** Default Destructor. **/
    virtual ~Topology();

    /** Reads data from a specific reader (readAmberInput for now) object **/
    void loadAtomAndBondInfoFromReader(readAmberInput *amberReader);

    /** Adds force field parameters read by the inputReader to DuMM **/
    void bAddAllParams(
        std::string resName
        , readAmberInput *amberReader
        , SimTK::DuMMForceFieldSubsystem& dumm
    );

    // Interface:
    /** Get the name of this molecule **/
    std::string getname(void){return this->name;}

    /** Set the name of this molecule **/
    void setName(std::string nameOfThisMolecule){this->name = nameOfThisMolecule;}

    /** Get the number of atoms. **/
    int getNAtoms(void) const;

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

    /** Builds the Compound's tree, closes the rings, matches the configuration
    on the graph using using Molmodels matchDefaultConfiguration and sets the
    general flexibility of the molecule. **/
    void build(
        SimTK::DuMMForceFieldSubsystem &dumm
        , std::string flexFN="ligand.flex"
        , std::string regimenSpec="IC"
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
    //RE bSpecificAtom *bAtomList;
    std::vector<bSpecificAtom > bAtomList; // RE

    int nbonds;
    //bBond *bonds; 
    std::vector<bBond > bonds; // RE

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

    /** Biotype is a Molmodel hook that is usually used to look up molecular
    force field specific parameters for an atom type. Gmolmodel defines a 
    new Biotype for each atom. The only thing that is specified is the element
    with info about name, atomic number, valence and mass. **/
    void bAddBiotypes(
        std::string resName
        , readAmberInput *amberReader
        , SimTK::DuMMForceFieldSubsystem& dumm
    );

    /** It calls DuMMs defineAtomClass, defineChargedAtomTye and
    setBiotypeChargedAtomType for every atom. These Molmodel functions contain
    information regarding the force field parameters. **/
    void bAddAtomClasses(
        std::string resName
        , readAmberInput *amberReader
        , SimTK::DuMMForceFieldSubsystem& dumm
    );

    /** Calls DuMM defineBondStretch. **/
    void bAddBondParams(
        std::string resName
        , readAmberInput *amberReader
        , SimTK::DuMMForceFieldSubsystem& dumm
    );

    /** Calls DuMM defineBondBend. **/
    void bAddAngleParams(
        std::string resName
        , readAmberInput *amberReader
        , SimTK::DuMMForceFieldSubsystem& dumm
    );

    /** Calls DuMM defineBondTorsion for 1, 2 and 3 periodicities **/
    void bAddTorsionParams(
        std::string resName
        , readAmberInput *amberReader
        , SimTK::DuMMForceFieldSubsystem& dumm
    );

    /** Graph building functions
    // Process a graph node **/
    //void process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode, int nofProcesses, int baseAtomNumber);
    void process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode);

    // Construct the molecule topology
    void walkGraph(bSpecificAtom *root);

};




#endif //TOPOLOGY_H_
