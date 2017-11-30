#include "Topology.hpp"

using namespace SimTK;

Topology::Topology(){}
Topology::Topology(std::string argName){setName(argName);}


/* ==================================================
 *    BUILD AN UNLINKED MolModel NESTED IN THIS CLASS 
 * ================================================== */
void Topology::setMolModel(void){
    bMolAtomList = new MolAtom[natms];
    mol_MolModelCreate ("MainModel", &model);

    for(int k=0; k<natms; k++){
        bAtomAssign(&bMolAtomList[k], &bAtomList[k]);
        bMolAtomList[k].id = k;
        bMolAtomList[k].insertion_code = ' ';
        bMolAtomList[k].het = 1;
        bMolAtomList[k].chain_id = 'A';

       mol_MolModelCurrStrucGet (model, &struc);
       mol_StructureAtomAdd (struc, MOL_FALSE, &bMolAtomList[k]);
    }
    #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
    for(int tz=0; tz<natms; tz++){
        std::cout<<"bMainRes SpecAtom tz name bonds freebonds: "
        <<tz<<' '<<bAtomList[tz].name<<' '<<bAtomList[tz].nbonds<<' '<<bAtomList[tz].freebonds<<std::endl;
    }
    #endif

    mol_StructureChainsBuild(struc, 1);
    struc[0].chain_names = new char*[1];
    struc[0].chain_names[0] = &bMolAtomList[0].chain_id;
}

/* ==================================================
 *    Scale all DuMM force field terms by scale_factor
 * ================================================== */
void Topology::setDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm, SimTK::Real scale_factor){    
    dumm.setBondStretchGlobalScaleFactor(scale_factor);
    dumm.setBondBendGlobalScaleFactor(scale_factor);
    dumm.setBondTorsionGlobalScaleFactor(scale_factor);
    dumm.setAmberImproperTorsionGlobalScaleFactor(scale_factor);
    dumm.setVdw12ScaleFactor(scale_factor);
    dumm.setVdw13ScaleFactor(scale_factor);
    dumm.setVdw14ScaleFactor(scale_factor);
    dumm.setVdw15ScaleFactor(scale_factor);
    dumm.setVdwGlobalScaleFactor(scale_factor);
    dumm.setCoulomb12ScaleFactor(scale_factor);
    dumm.setCoulomb13ScaleFactor(scale_factor);
    dumm.setCoulomb14ScaleFactor(scale_factor);
    dumm.setCoulomb15ScaleFactor(scale_factor);
    dumm.setCoulombGlobalScaleFactor(scale_factor);
    dumm.setGbsaGlobalScaleFactor(scale_factor);
}

/* ==================================================
 *    Scale DuMM force field terms by scale_factor
 * ================================================== */
void Topology::setSpecificDuMMScaleFactor(SimTK::DuMMForceFieldSubsystem &dumm){    

    dumm.setBondStretchGlobalScaleFactor(0.0);

    dumm.setBondBendGlobalScaleFactor(0.0);

    dumm.setBondTorsionGlobalScaleFactor(0.0);

    dumm.setAmberImproperTorsionGlobalScaleFactor(0.0);

    dumm.setVdw12ScaleFactor(0.0);
    dumm.setVdw13ScaleFactor(0.0);
    dumm.setVdw14ScaleFactor(0.5);
    dumm.setVdw15ScaleFactor(1.0);
    //dumm.setVdwGlobalScaleFactor(1.0);

    dumm.setCoulomb12ScaleFactor(0.0);
    dumm.setCoulomb13ScaleFactor(0.0);
    dumm.setCoulomb14ScaleFactor(0.0);
    dumm.setCoulomb15ScaleFactor(0.0);
    dumm.setCoulombGlobalScaleFactor(0.0);

    dumm.setGbsaGlobalScaleFactor(0.0);
}


  Topology::~Topology(){
  }

  // Topology interface
  // Set graph

  void Topology::insertAtom(bSpecificAtom *){}
  void Topology::insertBond(int, int, int bondOrder){}

  // Parameters
 /*
  void Topology::setDuMMAtomParam(int, SimTK::Real vdw, SimTK::Real well){assert(!"Not implemented.");}
  void Topology::setDuMMBondParam(int, int, SimTK::Real k, SimTK::Real equil){assert(!"Not implemented.");}
  void Topology::setDuMMAngleParam(int, int, int, SimTK::Real k, SimTK::Real equil){assert(!"Not implemented.");}

  void Topology::setDuMMDihedralParam(int, int, int, int,
      int periodicity, SimTK::Real ampInKJ, SimTK::Real phaseInDegrees
  ){assert(!"Not implemented.");}
  void Topology::setDuMMDihedralParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2
  ){assert(!"Not implemented.");}
  void Topology::setDuMMDihedralParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2,
      int periodicity3, SimTK::Real ampInKJ3, SimTK::Real phaseInDegrees3
  ){assert(!"Not implemented.");}

  void Topology::setDuMMImproperParam(int, int, int, int,
      int periodicity, SimTK::Real ampInKJ, SimTK::Real phaseInDegrees
  ){assert(!"Not implemented.");}
  void Topology::setDuMMImproperParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2
  ){assert(!"Not implemented.");}
  void Topology::setDuMMImproperParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2,
      int periodicity3, SimTK::Real ampInKJ3, SimTK::Real phaseInDegrees3
  ){assert(!"Not implemented.");}
  */
  // Get

  int Topology::getNAtoms(void) const{assert(!"Not implemented.");}
  int Topology::getNBonds(void) const{assert(!"Not implemented.");}

  bSpecificAtom * Topology::getAtomByNumber(int number) const{assert(!"Not implemented.");}
  bSpecificAtom * Topology::getAtomByAtomIx(int aIx) const{assert(!"Not implemented.");}
  bSpecificAtom * Topology::getAtomByName(std::string name) const{assert(!"Not implemented.");}

  std::vector<bSpecificAtom> Topology::getNeighbours(int) const{assert(!"Not implemented.");}
  bBond * Topology::getBond(int, int) const{assert(!"Not implemented.");}
  int Topology::getBondOrder(int, int) const{assert(!"Not implemented.");}

///////////////////////
// Other functions
///////////////////////

// Process a graph node
//void Topology::process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode, int nofProcesses, int baseAtomNumber)
void Topology::process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode)
{
    baseSetFlag = 0;
    std::cout << " switch to " << node->number << std::endl;

    if (node->visited == CurrentGeneration) {
        return;
    }

    ++nofProcesses;

    std::cout << " update generation " << std::endl;
    std::cout << " left bond " << previousNode->number << ' ' << node->number << std::endl;    

    // Mark Gmolmodel bond and create bond in Molmodel
    for(std::vector<bBond *>::iterator it = (node->bondsInvolved).begin();
    it != (node->bondsInvolved).end(); ++it){
        if ((*it)->isThisBond(node->number, previousNode->number) ){
            (*it)->setVisited(CurrentGeneration);

            // Create bond in Molmodel
            if(nofProcesses == 1){
                std::cout << "Skipped the first step" << std::endl;
            }else if(nofProcesses == 2){
                if( baseSetFlag == 0 ){
                    std::cout << "Set base atom" << std::endl;
                    this->setBaseAtom( *(previousNode->bAtomType) );
                    //this->setAtomBiotype(previousNode->name, (this->name), previousNode->biotype);
                    this->setAtomBiotype(previousNode->name, (this->name), previousNode->getName());
                    this->convertInboardBondCenterToOutboard();
                    baseSetFlag = 1;
                }

                std::stringstream sbuff;
                if(previousNode->number == baseAtomNumber){
                    sbuff << previousNode->name << "/bond" << previousNode->freebonds;
                }else{
                    sbuff << previousNode->name << "/bond" << (previousNode->nbonds - previousNode->freebonds + 1);
                }

                std::cout << "Trying to connect " << node->name << "(" << node->getInName() << ") " 
                    << node->number << " to " << previousNode->number << "(" << previousNode->getInName() << ") "
                    << (sbuff.str()).c_str() << " ... " << std::flush;

                std::cout << "Node inName: " << *(node->inName) << std::endl << std::flush;

                this->bondAtom( *(node->bAtomType), (sbuff.str()).c_str(), 0.149, 0); // (Compound::SingleAtom&, BondCenterPathName, Length, Angle
                //this->setAtomBiotype(node->name, (this->name), node->biotype);

                std::cout << "Topology::process_node: setAtomBiotype: " 
                    << " node->name " << node->name << " (this->name).c_str() " << (this->name).c_str() << " node->getName() " << node->getName() 
                    << std::endl;

                this->setAtomBiotype(node->name, (this->name).c_str(), node->getName());

                // Set bSpecificAtom atomIndex to the last atom added to bond
                node->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ; // Set bSpecificAtom atomIndex to the last atom added to bond
                previousNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 0) ; // The only time we have to set atomIndex to the previous node

                --previousNode->freebonds;
                --node->freebonds;

                std::cout << "done." << std::endl << std::flush;

            }else if(nofProcesses > 2){
                std::stringstream sbuff;
                if(previousNode->number == baseAtomNumber){
                    sbuff << previousNode->name << "/bond" << previousNode->freebonds;
                }else{
                    sbuff << previousNode->name << "/bond" << (previousNode->nbonds - previousNode->freebonds + 1);
                }


                std::cout << "Trying to connect " << node->name << "(" << node->getInName() << ") " 
                    << node->number << " to " << previousNode->number << "(" << previousNode->getInName() << ") "
                    << (sbuff.str()).c_str() << " ... " << std::flush;

                this->bondAtom( *(node->bAtomType), (sbuff.str()).c_str(), 0.149, 0);
                //this->setAtomBiotype(node->name, (this->name), node->biotype);
                std::cout << "Topology::process_node: setAtomBiotype: " 
                    << " node->name " << node->name << " (this->name).c_str() " << (this->name).c_str() << " node->getName() " << node->getName() 
                    << std::endl;

                this->setAtomBiotype(node->name, (this->name), node->getName());

                // Set bSpecificAtom atomIndex to the last atom added to bond
                node->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ; // Set bSpecificAtom atomIndex to the last atom added to bond

                --previousNode->freebonds;
                --node->freebonds;

                std::cout << "done." << std::endl << std::flush;

            }
            break;
        }
    }
    // ========
    previousNode = node;

    node->visited = CurrentGeneration;

    std::cout << " start checking neighbors " << std::endl;
    unsigned int i;
    for(i = 0; i < (node->neighbors).size(); i++) {
        //process_node( (node->neighbors)[i], CurrentGeneration, previousNode, nofProcesses, baseAtomNumber );
        process_node( (node->neighbors)[i], CurrentGeneration, previousNode);
    }

    // At the end of the graph walking
    // Compound can no longer be a child to the geometry of another compound
    if(node->number == 0){
        //this->convertInboardBondCenterToOutboard(); 
    }
    std::cout << " end processing " << node->number << std::endl;
}

// Construct the molecule topology
//void Topology::walkGraph(bSpecificAtom *root, int baseAtomNumber)
void Topology::walkGraph(bSpecificAtom *root)
{
    nofProcesses = 0;
    int CurrentGeneration = 0;

    CurrentGeneration += 1;
    bSpecificAtom *previousNode = root;

    baseSetFlag = 0;
    PrintStaticVars();
    process_node(root, CurrentGeneration, previousNode);
    std::cout << std::endl;
}

void Topology::build(
    SimTK::DuMMForceFieldSubsystem &dumm,
    int natoms,
    bSpecificAtom *bAtomList,
    int nbonds,
    bBond *bonds, 
    std::string flexFN,
    std::string ictdF
)
{
    std::cout << "TOPOLOGY BUILD" << std::endl;

    this->natms = natoms;
    this->bAtomList = bAtomList;
    this->nbnds = nbonds;
    this->bonds = bonds;
    this->ictdF = ictdF;

    assert(bAtomList != NULL);
    this->setCompoundName((this->name));
    for(int i=0; i<natoms; i++){
        bAtomList[i].visited = 0;
    }
    for(int i=0; i<nbonds; i++){
        bonds[i].setVisited(0);
    }

    // Walk graph
    std::cout << "Walk the graph" << std::endl;
    bSpecificAtom *root = &(bAtomList[0]);
    baseAtomNumber = root->number;
    walkGraph( &(bAtomList[0]));

    // Add ring closing bonds
    for(int i=0; i<nbonds; i++){
        if(bonds[i].isVisited() == 0){
            
            bSpecificAtom *leftNode  =  &(bAtomList[bonds[i].i]);
            bSpecificAtom *rightNode =  &(bAtomList[bonds[i].j]);

            std::stringstream sbuff;
            if(leftNode->number == baseAtomNumber){
                sbuff << leftNode->name << "/bond" << leftNode->freebonds;
            }else{
                sbuff << leftNode->name << "/bond" << (leftNode->nbonds - leftNode->freebonds + 1);
            }
    
            std::stringstream otsbuff;
            if(rightNode->number == baseAtomNumber){
                otsbuff << rightNode->name << "/bond" << rightNode->freebonds;
            }else{
                otsbuff << rightNode->name << "/bond" << (rightNode->nbonds - rightNode->freebonds + 1);
            }
                  
            this->addRingClosingBond(
                (sbuff.str()).c_str(),
                (otsbuff.str()).c_str(),
                0.14,
                109*Deg2Rad,
                BondMobility::Rigid // CHANGE
            );
    
            //this->setAtomBiotype(leftNode->name, (this->name), leftNode->biotype);
            //this->setAtomBiotype(rightNode->name, (this->name), rightNode->biotype);
            this->setAtomBiotype(leftNode->name, (this->name), leftNode->getName());
            this->setAtomBiotype(rightNode->name, (this->name), rightNode->getName());
    
            //leftNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 0) ; // Not necessary
            //rightNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1) ; // Not necessary

            --leftNode->freebonds;
            --rightNode->freebonds;

            std::cout << "Closed ring " 
                << leftNode->name << "(" << leftNode->getInName() 
                << ") " << leftNode->number << " " << (sbuff.str()).c_str() << " to " 
                << rightNode->name << "(" << rightNode->getInName() 
                << ") " << rightNode->number << " " << (otsbuff.str()).c_str() 
                << " ... " << std::flush;

            std::cout << "done." << std::endl << std::flush;
        }
    }

    std::cout << "Topology: name inName atomIndexi:" << std::endl;
    for(int ix = 0; ix < getNumAtoms(); ++ix){
        std::cout << bAtomList[ix].name << " " << bAtomList[ix].inName << " " << bAtomList[ix].atomIndex << std::endl;
    }
    std::cout << std::endl << std::flush;

    // Assign Compound coordinates by matching bAtomList coordinates
    std::cout << "atomTargets from passed coords array: " << std::endl;
    std::map<AtomIndex, Vec3> atomTargets;
    for(int ix = 0; ix < getNumAtoms(); ++ix){
        std::cout << bAtomList[ix].atomIndex << " " 
            << this->getAtomName(bAtomList[ix].atomIndex) << " " 
            << bAtomList[ix].getX() << " " << bAtomList[ix].getY() << " " << bAtomList[ix].getZ()
            << std::endl;
        Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
        atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
    }

    std::cout << "Trying matchDefaultConfiguration Match_Exact ... " << std::flush;
    matchDefaultConfiguration(atomTargets, Match_Exact, true, 150.0); //Compound::Match_Idealized
    std::cout << "done. " << std::endl << std::flush;

    PdbStructure  pdb(*this);
    std::ostringstream sstream;
    sstream << "pdbs/sb_" << this->name <<"_ini"<<".pdb";
    std::string ofilename = sstream.str();
    std::filebuf fb;
    std::cout<<"Writing pdb file: "<<ofilename<<std::endl;
    fb.open(ofilename.c_str(), std::ios::out);
    std::ostream os(&fb);
    pdb.write(os); // automatically multiplies by ten (nm to A)
    fb.close();

    /* Just checking *////////
    std::cout << "Checking after Topology" << std::endl;
    for(int i=0; i<getNumAtoms();i++){
        bAtomList[i].Print();
    }
    for(int i=0; i<getNumBonds(); i++){ // EU
        std::cout<<"bond: "<<bonds[i].i<<" "<<bonds[i].j<<std::endl;
        fflush(stdout);
    }
    ///////////////////////////

    setRegimen(ictdF, flexFN);

}
  
// Set regimen
void Topology::setRegimen(std::string argRegimen, std::string flexFN){
    
    if(argRegimen == "IC"){
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Free, Compound::BondIndex(r));
        }
        std::cout << "Changed regimen to: " << "IC" << std::endl;
    }else if(argRegimen == "TD"){
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
        }
        std::cout << "Changed regimen to: " << "TD" << std::endl;
    }else if(argRegimen == "RB"){

        // Get flexible bonds from file and put it in PrmFlexBonds
        std::string line;
        std::ifstream F(flexFN);
        std::vector<std::pair<int, int>> PrmFlexBonds;
        std::vector<std::pair<int, int>>::iterator PrmFlexBondsIt;
        std::vector<std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>> MolmodelFlexBonds;
        std::vector<std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>>::iterator MolmodelFlexBondsIt;
        int line_i = -1;
        while(F.good()){
            line_i++;
            std::getline(F, line);
            std::istringstream iss(line);
            std::string word;
            std::vector<std::string> LineWords;
          
            int word_i = -1;
            while(iss >> word){
                if(word[0] == '#'){
                    break;
                }
                word_i++;
                LineWords.push_back(std::move(word));
            }
            if(word_i > 0){
                assert((word_i >= 1) && "2 indeces needed on each line of ligand.flex.");
                PrmFlexBonds.push_back(std::pair<int, int>( std::stod(LineWords[0]), std::stod(LineWords[1]) ));
            }
        }

        // Set all bonds to rigid first
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Rigid, SimTK::Compound::BondIndex(r));
        }

        // Iterate through prmtop flexible bonds, get Molmodel Compound::AtomIndeces and put it in MolmodelFlexBonds
        for ( PrmFlexBondsIt = PrmFlexBonds.begin(); PrmFlexBondsIt != PrmFlexBonds.end(); ++PrmFlexBondsIt){
            MolmodelFlexBonds.push_back(std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>(
                SimTK::Compound::AtomIndex(bAtomList[(*PrmFlexBondsIt).first].atomIndex), 
                SimTK::Compound::AtomIndex(bAtomList[(*PrmFlexBondsIt).second].atomIndex)) );
        }

        // Print
        for ( PrmFlexBondsIt = PrmFlexBonds.begin(); PrmFlexBondsIt != PrmFlexBonds.end(); ++PrmFlexBondsIt){
            std::cout << "Topology::setRegimen: FlexBond Prm Indeces " << (*PrmFlexBondsIt).first << " " << (*PrmFlexBondsIt).second << std::endl;
        }
        for ( MolmodelFlexBondsIt = MolmodelFlexBonds.begin(); MolmodelFlexBondsIt != MolmodelFlexBonds.end(); ++MolmodelFlexBondsIt){
            std::cout << "Topology::setRegimen: FlexBond AtomIndeces " << (*MolmodelFlexBondsIt).first << " " << (*MolmodelFlexBondsIt).second << std::endl;
        }
        for (unsigned int j = 0 ; j < getNumBonds(); j++){
            std::cout << "Topology::setRegimen: Compound Bonds AtomIndeces " << getBondAtomIndex(Compound::BondIndex(j), 0) << " " << getBondAtomIndex(Compound::BondIndex(j), 1) << std::endl;
        }

        // Iterate through  Molmodel Compound bonds and match with MolmodelFlexBonds
        for (unsigned int i = 0 ; i < MolmodelFlexBonds.size(); i++){
            for (unsigned int j = 0 ; j < getNumBonds(); j++){
                if( ((getBondAtomIndex(Compound::BondIndex(j), 0) == MolmodelFlexBonds[i].first) && (getBondAtomIndex(Compound::BondIndex(j), 1) == MolmodelFlexBonds[i].second)) || 
                    ((getBondAtomIndex(Compound::BondIndex(j), 1) == MolmodelFlexBonds[i].first) && (getBondAtomIndex(Compound::BondIndex(j), 0) == MolmodelFlexBonds[i].second)) ){
                    setBondMobility(BondMobility::Torsion, Compound::BondIndex(j));
                    std::cout << "Topology::setRegimen: Bond " << j << " set to torsion" << std::endl;
                    break;
                }
            }
        }
        std::cout << "Changed regimen to: " << "RB" << std::endl;
    }
    this->regimen = argRegimen;
}







