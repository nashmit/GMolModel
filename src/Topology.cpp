#include "Topology.hpp"

using namespace SimTK;

Topology::Topology(){}


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

    dumm.setBondStretchGlobalScaleFactor(1.0);

    dumm.setBondBendGlobalScaleFactor(1.0);

    dumm.setBondTorsionGlobalScaleFactor(1.0);

    dumm.setAmberImproperTorsionGlobalScaleFactor(1.0);

    dumm.setVdw12ScaleFactor(1.0);
    dumm.setVdw13ScaleFactor(1.0);
    dumm.setVdw14ScaleFactor(1.0);
    dumm.setVdw15ScaleFactor(1.0);
    dumm.setVdwGlobalScaleFactor(1.0);

    dumm.setCoulomb12ScaleFactor(1.0);
    dumm.setCoulomb13ScaleFactor(1.0);
    dumm.setCoulomb14ScaleFactor(1.0);
    dumm.setCoulomb15ScaleFactor(1.0);
    dumm.setCoulombGlobalScaleFactor(1.0);

    dumm.setGbsaGlobalScaleFactor(0.0);
}

/*Any kind of molecule*/
void Topology::init(
    DuMMForceFieldSubsystem &dumm,
    int natms,
    bSpecificAtom *bAtomList,
    int nbnds, // EU
    bBond *bonds, // EU
    bool first_time,
    std::string flexFN,
    std::string ictdF
)
{
    std::cout << "Topology::init START" << std::endl;
    this->natms = natms;
    this->bAtomList = bAtomList;
    this->nbnds = nbnds;
    this->bonds = bonds;
    this->ictdF = ictdF;
    assert(bAtomList != NULL);
  
    int noDummies = 0;
    std::stringstream sbuff;
    std::stringstream otsbuff;  // other sbuff
    std::string buff[nbnds]; // EU
    std::string otbuff;        // other buff
    int bondFirstAtom[nbnds], bondSecondAtom[nbnds]; // EU
    int k;

    this->setCompoundName("Topology");

    // Set up the connectivity
    // First atom
    int found = 0;
    int inter = -1;
    int firstBond = 0;  // First bond to be used
    std::vector<int> pushed;
    std::vector<int> crbonds;
    std::vector<int>::iterator pushedIt1;
    std::vector<int>::iterator pushedIt2;

    for (k=0; k < nbnds; k++){ // EU
      bondFirstAtom[k] = bonds[k].i;
      bondSecondAtom[k] = bonds[k].j ;
      #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
      std::cout<<"bMainRes: k bondFirstAtom[k] bondSecondAtom[k] "<<k<<' '<<bondFirstAtom[k]<<' '<<bondSecondAtom[k]<<std::endl;
      #endif
    }
    
    // Set base atom and link second one if dummies not present
    if(noDummies == 0){
      for (k=0; k<nbnds; k++){
        if((!bonds[k].isRingClosing()) && (bAtomList[bondFirstAtom[k]].freebonds != 1)){
          firstBond = k;
          break;
        }
      }
    }
    else{    // if dummies
      firstBond = nbnds - 2;
    }
    #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
    std::cout<<"bMainRes: firstBond "<<firstBond<<std::endl;
    #endif
    
    this->setBaseAtom( *(bAtomList[bondFirstAtom[firstBond]].bAtomType) );
    this->setAtomBiotype(bAtomList[bondFirstAtom[firstBond]].name, "mainRes", bAtomList[bondFirstAtom[firstBond]].biotype);
    this->convertInboardBondCenterToOutboard();

    sbuff.str("");
    sbuff<<bAtomList[bondFirstAtom[firstBond]].name<<"/bond"<<1;
    buff[firstBond] = sbuff.str();

    this->bondAtom(*bAtomList[bondSecondAtom[firstBond]].bAtomType, buff[firstBond].c_str());
    this->setAtomBiotype(bAtomList[bondSecondAtom[firstBond]].name, "mainRes", bAtomList[bondSecondAtom[firstBond]].biotype);
    bAtomList[bondFirstAtom[firstBond]].freebonds = 2;    // linked to Ground and bj (++ in loop) // WATCHOUT freebonds hardcoded
    #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
    std::cout<<"MainRes: bond "<<bAtomList[bondSecondAtom[firstBond]].name<<" to "<<bAtomList[bondFirstAtom[firstBond]].name<<' '
      <<buff[firstBond].c_str()<<std::endl;
    #endif
 
    pushed.push_back(bondFirstAtom[firstBond]);
    pushed.push_back(bondSecondAtom[firstBond]);
    
    // Rearrange and connect in Molmodel Compound
    assert(nbnds > 0);
    while(pushed.size() < (unsigned int)(2*nbnds)){ // EU
        for(int m=0; m<nbnds; m++){ // EU
            if(!bonds[m].isRingClosing()){
                pushedIt1 = find(pushed.begin(), pushed.end(), bondSecondAtom[m]);
                pushedIt2 = find(pushed.begin(), pushed.end(), bondFirstAtom[m]);
                found = 0;

                #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
                bool boolI = false, boolJ = false;
                if(pushedIt1 != pushed.end()){boolI=true;}
                if(pushedIt2 != pushed.end()){boolJ=true;}
                std::cout<<"bMainRes: conn bi bj boolI boolJ "
                    <<bondFirstAtom[m]<<' '<<bondSecondAtom[m]<<' '<<boolI<<' '<<boolJ<<std::endl;
                #endif
                
                if((pushedIt1 == pushed.end()) && (pushedIt2 != pushed.end())){
                // bj not found, bi found
                    found = 1;
                }
                if((pushedIt1 != pushed.end()) && (pushedIt2 == pushed.end())){
                // bj found, bi not found => swap
                    found = 1;
                    inter = bondFirstAtom[m]; bondFirstAtom[m] = bondSecondAtom[m]; bondSecondAtom[m] = inter; // swap
                }
                if(found == 1){
                    if(m != firstBond){

                        /*
                        Compound & setDefaultBondLength (mdunits::Length length,
                            const AtomPathName &atom1,
                            const AtomPathName &atom2);
                        Compound & setDefaultBondAngle (Angle angle, 
                            const AtomPathName &atom1,
                            const AtomPathName &atom2,
                            const AtomPathName &atom3);
                        Compound & setDefaultDihedralAngle (Angle angle, 
                            Compound::AtomIndex atom1, 
                            Compound::AtomIndex atom2,
                            Compound::AtomIndex atom3,
                            Compound::AtomIndex atom4);
                        */

                        sbuff.str("");
                        sbuff<<bAtomList[bondFirstAtom[m]].name<<"/bond"<<bAtomList[bondFirstAtom[m]].freebonds;
                        buff[m] = sbuff.str();
                        this->bondAtom(*bAtomList[bondSecondAtom[m]].bAtomType, buff[m].c_str());
                        this->setAtomBiotype(bAtomList[bondSecondAtom[m]].name, "mainRes", bAtomList[bondSecondAtom[m]].biotype);
                        if(bondFirstAtom[m] == bondFirstAtom[firstBond]){
                            ++bAtomList[bondFirstAtom[m]].freebonds; // The first atom
                        }
                        else{
                            --bAtomList[bondFirstAtom[m]].freebonds;
                        }
                        #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
                        std::cout<<"MainRes: bond "<<bAtomList[bondSecondAtom[m]].name<<" to "<<bAtomList[bondFirstAtom[m]].name<<' '
                            <<buff[m].c_str()<<std::endl;
                        #endif
                    }
                    pushed.push_back(bondFirstAtom[m]);
                    pushed.push_back(bondSecondAtom[m]);
                    break;
                }
            }
        }
        if(found == 0){
            break;
        }
    }
    // Close the rings in Molmodel Compound
      for(int m=0; m<nbnds; m++){ // EU
        if(bonds[m].isRingClosing()){
          #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
          std::cout<<"MainRes: RingClosing bonds found"<<std::endl;
          #endif
          found = 0;
          pushedIt1 = find(crbonds.begin(), crbonds.end(), m);
          if(pushedIt1 == crbonds.end()){
            found = 1;
          }
          if(found == 1){
              #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
              std::cout<<"MainRes: RingClosing bonds found set to 1"<<std::endl;
              #endif
              sbuff.str("");
              sbuff<<bAtomList[bondFirstAtom[m]].name
                <<"/bond"<<bAtomList[bondFirstAtom[m]].freebonds;
              buff[m] = sbuff.str();

              otsbuff.str("");
              otsbuff<<bAtomList[bondSecondAtom[m]].name              
                <<"/bond"<<bAtomList[bondSecondAtom[m]].freebonds;
              otbuff = otsbuff.str();
              
              #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
              std::cout<<"MainRes: attempt ring bond "<<bAtomList[bondFirstAtom[m]].name<<" to "<<bAtomList[bondSecondAtom[m]].name<<' '
                <<buff[m].c_str()<<' '<<otbuff.c_str()<<' '<<bondFirstAtom[m]<<' '<<bondSecondAtom[m]<<std::endl;
              #endif
              
              this->addRingClosingBond(
                buff[m].c_str(),
                otbuff.c_str(),
                0.14,
                109*Deg2Rad,
                BondMobility::Torsion // TODO
                );

              this->setAtomBiotype(bAtomList[bondFirstAtom[m]].name, "mainRes", bAtomList[bondFirstAtom[m]].biotype);
              this->setAtomBiotype(bAtomList[bondSecondAtom[m]].name, "mainRes", bAtomList[bondSecondAtom[m]].biotype);

              if(bondFirstAtom[m] == bondFirstAtom[firstBond]){
                ++bAtomList[bondFirstAtom[m]].freebonds; // The first atom
              }
              else{
                --bAtomList[bondFirstAtom[m]].freebonds;
              }
          }
          crbonds.push_back(m);
        }
      }
    
    std::cout << "Before connectivity:" << std::endl;
    for (int bondIndex = 0; bondIndex < nbnds; bondIndex++){
        bonds[bondIndex].Print();
    }
    for (int bondIndex = 0; bondIndex < nbnds; bondIndex++){
    }


    //Create charged atom types in DuMM
    //Must be called AFTER first mainRes is declared,
    //so Biotypes and atom classes will be defined
    std::string abuff;
    for(k=0; k<natms; k++){
      abuff =  "rob";
      abuff += bAtomList[k].biotype;


      DuMM::ChargedAtomTypeIndex tempChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
      std::cout << "defineChargedAtomType atomTypeIx "<<tempChargedAtomTypeIndex 
          << " atomTypeName " << abuff.c_str() << " atomClassIx " << bAtomList[k].getAtomClassIndex()
          << " partialChargeInE " << bAtomList[k].charge << std::endl;
      dumm.defineChargedAtomType(
        tempChargedAtomTypeIndex,
        abuff.c_str(),
        bAtomList[k].getAtomClassIndex(), // dumm.getAtomClassIndex(bAtomList[k].fftype), from MOL2
        bAtomList[k].charge
      );

      bAtomList[k].setChargedAtomTypeIndex(tempChargedAtomTypeIndex);

      // Associate a ChargedAtomTypeIndex with a Biotype index
      std::cout << "setBiotypeChargedAtomType bAtomList["<< k << "].getChargedAtomTypeIndex() "
          << bAtomList[k].getChargedAtomTypeIndex()
          << " bAtomList[" << k << "].biotype " << bAtomList[k].biotype 
          << std::endl << std::flush;
      dumm.setBiotypeChargedAtomType( 
        bAtomList[k].getChargedAtomTypeIndex(),
        Biotype::get("mainRes", bAtomList[k].biotype).getIndex()
      );

    }

    // Assign AtomIndex values to atoms in bAtomList[] by name REVISE CHANGE
    int ix = 0; // MINE change type
    int jx = 0; // MINE change type
    std::string cname, myname;
 
    // Set atomIndeces
    for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
        std::cout << "this->getAtomName(aIx)" << this->getAtomName(aIx)  << std::endl;
    }
    std::cout << std::endl;
    for(ix = 0; ix<natms; ix++){
        std::cout << "std::string(bAtomList[ix].name) " << std::string(bAtomList[ix].name)  << std::endl;
    }

    for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
      for(ix = 0; ix<natms; ix++){
        if(this->getAtomName(aIx) == std::string(bAtomList[ix].name)){ // compare SimTK::String with std::string
          bAtomList[ix].atomIndex = aIx;
          break;
        }
      }
    }
    #ifdef MAIN_RESIDUE_DEBUG_SPECIFIC
    std::cout<<"bAtomList[ix].atomIndeces assigned"<<std::endl<<std::flush;
    #endif

    // Assign bBond::BondIndex values in bonds[] from Molmodel::BondIndex
    int inumber, jnumber;
    Compound::AtomIndex iIx, jIx;
    for ( Compound::BondIndex bondIx(0); bondIx < this->getNumBonds(); ++bondIx){
        iIx = this->getBondAtomIndex(Compound::BondIndex(bondIx), 0);
        jIx = this->getBondAtomIndex(Compound::BondIndex(bondIx), 1);

        for(ix = 0; ix<natms; ix++){
            std::cout << "bAtomList[ix].getCompoundAtomIndex() " << bAtomList[ix].getCompoundAtomIndex()
                << " iIx " << iIx  << std::endl;
            if(bAtomList[ix].getCompoundAtomIndex() == iIx){
                inumber = bAtomList[ix].number;
            }
        }
        for(jx = 0; jx<natms; jx++){
            if(bAtomList[jx].getCompoundAtomIndex() == jIx){
                jnumber = bAtomList[jx].number;
            }
        }

        for(int m=0; m<nbnds; m++){ // EU
            if(((bonds[m].i == inumber) && (bonds[m].j == jnumber)) ||
              ((bonds[m].i == jnumber) && (bonds[m].j == inumber))){
                bonds[m].setBondIndex(Compound::BondIndex(bondIx));
            }
        }
    }
    #ifdef MAIN_RESIDUE_DEBUG_SPECIFIC
    std::cout<<"bonds[m].bondIndexes assigned"<<std::endl<<std::flush;
    #endif

    // Create atomTargets from passed coords array REVISE CHANGE
    std::cout<<"Create atomTargets from passed coords array"<<std::endl;
    std::map<AtomIndex, Vec3> atomTargets;  
    ix = 0;
    if(first_time == true){ // Take coordinates from main (crd / previously MMTK)
      std::cout<<"Take coordinates from main (mol2)"<<std::endl;
      for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
        Vec3 v(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
        atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, v));
        ix++;
      }
    }

    // Assign Compound coordinates by matching bAtomList coordinates
    matchDefaultTopLevelTransform(atomTargets);
    matchDefaultConfiguration(atomTargets, Match_Exact, true, 150.0); //Compound::Match_Idealized
    //fitDefaultConfiguration(atomTargets, 150.0);
    //matchDefaultConfiguration(atomTargets, Match_Idealized, true, 150.0); //Compound::Match_Idealized
    std::cout << "Configuration matched" << std::endl << std::flush;

    // Reinsert coordinates after matching
    PdbStructure  pdb(*this);
    /*
    for(unsigned int i=0; i<natms; i++){
      std::cout << "String(bAtomList[i].name " << String(bAtomList[i].name) << std::endl << std::flush;
      const PdbAtom& P = pdb.getAtom(String(bAtomList[i].name), PdbResidueId(1), String(" "));
      std::string s(P.getName());
      const Vec3& PC = P.getCoordinates();
      bAtomList[i].x = PC[0];
      bAtomList[i].y = PC[1];
      bAtomList[i].z = PC[2];
    }
    */

    // Set rigidity
    if(ictdF=="RB"){

        std::cout << "General regimen: " << "RB" << std::endl;
        std::cout << "Reading specs from " << flexFN << std::endl;

        std::ifstream flexIfStream(flexFN);

        int noFlexBonds;
        std::string line;
        std::getline(flexIfStream, line);
        std::istringstream iss(line);
        iss >> noFlexBonds;
        std::cout << "No of bonds " << noFlexBonds << std::endl;

        std::vector<int> flexBondsIxs (noFlexBonds, 0);
        std::getline(flexIfStream, line);
        std::istringstream iss2(line);
        std::cout << "Bonds indeces read";
        for(int i=0; i<noFlexBonds; i++){
            iss2 >> flexBondsIxs[i];
            std::cout << " " << flexBondsIxs[i];
        }
        std::cout << std::endl;
        flexIfStream.close();

        std::cout << "Bonds set mobility:" << std::endl;
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            if(std::find(flexBondsIxs.begin(), flexBondsIxs.end(), r) != flexBondsIxs.end()){
                setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
                std::cout << "Bond " << r << " set to torsion" << std::endl;
            }
            else{
                setBondMobility(BondMobility::Rigid, Compound::BondIndex(r));
                std::cout << "Bond " << r << " set to rigid" << std::endl;
            }
        }
        std::cout << "Bonds rigidity:" << std::endl;
        for ( Compound::BondIndex bondIx(0); bondIx < this->getNumBonds(); ++bondIx){
            
        }

    }

    else if(ictdF=="IC"){
        std::cout << "General regimen: " << "IC" << std::endl; 
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Free, Compound::BondIndex(r));
        }
    }
    else if(ictdF=="TD"){
        std::cout << "General regimen: " << "TD" << std::endl; 
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
        }
    }
    else if(ictdF=="RR"){ // Torsional dynamics with rigid rings
        std::cout << "General regimen: " << "RR" << std::endl; 
        for(int m=0; m<nbnds; m++){ // EU
            if(bonds[m].isInRing()){
                setBondMobility(BondMobility::Rigid, bonds[m].getBondIndex());
                std::cout<<"Bond "<<m<<"("<<bonds[m].getBondIndex()<<")"<<" rigidized ";
                for(ix = 0; ix < getNumAtoms(); ++ix){
                    if(bAtomList[ix].number == bonds[m].i || bAtomList[ix].number == bonds[m].j){
                        std::cout<<" inName "<<bAtomList[ix].inName<<" name  "<<bAtomList[ix].name;
                    }
                }
                std::cout<<std::endl;
            }
        }
    }
    else{
      fprintf(stderr, "Dynamics type unknown\n");
      exit(1);
    }

    std::ostringstream sstream;
    sstream<<"pdbs/sb_"<<"ini"<<".pdb";
    std::string ofilename = sstream.str();
    std::cout<<"Writing pdb file: "<<ofilename<<std::endl;
    std::filebuf fb;
    fb.open(ofilename.c_str(), std::ios::out);
    std::ostream os(&fb);
    pdb.write(os); // automatically multiplies by ten (nm to A)
    fb.close();

  }

  Topology::~Topology(){
  }

  // Topology interface
  // Set graph

  void Topology::insertAtom(bSpecificAtom *){}
  void Topology::insertBond(int, int, int bondOrder){}

  // Parameters

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
void Topology::process_node(bSpecificAtom *node, int CurrentGeneration, bSpecificAtom *previousNode, int nofProcesses, int baseAtomNumber)
{
    static int baseSetFlag = 0;
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
                    this->setAtomBiotype(previousNode->name, "mainRes", previousNode->biotype);
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

                this->bondAtom( *(node->bAtomType), (sbuff.str()).c_str(), 0.149, 0); // (Compound::SingleAtom&, BondCenterPathName, Length, Angle
                this->setAtomBiotype(node->name, "mainRes", node->biotype);

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
                this->setAtomBiotype(node->name, "mainRes", node->biotype);

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
        process_node( (node->neighbors)[i], CurrentGeneration, previousNode, nofProcesses, baseAtomNumber );
    }

    // At the end of the graph walking
    // Compound can no longer be a child to the geometry of another compound
    if(node->number == 0){
        //this->convertInboardBondCenterToOutboard(); 
    }
    std::cout << " end processing " << node->number << std::endl;
}

// Construct the molecule topology
void Topology::walkGraph(bSpecificAtom *root, int baseAtomNumber)
{
    static int nofProcesses = 0;
    int CurrentGeneration = 0;

    CurrentGeneration += 1;
    bSpecificAtom *previousNode = root;

    //++root->freebonds; // add Ground
    //root->freebonds = 1; // only the Ground
    
    process_node(root, CurrentGeneration, previousNode, nofProcesses, baseAtomNumber);
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
    this->natms = natoms;
    this->bAtomList = bAtomList;
    this->nbnds = nbonds;
    this->bonds = bonds;
    this->ictdF = ictdF;

    assert(bAtomList != NULL);
    this->setCompoundName("mainRes");

    // Print bonds
    std::cout << "Bonds before walk the graph" << std::endl;
    for(int i=0; i<nbonds; i++){
        bonds[i].Print();
    }

    // Search for appropriate atom/node to start walking the graph from
    /*
    for(int i=0; i<natoms; i++){
        //if(bAtomList[i].freebonds == 2){
        if(bAtomList[i].getInName() == std::string("OH")){
            // Walk graph
            std::cout << "Walk the graph" << std::endl;
            bSpecificAtom *root = &(bAtomList[i]);
            static int baseAtomNumber = root->number;
            walkGraph( root, baseAtomNumber);
            break;
        }
    }
    */

    // Walk graph
    std::cout << "Walk the graph" << std::endl;
    bSpecificAtom *root = &(bAtomList[0]);
    static int baseAtomNumber = root->number;
    walkGraph( &(bAtomList[0]), baseAtomNumber);

    // Print bonds
    std::cout << "Bonds after walk the graph" << std::endl;
    for(int i=0; i<nbonds; i++){
        bonds[i].Print();
    }

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
    
            this->setAtomBiotype(leftNode->name, "mainRes", leftNode->biotype);
            this->setAtomBiotype(rightNode->name, "mainRes", rightNode->biotype);
    
            --leftNode->freebonds;
            --rightNode->freebonds;

            std::cout << "Close ring " 
                << leftNode->name << "(" << leftNode->getInName() 
                << ") " << leftNode->number << " " << (sbuff.str()).c_str() << " to " 
                << rightNode->name << "(" << rightNode->getInName() 
                << ") " << rightNode->number << " " << (otsbuff.str()).c_str() 
                << " ... " << std::flush;

            std::cout << "done." << std::endl << std::flush;
        }
    }

    // Define DuMM charged atom types
    for(int k=0; k<natoms; k++){
      std::string abuff =  "rob";
      abuff += bAtomList[k].biotype;

      DuMM::ChargedAtomTypeIndex tempChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
      std::cout << "defineChargedAtomType atomTypeIx "<<tempChargedAtomTypeIndex 
          << " atomTypeName " << abuff.c_str() << " atomClassIx " << bAtomList[k].getAtomClassIndex()
          << " partialChargeInE " << bAtomList[k].charge << std::endl;
      dumm.defineChargedAtomType(
        tempChargedAtomTypeIndex,
        abuff.c_str(),
        bAtomList[k].getAtomClassIndex(), // dumm.getAtomClassIndex(bAtomList[k].fftype), from MOL2
        bAtomList[k].charge
      );

      bAtomList[k].setChargedAtomTypeIndex(tempChargedAtomTypeIndex);

      // Associate a ChargedAtomTypeIndex with a Biotype index
      std::cout << "setBiotypeChargedAtomType bAtomList["<< k << "].getChargedAtomTypeIndex() "
          << bAtomList[k].getChargedAtomTypeIndex()
          << " bAtomList[" << k << "].biotype " << bAtomList[k].biotype 
          << std::endl << std::flush;
      dumm.setBiotypeChargedAtomType( 
        bAtomList[k].getChargedAtomTypeIndex(),
        Biotype::get("mainRes", bAtomList[k].biotype).getIndex()
      );  

    }   

    // Assign AtomIndex values to atoms in bAtomList[] by name
    for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
      for(int ix = 0; ix<natoms; ix++){
        //std::cout << "Compound atom name " << this->getAtomName(aIx) << " vs "
        //    << "bAtomList name " << std::string(bAtomList[ix].name) << std::endl;
        if(this->getAtomName(aIx) == std::string(bAtomList[ix].name)){ // compare SimTK::String with std::string
          std::cout << "bAtomList[" << ix << "].atomIndex set to " << aIx << " "
              << std::string(bAtomList[ix].name) << " " << std::string(bAtomList[ix].inName) << std::endl;
          bAtomList[ix].atomIndex = aIx;
          break;
        }
      }
    }

    // Assign bBond::BondIndex values in bonds[] from Molmodel::BondIndex
    /*
    int inumber, jnumber;
    Compound::AtomIndex iIx, jIx;
    for ( Compound::BondIndex bondIx(0); bondIx < this->getNumBonds(); ++bondIx){
        iIx = this->getBondAtomIndex(Compound::BondIndex(bondIx), 0);
        jIx = this->getBondAtomIndex(Compound::BondIndex(bondIx), 1);

        for(int ix = 0; ix<natms; ix++){
            std::cout << "bAtomList[ix].getCompoundAtomIndex() " << bAtomList[ix].getCompoundAtomIndex()
                << " iIx " << iIx  << std::endl;
            if(bAtomList[ix].getCompoundAtomIndex() == iIx){
                inumber = bAtomList[ix].number;
            }
        }

        for(int jx = 0; jx<natms; jx++){
            if(bAtomList[jx].getCompoundAtomIndex() == jIx){
                jnumber = bAtomList[jx].number;
            }
        }

        for(int m=0; m<nbonds; m++){
            if( ((bonds[m].i == inumber) && (bonds[m].j == jnumber)) ||
                ((bonds[m].i == jnumber) && (bonds[m].j == inumber)) ){
                bonds[m].setBondIndex(Compound::BondIndex(bondIx));
            }
        }
    }
    */

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

    //std::cout << "Trying matchDefaultTopLevelTransform ... " << std::flush;
    //matchDefaultTopLevelTransform(atomTargets);
    //std::cout << "done. " << std::endl << std::flush;

    std::cout << "Trying matchDefaultConfiguration Match_Exact ... " << std::flush;
    matchDefaultConfiguration(atomTargets, Match_Exact, true, 150.0); //Compound::Match_Idealized
    std::cout << "done. " << std::endl << std::flush;

    //std::cout << "Trying matchDefaultConfigurationi Match_Idealized ... " << std::flush;
    //matchDefaultConfiguration(atomTargets, Match_Idealized, true, 150.0); //Compound::Match_Idealized
    //std::cout << "done. " << std::endl << std::flush;

    PdbStructure  pdb(*this);
    std::ostringstream sstream;
    sstream<<"pdbs/sb_"<<"ini"<<".pdb";
    std::string ofilename = sstream.str();
    std::filebuf fb;
    std::cout<<"Writing pdb file: "<<ofilename<<std::endl;
    fb.open(ofilename.c_str(), std::ios::out);
    std::ostream os(&fb);
    pdb.write(os); // automatically multiplies by ten (nm to A)
    fb.close();

    // Print bonds
    std::cout << "Bonds after closing the rings" << std::endl;
    for(int i=0; i<nbonds; i++){
        bonds[i].Print();
    }

    if(ictdF=="IC"){
        std::cout << "General regimen: " << "IC" << std::endl;
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Free, Compound::BondIndex(r));
        }
    }else if(ictdF=="TD"){
        std::cout << "General regimen: " << "TD" << std::endl;
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
        }
    }


}
  









