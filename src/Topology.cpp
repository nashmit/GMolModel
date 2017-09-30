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

/*Any kind of molecule*/
void Topology::init(
    DuMMForceFieldSubsystem &dumm,
    unsigned int natms,
    bSpecificAtom *bAtomList,
    //std::vector<bBond> bonds, // RESTORE
    unsigned int nbnds, // EU
    bBond *bonds, // EU
    TARGET_TYPE **coords,
    TARGET_TYPE **indexMap,
    TARGET_TYPE *PrmToAx_po,
    TARGET_TYPE *MMTkToPrm_po,
    bool first_time,
    std::string ictdF
)
{
    this->natms = natms;
    this->bAtomList = bAtomList;
    this->nbnds = nbnds;
    this->bonds = bonds;
    this->ictdF = ictdF;
    this->PrmToAx_po = PrmToAx_po;
    this->MMTkToPrm_po = MMTkToPrm_po;
    assert(bAtomList != NULL);
  
    int noDummies = 0;
    std::stringstream sbuff;
    std::stringstream otsbuff;  // other sbuff
    std::string buff[nbnds]; // EU
    std::string otbuff;        // other buff
    int bondFirstAtom[nbnds], bondSecondAtom[nbnds]; // EU
    unsigned int k;

    // Reinsert setMolModel here if needed
     /* ========================================
      *    BUILD THIS CLASS' MOLECULE TOPOLOGY
      * ========================================*/

    //setDuMMScaleFactor(dumm, 0.0);

    this->setCompoundName("Topology");

    /*First atom*/
    int found = 0;
    int inter = -1;
    unsigned int firstBond = 0;  // First bond to be used
    std::vector<int> pushed;
    std::vector<int> crbonds;
    std::vector<int>::iterator pushedIt1;
    std::vector<int>::iterator pushedIt2;

    /* Set up the connectivity */
    for (k=0; k < nbnds; k++){ // EU
      //bondFirstAtom[k] = bonds[k].i - 1; // ONLY FOR MOL2
      //bondSecondAtom[k] = bonds[k].j - 1; // ONLY FOR MOL2
      bondFirstAtom[k] = bonds[k].i;
      bondSecondAtom[k] = bonds[k].j ;
      #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
      std::cout<<"bMainRes: k bondFirstAtom[k] bondSecondAtom[k] "<<k<<' '<<bondFirstAtom[k]<<' '<<bondSecondAtom[k]<<std::endl;
      #endif
    }
    
    /*Set base atom and link second one if dummies not present*/
    if(noDummies == 0){
      for (k=0; k < nbnds; k++){ //EU
        if((!bonds[k].isRingClosing()) && (bAtomList[bondFirstAtom[k]].freebonds != 1)){
          firstBond = k;
          break;
        }
      }
    }
    else{    // if dummies
      firstBond = nbnds - 2;  // EU
    }
    #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
    std::cout<<"bMainRes: firstBond "<<firstBond<<std::endl;
    #endif
    
    this->setBaseAtom( *(bAtomList[bondFirstAtom[firstBond]].bAtomType) );
    //this->setAtomBiotype(bAtomList[bondFirstAtom[firstBond]].name, "mainRes", bAtomList[bondFirstAtom[firstBond]].biotype);
    this->convertInboardBondCenterToOutboard();

    sbuff.str("");
    sbuff<<bAtomList[bondFirstAtom[firstBond]].name<<"/bond"<<1;
    buff[firstBond] = sbuff.str();

    this->bondAtom(*bAtomList[bondSecondAtom[firstBond]].bAtomType, buff[firstBond].c_str());
    //this->setAtomBiotype(bAtomList[bondSecondAtom[firstBond]].name, "mainRes", bAtomList[bondSecondAtom[firstBond]].biotype);
    bAtomList[bondFirstAtom[firstBond]].freebonds = 2;    // linked to Ground and bj (++ in loop) // WATCHOUT freebonds hardcoded
    #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
    std::cout<<"MainRes: bond "<<bAtomList[bondSecondAtom[firstBond]].name<<" to "<<bAtomList[bondFirstAtom[firstBond]].name<<' '
      <<buff[firstBond].c_str()<<std::endl;
    #endif
 
    pushed.push_back(bondFirstAtom[firstBond]);
    pushed.push_back(bondSecondAtom[firstBond]);
    
    // Rearrange and connect in Molmodel Compound
    bool boolI, boolJ;
    while(pushed.size() < 2*nbnds){ // EU
      for(unsigned int m=0; m<nbnds; m++){ // EU
        if(!bonds[m].isRingClosing()){
          pushedIt1 = find(pushed.begin(), pushed.end(), bondSecondAtom[m]);
          pushedIt2 = find(pushed.begin(), pushed.end(), bondFirstAtom[m]);
          found = 0;

          boolI = boolJ = false;
          if(pushedIt1 != pushed.end()){boolI=true;}
          if(pushedIt2 != pushed.end()){boolJ=true;}
          #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
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
              sbuff.str("");
              sbuff<<bAtomList[bondFirstAtom[m]].name<<"/bond"<<bAtomList[bondFirstAtom[m]].freebonds;
              buff[m] = sbuff.str();
              this->bondAtom(*bAtomList[bondSecondAtom[m]].bAtomType, buff[m].c_str());
              //this->setAtomBiotype(bAtomList[bondSecondAtom[m]].name, "mainRes", bAtomList[bondSecondAtom[m]].biotype);
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
      //for(unsigned int m=0; m<bonds.size(); m++){ // RESTORE
      for(unsigned int m=0; m<nbnds; m++){ // EU
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
              
              //std::cout<<bondSecondAtom[m]+1<<" "<<bondFirstAtom[m]+1<<std::endl
              //  <<buff[m].c_str()<<' '<<otbuff.c_str()<<std::endl;
              #ifdef MAIN_RESIDUE_DEBUG_LEVEL02
              std::cout<<"MainRes: attempt ring bond "<<bAtomList[bondFirstAtom[m]].name<<" to "<<bAtomList[bondSecondAtom[m]].name<<' '
                <<buff[m].c_str()<<' '<<otbuff.c_str()<<' '<<bondFirstAtom[m]<<' '<<bondSecondAtom[m]<<std::endl;
              #endif
              
              this->addRingClosingBond(
                buff[m].c_str(),
                otbuff.c_str(),
                0.14,
                109*Deg2Rad,
                //BondMobility::Rigid
                BondMobility::Torsion // TODO
                );
              //this->setAtomBiotype(bAtomList[bondFirstAtom[m]].name, "mainRes", bAtomList[bondFirstAtom[m]].biotype);
              //this->setAtomBiotype(bAtomList[bondSecondAtom[m]].name, "mainRes", bAtomList[bondSecondAtom[m]].biotype);
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

  /*
   Create charged atom types in DuMM
   Must be called AFTER first mainRes is declared,
   so Biotypes and atom classes will be defined
  */
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

      // Why is this necessary ??? 
      dumm.setBiotypeChargedAtomType( 
        bAtomList[k].getChargedAtomTypeIndex(),
        Biotype::get("mainRes", bAtomList[k].biotype).getIndex()
      );

    }

    // Assign AtomIndex values to atoms in bAtomList[] by name REVISE CHANGE
    unsigned int ix = 0; // MINE change type
    unsigned int jx = 0; // MINE change type
    std::string cname, myname;
    for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
      for(ix = 0; ix<natms; ix++){
        if(this->getAtomName(aIx) == std::string(bAtomList[ix].name)){ // compare SimTK::String with std::string
          bAtomList[ix].atomIndex = aIx;
          break;
        }
      }
    }
    #ifdef MAIN_RESIDUE_DEBUG_SPECIFIC
    std::cout<<"bAtomList[ix].atomIndexs assigend"<<std::endl<<std::flush;
    #endif

    // Assign BondIndex values to bonds in bonds[] REVISE CHANGE
    int inumber, jnumber;
    Compound::AtomIndex iIx, jIx;
    for (unsigned int r=0 ; r<getNumBonds(); r++){
      iIx = getBondAtomIndex(Compound::BondIndex(r), 0);
      jIx = getBondAtomIndex(Compound::BondIndex(r), 1);

      for(ix = 0; ix<natms; ix++){
        if(bAtomList[ix].atomIndex == iIx){
          inumber = bAtomList[ix].number;
        }
      }
      for(jx = 0; jx<natms; jx++){
        if(bAtomList[jx].atomIndex == jIx){
          jnumber = bAtomList[jx].number;
        }
      }

      //for(unsigned int m=0; m<bonds.size(); m++){ // RESTORE
      for(unsigned int m=0; m<nbnds; m++){ // EU
        if(((bonds[m].i == inumber) && (bonds[m].j == jnumber)) ||
           ((bonds[m].i == jnumber) && (bonds[m].j == inumber))){
          bonds[m].setBondIndex(Compound::BondIndex(r));
        }
      }
    }
    #ifdef MAIN_RESIDUE_DEBUG_SPECIFIC
    std::cout<<"bonds[m].bondIndexes assigned"<<std::endl<<std::flush;
    #endif


   // Fill indexMap
   // ORDER
   for(ix = 0; ix<natms; ix++){
      indexMap[ix][0] = ix;
      indexMap[ix][1] = bAtomList[ix].atomIndex;
    }

    // Inverse indexMap
    int Inverse[natms];
    std::cout << " natms " << natms << std::endl;
    assert(PrmToAx_po);
    
    for(unsigned int i=0; i<natms; i++){ // invert
      Inverse[ (int)(indexMap[i][1]) ] = i;
    }
    for(unsigned int i=0; i<natms; i++){ // copy
      PrmToAx_po[i] = TARGET_TYPE(Inverse[i]);
    }

    for(unsigned int i=0; i<natms; i++){ // invert
      Inverse[ (int)(indexMap[i][2]) ] = i;
    }
    for(unsigned int i=0; i<natms; i++){ // copy
      MMTkToPrm_po[i] = TARGET_TYPE(Inverse[i]);
    }

    #ifdef MAIN_RESIDUE_DEBUG_SPECIFIC
    std::cout<<"PrmToAx_po filled "<<std::endl<<std::flush;
    for(unsigned int k=0; k<natms; k++){
        std::cout<<"PrmToAx_po["<<k<<"]: "<<PrmToAx_po[k]<<std::endl;
    }
    std::cout<<"MMTkToPrm_po filled"<<std::endl<<std::flush;
    for(unsigned int k=0; k<natms; k++){
        std::cout<<"MMTkToPrm_po["<<k<<"]: "<<MMTkToPrm_po[k]<<std::endl;
    }
    std::cout<<"indexMap filled "<<std::endl<<std::flush;
    for(int i = 0; i<natms; i++){
      std::cout<<"indexMap["<<i<<"]: "<<indexMap[i][0]<<" "<<indexMap[i][1]<<" "<<indexMap[i][2]<<std::endl;
    }
    #endif

    // Create atomTargets from passed coords array REVISE CHANGE
    std::cout<<"Create atomTargets from passed coords array"<<std::endl;
    std::map<AtomIndex, Vec3> atomTargets;  
    ix = 0;
    if(first_time == true){ // Take coordinates from main (mol2 / previously MMTK)
      std::cout<<"Take coordinates from main (mol2)"<<std::endl;
      for (Compound::AtomIndex aIx(0); aIx < getNumAtoms(); ++aIx){
       //Vec3 v(coords[ ix ][0],
       //       coords[ ix ][1],
       //       coords[ ix ][2]);
        Vec3 v(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
        atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, v));
        ix++;
      }
    }
    else{ // Take coordinates from REVISE CHANGE
      std::cout<<"Take coordinates from MMTK"<<std::endl;
      int ixi, prmtopi_from_SimTK, prmtopi_from_MMTK;
      for(ix = 0; ix < getNumAtoms(); ++ix){
        ixi                 = indexMap[ix][0];
        prmtopi_from_SimTK  = indexMap[ix][1];
        prmtopi_from_MMTK   = indexMap[ix][2];
        Vec3 v(coords[prmtopi_from_MMTK][0]/10, coords[prmtopi_from_MMTK][1]/10, coords[prmtopi_from_MMTK][2]/10);

        atomTargets.insert(pair<AtomIndex, Vec3>
          (bAtomList[ix].atomIndex, v)
        );
      }
    }

    // Assign Compound coordinates by matching bAtomList coordinates
    matchDefaultTopLevelTransform(atomTargets);
    matchDefaultConfiguration(atomTargets, Match_Exact, true, 150.0); //Compound::Match_Idealized
    //matchDefaultConfiguration(atomTargets, Match_Idealized, true, 150.0); //Compound::Match_Idealized

    // Reinsert coordinates after matching
    PdbStructure  pdb(*this);
    for(unsigned int i=0; i<natms; i++){
      const PdbAtom& P = pdb.getAtom(String(bAtomList[i].name), PdbResidueId(1), String(" "));
      std::string s(P.getName());
      const Vec3& PC = P.getCoordinates();
      bAtomList[i].x = PC[0];
      bAtomList[i].y = PC[1];
      bAtomList[i].z = PC[2];
    }

    // Set rigidity
    if(ictdF=="IC"){
      for (int r=0 ; r<getNumBonds(); r++){
        setBondMobility(BondMobility::Free, Compound::BondIndex(r));
      }
    }
    else if(ictdF=="TD"){
      for (int r=0 ; r<getNumBonds(); r++){
        setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
      }
    }
    else if(ictdF=="RR"){ // Torsional dynamics with rigid rings
      for(unsigned int m=0; m<nbnds; m++){ // EU
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

  void Topology::setDuMMAtomParam(int, SimTK::Real vdw, SimTK::Real well){}
  void Topology::setDuMMBondParam(int, int, SimTK::Real k, SimTK::Real equil){}
  void Topology::setDuMMAngleParam(int, int, int, SimTK::Real k, SimTK::Real equil){}

  void Topology::setDuMMDihedralParam(int, int, int, int,
      int periodicity, SimTK::Real ampInKJ, SimTK::Real phaseInDegrees
  ){}
  void Topology::setDuMMDihedralParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2
  ){}
  void Topology::setDuMMDihedralParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2,
      int periodicity3, SimTK::Real ampInKJ3, SimTK::Real phaseInDegrees3
  ){}

  void Topology::setDuMMImproperParam(int, int, int, int,
      int periodicity, SimTK::Real ampInKJ, SimTK::Real phaseInDegrees
  ){}
  void Topology::setDuMMImproperParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2
  ){}
  void Topology::setDuMMImproperParam(int, int, int, int,
      int periodicity1, SimTK::Real ampInKJ1, SimTK::Real phaseInDegrees1,
      int periodicity2, SimTK::Real ampInKJ2, SimTK::Real phaseInDegrees2,
      int periodicity3, SimTK::Real ampInKJ3, SimTK::Real phaseInDegrees3
  ){}

  // Get

  int Topology::getNAtoms(void) const{}
  int Topology::getNBonds(void) const{}

  bSpecificAtom * Topology::getAtomByNumber(int number) const{}
  bSpecificAtom * Topology::getAtomByAtomIx(int aIx) const{}
  bSpecificAtom * Topology::getAtomByName(std::string name) const{}

  std::vector<bSpecificAtom> Topology::getNeighbours(int) const{}
  bBond * Topology::getBond(int, int) const{}
  int Topology::getBondOrder(int, int) const{}


  

