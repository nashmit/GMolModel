#include "bAddParams.hpp"

// Add parameters from prmtop PRMTOP
void bAddBiotypes(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
)
{
    // Add Biotypes
    for(int i=0; i<amberReader->getNumberAtoms(); i++){

       if (! SimTK::Biotype::exists(resName.c_str(), bAtomList[i].name, SimTK::Ordinality::Any) ){
            SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
                  SimTK::Element(
                      bAtomList[i].getAtomicNumber()
                      , (std::to_string(bAtomList[i].getElem())).c_str()
                      , (std::to_string(bAtomList[i].getElem())).c_str()
                      , bAtomList[i].getMass())
                , bAtomList[i].getNBonds()
                , resName.c_str()
                , (bAtomList[i].getName()).c_str()
                , SimTK::Ordinality::Any
            );

            bAtomList[i].setBiotypeIndex(biotypeIndex);
            std::cout << " bAddGaffParams: Defined Biotype: " 
                << resName << bAtomList[i].name << " " << SimTK::Ordinality::Any << "|" 
                << " with BiotypeIndex " << bAtomList[i].getBiotypeIndex() << std::endl;
        }else{
            std::cout << " bAddGaffParams: Biotype already set: " 
                << resName << bAtomList[i].name << " " << SimTK::Ordinality::Any << "|" 
                << " with BiotypeIndex " << bAtomList[i].getBiotypeIndex() << std::endl;
        }

    }
}

void bAddAtomClasses(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
)
{
    // Add atom classes
    SimTK::DuMM::AtomClassIndex aIx;
    for(int i=0; i<amberReader->getNumberAtoms(); i++){
        //if(int(bAtomList[i].getAtomClassIndex()) < 0){
        ////if( dumm.isValidAtomClass(bAtomList[i].getAtomClassIndex()) ){ // this is protected inherited ...
            aIx = dumm.getNextUnusedAtomClassIndex();
            bAtomList[i].setAtomClassIndex(aIx);
            dumm.defineAtomClass(
                (SimTK::DuMM::AtomClassIndex)aIx,
	        ( std::string("top") + resName + bAtomList[i].getFftype() 
                    + std::string("_") + std::to_string(bAtomList[i].getNumber()) ).c_str(),
                bAtomList[i].getAtomicNumber(), // int atomicNumber
                bAtomList[i].getNBonds(),
                bAtomList[i].getVdwRadius() / 10.0, // nm
                //bAtomList[i].getVdwRadius(), // A
                bAtomList[i].getLJWellDepth() * 4.184 // kcal to kJ
            );
            std::cout << "bAddGaffParams: Defined AtomClass: " << bAtomList[i].getAtomClassIndex() << std::endl;
        //}else{
        //    std::cout << "bAddGaffParams: AtomClass already set: " << bAtomList[i].getAtomClassIndex() << std::endl;
        //}


    }

    // Define DuMM charged atom types
    for(int k=0; k<amberReader->getNumberAtoms(); k++){
      std::string abuff =  resName;
      abuff += bAtomList[k].biotype;

      SimTK::DuMM::ChargedAtomTypeIndex tempChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
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
          //<< " Biotype::get(\"(this->name).c_str()\", bAtomList[k].biotype).getIndex() "
          //<< Biotype::get((this->name).c_str(), bAtomList[k].biotype).getIndex()
          << std::endl << std::flush;
      dumm.setBiotypeChargedAtomType(
        bAtomList[k].getChargedAtomTypeIndex(),
        //Biotype::get((this->name), bAtomList[k].biotype).getIndex()
        bAtomList[k].getBiotypeIndex()
      );

    }
    // #


}

void bAddBondParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
)
{
    // Suppose or try to have the same order as the reader
    for(int t=0; t<amberReader->getNumberBonds(); t++){
        dumm.defineBondStretch_KA(
            (bAtomList[bonds[t].i]).getAtomClassIndex(), 
            (bAtomList[bonds[t].j]).getAtomClassIndex(),
            amberReader->getBondsForceK(t),  //k1
            amberReader->getBondsEqval(t)   //equil1
        );

    }
}


void bAddAngleParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
)
{
    for(int t=0; t<amberReader->getNumberAngles(); t++){
      dumm.defineBondBend_KA(
        bAtomList[amberReader->getAnglesAtomsIndex1(t)].getAtomClassIndex(),
        bAtomList[amberReader->getAnglesAtomsIndex2(t)].getAtomClassIndex(),
        bAtomList[amberReader->getAnglesAtomsIndex3(t)].getAtomClassIndex(),
        amberReader->getAnglesForceK(t),
        ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getAnglesEqval(t))
      );
    }
}

void bAddTorsionParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
)
{
    std::vector<std::pair<int, int>> pairStartAndLens = amberReader->getPairStartAndLen(); 

    for(unsigned int index=0; index<pairStartAndLens.size(); index++){

        int first    = pairStartAndLens[index].first;
        int numberOf = pairStartAndLens[index].second;

        for(int t=first; t<(first+numberOf); t++){
            if(numberOf == 1){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t),   
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))
                );
            }
            else if(numberOf == 2){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t), 
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)),
                    amberReader->getDihedralsPeriod(t+1), amberReader->getDihedralsForceK(t+1), 
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1))
                );
            }
            else if(numberOf == 3){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)),
                    amberReader->getDihedralsPeriod(t+1), amberReader->getDihedralsForceK(t+1),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)),
                    amberReader->getDihedralsPeriod(t+2), amberReader->getDihedralsForceK(t+2),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2))
                );
            }
        }
    }
}


// Add all parameters
void bAddAllParams(
    std::string resName,
    readAmberInput *amberReader,
    SimTK::DuMMForceFieldSubsystem& dumm,
    bSpecificAtom *bAtomList,
    bBond *bonds
)
{
    bAddBiotypes(resName, amberReader, dumm, bAtomList, bonds);
    bAddAtomClasses(resName, amberReader, dumm, bAtomList, bonds);
    bAddBondParams(resName, amberReader, dumm, bAtomList, bonds);
    bAddAngleParams(resName, amberReader, dumm, bAtomList, bonds);
    bAddTorsionParams(resName, amberReader, dumm, bAtomList, bonds);
}





