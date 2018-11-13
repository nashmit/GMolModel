loadamberparams frcmod.ff14SB
dial = loadmol2 ligand.mol2
savepdb dial ligand.pdb
saveamberparm dial ligand.prmtop ligand.inpcrd

or 

dial = = sequence { ACE ALA NME }
