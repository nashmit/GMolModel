import sys, os, subprocess
#sys.path.insert(1, '/usr/lib64/python2.6/site-packages') # MMTK installation
#sys.path.insert(1, '/home/laurentiu/Downloads/MMTK-2.7.9/install_dir/lib64/python2.6/site-packages')
import copy
import numpy as np
import random

from MMTK import Database
from MMTK import Molecule
from MMTK.ParticleProperties import Configuration
from MMTK.ForceFields import Amber12SBForceField
from MMTK import Universe, InfiniteUniverse
from MMTK import Units

import GCHMC

mol_name = "ligand"

R = 8.3144621*Units.J/Units.mol/Units.K

gaff_FN = 'gaff.dat'
parm_dir = '/share/apps/amber/16/dat/leap/parm/' # Amber parameters directory
mol_dir = '2but'
newmol = Molecule('2but')
mol_frcmod = os.path.join(mol_dir, 'ligand.frcmod')
universe = InfiniteUniverse(Amber12SBForceField(parameter_file=gaff_FN, mod_files=['frcmod.ff12SB', mol_frcmod]))
universe.addObject(newmol)

universe.configuration();
configuration =  universe.configuration().array

GCintegrator = GCHMC.GCHMCIntegrator(universe, mol_dir, parm_dir + gaff_FN)
(confs, Es_MM, acc, ntrials, dt) = GCintegrator.Call(30, 10, 300, 0.0015, random.randint(1,300), 0, 1, 0.5)
print "GC: ", Es_MM
GCintegrator.Clear()






