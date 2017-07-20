import sys, os, subprocess
sys.path.insert(1, 'build-debug/lib/') #  GCHMC.so
import numpy as np

from MMTK import Database
from MMTK import Molecule
from MMTK.ParticleProperties import Configuration
from MMTK.ForceFields import Amber12SBForceField
from MMTK import Universe, InfiniteUniverse
from MMTK import Units

from VelocityVerlet import VelocityVerletIntegrator
import GCHMC

R = 8.3144621*Units.J/Units.mol/Units.K

mol_name = "ligand"
gaff_FN = 'gaff.dat'
parm_dir = './' 
mol_dir = '2but'
newmol = Molecule('2but')
mol_frcmod = os.path.join(mol_dir, 'ligand.frcmod')
universe = InfiniteUniverse(Amber12SBForceField(parameter_file=gaff_FN, mod_files=['frcmod.ff12SB', mol_frcmod]))
universe.addObject(newmol)

universe.configuration();
configuration =  universe.configuration().array

VVintegrator = VelocityVerletIntegrator(universe)
print "VV constructed"
GCintegrator = GCHMC.GCHMCIntegrator(universe, mol_dir, parm_dir + gaff_FN)
print "GCHMC constructed"

(confs, KEs_MM, Es_MM, acc, ntrials, dt) = VVintegrator(steps=200, steps_per_trial=50, T=300.0, delta_t=0.0015, \
     normalize=False, adapt=False, random_seed=random.randint(1,300))
print "Velocity Verlet finished"

(confs, Es_MM, acc, ntrials, dt) = GCintegrator.Call(30, 10, 300, 0.0015, random.randint(1,300), 0, 1, 0.5)
print "GC: ", Es_MM
GCintegrator.Clear()

