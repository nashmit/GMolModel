# Imports
import sys, os, glob
import errno
import numpy as np
import argparse
import subprocess
import time
import math
from random import seed
from random import randint

def human_readable_size(size):
    for unit in ['','k','M','G','T']:
        if size < 1024.0:
            break
        size /= 1024.0
    return "%.0f%s" % (size, unit)

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inputDir', default='2but',
  help='Molecular input files directory.')
parser.add_argument('--outputDir', default='2but.heat',
  help='Directory to store results')
parser.add_argument('--mols', default=["2but"], nargs='+',
  help='Molecules.')
parser.add_argument('--restartDir', default='2but.heat/restart',
  help='Restart / inpcrd files.')
parser.add_argument('--restartFNRoot', default='ligand.inpcrd',
  help='Restart / inpcrd files.')
parser.add_argument('--regimen', default="IC",
  help='Regimen.')
parser.add_argument('--sampler', default="IC",
  help='Sampler: IC, TD or MIX.')
parser.add_argument('--fpotential', default="FP",
  help='Fixman potential indicator string.')
parser.add_argument('--ftorque', default="NFT",
  help='Fixman torque indicator string.')
parser.add_argument('--MDns', default=[5], type=int, nargs='+',
  help='Number of molecular dynamics steps.')
parser.add_argument('--taus', default=[0.001], type=float, nargs='+',
  help='Timesteps.')
parser.add_argument('--Ts', default=[300], type=float, nargs='+',
  help='Temperatures.')
parser.add_argument('--sas', default=[1], type=int, nargs='+',
  help='Samples per round.')
parser.add_argument('--nrounds', default=10, type=int,
  help='Number of rounds.')
parser.add_argument('--printFreq', default=1, type=int,
  help='Output print frequency.')
parser.add_argument('--pdbFreq', default=0, type=int,
  help='Write pdb frequency.')
parser.add_argument('--nsims', default=1, type=int,
  help='Number of repeats.')
parser.add_argument('--attachPrevSeeds', action='store_true', default=False,
  help='Attach restart files seeds if found.')
parser.add_argument('--dry', action='store_true', default=False,
  help='Dry run.')
args = parser.parse_args()

# Get restart files
FNlist = glob.glob(os.path.join(args.restartDir, args.restartFNRoot + '*'))
nfiles = len(FNlist)

# Get previous seeds
prevseeds = nfiles * ['']
if args.attachPrevSeeds:
    for li in range(nfiles):
        FNWords = FNlist[li].split('.')
        for i in range(len(FNWords) - 1, 0, -1):
            if FNWords[i].find(args.mols[0]) == 0:
                prevseeds[li] = FNWords[i][len(args.mols[0]):] + '.'
                break

inpFilenames = (args.nsims * nfiles) * [None]
outFilenames = (args.nsims * nfiles) * [None]

for li in range(nfiles):
    # Get a set of new seeds
    seeds = set()
    nofTries = args.nsims * 10
    for i in range(nofTries):
        seeds.add(randint(1, 100000))
        if(len(seeds) == (args.nsims * 2)):
            break
    # Put them in aNumpy vector
    vseeds = np.zeros((args.nsims), dtype=int)
    vseeds2 = np.zeros((args.nsims), dtype=int)
    
    vIx = -1
    for s in seeds:
        vIx = vIx + 1
        if vIx < vseeds.size:
            vseeds[vIx] = s
        else:
            vseeds2[vIx - vseeds.size] = s

    # Generate filenames
    for repeat in range(args.nsims):
        hseed = "%.3E" % vseeds[repeat]
        hnrounds = human_readable_size(args.nrounds)
        
        if args.regimen == "MIX":
            suffix = args.mols[0] + '.' \
                + args.regimen + '.' \
                + args.sampler + '.' \
                + args.fpotential + args.ftorque + '.' \
                + 'MD' + str(args.MDns[0]) + '_' + str(args.MDns[1]) + '.' \
                + 'ts' + "%.0f" % (args.taus[0] * 1000.0) + '_' + "%.0f" % (args.taus[1] * 1000.0) + '.' \
                + 'T' + str(args.Ts[0]) + '.' \
                + 'sa' + str(args.sas[0]) + '_' + str(args.sas[1]) + '.' \
                + hnrounds + '.' \
                + prevseeds[li] + str(vseeds[repeat]) + '.' + str(repeat)
        else:
            suffix = args.mols[0] + '.' \
                + args.regimen + '.' \
                + args.sampler + '.' \
                + args.fpotential + args.ftorque + '.' \
                + 'MD' + str(args.MDns[0]) + '.' \
                + 'ts' + "%.0f" % (args.taus[0] * 1000.0) + '.' \
                + 'T' + str(args.Ts[0]) + '.' \
                + 'sa' + str(args.sas[0]) + '.' \
                + hnrounds + '.' \
                + prevseeds[li] + str(vseeds[repeat]) + '.' + str(repeat)
    
        inpFilenames[repeat] = 'inp.' + suffix
        outFilenames[repeat] = 'out.' + suffix
    
    print "Input filenames:"
    for repeat in range(args.nsims):
        print inpFilenames[repeat]
    #print "Output filenames:"
    #for repeat in range(args.nsims):
    #    print outFilenames[repeat]
    
    # Input files content
    for repeat in range(args.nsims):
        print "Generate input file", inpFilenames[repeat], "...",
        inpFile = open(args.inputDir + "/" + inpFilenames[repeat], "w+")
        fpotential = ''
        ftorque = ''
        if args.fpotential == "FP":
            fpotential = 'TRUE'
        else:
            fpotential = 'FALSE'
    
        if args.ftorque == "FT":
            ftorque = 'TRUE'
        else:
            ftorque = 'FALSE'
    
        content = """# Molecule information
MOLECULES %s                  # Molecules to load
PRMTOP ligand.prmtop        # Parameter file
INPCRD %s        # Coordinate / Restart file
RBFILE ligand.rb            # Rigid bodies definition file
FLEXFILE ligand.flex        # Flexible bonds definition file
OUTPUT_DIR %s

# Simulation
WORLDS IC TD                       # Regimen (IC, TD, MIX, RB, RBMIX)
SAMPLER %s %s            # HMC or MC
ROUNDS %d                   # MC steps (Acception-rejection steps)
TIMESTEPS %.4f %.4f                # Timesteps to be used with regimens
MDSTEPS  %d %d   # Number of MD trial steps
SAMPLES_PER_ROUND %d %d  # Number of MC trials within a mixing round
REPRODUCIBLE TRUE
SEED %d %d

# Thermodynamics
THERMOSTAT None None        # Thermostat (HMC + Andersen = MD)
TEMPERATURE_INI  %.2f %.2f       # Temperature for constant temperature simulations
TEMPERATURE_FIN  %.2f %.2f     # Temperature for constant temperature simulations
FFSCALE AMBER AMBER        # Force field scale factors
GBSA 0.0 0.0              # GBSA scale factor

# Generalized coordinates related
FIXMAN_POTENTIAL FALSE %s # Use Fixman potential
FIXMAN_TORQUE FALSE %s          # Use Fixman torque

# Output
VISUAL FALSE FALSE                 # Use the visualizer
PRINT_FREQ %d
WRITEPDBS %d %d      # Write pdbs
GEOMETRY TRUE TRUE                # Calculate geometric features
DIHEDRAL 354 356 358 364 356 358 364 366
DISTANCE 0 718

# Software specs
THREADS 1 1
OPENMM TRUE
""" % (args.mols[0] \
        , '../' + FNlist[li] \
        , args.outputDir \
        , args.sampler, args.sampler \
        , args.nrounds \
        , args.taus[0] , args.taus[1] \
        , args.MDns[0], args.MDns[1] \
        , args.sas[0], args.sas[1]  \
        , vseeds[repeat], vseeds2[repeat]
        , args.Ts[0], args.Ts[0], args.Ts[1], args.Ts[1] \
        , fpotential, ftorque
        , args.printFreq \
        , args.pdbFreq, args.pdbFreq \
        )
        inpFile.write(content)
        inpFile.close()
        print "Done."
    
    if args.dry:
        for repeat in range(args.nsims):
            #print "Starting simulation ", inpFilenames[repeat], "..."
            print "./tests/Robosample " + args.inputDir + '/' + inpFilenames[repeat] + " > " + args.outputDir + '/' + outFilenames[repeat]
            #print "Done."
    else:
        for repeat in range(args.nsims):
            print "Starting simulation ", inpFilenames[repeat], "..."
            os.system( "nohup ./tests/Robosample " + args.inputDir + '/' + inpFilenames[repeat] + " > " + args.outputDir + '/' + outFilenames[repeat] + " &")
            print "Done."






