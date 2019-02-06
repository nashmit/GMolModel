# Imports
import sys, os, glob
import numpy as np
import argparse
import subprocess
import time
import errno

# Arguments
parser = argparse.ArgumentParser()

parser.add_argument('--pdbsDir', default=None,
  help='Directory with input data files')
parser.add_argument('--inpcrdsDir', default=None,
  help='Directory with input data files')
parser.add_argument('--inFNRoots', default=None, nargs='+',
  help='Name root for first input data files')
parser.add_argument('--skip_header', default=0, type=int,
  help='# of lines to be skipped by numpy.loadtxt from the beggining')
parser.add_argument('--skip_footer', default=0, type=int,
  help='# of lines to be skipped by numpy.loadtxt at the end')

args = parser.parse_args()

# Useful variables
pdbsDir = args.pdbsDir
inpcrdsDir = args.inpcrdsDir

# Main loop
try:
    os.makedirs(inpcrdsDir)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise
    else:
        print "\nDirectory %s already exists!" % inpcrdsDir

# Iterate through roots
nroots = len(args.inFNRoots)
natoms = 0

for ri in range(nroots):
    FNlist = glob.glob(os.path.join(pdbsDir, args.inFNRoots[ri] + '*'))
    nfiles = len(FNlist)
    for li in range(nfiles):

        FNlist[li] = FNlist[li].split("/")[-1]
        FNwords = FNlist[li].split(".")
        inpcrdFN = ''
        for wordIx in range(len(FNwords) - 1):
            inpcrdFN = inpcrdFN + FNwords[wordIx] + '.'
        inpcrdFN = inpcrdFN + 'inpcrd'
        print inpcrdFN

        file_length = len(open(pdbsDir + "/" + FNlist[li],'r').read().split('\n'))
        if file_length > 1: # Check
            with open(pdbsDir + "/" + FNlist[li], 'r') as in_FN1:
                print 'Reading', FNlist[li], ' ',
                lines = in_FN1.readlines()
                natoms = (len(lines))
                xyz = np.zeros((natoms, 3))
                for lineIx in range(len(lines)):
                    xyz[lineIx] = [float(lines[lineIx][30:38]), float(lines[lineIx][38:46]), float(lines[lineIx][46:54])]
                print " done."

            # Write inpcrd file
            print "Writing file ", inpcrdFN, "... ",
            if (os.path.isfile(inpcrdsDir + "/" + inpcrdFN) ):
                answer = raw_input(inpcrdsDir + "/" + inpcrdFN \
                    + "  already exists! Do you want to overwrite it? (y/n)")
            else:
                answer = 'y'
            if(answer == 'y'):
                inpcrdF = open(inpcrdsDir + "/" + inpcrdFN, "w+")
                inpcrdF.write( 'ACE\n')
                inpcrdF.write( '  ' + str(natoms) + '\n' )
                for atomIx in range(natoms):
                    line = ("%12.7f%12.7f%12.7f" % (xyz[atomIx][0], xyz[atomIx][1], xyz[atomIx][2]))
                    inpcrdF.write(line)
                    if( ((atomIx + 1)% 2) == 0 ):
                        inpcrdF.write( '\n' )
                inpcrdF.close()
                print " done."
            else:
                print 




