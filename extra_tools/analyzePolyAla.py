import sys, os, glob
import numpy as np
import copy
import argparse
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# Return a tuple indicating if Fixman potential or torque is present
def FPFTLabel(inpString):
    try:
        if inpString.find('NFPNFT') > 0:
            return (False, False)
        elif inpString.find('FPNFT') > 0:
            return (True, False)
        elif inpString.find('NFPFT') > 0:
            return (False, True)
        elif inpString.find('FPFT') > 0:
            return (True, True)
        else:
            raise Exception('FPFTLabel: string not found')
    except Exception as exception:
        print exception.args
#

#
def FPFTLabel2String(FPFTTuple):
    try:
        if FPFTTuple == (False, False):
            return 'NFPNFT'
        elif FPFTTuple == (True, False):
            return 'FPNFT'
        elif FPFTTuple == (False, True):
            return 'NFPFT'
        elif FPFTTuple == (True, True):
            return 'FPFT'
        else:
            raise TypeError
    except TypeError:
        print "FPFTLabel2String: need a tuple of two bools"
#

# Print array nicely
def PrintArray(Arr):
    for a in Arr:
        print a
#

# Print transition matrix nicely
def Print2DArray(M):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            print M[i][j],
        print
#

# Print transition matrix nicely
def PrintMatrix(M):
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            print M[i, j],
        print
#

# What interval is the value in
def which_basin1d(basins, value):
    for i in range(basins.shape[0] - 1):
        if (value >= basins[i]) and (value < basins[i + 1]):
            return i
#

# What interval is the value in
def which_basin2d(basins, phi, psi):
    for i in range(basins.shape[0]):
        if (phi >= basins[i][0, 0]) and (phi < basins[i][0, 1]):
            if (psi >= basins[i][1, 0]) and (psi < basins[i][1, 1]):
                return i
    print 'Values', phi, psi, 'out of basins.'
    
#

parser = argparse.ArgumentParser()
parser.add_argument('--dir', default=None,
  help='Directory with input data files')
parser.add_argument('--inFNRoots', default=None, nargs='+',
  help='Name root for first input data files')
parser.add_argument('--skip_header', default=0, type=int,
  help='# of lines to be skipped by numpy.loadtxt from the beggining')
parser.add_argument('--skip_footer', default=0, type=int,
  help='# of lines to be skipped by numpy.loadtxt at the end')
parser.add_argument('--usecols', default=-1, type=int, nargs='+',
  help='Columns to be read by numpy.loadtxt')
parser.add_argument('--T', default=0.0, type=float,
  help='Simulation temperature')
parser.add_argument('--nbins', default=10, type=int,
  help='Number of bins for histograms. Default for 7.2 spacing')
parser.add_argument('--makeplots', default=[], nargs='+',
  help='Which plots to make')
parser.add_argument('--detailedplots', action='store_true', default=False,
  help='Make additional plots')
parser.add_argument('--allplots', action='store_true', default=False,
  help='Make all plots')
parser.add_argument('--savefigs', action='store_true', default=False,
  help='Save the plot into a file')
args = parser.parse_args()

# General plot parameters
if args.makeplots:
    matplotlib.rc('font',**{\
      'family':'serif',\
      'serif':['Computer Modern Roman'], \
      'size':8}) # EU mod
    matplotlib.rc('text', usetex=True)
    fig_width = 3.36  # One-column figure
    fig_width2 = 6.68 # Two-column figure
    fig_height = 3.0
    fignum = 1
    #plt.figure(num=fignum, figsize=(fig_width, fig_height), linewidth=3.0)
    plt.subplots_adjust(left=0.20, right=0.95, bottom=0.2, top=0.95)


# Thermodynamics
R = 0.008314472 # kJ / (mol K)
RT = R * args.T

# Colors
BasicColors = ['black', 'red', 'green', 'blue', 'cyan', 'magenta', 'yellow']

# Values for plots
figno = 0
axno = 0
nofplots = 2
nofdetailedplots = 2
nofallplots = 3
figs = []

# Alanine dipeptide basin definitions
dialBasins = np.zeros((8,2,2))
tol = 0.001
mPI = -1.0 * np.pi - tol
PI = np.pi + tol

dialBasins[0] = np.array([[mPI, 0.175], [ 2.15,    PI]]) # C5
dialBasins[1] = np.array([[0.175,  PI], [ 2.15,    PI]]) # m1u

dialBasins[2] = np.array([[mPI, 0.175], [ 0.32,  2.15]]) # C7
dialBasins[3] = np.array([[0.175,  PI], [ 0.32,  2.15]]) # aL

dialBasins[4] = np.array([[mPI, 0.175], [-2.00,  0.32]]) # aR
dialBasins[5] = np.array([[0.175,  PI], [-2.00,  0.32]]) # C7

dialBasins[6] = np.array([[mPI, 0.175], [  mPI, -2.00]]) # C5d
dialBasins[7] = np.array([[0.175,  PI], [  mPI, -2.00]]) # m1d

# Iterate through roots
nroots = len(args.inFNRoots)
printFreq = 1
for ri in range(nroots):
    FNlist = glob.glob(os.path.join(args.dir, args.inFNRoots[ri] + '*'))

    # Iterate through files
    nfiles = len(FNlist)
    anglePhi = nfiles * [None]
    anglePhiTD = nfiles * [None]
    anglePsi = nfiles * [None]
    anglePsiTD = nfiles * [None]
    PEs = nfiles * [None]
    PEsIC = nfiles * [None]
    PEsTD = nfiles * [None]
    FixPs = nfiles * [None]
    accFixPs = nfiles * [None]
    rejFixPs = nfiles * [None]
    meanFixPs = nfiles * [None]
    meanAccFixPs = nfiles * [None]
    meanRejFixPs = nfiles * [None]
    stdFixPs = nfiles * [None]
    stdAccFixPs = nfiles * [None]
    stdRejFixPs = nfiles * [None]

    accs = nfiles * [None]
    accsIC = nfiles * [None]
    accsTD = nfiles * [None]
    acc = nfiles * [None]
    accIC = nfiles * [None]
    accTD = nfiles * [None]

    _PE_hists = nfiles * [None]
    PE_hists = nfiles * [None]
    _PETIC_hists = nfiles * [None]
    PEIC_hists = nfiles * [None]
    _PETD_hists = nfiles * [None]
    PETD_hists = nfiles * [None]
    _FixP_hists = nfiles * [None]
    FixP_hists = nfiles * [None]

    acc_FixP_hists = nfiles * [None]
    accFixP_hists = nfiles * [None]
    rej_FixP_hists = nfiles * [None]
    rejFixP_hists = nfiles * [None]

    _anglePhi_hists = nfiles * [None]
    anglePhi_hists = nfiles * [None]
    _anglePsi_hists = nfiles * [None]
    anglePsi_hists = nfiles * [None]
    _angles_hists = nfiles * [None]
    angles_hists = nfiles * [None]

    Psyms = nfiles * [None]

    for li in range(nfiles):
        with open(FNlist[li], 'r') as in_FN1:
            print 'Reading', FNlist[li]
            if args.usecols == -1:
                alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer)
            else:
                alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer, usecols=args.usecols)

            # Formats:
            # 0    1        2    3      4       5       6     7       8       9 10  11
            # dofs accsteps pe_o pe_set ke_last ke_prop fix_o fix_set fix_n   d phi psi
            # dofs accsteps pe_o pe_set ke_prop ke_n    fix_o fix_n   fix_set d phi psi 

            # Get constrained dynamics # of DOFs
            minndofs = np.minimum(alldata[0][0], alldata[1][0])
            
            # Get all potential energies
            PEs[li] = np.array([a[3] for a in alldata])
            PEsIC[li] = np.array([a[3] for a in alldata if a[0] != minndofs])
            PEsTD[li] = np.array([a[3] for a in alldata if a[0] == minndofs])
            print "Mean Potential Energy", np.mean(PEs[li])
            print "Stdev Potential Energy", np.std(PEs[li])

            # Get Fixman potentials
            FixPs[li] = np.array([a[8] for a in alldata if a[0] == minndofs])
            accFixPs[li] = np.array([a[7] for a in alldata if ((a[0] == minndofs) and (a[3] != a[2]))])
            rejFixPs[li] = np.array([a[7] for a in alldata if ((a[0] == minndofs) and (a[3] == a[2]))])
           
            meanFixPs[li] = np.mean(FixPs[li])
            meanAccFixPs[li] = np.mean(accFixPs[li])
            meanRejFixPs[li] = np.mean(rejFixPs[li])
            stdFixPs[li] = np.std(FixPs[li])
            stdAccFixPs[li] = np.std(accFixPs[li])
            stdRejFixPs[li] = np.std(rejFixPs[li])
            print "Simulation", li, "Fixman Potential mean and std", meanFixPs[li], stdFixPs[li]
            print "Simulation", li, "Fixman Potential mean and std", meanAccFixPs[li], stdAccFixPs[li]
            print "Simulation", li, "Fixman Potential mean and std", meanRejFixPs[li], stdRejFixPs[li]

            # Get angles
            anglePhi[li] = np.array([a[10] for a in alldata])
            anglePhiTD[li] = np.array([a[10] for a in alldata if a[0] == minndofs])
            anglePsi[li] = np.array([a[11] for a in alldata])
            anglePsiTD[li] = np.array([a[11] for a in alldata if a[0] == minndofs])

            # Change form -pi pi to 0 2pi
            #for i in range(anglePhi[li].shape[0]):
                #if anglePhi[li][i] < 0:
                #    anglePhi[li][i] += (2.0*np.pi)
                #if anglePsi[li][i] < 0:
                #    anglePsi[li][i] += (2.0*np.pi)
            for i in range(anglePhiTD[li].shape[0]):
                if anglePhiTD[li][i] < 0:
                    anglePhiTD[li][i] += (2.0*np.pi)
                if anglePsiTD[li][i] < 0:
                    anglePsiTD[li][i] += (2.0*np.pi)

            # Get acceptance rate
            accs[li] = np.array([a[1] for a in alldata])
            accsIC[li] = np.array([a[1] for a in alldata if a[0] != minndofs])
            accsTD[li] = np.array([a[1] for a in alldata if a[0] == minndofs])
            acc[li] = accs[li][-1] / (accs[li].size * printFreq)
            accIC[li] = accsIC[li][-1] / (accsIC[li].size * printFreq)
            accTD[li] = accsTD[li][-1] / (accsTD[li].size * printFreq)
            print "Acceptance", acc[li]
            print "Acceptance in IC", accIC[li]
            print "Acceptance in TD", accTD[li]
        
            # Free energy
            _anglePhi_hists[li] = np.histogram(anglePhi[li], bins=args.nbins, range=(-1.0*np.pi, np.pi), density=True)
            anglePhi_hists[li] = _anglePhi_hists[li][0]
            anglePhiPrs = _anglePhi_hists[li][0] * (_anglePhi_hists[li][1][1] - _anglePhi_hists[li][1][0])
            _anglePsi_hists[li] = np.histogram(anglePsi[li], bins=args.nbins, range=(-1.0*np.pi, np.pi), density=True)
            anglePsi_hists[li] = _anglePsi_hists[li][0]
            anglePsiPrs = _anglePsi_hists[li][0] * (_anglePsi_hists[li][1][1] - _anglePsi_hists[li][1][0])

            #_angles_hists[li] = np.histogram2d(anglePhi[li], anglePsi[li], bins=args.nbins, range=np.array([[-1.0*np.pi, np.pi], [-1.0*np.pi, np.pi]]) )
            _angles_hists[li] = np.histogram2d(anglePhi[li], anglePsi[li], bins = [args.nbins, args.nbins] )
            angles_hists[li] = _angles_hists[li][0]

            Sum = 0
            for a in angles_hists[li]:
                Sum = Sum + a

            angles_hists[li] = angles_hists[li] / Sum
            anglesPrs = angles_hists[li][0] * (_angles_hists[li][1][1] - _angles_hists[li][1][0])

            _PE_hists[li] = np.histogram(PEs[li], bins=args.nbins, density=True)
            PE_hists[li] = _PE_hists[li][0]
            PE_Prs = _PE_hists[li][0] * (_PE_hists[li][1][1] - _PE_hists[li][1][0])

            _PETD_hists[li] = np.histogram(PEsTD[li], bins=args.nbins, density=True)
            PETD_hists[li] = _PETD_hists[li][0]
            PETD_Prs = _PETD_hists[li][0] * (_PETD_hists[li][1][1] - _PETD_hists[li][1][0])

            _FixP_hists[li] = np.histogram(FixPs[li], bins=args.nbins, density=True)
            FixP_hists[li] = _FixP_hists[li][0]
            FixP_Prs = _FixP_hists[li][0] * (_FixP_hists[li][1][1] - _FixP_hists[li][1][0])

            acc_FixP_hists[li] = np.histogram(accFixPs[li], bins=args.nbins, density=True)
            accFixP_hists[li] = acc_FixP_hists[li][0]

            rej_FixP_hists[li] = np.histogram(rejFixPs[li], bins=args.nbins, density=True)
            rejFixP_hists[li] = rej_FixP_hists[li][0]

            # MFPT
            # Stationary probability vector
            _pi_ = np.zeros((dialBasins.shape[0]))
            P = np.zeros((dialBasins.shape[0], dialBasins.shape[0]), dtype=np.float64)
            W = np.zeros((dialBasins.shape[0], dialBasins.shape[0]), dtype=np.float64)
            Z = np.zeros((dialBasins.shape[0], dialBasins.shape[0]), dtype=np.float64)
            IdimP = np.zeros((dialBasins.shape[0], dialBasins.shape[0]), dtype=np.float64)
            for i in range(IdimP.shape[0]):
                IdimP[i, i] = 1.0
            MFPT = np.zeros((dialBasins.shape[0], dialBasins.shape[0]), dtype=np.float64)
            

            # Count
            for ai in range(1, anglePhi[li].shape[0]):
                i = which_basin2d(dialBasins, anglePhi[li][ai - 1], anglePsi[li][ai - 1])
                j = which_basin2d(dialBasins, anglePhi[li][ai],     anglePsi[li][ai])
                P[i][j] = P[i][j] + 1

                # Stationary probability vector counts accumulation
                _pi_[i] = _pi_[i] + 1
            _pi_[j] = _pi_[j] + 1 # last count

            #PrintArray(_pi_)
            #print 'Tranz counts'
            #Print2DArray(P)

            # Stationry counts to probabilities
            Sum = 0.0
            for i in range(_pi_.shape[0]):
                Sum = Sum + _pi_[i]
            _pi_ = _pi_ / Sum
            #print "Stationary vector"
            #print _pi_

            # Build W matrix
            for i in range(W.shape[0]):
                for j in range(W.shape[1]):
                    W[i][j] = _pi_[j]
            #print 'W matrix'
            #Print2DArray(W)

            # Transition probabilties
            rowsum = np.zeros((P.shape[0]))
            for i in range(P.shape[0]):
                for j in range(P.shape[1]):
                    rowsum[i] = rowsum[i] + P[i][j]

            for i in range(P.shape[0]):
                for j in range(P.shape[1]):
                    P[i][j] = P[i][j] / rowsum[i]

            #print 'Tranz matrix'
            #Print2DArray(P)
            #print P
            #print 'P determinant', np.linalg.det(P)

            # Swith to numpy matrices
            _pi_Mat = np.matrix(_pi_, copy = True)
            PMat = np.matrix(P, copy = True)
            WMat = np.matrix(W, copy = True)
            IdimPMat = np.matrix(IdimP, copy = True)
            ZMat = np.matrix(Z, copy = True)
            MFPTMat = np.matrix(MFPT, copy = True)

            print 'pi mat'
            PrintMatrix(_pi_Mat)
            #print 'pi mat', _pi_Mat
            #print 'piT * P', _pi_Mat.dot(PMat)

            print 'P mat'
            PrintMatrix(PMat)
            #print 'W mat'
            #PrintMatrix(WMat)
            #print 'I mat'
            #PrintMatrix(IdimPMat)

            ZMat = (IdimPMat - PMat + WMat).getI()
            #print 'Fundamental matrix Z mat'
            #PrintMatrix(ZMat)


            for i in range(MFPTMat.shape[0]):
                for j in range(MFPTMat.shape[1]):
                    MFPTMat[i, j] = (ZMat[j, j] - ZMat[i, j]) / _pi_Mat[0, j]

            print 'MFPT mat'
            PrintMatrix(MFPTMat)

    # Print acceptance
    print "Acceptance mean and std", np.mean(acc), np.std(acc)
    print "Acceptance in IC mean and std", np.mean(accIC), np.std(accIC)
    print "Acceptance in TD mean and std", np.mean(accTD), np.std(accTD)

    print "Fixman potentials mean of means and stds", np.mean(meanFixPs), np.mean(stdFixPs)
    print "Accepted Fixman potentials mean of means and stds", np.mean(meanAccFixPs), np.mean(stdAccFixPs)
    print "Rejected Fixman potentials in TD mean of means and stds", np.mean(meanRejFixPs), np.mean(stdRejFixPs)

    print "Fixman potentials stds of means and stds", np.std(meanFixPs), np.std(stdFixPs)
    print "Accepted Fixman potentials stds of means and stds", np.std(meanAccFixPs), np.std(stdAccFixPs)
    print "Rejected Fixman potentials in TD stds of means and stds", np.std(meanRejFixPs), np.std(stdRejFixPs)

    # Plot average angle distribution with error bars
    if args.makeplots:
        anglePhi_hists = np.array(anglePhi_hists)
        anglePsi_hists = np.array(anglePsi_hists)

        PE_hists = np.array(PE_hists)
        PETD_hists = np.array(PETD_hists)

        FixP_hists = np.array(FixP_hists)
        accFixP_hists = np.array(accFixP_hists)
        rejFixP_hists = np.array(rejFixP_hists)

        avg_anglePhi_hist = np.zeros((anglePhi_hists.shape[1]))
        std_anglePhi_hist = np.zeros((anglePhi_hists.shape[1]))
        avg_anglePsi_hist = np.zeros((anglePsi_hists.shape[1]))
        std_anglePsi_hist = np.zeros((anglePsi_hists.shape[1]))
        avg_PE_hist = np.zeros((PE_hists.shape[1]))
        std_PE_hist = np.zeros((PE_hists.shape[1]))
        avg_PETD_hist = np.zeros((PETD_hists.shape[1]))
        std_PETD_hist = np.zeros((PETD_hists.shape[1]))

        avg_FixP_hist = np.zeros((FixP_hists.shape[1]))
        std_FixP_hist = np.zeros((FixP_hists.shape[1]))
        avg_accFixP_hist = np.zeros((accFixP_hists.shape[1]))
        std_accFixP_hist = np.zeros((accFixP_hists.shape[1]))
        avg_rejFixP_hist = np.zeros((rejFixP_hists.shape[1]))
        std_rejFixP_hist = np.zeros((rejFixP_hists.shape[1]))


        for i in range(avg_PE_hist.shape[0]):
            avg_PE_hist[i] = np.mean(PE_hists[:, i])
            std_PE_hist[i] = np.std(PE_hists[:, i])
            avg_PETD_hist[i] = np.mean(PETD_hists[:, i])
            std_PETD_hist[i] = np.std(PETD_hists[:, i])

        for i in range(avg_anglePhi_hist.shape[0]):
            avg_anglePhi_hist[i] = np.mean(anglePhi_hists[:, i])
            std_anglePhi_hist[i] = np.std(anglePhi_hists[:, i])
            avg_anglePsi_hist[i] = np.mean(anglePsi_hists[:, i])
            std_anglePsi_hist[i] = np.std(anglePsi_hists[:, i])

        for i in range(avg_FixP_hist.shape[0]):
            avg_FixP_hist[i] = np.mean(FixP_hists[:, i])
            std_FixP_hist[i] = np.std(FixP_hists[:, i])

        for i in range(avg_accFixP_hist.shape[0]):
            avg_accFixP_hist[i] = np.mean(accFixP_hists[:, i])
            std_accFixP_hist[i] = np.std(accFixP_hists[:, i])

        for i in range(avg_rejFixP_hist.shape[0]):
            avg_rejFixP_hist[i] = np.mean(rejFixP_hists[:, i])
            std_rejFixP_hist[i] = np.std(rejFixP_hists[:, i])

        dx = (_anglePhi_hists[0][1][1] - _anglePhi_hists[0][1][0])

    # Detailed plots 
    for li in range(nfiles):
        if 'detailanglesmap' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            #figs[-1].suptitle("Angles", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
    
            #X = _angles_hists[li][1] + ((_angles_hists[li][1][1] - _angles_hists[li][1][0]) / 2.0)
            X = np.linspace(_angles_hists[li][1][0], _angles_hists[li][1][-1], args.nbins)
            Y = copy.deepcopy(X)
            Z = -1.0 * RT * np.transpose((angles_hists[li]) / (np.max(angles_hists[li])))

            #cset1 = ax.contourf(X , Y, np.transpose(angles_hists[li]), 40, cmap = plt.cm.RdYlBu)
            cset1 = ax.contourf(X , Y, Z, 40, cmap = plt.cm.bwr)
            plt.colorbar(cset1)
            #ax.contour(X , Y, np.transpose(angles_hists[li]), 40)
            #ax.scatter(anglePhi[li], anglePsi[li], label = currLabel, color='black', s=0.3)

            #ax.legend()
            ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=20, fontweight  = 'bold')
            ax.set_ylabel(r'$\mathrm{\psi}$', fontsize=20, fontweight  = 'bold')
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=20, fontweight  = 'bold')
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=20, fontweight  = 'bold')
            #plt.xlim(0.0, 2*np.pi)
            #plt.ylim(0.0, 2*np.pi)
            plt.xlim(-1.0*np.pi, np.pi)
            plt.ylim(-1.0*np.pi, np.pi)

            if args.savefigs:
                figFN = "temp." + "angmap" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

        if 'detailangles' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            #figs[-1].suptitle("Angles", fontsize=20) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
    
            ax.scatter(anglePhi[li], anglePsi[li], label = currLabel, color='black', s=0.3)

            #ax.legend()
            ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=20, fontweight  = 'bold')
            ax.set_ylabel(r'$\mathrm{\psi}$', fontsize=20, fontweight  = 'bold')
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=20, fontweight  = 'bold')
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=20, fontweight  = 'bold')
            #plt.xlim(0.0, 2*np.pi)
            #plt.ylim(0.0, 2*np.pi)
            plt.xlim(-1.0*np.pi, np.pi)
            plt.ylim(-1.0*np.pi, np.pi)

            if args.savefigs:
                figFN = "temp." + "angscat" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

        if 'detailtranz' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Angles", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
    
            X = _angles_hists[li][1] + ((_angles_hists[li][1][1] - _angles_hists[li][1][0]) / 2.0)
            #X = (-1.0) * X[:-1]
            X = X[:-1]
            Y = copy.deepcopy(X)
            ax.scatter(anglePhi[li], anglePsi[li], label = currLabel, color='black', s=0.5)
            ax.plot(anglePhi[li], anglePsi[li], label = currLabel, color='black', linewidth = 0.5)
            #ax.contour(X, Y, np.transpose(angles_hists[li]), 40)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\psi}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)
            #plt.xlim(0.0, 2*np.pi)
            #plt.ylim(0.0, 2*np.pi)
            plt.xlim(-1.0*np.pi, np.pi)
            plt.ylim(-1.0*np.pi, np.pi)

            if args.savefigs:
                figFN = "temp." + "angtranz" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')
        if 'detailPE' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Potential energy TD", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            #ax.scatter(anglePhi[li], PEs, label = currLabel, color=BasicColors[ri], s=0.5)
            ax.scatter(anglePhiTD[li], PEsTD[li], label = currLabel, color='red', s=0.5)
            #plt.plot(anglePhi[li], PEs[0], linewidth=0.5)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{angle}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{U}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "detPE" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

        if 'detailFixman' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Angles Fixman Potential", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            ax.scatter(anglePhiTD[li], FixPs, label = currLabel, color=BasicColors[ri], s=0.5)
            #plt.plot(anglePhi[li], PEs[0], linewidth=0.5)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{U_f}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "angPhiFix" + ".pdf"
           
        # Other functions
        if 'detailFuncs' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fix Funcs", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
   
            normPEsTD = (PEsTD[li] - np.min(PEsTD[li])) / (np.max(PEsTD[li]) - np.min(PEsTD[li])) 
            normFixPs = (FixPs[li] - np.min(FixPs[li])) / (np.max(FixPs[li]) - np.min(FixPs[li])) 
            #ax.scatter(anglePhiTD[li][::3], normPEsTD[::3], label = currLabel, color=BasicColors[ri], s=0.3)
            #ax.scatter(anglePhiTD[li][::3], normFixPs[::3], label = currLabel, color='red', s=0.3)
            ax.scatter(anglePhiTD[li][::3], normFixPs[::3] - normPEsTD[::3], label = currLabel, color='red', s=0.3)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{anglePhi}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{f()}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "someFunc" + ".pdf"
    # END detailed plots

    if ('Fixman' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fixman Potential Distribution", fontsize=16) 
            ax = plt.subplot(1,1,1)         
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            (_, caps, _) = plt.errorbar(_FixP_hists[0][1][:-1], avg_FixP_hist, yerr=std_FixP_hist, color = BasicColors[ri], \
                label = currLabel, capsize=2)


    if ('angles' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Angles Distributions", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            (_, caps, _) = plt.errorbar(_anglePhi_hists[0][1][:-1], avg_anglePhi_hist, yerr=std_anglePhi_hist, color = BasicColors[ri], \
                label = currLabel, capsize=2)
            (_, caps, _) = plt.errorbar(_anglePsi_hists[0][1][:-1], avg_anglePsi_hist, yerr=std_anglePsi_hist, color = BasicColors[ri], \
                label = currLabel, capsize=2)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\rho}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "ang" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

    if ('FixmanAll' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fixman Potential Distributions", fontsize=16) 
            ax = plt.subplot(1,1,1)         
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            #(_, caps, _) = plt.errorbar(_FixP_hists[0][1][:-1], avg_FixP_hist, yerr=std_FixP_hist, color = BasicColors[ri], \
            #    label = currLabel, capsize=2)

            for hist in _FixP_hists:
                plt.plot(_FixP_hists[0][1][:-1], hist[0])

            ax.legend()
            ax.set_xlabel(r'$\mathrm{U_f}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\rho}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)

            if args.savefigs:
                figFN = "temp." + "Fix" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

    if ('FixmanDiff' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            print 'FixmanDiff figure number', figno
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fixman Potential Distributions", fontsize=16) 
            ax = plt.subplot(1,1,1)         
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            (_, caps, _) = plt.errorbar(_FixP_hists[0][1][:-1], avg_FixP_hist, yerr=std_FixP_hist, color = BasicColors[0], \
                label = currLabel, capsize=2)
            plt.errorbar(acc_FixP_hists[0][1][:-1], avg_accFixP_hist, yerr=std_FixP_hist, color = BasicColors[1], \
                label = currLabel + 'acc', capsize=2)
            plt.errorbar(rej_FixP_hists[0][1][:-1], avg_rejFixP_hist, yerr=std_FixP_hist, color = BasicColors[2], \
                label = currLabel + 'rej', capsize=2)
            #plt.plot(_FixP_hists[0][1][:-1], std_FixP_hist)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{U_f}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\rho}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)

            if args.savefigs:
                figFN = "temp." + "FixDiff" + str(ri) + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

    if ('FixmanStds' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fixman Potential Stds", fontsize=16) 
            ax = plt.subplot(1,1,1)         
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            plt.plot(_FixP_hists[0][1][:-1], std_FixP_hist)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{stdU_f}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\rho}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)

            if args.savefigs:
                figFN = "temp." + "FixStds" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

    if ('PFixP' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Potential vs Fixman Potential", fontsize=16) 
            ax = plt.subplot(1,1,1)         
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))

            plt.errorbar(range(args.nbins), PE_Prs, yerr=std_PE_hist, color = BasicColors[0], \
                label = currLabel, capsize=2)
            plt.errorbar(range(args.nbins), FixP_Prs, yerr=std_FixP_hist, color = BasicColors[1], \
                label = currLabel, capsize=2)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{U}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\rho}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)

            if args.savefigs:
                figFN = "temp." + "PFixP" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

    if ('Fixtimeseries' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fixman Potential Timeseries", fontsize=16) 
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))

            for li in range(nfiles):
                ax = plt.subplot(nfiles, 1, li + 1)         
                #ax.plot(PEsTD[li], label = currLabel, color = BasicColors[ri])
                ax.plot(FixPs[li], label = currLabel, color = 'red')

            plt.xticks(np.arange(0, FixPs[li].size, 100))
            ax.legend()
            ax.set_xlabel(r'$\mathrm{t}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{U_f}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)

            if args.savefigs:
                figFN = "temp." + "Fixtime" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

    if ('PEtimeseries' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fixman Potential Timeseries", fontsize=16) 
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))

            for li in range(nfiles):
                ax = plt.subplot(nfiles, 1, li + 1)         
                #ax.plot(PEsTD[li], label = currLabel, color = BasicColors[ri])
                ax.plot(PEs[li], label = currLabel, color = 'black')

            plt.xticks(np.arange(0, FixPs[li].size, 100))
            ax.legend()
            ax.set_xlabel(r'$\mathrm{t}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{U_f}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)

            if args.savefigs:
                figFN = "temp." + "PEtime" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')


# Show
if args.makeplots:
  #if args.savefigs:
  #  figFN = "temp." + str(args.makeplots[figno]) + ".pdf"
  #  plt.savefig(figFN, dpi=600, format='pdf')
  plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
