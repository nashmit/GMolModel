import sys, os, glob
import numpy as np
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
            print M[i][j],
        print
#

# What interval is the value in
def which_basin(basins, value):
    for i in range(basins.shape[0] - 1):
        if (value >= basins[i]) and (value < basins[i + 1]):
            return i
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
    #plt.subplots_adjust(left=0.20, right=0.95, bottom=0.2, top=0.95)


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

# Iterate through roots
nroots = len(args.inFNRoots)
for ri in range(nroots):
    FNlist = glob.glob(os.path.join(args.dir, args.inFNRoots[ri] + '*'))

    # Iterate through files
    nfiles = len(FNlist)
    angles = nfiles * [None]
    anglesTD = nfiles * [None]
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

    _angle_hists = nfiles * [None]
    angle_hists = nfiles * [None]

    freeEs = nfiles * [None]
    Psyms = nfiles * [None]

    # ndofs accs pe_o pe_set ke_o ke_n fix_o fix_set fix_n geometry
    printFreq = 100
    for li in range(nfiles):
        with open(FNlist[li], 'r') as in_FN1:
            print 'Reading', FNlist[li]
            if args.usecols == -1:
                alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer)
            else:
                alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer, usecols=args.usecols)

            # Get constrained dynamics # of DOFs
            minndofs = np.minimum(alldata[0][0], alldata[1][0])
            
            # Get all potential energies
            PEs[li] = np.array([a[3] for a in alldata])
            PEsIC[li] = np.array([a[3] for a in alldata if a[0] != minndofs])
            PEsTD[li] = np.array([a[3] for a in alldata if a[0] == minndofs])
            print "Mean Potential Energy", np.mean(PEs[li])
            print "Stdev Potential Energy", np.std(PEs[li])

            # Get Fixman potentials
            FixPs[li] = np.array([a[7] for a in alldata if a[0] == minndofs])
            accFixPs[li] = np.array([a[7] for a in alldata if ((a[0] == minndofs) and (a[3] != a[2]))])
            rejFixPs[li] = np.array([a[8] for a in alldata if ((a[0] == minndofs) and (a[3] == a[2]))])
           
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
            angles[li] = np.array([a[9] for a in alldata])
            anglesTD[li] = np.array([a[9] for a in alldata if a[0] == minndofs])
            # Change form -pi pi to 0 2pi
            for i in range(angles[li].shape[0]):
                if angles[li][i] < 0:
                    angles[li][i] += (2.0*np.pi)
            for i in range(anglesTD[li].shape[0]):
                if anglesTD[li][i] < 0:
                    anglesTD[li][i] += (2.0*np.pi)

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
            _angle_hists[li] = np.histogram(angles[li], bins=args.nbins, range=(0.0, 2 * np.pi), density=True)
            angle_hists[li] = _angle_hists[li][0]
            anglePrs = _angle_hists[li][0] * (_angle_hists[li][1][1] - _angle_hists[li][1][0])

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


            freeEs[li] = -1.0 * RT * np.log(anglePrs)

            # MFPT (assume equaly spaced intervals for angles
            # Stationary probabilities
            nintervals = 24
            stationaries = np.histogram(angles[li], bins=nintervals, range=(0.0, 2 * np.pi), density=True)
            basins = stationaries[1]

            # Transition matrix
            # Count
            jumps = []
            P = np.zeros((nintervals, nintervals), dtype=np.float64)
            previ = 0
            tauCount = np.zeros((nintervals), dtype=np.float64)
            tauStops = np.zeros((nintervals), dtype=np.float64)
            for ai in range(1, alldata.shape[0]):
                i = which_basin(basins, angles[li][ai - 1])
                j = which_basin(basins, angles[li][ai])
                P[i][j] += 1

                # Tau counts
                if (i == previ) and (j == previ):
                    tauCount[i] = tauCount[i] + 1
                else:
                    tauStops[i] = tauStops[i] + 1

                if i == j:
                    previ = i

                if i != j:
                    jumps.append(alldata[ai])
            jumps = np.array(jumps)

            print 'Basins', basins
            print 'Counts:'
            for i in range(P.shape[0]):
                for j in range(P.shape[1]):
                    print P[i][j], 
                print
            
            print 'Diagonal:'
            for i in range(P.shape[0]):
               print P[i][i], 
            print

            print 'tau: '
            for i in range(tauCount.shape[0]):
                print tauCount[i] / tauStops[i],
            print


            # Normalize transition matrix
            # TODO: Numpy sum is incorect ??
            #for i in range(P.shape[0]):
            #    for j in range(P.shape[1]):
            #        print P[i][j],
            #    print
            rowsums = np.zeros((nintervals), dtype=np.float64)
            for i in range(P.shape[0]):
                rowsums[i] = np.sum(P[:, i], dtype=np.float64)
            for i in range(P.shape[0]):
                for j in range(P.shape[1]):
                    P[i][j] = P[i][j] / rowsums[i] 

            # Symetrize transition matrix
            Psyms[li] = 0.5 * (P + P.transpose())
            #print "Symmetrized transition matrix:"
            #Print2DArray(Psyms[li])

    # Print acceptance
    print "Acceptance mean and std", np.mean(acc), np.std(acc)
    print "Acceptance in IC mean and std", np.mean(accIC), np.std(accIC)
    print "Acceptance in TD mean and std", np.mean(accTD), np.std(accTD)

    print "Fixman potentials mean of means and stds", np.mean(meanFixPs), np.mean(stdFixPs)
    print "Accepted Fixman potentials mean of means and stds", np.mean(meanAccFixPs), np.mean(stdAccFixPs)
    print "Rejected Fixman potentials in TD mean of means and stds", np.mean(meanRejFixPs), np.mean(stdRejFixPs)

    # Plot average angle distribution with error bars
    if args.makeplots:
        angle_hists = np.array(angle_hists)

        PE_hists = np.array(PE_hists)
        PETD_hists = np.array(PETD_hists)

        FixP_hists = np.array(FixP_hists)
        accFixP_hists = np.array(accFixP_hists)
        rejFixP_hists = np.array(rejFixP_hists)

        freeEs = np.array(freeEs)

        avg_angle_hist = np.zeros((angle_hists.shape[1]))
        std_angle_hist = np.zeros((angle_hists.shape[1]))
        avg_PE_hist = np.zeros((PE_hists.shape[1]))
        std_PE_hist = np.zeros((PE_hists.shape[1]))
        avg_PETD_hist = np.zeros((PETD_hists.shape[1]))
        std_PETD_hist = np.zeros((PETD_hists.shape[1]))

        avg_freeEs = np.zeros((angle_hists.shape[1]))
        std_freeEs = np.zeros((angle_hists.shape[1]))

        avg_FixP_hist = np.zeros((FixP_hists.shape[1]))
        std_FixP_hist = np.zeros((FixP_hists.shape[1]))
        avg_accFixP_hist = np.zeros((accFixP_hists.shape[1]))
        std_accFixP_hist = np.zeros((accFixP_hists.shape[1]))
        avg_rejFixP_hist = np.zeros((rejFixP_hists.shape[1]))
        std_rejFixP_hist = np.zeros((rejFixP_hists.shape[1]))

        Psym = np.zeros((angle_hists.shape[0], angle_hists.shape[1]))
        Psym = np.mean(Psyms, axis = 0)
        #print 'Average symmetrized transition matrix P:'
        #Print2DArray(Psym)

        for i in range(avg_PE_hist.shape[0]):
            avg_PE_hist[i] = np.mean(PE_hists[:, i])
            std_PE_hist[i] = np.std(PE_hists[:, i])
            avg_PETD_hist[i] = np.mean(PETD_hists[:, i])
            std_PETD_hist[i] = np.std(PETD_hists[:, i])

        for i in range(avg_angle_hist.shape[0]):
            avg_angle_hist[i] = np.mean(angle_hists[:, i])
            std_angle_hist[i] = np.std(angle_hists[:, i])

            avg_freeEs[i] = np.mean(freeEs[:, i])
            std_freeEs[i] = np.std(freeEs[:, i])
        #print "Angles histogram", avg_angle_hist

        for i in range(avg_FixP_hist.shape[0]):
            avg_FixP_hist[i] = np.mean(FixP_hists[:, i])
            std_FixP_hist[i] = np.std(FixP_hists[:, i])

        for i in range(avg_accFixP_hist.shape[0]):
            avg_accFixP_hist[i] = np.mean(accFixP_hists[:, i])
            std_accFixP_hist[i] = np.std(accFixP_hists[:, i])

        for i in range(avg_rejFixP_hist.shape[0]):
            avg_rejFixP_hist[i] = np.mean(rejFixP_hists[:, i])
            std_rejFixP_hist[i] = np.std(rejFixP_hists[:, i])

        dx = (_angle_hists[0][1][1] - _angle_hists[0][1][0])
        avg_freeE = -1.0 * RT * np.log(avg_angle_hist) * dx
        #print "Free energy", avg_freeE
        #for a in avg_freeE:
        #    print a,
        #print

    print args.makeplots, 'makeplots list'   
    # Detailed plots 
    for li in range(nfiles):
        if 'detailPE' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Potential energy TD", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            #ax.scatter(angles[li], PEs, label = currLabel, color=BasicColors[ri], s=0.5)
            ax.scatter(anglesTD[li], PEsTD[li], label = currLabel, color='black', s=0.5)
            #plt.plot(angles[li], PEs[0], linewidth=0.5)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{angle}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{U}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=18)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=18)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "detPE" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

        if 'detailangle' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Phase Space Transitions", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            plt.plot(angles[li], linewidth=0.5, color='black')

            ax.legend()
            ax.set_xlabel(r'$\mathrm{step}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\phi}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=18)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=18)
            plt.ylim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "detang" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

        if 'detailFixman' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Angles Fixman Potential", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            ax.scatter(anglesTD[li], FixPs, label = currLabel, color=BasicColors[ri], s=0.5)
            #plt.plot(angles[li], PEs[0], linewidth=0.5)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{angle}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{U_f}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=12)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=12)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "angFix" + ".pdf"
           
        # Other functions
        if 'detailFuncs' in args.makeplots:
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fix Funcs", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
   
            normPEsTD = (PEsTD[li] - np.min(PEsTD[li])) / (np.max(PEsTD[li]) - np.min(PEsTD[li])) 
            normFixPs = (FixPs[li] - np.min(FixPs[li])) / (np.max(FixPs[li]) - np.min(FixPs[li])) 
            #ax.scatter(anglesTD[li][::3], normPEsTD[::3], label = currLabel, color=BasicColors[ri], s=0.3)
            #ax.scatter(anglesTD[li][::3], normFixPs[::3], label = currLabel, color='red', s=0.3)
            ax.scatter(anglesTD[li][::3], normFixPs[::3] - normPEsTD[::3], label = currLabel, color='red', s=0.3)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{angle}$', fontsize=14)
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
            figs[-1].suptitle("Middle Angle Distribution", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            (_, caps, _) = plt.errorbar(_angle_hists[0][1][:-1], avg_angle_hist, yerr=std_angle_hist, color = BasicColors[ri], \
                label = currLabel, capsize=2)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{\rho}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=18)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=18)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "ang" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Free Energy Surface", fontsize=16) 
            ax = plt.subplot(1,1,1)
            currLabel = FPFTLabel2String(FPFTLabel(args.inFNRoots[ri]))
     
            (_, caps, _) = plt.errorbar(_angle_hists[0][1][:-1], avg_freeE, yerr=std_angle_hist, color = BasicColors[ri], \
                label = currLabel, capsize=2)

            ax.legend()
            ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=14)
            ax.set_ylabel(r'$\mathrm{F}$', fontsize=14)
            plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=32)
            plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=32)
            plt.xlim(0.0, 2*np.pi)

            if args.savefigs:
                figFN = "temp." + "FreeE" + ".pdf"
                plt.savefig(figFN, dpi=600, format='pdf')

    if ('Fixman' in args.makeplots) or ('all' in args.makeplots):
            figno = figno + 1
            figs.append(plt.figure(figno))
            figs[-1].suptitle("Fixman Potential Distribution", fontsize=16) 
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
            plt.errorbar(range(args.nbins), FixP_Prs, yerr=std_angle_hist, color = BasicColors[1], \
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


# Show
if args.makeplots:
  plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
