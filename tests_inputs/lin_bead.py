import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter



parser = argparse.ArgumentParser()
parser.add_argument('--dir', default=None,
  help='Directory with input data files')
parser.add_argument('--inFNRoots', default=None, nargs='+',
  help='Name root for first input data files')
parser.add_argument('--skiprows', default=300, type=int,
  help='# of lines to be skipped by numpy.loadtxt')
parser.add_argument('--nbins', default=50, type=int,
  help='Number of bins for histograms. Default for 7.2 spacing')
parser.add_argument('--xmin', default=-180, type=float,
  help='Minimum x value for histograms')
parser.add_argument('--xmax', default=180, type=float,
  help='Maximum x value for histograms')
parser.add_argument('--makehist', action='store_true', default=False,
  help='Make histogram')
parser.add_argument('--printhist', action='store_true', default=False,
  help='Make histogram')
parser.add_argument('--makeplot', action='store_true', default=False,
  help='Make histogram plot')
parser.add_argument('--minmax', action='store_true', default=False,
  help='Print min and max')
parser.add_argument('--moments', action='store_true', default=False,
  help='Print average and standard deviation')
args = parser.parse_args()

import genfuncs

# Additional data
flex_expect_rho = np.array([(1.0/360.0) for i in range(args.nbins)])

# General plot parameters
if args.makeplot:
  matplotlib.rc('font',**{\
    'family':'serif',\
    'serif':['Computer Modern Roman'], \
    'size':8}) # EU mod
  matplotlib.rc('text', usetex=True)
  fig_width = 3.36  # One-column figure
  fig_width2 = 6.68 # Two-column figure
  fig_height = 3.0
  fignum = 1
  plt.figure(num=fignum, figsize=(fig_width, fig_height), linewidth=3.0)
  #plt.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.95)
  plt.subplots_adjust(left=0.20, right=0.95, bottom=0.2, top=0.95)

subtitle = 3 * [None]
subtitle[0] = 'NO FIXMAN'
subtitle[1] = 'CDHMC'
subtitle[2] = 'MIXED'

for ri in range(len(args.inFNRoots)): # Iterate through roots
  FNlist = glob.glob(os.path.join(args.dir, args.inFNRoots[ri] + '*'))
  nrounds = len(FNlist)
  num_lines = sum(1 for line in open(FNlist[0]))
  #hists = []
  hists = np.zeros((nrounds, 2, args.nbins))
  relHists = np.zeros((nrounds, 2, args.nbins))
  print FNlist
  
  for li in range(nrounds):
    with open(FNlist[li], 'r') as in_FN1:
      alldata = np.loadtxt(in_FN1, skiprows=args.skiprows)
      print "alldata", alldata
  
    #col_size = int(np.ceil(alldata.shape[0]/2))
    col_size = alldata.shape[0]
    col_size_1 = col_size - 1
    col_size_2 = col_size - 2
  
    # Get col
    col = np.zeros((col_size))
    for ci in range(col_size):
      col[ci] = alldata[ci]
      #if alldata[ci] == 1: # Move accepted
      #  col[ci] = alldata[ci]
      #else:                   # Move not accepted
      #  col[ci] = alldata[ci]
    if args.minmax:
      print col.min(), col.max()
   
    if args.moments: 
      mean = np.mean(col)
      var = np.var(col)
      print mean, var, scipy.stats.skew(col) , scipy.stats.kurtosis(col)
   
    # Histogram
    if args.makehist:
      h = np.histogram(col, bins=args.nbins, range=(args.xmin, args.xmax), density=True)
      hists[li][0] = h[0]
      hists[li][1] = h[1][:-1]
      #hists.append(np.histogram(col, bins=args.nbins, range=(args.xmin, args.xmax), density=True))
        #print hists[li][1][i] + ((hists[li][1][i+1] - hists[li][1][i])/2), hists[li][0][i]
      #print "0", np.sum(hists[li][0])
      #print "X Cnt"

  print 'hists', hists
  dx  = (hists[0][1][1] - hists[0][1][0])
  dx2 = (hists[0][1][1] - hists[0][1][0])/2
  
  if args.makehist:
    # Get probabilties
    meanHist = np.zeros((args.nbins))
    stdHist  = np.zeros((args.nbins))
  
    # Get mean of probabilities  
    #if args.printhist:
    #  print "X <P(X)>"
    #xticks = np.zeros((args.nbins))
    #for li in range(nrounds):
    #  for i in range(args.nbins):
    #    meanHist[i] += hists[li][0][i]
    #for i in range(args.nbins):
    #  meanHist[i] /= nrounds
    #  #xticks[i] = hists[0][1][i] + ((hists[0][1][i+1] - hists[0][1][i])/2)
    #  xticks[i] = hists[0][1][i] + dx2
    #  if args.printhist:
    #    print xticks[i], meanHist[i]
    
    xticks = np.zeros((args.nbins))
    for i in range(args.nbins):
      print 'hists[:, [0], [%d]]' % (i) , np.ravel(hists[:, [0], [i]])
      meanHist[i] = np.mean(np.ravel(hists[:, [0], [i]]))
      stdHist[i]  =  np.std(np.ravel(hists[:, [0], [i]]))
      xticks[i] = hists[0][1][i] + dx2
    print 'xticks', xticks
    print 'meanHist', meanHist
    print 'stdHist', stdHist
 
    # Plot
    if args.makeplot:

      ax = plt.subplot(len(args.inFNRoots), 1, ri+1)
      meanL = ax.errorbar(xticks, meanHist, color='#000000', yerr=stdHist, linewidth=1.0, \
        capthick=0.5, capsize=1.0, elinewidth=0.5)
      flexL, = ax.plot(xticks, flex_expect_rho, 'r--', linewidth=1.0, color='#009900')
      #ax.legend([meanL, flexL], [args.inFNRoots[ri][0:5], 'Expected'], \
      #  loc=0, fontsize=6)

      plt.ylabel(r'$\mathrm{\rho(\phi)}$', fontsize=8)
      plt.xlim((-180, 180))
      plt.ylim((0.0020, 0.0035))

      if ri != 2:
        plt.xticks([])
      else:
        plt.xticks(np.arange(-180, 181, 90))
        ax.get_xaxis().set_ticklabels(np.arange(-180, 181, 90))
        plt.setp(ax.get_xticklabels(), rotation='horizontal', fontsize=8)
        plt.xlabel(r'$\mathrm{\phi(degrees)}$', fontsize=8)

      plt.yticks(np.arange(0.0020, 0.0040, 0.0005))
      yticklabels = ["%.4f"%(i) for i in np.arange(0.0020, 0.0040, 0.0005)]

      ax.get_yaxis().set_ticklabels(yticklabels)
      plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=8)

      textrelhpos = 0.5
      textrelvpos = 0.9
      ax.text(textrelhpos, textrelvpos, subtitle[ri], transform=ax.transAxes, \
        horizontalalignment='center', verticalalignment='top', \
        fontsize=8)

if args.makeplot:
  figFN = 'temp.pdf'
  plt.savefig(figFN, dpi=600, format='pdf')
  
  
  
  
  
  
  
  
  
  
  
  
