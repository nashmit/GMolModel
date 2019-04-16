import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d

from matplotlib import cm
from extractLabels import *


parser = argparse.ArgumentParser()
parser.add_argument('--dir', default=None,
  help='Directory with input data files')
parser.add_argument('--inFNRoots', default=None, nargs='+',
  help='Name root for first input data files')
parser.add_argument('--skip_header', default=0, type=int,
  help='# of lines to be skipped by numpy.loadtxt from the beggining')
parser.add_argument('--skip_footer', default=0, type=int,
  help='# of lines to be skipped by numpy.loadtxt at the end')
parser.add_argument('--usecols', default=0, type=int, nargs='+',
  help='Columns to be read by numpy.loadtxt')
parser.add_argument('--stride', default=1, type=int,
  help='Stride for the read lines.')
parser.add_argument('--nbins', default=50, type=int,
  help='Number of bins for histograms. Default for 7.2 spacing')
parser.add_argument('--Ms', default=[100], type=int, nargs='+', 
  help='Number of steps to sum the correlation function on. One M per column')
parser.add_argument('--fft', action='store_true', default=False,
  help='Compute Fourier transform.')
parser.add_argument('--fitacf', action='store_true', default=False,
  help='Fit autocorrelation function to an exponential')
parser.add_argument('--acf', action='store_true', default=False,
  help='Compute autocorrelation function.')
parser.add_argument('--bse', action='store_true', default=False,
  help='Use block averaging')
parser.add_argument('--nofAddMethods', default=0, type=int,
  help='Number of additional methods wrote by James on StackOverflow.')
parser.add_argument('--makeplots', default='all', nargs='+', 
  help='Make plots')
parser.add_argument('--dataxlim', default=[0.0, 6.29], type=float, nargs='+',
  help='X limits for data plot.')
parser.add_argument('--dataylim', default=[0.0, 6.29], type=float, nargs='+',
  help='Y limits for data plot.')
parser.add_argument('--MD', default=[10, 10], type=int, nargs='+',
  help='Molecular dynamics stepsfor IC and TD respectively.')
parser.add_argument('--savefig', action='store_true', default=False,
  help='Save the plot into a file')
args = parser.parse_args()

import genfuncs

# General plot parameters
if args.makeplots:
  matplotlib.rc('font',**{\
    'family':'serif',\
    'serif':['Computer Modern Roman'], \
    'size':8}) # EU mod
  matplotlib.rc('text', usetex=True)
#  fig_width = 3.36  # One-column figure
#  fig_width2 = 6.68 # Two-column figure
#  fig_height = 3.0
#  fignum = 1
#  plt.figure(num=fignum, figsize=(fig_width, fig_height), linewidth=3.0)
#  plt.subplots_adjust(left=0.20, right=0.95, bottom=0.2, top=0.95)

ncols = len(args.usecols)
variances = np.zeros((ncols))
stds = np.zeros((ncols))
means = np.zeros((ncols))

tiny = 0.0000001
colors = ['black', 'red']

for ri in range(len(args.inFNRoots)): # Iterate through roots
  FNlist = glob.glob(os.path.join(args.dir, args.inFNRoots[ri] + '*'))
  nfiles = len(FNlist)
  hists = np.zeros((nfiles, 2, args.nbins))
  relHists = np.zeros((nfiles, 2, args.nbins))
  #print FNlist
  
  for li in range(nfiles):
    with open(FNlist[li], 'r') as in_FN1:
      alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer, usecols=args.usecols)
      alldata = alldata + (np.pi)
      #print 'alldata', alldata

      # Reshape data (transpose)
      if np.size(alldata.shape) == 1:
          cols = np.zeros((1, alldata.shape[0]))
          cols[0] = alldata
      else:
          shape_list = list(alldata.transpose().shape)
          #shape_list[0] = shape_list[0] + 1
          shape_tuple = tuple(shape_list)
          cols = np.zeros(shape_tuple)
          for coli in range(ncols):
              cols[coli] = alldata.transpose()[coli]
      #print 'cols', cols

      # 2D Plot
      if(('data' in args.makeplots) or ('all' in args.makeplots)):
          fig1= plt.figure(1)
          fig1.suptitle("Raw data")
          ax = plt.subplot(1, 1, 1)
          MDtot = args.MD[0] + args.MD[1]
          for i in range(0, cols.shape[-1], MDtot):
              j = i + args.MD[0]
              k = i + MDtot
              if(i == 0): # First set labels
                  ax.scatter(cols[1][i : j], cols[2][i : j], s = 1, c = 'blue', label = cols[0][i])
                  ax.scatter(cols[1][j : k], cols[2][j : k], s = 1, c = 'red', label = cols[0][j])
                  ax.plot(cols[1][i : j], cols[2][i : j], color = 'blue', linewidth = 1)
                  ax.plot(cols[1][j : k], cols[2][j : k], color = 'red', linewidth = 1)
              ax.scatter(cols[1][i : j], cols[2][i : j], s = 1, c = 'blue')
              ax.scatter(cols[1][j : k], cols[2][j : k], s = 1, c = 'red')
              ax.plot(cols[1][i : j], cols[2][i : j], color = 'blue', linewidth = 1)
              ax.plot(cols[1][j : k], cols[2][j : k], color = 'red', linewidth = 1)
          ax.scatter(cols[1][0], cols[2][0], s = 25, c = 'magenta')
          ax.legend()
  
          ax.set_xlabel(r'$\mathrm{\phi}$', fontsize=8)
          ax.set_ylabel(r'$\mathrm{\psi}$', fontsize=8)
          ax.set_xlim(args.dataxlim[0], args.dataxlim[1])
          ax.set_ylim(args.dataylim[0], args.dataylim[1])

          if args.savefig:
                figFN = 'temp.data.pdf'
                plt.savefig(figFN, dpi=600, format='pdf')
  
      # PLot a torus
      if(('torus' in args.makeplots) or ('all' in args.makeplots)):
          n = 100
          phi = np.linspace(0, 2.*np.pi, n)
          psi = np.linspace(0, 2.*np.pi, n)
          phi, psi = np.meshgrid(phi, psi)
          #print 'phi meshgrid', phi
          c = np.linspace(1, 1, n)
          a = np.linspace(1, 1, n)
          a[ 0 : n/2] = np.linspace(0, 1, n / 2)
          a[ n/2 : n] = np.linspace(1, 0, n / 2)
          c, a = np.meshgrid(c, a)
          #c, a = 1, 1
          x = (c + a*np.cos(phi)) * np.cos(psi)
          y = (c + a*np.cos(phi)) * np.sin(psi)
          z = a * np.sin(phi)

          aoffset = 0.0
          c = np.ones((cols.shape[-1]))
          a = np.ones((cols.shape[-1]))
          X = (c + (a + aoffset)*np.cos(cols[1])) * np.cos(cols[2])
          Y = (c + (a + aoffset)*np.cos(cols[1])) * np.sin(cols[2])
          Z = (a + aoffset) * np.sin(cols[1])

          fig2 = plt.figure(2)
          ax1 = fig2.add_subplot(111, projection='3d')
          #ax1.set_xlim((-1.0 * a) - 0.1, (1.0 * a) + 0.1)
          #ax1.set_ylim((-2.0 * a) - 0.1, (2.0 * a) + 0.1)
          #ax1.set_zlim((-1.0 * a) - 0.1, (1.0 * a) + 0.1)
          ax1.contour(x, y, z, 50, linewidth = 1, cmap='binary', edgecolor = 'grey')
          #ax1.plot_surface(x, y, z, rstride=5, cstride=5, color='k', edgecolors='w')
          #ax1.view_init(0, 0)

#          for i in range(0, cols.shape[-1], MDtot):
#              j = i + args.MD[0]
#              k = i + MDtot
#              if(i == 0): # First set labels
#                  ax1.scatter(X[i : j], Y[i : j], Z[i : j], s = 4, c = 'blue', label = cols[0][i])
#                  ax1.scatter(X[j : k], Y[j : k], Z[j : k], s = 4, c = 'red', label = cols[0][j])
#                  ax1.plot(X[i : j], Y[i : j], Z[i : j], color = 'blue', linewidth = 2)
#                  ax1.plot(X[j : k], Y[j : k], Z[j : k], color = 'red', linewidth = 2)
#              ax1.scatter(X[i : j], Y[i : j], Z[i : j], s = 4, c = 'blue')
#              ax1.scatter(X[j : k], Y[j : k], Z[j : k], s = 4, c = 'red')
#              ax1.plot(X[i : j], Y[i : j], Z[i : j], color = 'blue', linewidth = 2)
#              ax1.plot(X[j : k], Y[j : k], Z[j : k], color = 'red', linewidth = 2)
#          ax1.scatter(X[0], Y[0], Z[0], s = 25, c = 'magenta')
#          ax1.legend()

# Show plot
if args.makeplots[0] != '' :
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
