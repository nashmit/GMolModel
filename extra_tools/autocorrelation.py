import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import numpy
import matplotlib.pyplot as plt

# Functions from https://stackoverflow.com/users/2005415/jason
def autocorr1(x,lags):
    '''numpy.corrcoef, partial'''

    corr=[1. if l==0 else numpy.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
    return numpy.array(corr)

def autocorr2(x,lags):
    '''manualy compute, non partial'''

    mean=numpy.mean(x)
    var=numpy.var(x)
    xp=x-mean
    corr=[1. if l==0 else numpy.sum(xp[l:]*xp[:-l])/len(x)/var for l in lags]

    return numpy.array(corr)

def autocorr3(x,lags):
    '''fft, pad 0s, non partial'''

    n=len(x)
    # pad 0s to 2n-1
    ext_size=2*n-1
    # nearest power of 2
    fsize=2**numpy.ceil(numpy.log2(ext_size)).astype('int')

    xp=x-numpy.mean(x)
    var=numpy.var(x)

    # do fft and ifft
    cf=numpy.fft.fft(xp,fsize)
    sf=cf.conjugate()*cf
    corr=numpy.fft.ifft(sf).real
    corr=corr/var/n

    return corr[:len(lags)]

def autocorr4(x,lags):
    '''fft, don't pad 0s, non partial'''
    mean=x.mean()
    var=numpy.var(x)
    xp=x-mean

    cf=numpy.fft.fft(xp)
    sf=cf.conjugate()*cf
    corr=numpy.fft.ifft(sf).real/var/len(x)

    return corr[:len(lags)]

def autocorr5(x,lags):
    '''numpy.correlate, non partial'''
    mean=x.mean()
    var=numpy.var(x)
    xp=x-mean
    corr=numpy.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)

    return corr[:len(lags)]
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
parser.add_argument('--usecols', default=0, type=int, nargs='+',
  help='Columns to be read by numpy.loadtxt')
parser.add_argument('--stride', default=1, type=int,
  help='Stride for the read lines.')
parser.add_argument('--nbins', default=50, type=int,
  help='Number of bins for histograms. Default for 7.2 spacing')
parser.add_argument('--Ms', default=[100], type=int, nargs='+', 
  help='Number of steps to sum the correlation function on. One M per column')
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
  fig_width = 3.36  # One-column figure
  fig_width2 = 6.68 # Two-column figure
  fig_height = 3.0
  fignum = 1
  plt.figure(num=fignum, figsize=(fig_width, fig_height), linewidth=3.0)
  plt.subplots_adjust(left=0.20, right=0.95, bottom=0.2, top=0.95)

ncols = len(args.usecols)
variances = np.zeros((ncols))
stds = np.zeros((ncols))
means = np.zeros((ncols))

tiny = 0.0000001

for ri in range(len(args.inFNRoots)): # Iterate through roots
  FNlist = glob.glob(os.path.join(args.dir, args.inFNRoots[ri] + '*'))
  nfiles = len(FNlist)
  hists = np.zeros((nfiles, 2, args.nbins))
  relHists = np.zeros((nfiles, 2, args.nbins))
  #print FNlist
  
  for li in range(nfiles):
    with open(FNlist[li], 'r') as in_FN1:
      alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer, usecols=args.usecols)
      print alldata
      alldata = alldata[::args.stride]

      print "alldata moments", np.mean(alldata), np.var(alldata), scipy.stats.skew(alldata) , scipy.stats.kurtosis(alldata)

      # Reshape data (transpose)
      if np.size(alldata.shape) == 1:
          cols = np.zeros((1, alldata.shape[0]))
          cols[0] = alldata
      else:
          shape_list = list(alldata.transpose().shape)
          shape_list[0] = shape_list[0] + 1
          shape_tuple = tuple(shape_list)
          cols = np.zeros(shape_tuple)
          for coli in range(ncols):
              cols[coli] = alldata.transpose()[coli]

      for coli in range(ncols):
          print 'col', coli, 'shape', cols[coli].shape

      # Do statistics on each columns
      for coli in range(ncols):
          print "Column", coli
          # Moments 
          means[coli] = np.mean(cols[coli])
          variances[coli] = np.var(cols[coli])
          stds[coli] = np.std(cols[coli])
          #print 'mean', means[coli]
          #print 'variance', variances[coli]

          if args.acf:
              # Autocorrelation time estimate (2009 Grossfield)
              # Compute correlation function vor M values or until it reaches 0
              n = cols[coli].size # Sample size
              M = args.Ms[coli]
              if M == 0: M = n
    
              corr = np.zeros((7, M))
              lags = np.arange(0, M)
            
              cut = 0
    
              CestGrossfield = np.zeros((M))
              for lag in lags:
                  # Take the product (f(t))(f(t+lag)) and average over (n - lag) points
                  # Normalize. The variance is also the correlation at time 0 and
                  # should be highest
                  CestGrossfield[lag] = np.mean( \
                      (cols[coli][0 : (n - lag)] - means[coli]) * \
                      (cols[coli][lag : n]       - means[coli]) ) \
                      / variances[coli]
                  #print CestGrossfield[lag]
    
                  # If the correlation function reaches 0 the rest is considered noise
                  if CestGrossfield[lag] < tiny:
                      cut = lag 
                      break
    
              corr[0] = CestGrossfield
              
              # Integrated autocorrelation time
              Iac = np.sum(corr[0])
              print "Integrated autocorrelation time", Iac 
              ESS = n / Iac
              print "Independent samples", ESS
    
              # Fit to an exponential Aexp(Bx)
              if args.fitacf:
                  y = CestGrossfield 
                  x = np.arange(0,n)
                  if cut > 1:
                      AB = np.polyfit( x[1:cut], np.log(y[1:cut]), 1, w=(y[1:cut]**2) )
                      corr[1][1:] = np.exp(AB[1]) * np.exp(AB[0] * x[1:])
                  else:
                      corr[1][1:] = 0.0
                  corr[1][0] = 1.0
                  print 'y = ', np.exp(AB[1]), ' * exp(', AB[0], ' * x)'
    
                  # Fit Integrated autocorrelation time
                  Iac = np.sum(corr[1])
                  print "Fit Integrated autocorrelation time", Iac 
                  ESS = n / Iac
                  print "Fit Independent samples", ESS
    
    
              # Additional methods
              add_funcs = [autocorr1, autocorr2, autocorr3, autocorr4, autocorr5]
              for i in range(2, args.nofAddMethods):
                  corr[i] = add_funcs[i](cols[coli], lags)
    
    
              # Standard error based on autocorrelation time fitting
              SE = stds[coli] * np.sqrt(Iac / float(n))
              print "Standard error", SE
    
              # Plots
              if(('acf' in args.makeplots) or ('all' in args.makeplots)):
                  fig1 = plt.figure(1)
                  fig1.suptitle("Autocorrelation time function") 
                  labels = ['Grossfield fit', 'Grossfield', 'np.corrcoef, partial', 'manual, non-partial',
                  'fft, pad 0s, non-partial', 'fft, no padding, non-partial',
                  'np.correlate, non-partial']
                  colors = ['red', 'orange', 'blue', 'magenta', 'cyan', 'black', 'green']
                  ax = plt.subplot(ncols, 1, coli+1)
                  line, = ax.plot(corr[0], label='Col ' +  str(coli) + ' ' + ' real', color='black')
                  line.set_dashes([2,2,10,2])
                  ax.plot(corr[1], label='Col ' +  str(coli) + ' ' + ' fitted', color='red')
                  ax.legend()
      
                  ax.set_xlabel(r'$\mathrm{\tau}$', fontsize=8)
                  ax.set_ylabel(r'$\mathrm{C(\tau)}$', fontsize=8)
            
                  plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=8)
                  plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=8)
            
                  textrelhpos = 0.5
                  textrelvpos = 0.9
                  #ax.text(textrelhpos, textrelvpos, FNlist[li], transform=ax.transAxes, \
                  #    horizontalalignment='center', verticalalignment='top', \
                  #    fontsize=8)
    
                  if args.savefig:
                        figFN = 'temp.acf.pdf'
                        plt.savefig(figFN, dpi=600, format='pdf')
    
          if args.bse:
              # Block averaging (Grossfield)
              stdsOfN = np.zeros((n/2 - 1))
              BSEOfN = np.zeros((n/2 - 1))
              for blockSize in range(1, n/2):
                  nofBlocks = int(np.floor(float(n) / float(blockSize)))
                  meansOfN = np.zeros((nofBlocks))
                  #print 'n, blockSize, nofBlocks', n, blockSize, nofBlocks
                  for i in range(0, nofBlocks):
                      pass
                      meansOfN[i] = np.mean(cols[coli][ (blockSize*i) : (blockSize*(i+1)) ])
                      #print (blockSize*i), (blockSize*(i+1)), ' = ',
                      #print meansOfN[i], '|',
                  #print 
                  stdsOfN[blockSize - 1] = np.std( meansOfN )
                  BSEOfN[blockSize - 1] = stdsOfN[blockSize - 1] / np.sqrt(nofBlocks)
              #print 'BSEOfN', BSEOfN
    
              # Fit to a log curve
              BSEOfNFit = np.zeros((n/2 - 1))
              x = np.arange(0, n/2-1)
              y = BSEOfN
              degree = 2
              AB = np.polyfit( x, y, degree, w=(np.exp(-1.*y)))
              #AB = np.polyfit( x, np.exp(y), 1, w=(1./y) )
              for i in range(BSEOfNFit.size):
                  for j in range(degree, -1, -1):
                      BSEOfNFit[i] += AB[degree - j] * np.power(x[i], j)
              ABSE = BSEOfNFit[n/2 - 2]
              print 'Asymptotic standard error', ABSE
    
              # Asymptotic autocorrelation time
              AESS = (stds[coli] / ABSE)**2
              print 'Asymptotic independent samples', AESS
    
              AIac = n / AESS
              print 'Asymptotic autocorrelation time', AIac
    
              if(('bse' in args.makeplots) or ('all' in args.makeplots)):
                  fig2 = plt.figure(2)
                  fig2.suptitle("BSE")
                  ax = plt.subplot(ncols, 1, coli+1)
                  ax.plot(BSEOfN, label='Col ' +  str(coli) + ' real', color='black')
                  line, = ax.plot(BSEOfNFit, label='Col ' +  str(coli) + ' fit', color = 'red')
                  line.set_dashes([2,2,10,2])
                  ax.legend()
      
                  ax.set_xlabel(r'$\mathrm{n(block size)}$', fontsize=8)
                  ax.set_ylabel(r'$\mathrm{BSE(n)}$', fontsize=8)
       
                  if args.savefig:
                        figFN = 'temp.bse.pdf'
                        plt.savefig(figFN, dpi=600, format='pdf')
      
          if(('data' in args.makeplots) or ('all' in args.makeplots)):
              fig3 = plt.figure(3)
              fig3.suptitle("Raw data")
              ax = plt.subplot(ncols, 1, coli+1)
              ax.plot(cols[coli], label='Col ' +  str(coli) + ' real', color='black')
              ax.legend()
  
              ax.set_xlabel(r'$\mathrm{t}$', fontsize=8)
              ax.set_ylabel(r'$\mathrm{f(t)}$', fontsize=8)

              if args.savefig:
                    figFN = 'temp.data.pdf'
                    plt.savefig(figFN, dpi=600, format='pdf')
  
if args.makeplots[0] != '' :
    plt.show()
  
  
  
  
  
  
  
  
  
  
  
  
