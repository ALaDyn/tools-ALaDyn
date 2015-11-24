#!/usr/bin/python
######################################################################
# Name:         PWFA_bunch_values
# Author:       A. Marocchino
# Date:			24-11-2015
# Purpose:      matching condition for PWFA
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy, pylab
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#from mpl_toolkits.mplot3d import Axes3D
# - #
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/PlasmaParameters'))
from plasma_basic_parameters import *
# --- #

### --- ### shell inputs
if(len(sys.argv)<2):
	print 'Input  [1]: n_0: plasma background density (cm-3)'
	print 'Input  [2]: Lorentz-gamma'
	print 'Input  [3]: normalized transverse emittance [mm-mrad]'
	exit(0)


#--- read input ---#
n0 	  = float(	sys.argv[1])
gamma = float(  sys.argv[2])
epsn  = float(  sys.argv[3]) 
### --- ###

#- -#
n0_m3=n0*1e6
epsn_mmmrad = epsn*1e-6


k0 = electron_plasma_wavenumber(n0_m3)
sigma_matching = np.sqrt(np.sqrt(2./gamma)) * np.sqrt(epsn_mmmrad/k0)

print('\n')
print('matched sigma: %8.4f [um]' % (sigma_matching*1e6) )
print('\n')


