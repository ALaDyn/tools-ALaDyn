#!/usr/bin/python
######################################################################
# Name:         beam_loading
# Author:       A. Marocchino
# Date:			2015 12 18
# Purpose:      beam loading condition for PWFA
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
from blessings import Terminal
from termcolor import colored
#from mpl_toolkits.mplot3d import Axes3D
# - #
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/PlasmaParameters'))
from plasma_basic_parameters import *
# --- #

### --- ### shell inputs
if(len(sys.argv)<2):
    print 'Input  [1]: n_0: plasma background density (cm-3)'
    print 'Input  [2]: sigma_x,y driver bunch [um]'
    print 'Input  [3]: driver-alpha'
    print 'Input  [4]: desired E-field within bunch [GV/m]'
    exit(0)


#--- read input ---#
n0 	    = float( sys.argv[1])
sigma_x = float( sys.argv[2])
alpha   = float( sys.argv[3])
Ein     = float( sys.argv[4])
### --- ###

#- -#
n0_m3=n0*1e6
sigma_x_m = sigma_x*1e-6
Ein_GVm = Ein*1e9

k0 = electron_plasma_wavenumber(n0_m3)
E0 = electron_mass*c/electron_charge * electron_plasma_frequency(n0_m3)
Rb = 2.*np.sqrt(alpha)*sigma_x_m

Qw = 0.047*np.sqrt(1e16*1e6/n0_m3)*(k0*Rb)**4 * (E0/Ein_GVm)


print('\n')
print colored('matched charge: %8.4f [pC]' % (Qw*1e3) , 'cyan')
print('\n')
