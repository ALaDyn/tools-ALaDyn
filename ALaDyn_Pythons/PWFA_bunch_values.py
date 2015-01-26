#!/usr/bin/python
######################################################################
# Name:         PWFA_bunch_values
# Author:       A. Marocchino
# Date:			18-02-2014
# Purpose:      from ALaDyn input -> to some Relevant physical quantities
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy, pylab
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as plt
from pylab import *
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


### --- ### shell inputs
if(len(sys.argv)<2):
	print 'Input  [1]: bch (n_bunch/n_0)'
	print 'Input  [2]: sx    [cm]'
	print 'Input  [3]: sy,sz [cm]'
	print 'Input  [4]: n_0   [cm^-3]'
	exit(0)

try:
	switch=int(sys.argv[1])
except:
	switch=float(sys.argv[1])

if switch == -1:
	bch= 2.6
	sx= 50 * 1e-4 #in [cm]
	sy= 6 * 1e-4 #in [cm] - sy == sz
	n0 		= 0.01*1e18	 #[cm-3]
else:
	bch	= float(	sys.argv[1])
	sx	= float(	sys.argv[2])
	sy	= float(	sys.argv[3])
	n0 	= float(	sys.argv[4])
### --- ###


# - Phys. Constants - #
e  		= 1.6e-19 #[C]
e_cgs 	= 4.8e-10 #[statC]
me 		= 9.1093e-31	#[kg]
me_cgs 	= 9.1093e-28 #[g]
c		= 3.0e8 #[m/s]
c_cgs	= 3.0e10 #[cm/s]
# - #



# - #
volume_N3		= (2.*np.pi)**(3./2.) * (sx*sy*sy)
wp 				= ( 4.*np.pi* n0 *e_cgs**2 / me_cgs )**(1./2.)
kp				= wp/c_cgs
lp				= 2.*np.pi/kp
Qtilde			= bch * volume_N3 * kp**3	# reduced charge
total_Charge	= (bch*n0) * volume_N3 * e / 1e-12 	#total charge in [pC]
# - #

#print 'Qtilde >>>',Qtilde
print('\n')
print('lambda_p: %e [um] - omega_p: %e [Hz]' % (lp*1e4,wp) )
print('total bunch Charge: %i [pC]' % (total_Charge) )
print('Qtilde - normalized charge: %5.3e' % (Qtilde) )
print('\n')


