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

#-- Inputs
bch= 1.4
sx= 50 * 1e-4 #in [cm]
sy= 8  * 1e-4 #in [cm] - sy == sz
n0 		= 0.01*1e18	 #[cm-3]
#-QFLUID
# n0		= 1.0e16 #[cm-3]
# sx 		= 50 * 1e-4
# sy      = 20 * 1e-4
# bch		= 0.05e-9 / (n0*1.6e-19* (2.*np.pi)**(3./2.) * (sx*sy*sy) )
# - #


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
total_Charge	= bch * volume_N3 * e / 1e-12 	#total charge in [pC]
# - #

#print 'Qtilde >>>',Qtilde
print('\n')
print('lambda_p: %e [um] - omega_p: %e [Hz]' % (lp*1e4,wp) )
print('total bunch Charge: %e [pC]' % (total_Charge) )
print('Qtilde - normalized charge: %e' % (Qtilde) )
print('\n')


