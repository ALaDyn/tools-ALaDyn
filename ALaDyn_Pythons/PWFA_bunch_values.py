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
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as plt
from pylab import *
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#-- Inputs
bch= 2
sx= 50 #in [um]
sy= 8  #in [um]
# - #


# - Phys. Quantities - #
n0 		= 1e18	 #[cm-3]
e  		= 1.6e-19 #[C]
e_cgs 	= 4.8e-10 #[statC]
me 		= 9.1093e-31	#[kg]
me_cgs 	= 9.1093e-28 #[g]
# - #

# - #
volume 			= (2.*np.pi)**(3./2.) * (sx*sy*sy)**(3./2.)
total_charge 	= volume * bch*n0*e
wp 				= ( 4.*np.pi* (bch*n0) *e_cgs**2 / me_cgs )**(1./2.)
kp				= c_cgs/wp
reduced_charge 	= bch * volume * kp**3
# - #













