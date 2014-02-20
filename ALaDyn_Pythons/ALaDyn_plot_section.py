#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:       
# Date:			2014-02-18
# Purpose:      it nests into 'ALaDyn_read_binary' to plot sections
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy, numpy, pylab, matplotlib
import struct
from scipy import *
import numpy as np
from matplotlib import *
from pylab import *
import matplotlib as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
### --- ###
from ALaDyn_read_binary import *
### --- ###


#E-fields in GV/m

# - #
path      = '/Users/alberto/sims/ALaDyn_sims/00001' 
file_name = 'Bdenout04.bin'
matrix = read_ALaDyn_bin(path,file_name)
file_name = 'Edenout04.bin'
matrix2 = read_ALaDyn_bin(path,file_name)
matrix = matrix2 # + matrix2


# file_name = 'Ezbout00.bin'
# matrix = read_ALaDyn_bin(path,file_name)


xlist = linspace(0.,1.,256.) 
ylist = linspace(0.,1.,128.) 
X, Y = meshgrid (ylist, xlist)
ax  = matplotlib.pyplot.subplot(111) #matplotlib.pyplot.subplot(111)
pyplot.imshow(-matrix[:,:,64].T)
pyplot.colorbar()
#pyplot.plot(-matrix[:,64,64])
show()


# 	#np.savetxt('test.txt', r[:,:,64])

















