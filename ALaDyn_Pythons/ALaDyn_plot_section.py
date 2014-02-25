#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:       
# Date:			2014-02-18
# Purpose:      it nests into 'ALaDyn_read_binary' to plot sections
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy
import struct
from scipy import *
import numpy as np
from pylab import *
import matplotlib as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
###>>>
#sys.path.append(os.path.split(os.getcwd())[0])
###>>>
### --- ###
from read_ALaDyn_bin import *
from ALaDyn_from_XYZ_to_surface import *
### --- ###


#E-fields in GV/m

# - #
path      = '/Users/alberto/sims/ALaDyn_sims/00002' 
file_name = 'Bdenout15.bin'
matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Edenout15.bin'
matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
matrix = matrix+matrix2

fig = figure()
ax  = fig.add_subplot(111) #, aspect='equal')
ax.contourf(x,y,-matrix[:,:,64].T,100) #, 15, linewidths = 0.5, colors = 'k')
show()


# fig = figure()
# ax  = fig.add_subplot(111, aspect='equal')
# ax.imshow(-matrix[:,:,64])
# show()




