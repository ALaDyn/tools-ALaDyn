#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:       
# Date:			2014-02-18
# Purpose:      it nests into 'ALaDyn_read_binary' to plot sections
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil
###>>>
# home_path = os.path.expanduser('~')
# sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
### --- ###
from read_ALaDyn_bin import *
from ALaDyn_plot_utilities_1 import *
from ALaDyn_plot_utilities_density import *
from ALaDyn_plot_utilities_Efield import *
from ALaDyn_plot_utilities_Bfield import *
### --- ###




#--- *** ---#
file_name = 'Bdenout10.bin'
matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Edenout10.bin'
matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

p = matrix.shape
x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;


fig = figure()
ax  = fig.add_subplot(111) #, aspect='equal')
ax.contourf(x,y,-matrix[:,:,z2].T - matrix2[:,:,z2].T,100, linewidths = 0.0001)
axis('tight')
name_output = 'rho_tot_XY_'+s+'.png'
savefig( os.path.join(path,'plots','rho',name_output) )













