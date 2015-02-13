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
from matplotlib import colors, ticker, cm
###>>>
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
### --- ###
from read_ALaDyn_bin import *
from ALaDyn_plot_utilities_1 import *
from ALaDyn_plot_utilities_density import *
from ALaDyn_plot_utilities_Efield import *
from ALaDyn_plot_utilities_Bfield import *
### --- ###

#--- *** ---#
path = os.getcwd()

# --- #
# FRAME NUMBER
frame = 6


# --- #
magnification_fig = 3.0

rho_min = 10
rho_max = 0.001




#--- *** ---#
file_name = 'Bdenout'+('%2.2i'%frame)+'.bin'
matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Edenout'+('%2.2i'%frame)+'.bin'
matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')


#--- *** ---#
file_name = 'Exbout'+('%2.2i'%frame)+'.bin'
Ex_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Eybout'+('%2.2i'%frame)+'.bin'
Ey_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Ezbout'+('%2.2i'%frame)+'.bin'
Ez_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

file_name = 'Exfout'+('%2.2i'%frame)+'.bin'
Ex_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Eyfout'+('%2.2i'%frame)+'.bin'
Ey_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Ezfout'+('%2.2i'%frame)+'.bin'
Ez_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')





#- cut
matrix = - matrix
matrix2 = - matrix2
matrix[ (matrix<rho_min) ] = rho_min
matrix[ (matrix>rho_max) ] = rho_max
matrix2[ (matrix2<rho_min) ] = rho_min
matrix2[ (matrix2>rho_max) ] = rho_max

p = matrix.shape
x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;


size1, size2 = figure_dimension_inch(x,y,z,magnification_fig)


# --- #
fig = figure(1, figsize=(size1, size2))
levs = np.linspace(rho_min,rho_max,50)
contourf(x,y,matrix[:,:,z2].T + matrix2[:,:,z2].T, levs, linewidths = 0.0001)
levs = np.logspace(log10(rho_min),log10(rho_max),50)
contourf(x,y,matrix[:,:,z2].T + matrix2[:,:,z2].T, levs, norm=colors.LogNorm(), linewidths = 0.0001)
#axis('tight')
##name_output = 'rho_tot_XY_'+s+'.png'
##savefig( os.path.join(path,'plots','rho',name_output) )
show()

# --- #
fig = figure(1, figsize=(size1, size2))
levs = np.linspace(-1.0,1.0,100)
#contourf(x,y,-matrix[:,:,z2].T,100, linewidths = 0.0001)
contourf(x,y,Ex_w[:,:,z2].T,levs, linewidths = 0.0001)
#clim(0.001,30.0)
show()

# --- #
fig = figure()
plot(x,Ez_w[:,y2,z2].T)
#plot(x,-matrix[:,y2,z2].T - matrix2[:,y2,z2].T)
#plot(x,-matrix[:,y2+1,z2+1].T - matrix2[:,y2+1,z2+1].T)
#plot(x,-matrix[:,y2+2,z2+2].T - matrix2[:,y2+2,z2+2].T)
show()

# --- #
Matrix = np.column_stack( (x,matrix[:,y2,z2].T) )
np.savetxt( 'section.dat' ,Matrix,fmt='%15.14e')
np.savetxt( 'section.dat' ,matrix[:,:,z2].T,fmt='%15.14e')









