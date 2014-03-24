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

#--- *** ---#
file_name = 'Bdenout47.bin'
matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Edenout47.bin'
matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

p = matrix.shape
x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;


size1, size2 = figure_dimension_inch(x,y,z,2.0)


fig = figure(1, figsize=(size1, size2))	
#ax  = fig.add_subplot(111) #, aspect='equal')
contourf(x,y,-matrix[:,:,z2].T - matrix2[:,:,z2].T,100, linewidths = 0.0001)
clim(0.001,30.0)
#axis('tight')
##name_output = 'rho_tot_XY_'+s+'.png'
##savefig( os.path.join(path,'plots','rho',name_output) )
show()


fig = figure(1, figsize=(size1, size2))
contourf(x,y,-matrix[:,:,z2].T,100, linewidths = 0.0001)
clim(0.001,30.0)
show()




fig = figure()
plot(x,-matrix[:,y2,z2].T - matrix2[:,y2,z2].T)
plot(x,-matrix[:,y2+1,z2+1].T - matrix2[:,y2+1,z2+1].T)
plot(x,-matrix[:,y2+2,z2+2].T - matrix2[:,y2+2,z2+2].T)
show()





f = open(full_file_name,'w+')
f.close()
Matrix = np.column_stack( (x,matrix[:,:,z2].T) )
np.savetxt( 'section.dat' ,Matrix,fmt='%15.14e')
np.savetxt( 'section.dat' ,matrix[:,:,z2].T,fmt='%15.14e')









