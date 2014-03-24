#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_utilities_Efield.py
# Author:       A. Marocchino
# Date:			2014-02-18
# Purpose:      it is a module of: ALaDyn_plot_sections - plots densities
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime
import numpy as np
from pylab import *
from matplotlib import colors, ticker, cm
###>>>
# home_path = os.path.expanduser('~')
# sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
from read_ALaDyn_bin import *
from ALaDyn_plot_utilities_1 import *
### --- ###






#- plot Sections
def plot_density_sections(path,frame,rho_min,rho_max,magnification_fig):
	s='%2.2i'%frame 				#conversion to 2-character-long-string

	
	file_name = 'Bdenout'+s+'.bin'
	matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Edenout'+s+'.bin'
	matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	
	#- cut & sign
	matrix = - matrix
	matrix2 = - matrix2
	matrix[ (matrix<rho_min) ] = rho_min
	matrix[ (matrix>rho_max) ] = rho_max
	matrix2[ (matrix2<rho_min) ] = rho_min
	matrix2[ (matrix2>rho_max) ] = rho_max


	p = matrix.shape
	x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;
	
	sizeX, sizeZ = figure_dimension_inch(x,y,z,magnification_fig)

	levs = np.logspace(log10(rho_min),log10(rho_max),40)


	#--------------------#
	#--- Linear plots ---#
	#--------------------#
	#- Plot Bdenout -#
	fig = figure(1, figsize=(sizeX, sizeZ))
	contourf(x,y,matrix[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'rho_Bunch_XY_lin_'+s+'.png'
	savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,matrix[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'rho_Bunch_XZ_lin_'+s+'.png'
	fig.savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,matrix[x2,:,:].T,100, linewidths = 0.00001)
#	axis('tight')
# 	name_output = 'rho_Bunch_YZ_lin_'+s+'.png'
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )
#	close(fig)



	#- Plot Edenout -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,matrix2[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'rho_Background_XY_lin_'+s+'.png'
	savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,matrix2[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'rho_Background_XZ_lin_'+s+'.png'
	fig.savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,matrix2[x2,:,:].T,100, linewidths = 0.00001)
#	axis('tight')
# 	name_output = 'rho_Background_YZ_lin_'+s+'.png'
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )
#	close(fig)



	#- Plot Bdenout+Edenout -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,matrix[:,:,z2].T + matrix2[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'rho_tot_XY_lin_'+s+'.png'
	savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,matrix[:,y2,:].T + matrix2[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'rho_tot_XZ_lin_'+s+'.png'
	fig.savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

# 	fig = figure(1, figsize=(sizeX, sizeX))	
# 	contourf(y,z,matrix[x2,:,:].T - matrix2[x2,:,:].T,100, linewidths = 0.00001)
#	axis('tight')
# 	name_output = 'rho_tot_YZ_lin_'+s+'.png'
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )
#	close(fig)










	#--------------------#
	#---  Log plots   ---#
	#--------------------#
	#- Plot Bdenout -#
	fig = figure(1, figsize=(sizeX, sizeZ))
#	contourf(x,y,matrix[:,:,z2].T,100, locator=ticker.LogLocator(), norm=colors.LogNorm())
	contourf(x,y,matrix[:,:,z2].T,100, levs, norm=colors.LogNorm())
	axis('tight')
	name_output = 'rho_Bunch_XY_log_'+s+'.png'
	savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,matrix[:,y2,:].T,100, locator=ticker.LogLocator(), norm=colors.LogNorm())
	axis('tight')
	name_output = 'rho_Bunch_XZ_log_'+s+'.png'
	fig.savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,matrix[x2,:,:].T,100, levs, norm=colors.LogNorm())
#	axis('tight')
# 	name_output = 'rho_Bunch_YZ_log_'+s+'.png'
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )
#	close(fig)



	#- Plot Edenout -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,matrix2[:,:,z2].T,100, locator=ticker.LogLocator(), norm=colors.LogNorm())
	axis('tight')
	name_output = 'rho_Background_XY_log_'+s+'.png'
	savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,matrix2[:,y2,:].T,100, locator=ticker.LogLocator(), norm=colors.LogNorm())
	axis('tight')
	name_output = 'rho_Background_XZ_log_'+s+'.png'
	fig.savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,matrix2[x2,:,:].T,100, levs, norm=colors.LogNorm())
#	axis('tight')
# 	name_output = 'rho_Background_YZ_log_'+s+'.png'
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )
#	close(fig)



	#- Plot Bdenout+Edenout -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,matrix[:,:,z2].T + matrix2[:,:,z2].T,100, locator=ticker.LogLocator(), norm=colors.LogNorm())
	axis('tight')
	name_output = 'rho_tot_XY_log_'+s+'.png'
	savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,matrix[:,y2,:].T + matrix2[:,y2,:].T,100, locator=ticker.LogLocator(), norm=colors.LogNorm())
	axis('tight')
	name_output = 'rho_tot_XZ_log_'+s+'.png'
	fig.savefig( os.path.join(path,'plots','rho',name_output) )
	close(fig)

# 	fig = figure(1, figsize=(sizeX, sizeX))	
# 	contourf(y,z,matrix[x2,:,:].T + matrix2[x2,:,:].T,100, levs, norm=colors.LogNorm())
#	axis('tight')
# 	name_output = 'rho_tot_YZ_log_'+s+'.png'
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )
#	close(fig)


