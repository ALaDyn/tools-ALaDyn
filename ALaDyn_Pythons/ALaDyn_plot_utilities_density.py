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
from pylab import *
###>>>
# home_path = os.path.expanduser('~')
# sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
from read_ALaDyn_bin import *
### --- ###






#- plot Sections
def plot_density_sections(path,frame):
	s='%2.2i'%frame 				#conversion to 2-character-long-string

	
	file_name = 'Bdenout'+s+'.bin'
	matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Edenout'+s+'.bin'
	matrix2,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	

	p = matrix.shape
	x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;


	#- Plot Bdenout -#
	fig = figure()
	ax  = fig.add_subplot(111, aspect='equal')
	ax.contourf(x,y,-matrix[:,:,z2].T,100) #, 15, linewidths = 0.5, colors = 'k')
	name_output = 'rho_Bunch_XY_'+s+'.png'
	ax.tight_layout()
	savefig( os.path.join(path,'plots','rho',name_output) )

	fig = figure()
	ax  = fig.add_subplot(111, aspect='equal')
	ax.contourf(x,z,-matrix[:,y2,:].T,100)
	name_output = 'rho_Bunch_XZ_'+s+'.png'
	ax.tight_layout()
	fig.savefig( os.path.join(path,'plots','rho',name_output) )

# 	fig = figure()
# 	ax  = fig.add_subplot(111, aspect='equal')
# 	ax.contourf(y,z,-matrix[x2,:,:].T,100)
# 	name_output = 'rho_Bunch_YZ_'+s+'.png'
#	ax.tight_layout()
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )



	#- Plot Edenout -#
	fig = figure()
	ax  = fig.add_subplot(111, aspect='equal') #, aspect='equal')
	ax.contourf(x,y,-matrix2[:,:,z2].T,100) #, 15, linewidths = 0.5, colors = 'k')
	name_output = 'rho_Background_XY_'+s+'.png'
	ax.tight_layout()
	savefig( os.path.join(path,'plots','rho',name_output) )

	fig = figure()
	ax  = fig.add_subplot(111, aspect='equal')
	ax.contourf(x,z,-matrix2[:,y2,:].T,100)
	name_output = 'rho_Background_XZ_'+s+'.png'
	ax.tight_layout()
	fig.savefig( os.path.join(path,'plots','rho',name_output) )

# 	fig = figure()
# 	ax  = fig.add_subplot(111, aspect='equal')
# 	ax.contourf(y,z,-matrix2[x2,:,:].T,100)
# 	name_output = 'rho_Background_YZ_'+s+'.png'
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )



	#- Plot Bdenout+Edenout -#
	fig = figure()
	ax  = fig.add_subplot(111, aspect='equal') #, aspect='equal')
	ax.contourf(x,y,-matrix[:,:,z2].T - matrix2[:,:,z2].T,100) #, 15, linewidths = 0.5, colors = 'k')
	name_output = 'rho_tot_XY_'+s+'.png'
	ax.tight_layout()
	savefig( os.path.join(path,'plots','rho',name_output) )

	fig = figure()
	ax  = fig.add_subplot(111, aspect='equal')
	ax.contourf(x,z,-matrix[:,y2,:].T - matrix2[:,y2,:].T,100)
	name_output = 'rho_tot_XZ_'+s+'.png'
	ax.tight_layout()
	fig.savefig( os.path.join(path,'plots','rho',name_output) )

# 	fig = figure()
# 	ax  = fig.add_subplot(111, aspect='equal')
# 	ax.contourf(y,z,-matrix[x2,:,:].T - matrix2[x2,:,:].T,100)
# 	name_output = 'rho_tot_YZ_'+s+'.png'
#	ax.tight_layout()
# 	fig.savefig( os.path.join(path,'plots','rho',name_output) )


