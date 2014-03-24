#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_utilities_Efield.py
# Author:       A. Marocchino
# Date:			2014-02-18
# Purpose:      it is a module of: ALaDyn_plot_sections - plot B-field
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime
import numpy as np
from pylab import *
###>>>
# home_path = os.path.expanduser('~')
# sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
from read_ALaDyn_bin import *
### --- ###






#- plot Sections
def plot_Bfield_sections(path,frame):
	s='%2.2i'%frame 				#conversion to 2-character-long-string

	
	#- background -#
	file_name = 'Bxfout'+s+'.bin'
	Bx,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Byfout'+s+'.bin'
	By,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Bzfout'+s+'.bin'
	Bz,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

	#- bunch(es) -#
	file_name = 'Bxbout'+s+'.bin'
	Bxb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Bybout'+s+'.bin'
	Byb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Bzbout'+s+'.bin'
	Bzb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	
	#-sum-#
	Bx = Bx+Bxb
	By = By+Byb
	Bz = Bz+Bzb

	p = Bx.shape
	x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;
	
	sizeX, sizeZ = figure_dimension_inch(x,y,z,scale_factor)



	#- Bx_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,Bx[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Bx_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )
	#- By_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,By[:,:,z2].T,100, linewidths = 0.00001) 
	axis('tight')
	name_output = 'By_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )
	#- Bz_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,Bz[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Bz_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )




	#- Bx_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,Bx[:,y2,:].T,100, linewidths = 0.00001) 
	axis('tight')
	name_output = 'Bx_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )
	#- By_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,By[:,y2,:].T,100, linewidths = 0.00001) 
	axis('tight')
	name_output = 'By_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )
	#- Bz_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,Bz[:,y2,:].T,100, linewidths = 0.00001) 
	axis('tight')
	name_output = 'Bz_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )



# 	#- Bx_YZ -#
# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,Bx[x2,:,:].T,100, linewidths = 0.00001) 
#	axis('tight')
# 	name_output = 'Bx_YZ_'+s+'.png'
# 	savefig( os.path.join(path,'plots','B_field',name_output) )
# 	#- By_YZ -#
# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,By[x2,:,:].T,100, linewidths = 0.00001) 
#	axis('tight')
# 	name_output = 'By_YZ_'+s+'.png'
# 	savefig( os.path.join(path,'plots','B_field',name_output) )
# 	#- Bz_YZ -#
# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,Bz[x2,:,:].T,100, linewidths = 0.00001)
#	axis('tight')
# 	name_output = 'Bz_YZ_'+s+'.png'
# 	savefig( os.path.join(path,'plots','B_field',name_output) )


	#- norm(B) -#
	norm_B = np.power( np.power(Bx,2)+np.power(By,2)+np.power(Bz,2) , 0.5 )


	#- norm_E_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,norm_B[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'B_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )
	#- norm_E_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,norm_B[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'B_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )
	#- norm_E_YZ -#
	fig = figure(1, figsize=(sizeZ, sizeZ))	
	contourf(y,z,norm_B[x2,:,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'B_YZ_'+s+'.png'
	savefig( os.path.join(path,'plots','B_field',name_output) )






