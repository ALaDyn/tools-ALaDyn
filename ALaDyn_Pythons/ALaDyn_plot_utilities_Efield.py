#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_utilities_Efield.py
# Author:       A. Marocchino
# Date:			2014-02-18
# Purpose:      it is a module of: ALaDyn_plot_sections - plot E-field
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
def plot_Efield_sections(path,frame):
	s='%2.2i'%frame 				#conversion to 2-character-long-string

	
	#- background -#
	file_name = 'Exfout'+s+'.bin'
	Ex,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eyfout'+s+'.bin'
	Ey,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezfout'+s+'.bin'
	Ez,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

	#- bunch(es) -#
	file_name = 'Exbout'+s+'.bin'
	Exb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Eybout'+s+'.bin'
	Eyb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	file_name = 'Ezbout'+s+'.bin'
	Ezb,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
	
	#-sum-#
	Ex = Ex+Exb
	Ey = Ey+Eyb
	Ez = Ez+Ezb

	p = Ex.shape
	x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;
	
	sizeX, sizeZ = figure_dimension_inch(x,y,z,scale_factor)



	#- Ex_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,Ex[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Ex_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)
	#- Ey_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,Ey[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Ey_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)
	#- Ez_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,Ez[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Ez_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)




	#- Ex_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,Ex[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Ex_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)
	#- Ey_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,Ey[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Ey_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)
	#- Ez_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,Ez[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'Ez_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)



# 	#- Ex_YZ -#
# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,Ex[x2,:,:].T,100, linewidths = 0.00001)
#	axis('tight')
# 	name_output = 'Ex_YZ_'+s+'.png'
# 	savefig( os.path.join(path,'plots','E_field',name_output) )
#	close(fig)
# 	#- Ey_YZ -#
# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,Ey[x2,:,:].T,100, linewidths = 0.00001)
#	axis('tight')
# 	name_output = 'Ey_YZ_'+s+'.png'
# 	savefig( os.path.join(path,'plots','E_field',name_output) )
#	close(fig)
# 	#- Ez_YZ -#
# 	fig = figure(1, figsize=(sizeZ, sizeZ))	
# 	contourf(y,z,Ez[x2,:,:].T,100, linewidths = 0.00001)
#	axis('tight')
# 	name_output = 'Ez_YZ_'+s+'.png'
# 	savefig( os.path.join(path,'plots','E_field',name_output) )
#	close(fig)


	#- norm(E) -#
	norm_E = np.power( np.power(Ex,2)+np.power(Ey,2)+np.power(Ez,2) , 0.5 )


	#- norm_E_XY -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,y,norm_E[:,:,z2].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'E_XY_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)
	#- norm_E_XZ -#
	fig = figure(1, figsize=(sizeX, sizeZ))	
	contourf(x,z,norm_E[:,y2,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'E_XZ_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)
	#- norm_E_YZ -#
	fig = figure(1, figsize=(sizeZ, sizeZ))	
	contourf(y,z,norm_E[x2,:,:].T,100, linewidths = 0.00001)
	axis('tight')
	name_output = 'E_YZ_'+s+'.png'
	savefig( os.path.join(path,'plots','E_field',name_output) )
	close(fig)






