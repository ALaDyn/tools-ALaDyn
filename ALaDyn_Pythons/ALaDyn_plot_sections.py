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



### --- ### shell inputs
if(len(sys.argv)<2):
	print 'Input [1]: frame_begin'
	print 'Input [2]: frame_end'
	print 'Input [3]: magnification_fig'
	print 'Input [4]: rho_min'
	print 'Input [5]: rho_max'
	print 'Input [6]: iso-lines'
	print 'Input [7]: cell to cut'

if sys.argv[1] == -1:
	frame_begin		  = 0
	frame_end         = last_output(os.getcwd())
	magnification_fig = 3.0
	rho_min 		  = 0.0001
	rho_max 		  = 20.
	isolines		  = 30
else:
	frame_begin 		= int(		sys.argv[1])	
	frame_end			= int(		sys.argv[2])
	magnification_fig 	= float(	sys.argv[3])
	rho_min 		  	= float(	sys.argv[4])
	rho_max 		  	= float(	sys.argv[5])
	isolines			= int(		sys.argv[6])
	celltocut			= int(		sys.argv[7])
### --- ###



#--- *** ---#
if __name__ == '__main__':
	
	#-path
	path = os.getcwd()
	
	#-folder output structure
	generate_folder_output_structure(path)


# 	N = last_output(path)
# 	print '-------------------'
# 	print 'last output number > ',N
# 	print '-------------------'
		
	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
		print '-------------------'
		if output_exists(path,'rho',i) == True:
			print 'rho --- frame >>> ',i
			plot_density_sections(path,i,rho_min,rho_max,isolines,celltocut,magnification_fig)

		if output_exists(path,'E',i) == True:
			print 'E --- frame >>> ',i
			plot_Efield_sections(path,i,magnification_fig)
			
		if output_exists(path,'B',i) == True:
			print 'B --- frame >>> ',i		
			plot_Bfield_sections(path,i,magnification_fig)




