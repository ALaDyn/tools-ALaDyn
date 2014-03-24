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



### --- ###
magnification_fig = 3.0



#--- *** ---#
if __name__ == '__main__':
	
	#-path
	path = os.getcwd()
	
	#-folder output structure
	generate_folder_output_structure(path)


	N = last_output(path)
	print '-------------------'
	print 'last output number > ',N
	print '-------------------'
		
	for i in range(0,N+1):
		print '-------------------'
		if output_exists(path,'rho',i) == True:
			print 'rho --- frame >>> ',i
			plot_density_sections(path,i,magnification_fig)

		if output_exists(path,'E',i) == True:
			print 'E --- frame >>> ',i
			plot_Efield_sections(path,i,magnification_fig)
			
		if output_exists(path,'B',i) == True:
			print 'B --- frame >>> ',i		
			plot_Bfield_sections(path,i,magnification_fig)




