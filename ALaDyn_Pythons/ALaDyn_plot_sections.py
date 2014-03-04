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
from ALaDyn_from_XYZ_to_surface import *
from ALaDyn_plot_utilities_1 import *
from ALaDyn_plot_utilities_density import *
### --- ###





#--- *** ---#
if __name__ == '__main__':
	
	#-path
	path = os.getcwd()
	
	#-folder output structure
	generate_folder_output_structure(path)

	N = last_output(path)
	print '-------------------'
	print 'N > ',N
	print '-------------------'
		
	for i in range(0,N):
		plot_density_sections(path,i)




