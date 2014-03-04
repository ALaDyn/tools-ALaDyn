#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:       
# Date:			2014-02-18
# Purpose:      it nests into 'ALaDyn_read_binary' to plot sections
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime
###>>>
# home_path = os.path.expanduser('~')
# sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
### --- ###




#- get last output number
def last_output(path):
	n_last_output = -1
	for root_dir, sub_dirs, files in os.walk(path):
		for file in files:
			if os.path.splitext(file)[1] == '.bin':
				print os.path.basename(file)[7:8]
				n_last_output = max( n_last_output, int(os.path.basename(file)[7:8]) )
	return n_last_output

#- folder structure for outputs -#
def	generate_folder_output_structure(path):
	directory = os.path.join(path,'plots')
	if not os.path.exists( directory ):
		os.makedirs(directory)
	
	directory_rho = os.path.join(directory,'rho')
	directory_E   = os.path.join(directory,'E_field')
	directory_B   = os.path.join(directory,'B_field')

	if not os.path.exists( directory_rho ):
		os.makedirs(directory_rho)
	if not os.path.exists( directory_E ):
		os.makedirs(directory_E)
	if not os.path.exists( directory_B ):
		os.makedirs(directory_B)

