#!/usr/bin/python
######################################################################
# Name:         clean_3D_frame
# Author:       A. Marocchino
# Date:			2014-12-10
# Purpose:      delete a 3D frame: as input - frame to delete
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, sys, shutil, time, string
### --- ###


#Path
path = os.getcwd()

#namelist files
nml_file = ['ALaDyn','aladyn','input.nml'];

#Exclusion List
Exc_list = ['.nml','nml_1D','nml_2D','.f','.f90','.c','.cpp','.py'];

#---Frame name
frame_name=['Bdenout','Edenout','Bxfout','Bybout','Byfout','Bzbout','Bzfout','Exbout','Exfout','Eybout','Eyfout','Ezbout','Ezfout','Jxbout'];
extensions=['.dat','.bin'];


#Search and Delete
for root_dir, sub_dirs, files in os.walk(path):
	for file in files:
		for left in frame_name:
			for end in extensions:
				if file == [left,'10',end]:
					print file
	
	
	
	
		
# 	nml_exists = 0
# 	for file in files:
# 		if file in nml_file:
# 			nml_exists += 1
# 		if file not in nml_file and os.path.splitext(file)[1] not in Exc_list:
# 			os.remove(os.path.join(root_dir,file))
# 			
# 	if nml_exists > 0 and os.path.exists(os.path.join(root_dir,'out')) == True:
#  		shutil.rmtree(os.path.join(root_dir,'out'))
# 	if nml_exists > 0 and os.path.exists(os.path.join(root_dir,'dump')) == True:
#  		shutil.rmtree(os.path.join(root_dir,'dump'))
# 	
# 		
	
				