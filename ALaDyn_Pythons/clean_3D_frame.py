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

if len(sys.argv)!=2:
    print "Usage path+clean_3D_frame.py frame_number"
    exit(0)
    pass
frame = int(sys.argv[1])

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
				if file == [left+'%2.2i'%frame+end][0]:
					print 'deleting:\t',file
					#print root_dir,file
					os.remove(os.path.join(root_dir,file))
	
	
	
	
		
	
				
