#!/usr/bin/python
######################################################################
# Name:         Particles_reader_utility.py
# Author:       A Marocchino      
# Date:         2014-06-18
# Purpose:      reads dued binary from output
# Source:       python
#####################################################################

### loading shell commands
import os,sys
import struct
import numpy as np
### --- ###



# - #
def read_particle_phasespace_bycomponent(dir_path,file_name,component):
	

	path     = os.path.join(dir_path,file_name)
	f        = open(path,'rb')
	cmp = []
	while True:
		try:
			Np=struct.unpack('i', f.read(4))
		except struct.error:
			#print "End of File"
			break
		vars=[]
		
		try:
			pass
			for i in range(0,Np[0]):
				vars=struct.unpack('ffffffff', f.read(4*8))
				if(component == 'X'): 	cmp.append(vars[0])
				if(component == 'Y'): 	cmp.append(vars[1])
				if(component == 'Z'): 	cmp.append(vars[2])
				if(component == 'Px'): 	cmp.append(vars[3])
				if(component == 'Py'): 	cmp.append(vars[4])
				if(component == 'Pz'): 	cmp.append(vars[5])
				if(component == 'W'): 	cmp.append(vars[6])
				if(component == 'Q'): 	cmp.append(vars[7])					
		except:
			pass	
	f.close()

	return cmp

#- folder structure for outputs -#
def	generate_folder_phasespace(path):
	directory = os.path.join(path,'plots','phasespace')
	if not os.path.exists( directory ):
		os.makedirs(directory)
	directory = os.path.join(path,'data','phasespace')
	if not os.path.exists( directory ):
		os.makedirs(directory)

#- prune dirs -#
def prune_dirs(dirs):
	if 'data' in dirs:
		dirs.remove('data')
	if 'plots' in dirs:
		dirs.remove('plots')
