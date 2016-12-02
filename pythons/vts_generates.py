#!/usr/bin/python
######################################################################
# Name:         ALaDyn_generates_vts.py
# Author:       A. Marocchino
# Date:			2014-02-18
# Purpose:      generates vts file for 3D plots
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, base64
import numpy as np
###>>>
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/ALaDyn_Code/tools-ALaDyn/pythons'))
### --- ###
from read_ALaDyn_bin import *
from utilities_1 import *
from vts_write import *
### --- ###



### --- ### shell inputs
if(len(sys.argv)<2):
	print 'Input [1]: frame_begin'
	print 'Input [2]: frame_end'
	print 'Input [3]: yz-shaving'
	print 'Input [4]: sliceposition_y'
	exit(0)

if sys.argv[1] == -1:
	frame_begin		  = 0
	frame_end         = last_output(os.getcwd())
	cell_cut		  = 0
	sliceposition_y	  = 0
else:
	frame_begin 		= int(		sys.argv[1])
	frame_end			= int(		sys.argv[2])
	cell_cut		 	= int(		sys.argv[3])
	sliceposition_y		= int(		sys.argv[4])
### --- ###


#--- *** ---#
if __name__ == '__main__':

	#-path
	path = os.getcwd()

	#-folder output structure
	generate_folder_vts(path)


#	write_vts(path,2)

	frame =0; s='%2.2i'%frame
 	rhobunch, X,Y,Z	= read_ALaDyn_bin(path,'Bdenout'+s+'.bin','grid')

	print '--- Generating VTS ---'
	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
# 		print '--- COMPLETE vts: frame >',i
#		write_vts(path,i,X,Y,Z,cell_cut)

		print '--- write vts-section longitudinal: frame >',i
		write_vts_section_longitudinal(path,i,X,Y,Z,cell_cut,sliceposition_y)

		#- in such a way it is much longer, but it is neater to write
		#- plus: it has to be done only once
		print '--- density only vts: frame >',i
		write_density_vts(path,i,X,Y,Z,cell_cut)

		print '--- E-field only vts: frame >',i
		write_E_vts(path,i,X,Y,Z,cell_cut)

		print '--- B-field only vts: frame >',i
		write_B_vts(path,i,X,Y,Z,cell_cut)
