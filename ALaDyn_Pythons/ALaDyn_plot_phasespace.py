#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_phasespace.py
# Author:       F. Mira
# Date:         2016-11-03
# Purpose:      reads ALaDyn binary from PDBunch* and plots phasespace
# Source:       python
#####################################################################

### loading shell commands
import os,sys
import struct
from scipy import *
import numpy as np
import pylab as pyl
import matplotlib.pyplot as plt
###>>>
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
### --- ###
from Particle_reader_utilities import *
from ALaDyn_plot_utilities_1 import *
### --- ###


# - #
### --- ### shell inputs
if(len(sys.argv)<4):
	print 'Usage:'
	print 'Phasespace Component1:'
	print 'Phasespace Component2:'
	print 'frame begin: '
	print 'frame end: '
	exit(0)
#---
component1  = str(sys.argv[1])
component2  = str(sys.argv[2])
frame_begin = int(sys.argv[3])
frame_end   = int(sys.argv[4])    

#--- *** ---#
if __name__ == '__main__':
	
	#-path
	path = os.getcwd()
	#-folder output structure
	generate_folder_phasespace(path)
	#-hunting and printing
	for i in range(frame_begin, frame_end + 1 ):
		print '-------------------'
		s='%2.2i'%i
		path_read  = os.path.join(path,'%4.4i'%i)
		path_write = path
		if output_exists(path_read,'phasespace',i) == True:
			
			print 'phasespace --- frame >>> ',i
			
			C1 = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin',component1 )
			C2 = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin',component2 )
			C1 = np.array(C1);    C2 = np.array(C2);
			
			plt.figure(1, figsize=(11.0,11.0) )
			#plt.plot(C1,C2,'ro')
			plt.scatter(C1[range(0,len(C1),10)],C2[range(0,len(C2),10)],s=.1,edgecolors='None')
		        plt.xlabel(component1); plt.ylabel(component2)	
			plt.axis('tight')
			name_output =  component1 +'_'+ component2 + '_' +('%2.2i'%i)+ '.png'
			plt.savefig( os.path.join(path_write,'data','phasespace',name_output) )
			plt.close()

