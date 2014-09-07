#!/usr/bin/python
######################################################################
# Name:         ALaDyn_Particles_hists.py
# Author:       A Marocchino      
# Date:         2014-06-18
# Purpose:      reads dued binary from PDBunch* and plots hists
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
### --- ###


# - #
### --- ### shell inputs
if(len(sys.argv)<2):
	print 'Usage:'
	print 'Input [1]: bin number'
	exit(0)
#---
nbins = int(sys.argv[1])


#--- *** ---#
if __name__ == '__main__':
	
	#-path
	path = os.getcwd()
	#-folder output structure
	generate_folder_phasespace(path)
	#-hunting and printing
	for root_dir, sub_dirs, files in os.walk(path):
		prune_dirs(sub_dirs)
		for file in files:
			if file[0:7] == 'PSBunch' and os.path.splitext(file)[1] == '.bin':
				#-operate here!
				#PSBunch2_21.bin
				print file
				X, Y, Z, Px, Py, Pz, W = read_particle_phasespace( os.path.join(root_dir,file) )
				X=np.array(X);    Y=np.array(Y);   Z=np.array(Z)
				Px=np.array(Px); Py=np.array(Py); Pz=np.array(Pz)
				
				H, xedges, pxedges = np.histogram2d(X, Px, bins=(nbins, nbins))
				K, yedges, pyedges = np.histogram2d(Y, Py, bins=(nbins, nbins))
				L, zedges, pzedges = np.histogram2d(Z, Pz, bins=(nbins, nbins))

				fig = plt.figure(1, figsize=(11.0,11.0) )
				ax1 = fig.add_subplot(2,3,1)
				ax1.contourf(0.5*(xedges[0:nbins]+xedges[1:nbins+1]), 0.5*(pxedges[0:nbins]+pxedges[1:nbins+1]), H.T,100,linewidths=.001)
				#ax1.xlabel('x'); ax1.ylabel('Px')
				ax2 = fig.add_subplot(2,3,2)
				ax2.contourf(0.5*(yedges[0:nbins]+yedges[1:nbins+1]), 0.5*(pyedges[0:nbins]+pyedges[1:nbins+1]), K.T,100,linewidths=.001)
				#ax2.xlabel('y'); ax2.ylabel('Py')
				ax3 = fig.add_subplot(2,3,3)
				ax3.contourf(0.5*(zedges[0:nbins]+zedges[1:nbins+1]), 0.5*(pzedges[0:nbins]+pzedges[1:nbins+1]), L.T,100,linewidths=.001)
				#ax3.xlabel('z'); ax3.ylabel('Pz')			
			
				ax4 = fig.add_subplot(2,3,4)
				ax4.scatter(X[range(0,len(X),10)],Px[range(0,len(Px),10)],s=.1,edgecolors='None')
				#ax4.xlabel('x'); ax4.ylabel('Px')
				ax5 = fig.add_subplot(2,3,5)
				ax5.scatter(Y[range(0,len(Y),10)],Py[range(0,len(Py),10)],s=.1,edgecolors='None')
				#ax5.xlabel('y'); ax5.ylabel('Py')
				ax6 = fig.add_subplot(2,3,6)
				ax6.scatter(Z[range(0,len(Z),10)],Pz[range(0,len(Pz),10)],s=.1,edgecolors='None')
				#ax6.xlabel('z'); ax6.ylabel('Pz')

			
				plt.axis('tight')
				name_output =  'phasespace_bunch_' + str(file[7]) + '_' + str(file[9:11])+ '.png'
				plt.savefig( os.path.join(path,'data','phasespace',name_output) )
				plt.close(fig)

