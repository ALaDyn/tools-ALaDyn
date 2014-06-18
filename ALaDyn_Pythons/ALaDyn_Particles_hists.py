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
nb = sys.argv[1]


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
		if file[0:6] == 'PSBunch':
			#-operate here!
			#PSBunch2_21.bin
			X, Y, Z, Px, Py, Pz, W = read_particle_phasespace(file)
			H, xedges, pxedges = np.histogram2d(X, Px, bins=(nb, nb))
			K, yedges, pyedges = np.histogram2d(Y, Py, bins=(nb, nb))
			L, zedges, pzedges = np.histogram2d(Z, Pz, bins=(nb, nb))

			fig = plt.figure(1, figsize(3.0,3.0) )
			ax1 = fig.add_subplot(3,2,1)
			ax1.imshow(H)
			ax1.xlabel('x'); ax1.ylabel('Px')
			ax2 = fig.add_subplot(3,2,2)
			ax2.imshow(K)
			ax2.xlabel('y'); ax2.ylabel('Py')
			ax3 = fig.add_subplot(3,2,3)
			ax3.imshow(L)
			ax3.xlabel('z'); ax3.ylabel('Pz')			
			
			ax4 = fig.add_subplot(3,2,4)
			ax4.scatter(X,Px)
			ax4.xlabel('x'); ax4.ylabel('Px')
			ax5 = fig.add_subplot(3,2,5)
			ax5.scatter(Y,Py)
			ax5.xlabel('y'); ax5.ylabel('Py')
			ax6 = fig.add_subplot(3,2,6)
			ax6.scatter(Z,Pz)
			ax6.xlabel('z'); ax6.ylabel('Pz')

			
			axis('tight')
			name_output = 'phasespace_bunch_'+file[7]+'_'+file[9:10]+'.png'
			savefig( os.path.joing(path,'plot','phasespace',name_output) )
			close(fig)

