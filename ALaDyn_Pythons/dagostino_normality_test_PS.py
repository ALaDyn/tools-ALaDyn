#!/usr/bin/python
######################################################################
# Name:         dagostino_normality_test_PS.py
# Author:       A Marocchino      
# Date:         2014-06-18
# Purpose:      evalute with D'agostino test normality of distribution for the whole Phase Space
# Source:       python
#####################################################################

### loading shell commands
import os,sys
import struct
import numpy as np
import pylab as pyl
import matplotlib.pyplot as plt
import scipy as sci
from scipy.stats import mstats
from scipy.stats import shapiro
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
	print 'bunch number'
	exit(0)
#---
nbunch = int(sys.argv[1])


#--- *** ---#
if __name__ == '__main__':
	
	#-path
	path = os.getcwd()
	#-Phase Space file hunting
	for root_dir, sub_dirs, files in os.walk(path):
		#prune_dirs(sub_dirs)
		for file in files:
			if file[0:7] == 'PSBunch' and int(file[7]) == nbunch and os.path.splitext(file)[1] == '.bin':
				
				X, Y, Z, Px, Py, Pz, W = read_particle_phasespace( os.path.join(root_dir,file) )
				X=np.array(X);    Y=np.array(Y);   Z=np.array(Z)
				Px=np.array(Px); Py=np.array(Py); Pz=np.array(Pz)

				nX = X.shape[0]; nY = Y.shape[0]; nZ = Z.shape[0]
				nPx = Px.shape[0]; nPy = Py.shape[0]; nPz = Pz.shape[0]
				
				#- Shapito test
				f=np.random.randint(nX,size=(4999,))
				KX,pX = shapiro(X[f])
				print KX,pX
# 				KX,pX = shapiro(X)
# 				KX,pX=mstats.normaltest(X)
				#t = mstats.chisquare(2).ppf( 0.95 )
				#print SX,SY,SZ,SPx,SPy,SPz
	
	
# 	#-hunting and printing
# 
# 			
# 				plt.axis('tight')
# 				name_output =  'phasespace_bunch_' + str(file[7]) + '_' + str(file[9:11])+ '.png'
# 				plt.savefig( os.path.join(path,'data','phasespace',name_output) )
# 				plt.close(fig)
# 
