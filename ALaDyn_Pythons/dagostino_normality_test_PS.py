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
#from scipy.stats import mstats
import scipy.stats as mstats
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
	print '			bunch number'
	print '			stat test: shapiro / dagostino'
	exit(0)
#---
nbunch = int(sys.argv[1])
test_kind = str( sys.argv[2] )

#---
test=np.zeros((7,50))
i=0
pX_shapiro=[]; pY_shapiro=[]; pZ_shapiro=[];
pPx_shapiro=[]; pPy_shapiro=[]; pPz_shapiro=[];
#---

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
				
				#- Shapito test: iterative
				if test_kind == 'shapiro':
	 				for i in range(0,1500):	
	 					f=np.random.randint(nX,size=(4999,))
						KX,pX = shapiro(X[f]); 			KY,pY = shapiro(Y[f]); 			KZ,pZ = shapiro(Z[f])
						KPx,pPx=mstats.shapiro(Px[f]); 	KPy,pPy=mstats.shapiro(Py[f]); 	KPz,pPz=mstats.shapiro(Pz[f])
						pX_shapiro.append(pX); 		pY_shapiro.append(pY);		pZ_shapiro.append(pZ);
						pPx_shapiro.append(pPx);	pPy_shapiro.append(pPy);	pPz_shapiro.append(pPz);
					pX=np.mean(pX_shapiro); 	pY=np.mean(pY_shapiro); 	pZ=np.mean(pZ_shapiro);
					pPx=np.mean(pPx_shapiro); 	pPy=np.mean(pPy_shapiro); 	pPz=np.mean(pPz_shapiro);

				#- normal test - full sample - d'Agostino test
				if test_kind == 'dagostino':
					KX,pX=mstats.normaltest(X); 	KY,pY=mstats.normaltest(Y); 	KZ,pZ=mstats.normaltest(Z)
					KPx,pPx=mstats.normaltest(Px); 	KPy,pPy=mstats.normaltest(Py); 	KPz,pPz=mstats.normaltest(Pz)

				#- screen output
				print '[',int(file[9:11]),',', pX,',', pY,',', pZ,',', pPx,',', pPy,',', pPz,'],'
# 				test[i][:] = [int(file[7]), pX, pY, pZ, pPx, pPy, pPz]; i+=1
	
	#-printing
# 	np.savetxt( os.path.join(path, 'dagostino_normality_test.txt' ), test[0:i-1,:] )
