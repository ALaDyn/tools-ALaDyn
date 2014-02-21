#!/usr/bin/python
######################################################################
# Name:         ALaDyn_from_XYZ_to_surface.py
# Author:       
# Date:			2014-02-18
# Purpose:      it nests into 'ALaDyn_read_binary' to plot sections
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy, numpy, pylab, matplotlib
import struct
from scipy import *
import numpy as np
from matplotlib import *
from pylab import *
import matplotlib as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
### --- ###


def ALaDyn_from_XYZ_to_surface(X,Y,Z):

# 	X=X.T; Y=Y.T; Z=Z.T;

	X = np.array(X).flatten()
	Y = np.array(Y).flatten()
	Z = np.array(Z).flatten()

	nx = 100
	
	xi = np.linspace(X.min(), X.max(), nx)
	yi = np.linspace(Y.min(), Y.max(), nx)
# 	for i in range(0,len(X)):
# 		print X[i],Y[i],Z[i]
# 	for i in range(0,256):
# 		for j in range(0,128):
# 			X[i]=i
# 			Y[j]=j
# 			print i,j
	zi = griddata(X, Y, Z, xi, yi)
	

	ax  = matplotlib.pyplot.subplot(111)
	ax.contour(xi, yi, zi) #, 15, linewidths = 0.5, colors = 'k')
	show()






# ax  = matplotlib.pyplot.subplot(111) #matplotlib.pyplot.subplot(111)
# pyplot.imshow(-matrix[:,:,64].T)
# pyplot.colorbar()
# #pyplot.plot(-matrix[:,64,64])
# show()
















