#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
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
from read_ALaDyn_bin import *
from ALaDyn_from_XYZ_to_surface import *
### --- ###


#E-fields in GV/m

# - #
path      = '/Users/alberto/sims/ALaDyn_sims/00001_test' 
file_name = 'Bdenout00.bin'
matrix,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# file_name = 'Edenout04.bin'
# matrix2 = read_ALaDyn_bin(path,file_name)
matrix = matrix #+ matrix2




ax  = matplotlib.pyplot.subplot(111)
# a=np.zeros((128,64)); b=a;
# for i in range(0,128):
# 	for j in range(0,64):
# 		a(i,j) = 
# 		b(i,j) = 
a=np.linspace(0,1,256);
b=np.linspace(0,1,128);
print '>>>',len(x),len(y),matrix.shape
ax.contourf(y,x,-matrix[:,:,64]) #, 15, linewidths = 0.5, colors = 'k')
ax.equal()
show()





#ALaDyn_from_XYZ_to_surface(x[:,:,64],y[:,:,64],-matrix[:,:,64])
#os.pause()



# 
# 
# xlist = linspace(0.,1.,256.) 
# ylist = linspace(0.,1.,128.) 
# X, Y = meshgrid (ylist, xlist)
ax  = matplotlib.pyplot.subplot(111) #matplotlib.pyplot.subplot(111)
pyplot.imshow(-matrix[:,:,64].T)
# pyplot.colorbar()
# #pyplot.plot(-matrix[:,64,64])
show()
# 

# 	#np.savetxt('test.txt', r[:,:,64])

















