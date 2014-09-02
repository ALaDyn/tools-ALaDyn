#!/usr/bin/python
######################################################################
# Name:         rm_duplicate_lines_afterestart.py
# Author:       A. Marocchino
# Date:			2014-07-25
# Purpose:      remove duplicate lines for lineout
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil
from matplotlib import colors, ticker
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
###>>>
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
from read_ALaDyn_bin import *
from ALaDyn_plot_utilities_1 import *
### --- ###

# - #
path = os.getcwd()

# - #
#read standard output
std_output = np.genfromtxt( os.path.join(path, 'bunch_1.dat') , unpack=True);

# - # find duplicates lines -- find jumps
s = std_output.shape
to_save=[]; to_delete=[]; ckc_point=-1.;
for i in range(1,s[1]-1):
	if std_output[0,i] > ckc_point:
		if std_output[0,i] < std_output[0,i-1]:
			to_delete.append( i )
			ckc_point = std_output[0,i-1]
		else:
			to_save.append( i )
	else:
		to_delete.append( i )
#	raw_input("Press enter to continue")

# - # plot
fig = plt.figure(1)
ax  = plt.subplot(111) #matplotlib.pyplot.subplot(111)
ax.plot( std_output[0,to_save], std_output[7,to_save], 'o') #, markersize=5.0)
ax.plot( std_output[0,to_delete], std_output[7,to_delete], 'o-')#, markersize=5.0)
ax.grid()
plt.show()

k_input = raw_input('what to do? 0: do nothing and exit --- 1: clean duplicate lines      :\t')
print k_input

if k_input is '1':
	file_to_clean = raw_input('what file do you want to clean? \t')
	to_clean      = np.genfromtxt( os.path.join(path, file_to_clean) , unpack=True)
	cleaned       = to_clean[:,to_save]
	print cleaned.shape
	np.savetxt( os.path.join(path, os.path.splitext(file_to_clean)[0]+'_cleaned'+os.path.splitext(file_to_clean)[1] ), cleaned )
	
# 	fig = plt.figure(1)
# 	ax  = plt.subplot(111) #matplotlib.pyplot.subplot(111)
# 	ax.plot( cleaned[0,:], cleaned[7,:], 'o') #, markersize=5.0)
# 	ax.grid()
# 	plt.show()



















