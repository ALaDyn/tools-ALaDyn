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


### --- ### shell inputs
if(len(sys.argv)<2):
        print 'Usage:'
        print 'Input [1]: frame number'
        exit(0)
#---
frame = int(sys.argv[1])



# - #
s='%2.2i'%frame
X1, Y1, Z1, Px1, Py1, Pz1, W1 = read_particle_phasespace( os.path.join(os.getcwd(),'PSBunch1_'+s+'.bin') )
X1=np.array(X1);    Y1=np.array(Y1);   Z1=np.array(Z1)
Px1=np.array(Px1); Py1=np.array(Py1); Pz1=np.array(Pz1)

X2, Y2, Z2, Px2, Py2, Pz2, W2 = read_particle_phasespace( os.path.join(os.getcwd(),'PSBunc\
h2_'+s+'.bin') )
X2=np.array(X2);    Y2=np.array(Y2);   Z2=np.array(Z2)
Px2=np.array(Px2); Py2=np.array(Py2); Pz2=np.array(Pz2)

#-overlapping particle counter
lb=min(X2)
rb=max(X2)
counter=0
for i in np.arange(0,len(X1)):
	if X1[i] > lb and X1[i] < rb:
		counter=counter+1


fig = plt.figure(1, figsize=(11.0,11.0) )
ax1 = fig.add_subplot(111)
ax1.scatter(X1,Px1,s=.1,color='blue',edgecolors='None')
ax1.scatter(X2,Px2,s=.1,color='magenta',edgecolors='None')
ax1.plot(max(X2),Px2[X2.argmax()],marker='o', color='green',markersize=7.5)
ax1.plot(min(X1),Px1[X1.argmin()],marker='o', color='black',markersize=7.5)
plt.xlabel('overlapping particles %d'%counter)
plt.savefig( os.path.join(os.getcwd(),'plots','phasespace','overlapping_X_'+s+'.eps') )
plt.savefig( os.path.join(os.getcwd(),'plots','phasespace','overlapping_X_'+s+'.png') )
plt.close(fig)

