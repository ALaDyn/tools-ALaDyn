#!/usr/bin/python
######################################################################
# Name:         read_ALaDyn_Particles.py
# Author:       C Gatti and A Marocchino      
# Date:            2014-06-11
# Purpose:      reads dued binary from output
# Source:       python
#####################################################################

### loading shell commands
import os,sys
import struct
from scipy import *
import numpy as np
# from matplotlib import *
# from pylab import *
import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
### --- ###



# - #


if len(sys.argv)!=2:
    print "Usage ./read_ALadyn_Particles.py  filename.bin"
    exit(0)
    pass

file_name = sys.argv[1]


f  = open(file_name,'rb')

X=[]  
Y=[] 
Z=[]
Px=[]
Py=[]
Pz=[]
W=[]

while True:
    try:
        Np=struct.unpack('i', f.read(4))
    except struct.error:
        print "End of File"
        break
    vars=[]
    for i in range(0,Np[0]):
        vars=struct.unpack('fffffff', f.read(4*7))
        X.append(vars[0])
        Y.append(vars[1])
        Z.append(vars[2])
        Px.append(vars[3])
        Py.append(vars[4])
        Pz.append(vars[5])
        W.append(vars[6])
        pass
    pass
f.close()


H, xedges, yedges = np.histogram2d(Y, Py, bins=(120, 120))

Xhist,Xbins=np.histogram(X,100)
print Xhist
print Xbins

#use interactive mode
plt.ion()
#plt.hist(X)
#plt.plot(Y,Py,'.')i
plt.imshow(H)
plt.draw()


while True:
    pass





