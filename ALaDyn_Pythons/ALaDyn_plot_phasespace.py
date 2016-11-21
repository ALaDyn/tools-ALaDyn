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
	print 'Usage: python ~/Codes/tools-ALaDyn/ALaDyn_Pythons/ALaDyn_plot_phasespace.py 0 9 30'
	print 'frame begin: '
	print 'frame end: '
	print 'gamma cut: only particle with gamma greater then'
	exit(0)
#---
frame_begin = int(sys.argv[1])
frame_end   = int(sys.argv[2])
gamma_threshold   = float(sys.argv[3])

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

			Px = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Px')
			gamma = 1.+Px**2
			del Px

			Py = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Py')
			gamma = gamma+Py**2
			del Py

			Pz = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Pz')
			gamma = gamma+Pz**2
			del Pz

			gamma=np.array(np.sqrt(gamma))
			gamma_selected = (gamma>gamma_threshold)
		
   			en_spread = np.std(gamma[gamma_selected])/np.mean(gamma[gamma_selected])

			Px = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Px')
			Px = Px[gamma_selected]

			Py = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Py')
			Py = Py[gamma_selected]

			Pz = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Pz')
			Pz = Pz[gamma_selected]

			X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
			X = X[gamma_selected]

			Y = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Y')
			Y = Y[gamma_selected]

			Z = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Z')
			Z = Z[gamma_selected]

			W = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','W')
			W = W[gamma_selected]
			
			en_spread = np.std(gamma[gamma_selected])/np.mean(gamma[gamma_selected])
			
			emittance_y = np.sqrt( np.std(Y)**2*np.std(Py)**2-np.cov(Y,Py)[0][1]**2);
			emittance_z = np.sqrt( np.std(Z)**2*np.std(Pz)**2-np.cov(X,Pz)[0][1]**2);

			print 'Energy spread: ', en_spread ,'%'
			print 'Normalized Emittance Y: ', ('%3.2e' % emittance_y), 'mm-mrad'
			print 'Normalized Emittance Z: ', ('%3.2e' % emittance_z), 'mm-mrad'
			print 'Charge:',('%3.2e' % (sum(W)*0.02*3.072e-17) )

			plt.figure(1, figsize=(11.0,11.0) )
			plt.subplot(311)
			plt.scatter(X,Px,s=.1,edgecolors='None')
			plt.xlabel(r'X $\mu m$');plt.ylabel('Px/mc^2');
			plt.subplot(312)
			plt.scatter(Y,Py,s=.1,edgecolors='None')
			plt.xlabel('Y');plt.ylabel('Py/mc^2');
			plt.subplot(313)
			plt.hist(gamma[gamma_selected]*0.511,50)						
	        	plt.xlabel('Energy');plt.ylabel('dN/dE');
			# plt.xlabel(component1); plt.ylabel(component2)
			plt.axis('tight')
			# name_output =  component1 +'_'+ component2 + '_' +('%2.2i'%i)+ '.png'
			# plt.savefig( os.path.join(path_write,'data','phasespace',name_output) )
			# plt.close()
			plt.show()
