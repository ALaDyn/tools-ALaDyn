#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_phasespace.py
# Author:       F. Mira
# Date:         2016-11-03
# Purpose:      reads ALaDyn binary from PDBunch* and plots phasespace
# Source:       python
#####################################################################
import os, sys, struct
from scipy import *
import numpy as np
import pylab as pyl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
###>>>
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/ALaDyn_Code/tools-ALaDyn/pythons'))
###>>>
### --- ###
from read_particle_phasespace import *
from utilities_1 import *
### --- ###


# --- # --- # --- #
### --- ### shell inputs
if(len(sys.argv)<5):
	print 'Usage: python ~/Codes/tools-ALaDyn/ALaDyn_Pythons/ALaDyn_plot_phasespace.py 0 9 30'
	print 'frame begin'
	print 'frame end'
	print 'selection via Weights: -1 for all particles'
	print 'gamma cut: only particle with gamma greater then: -1 for all particles'
	exit('---not enough arguments---')
#---
frame_begin = int(sys.argv[1])
frame_end   = int(sys.argv[2])
W_threshold = float(sys.argv[3])
gamma_threshold   = float(sys.argv[4])

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

		X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
		number_of_particles = len(X); del X
		gamma_selected = np.full((number_of_particles), True, dtype=bool)
		W_selected = np.full((number_of_particles), True, dtype=bool)

		if(gamma_threshold>-1.):
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
			del gamma

		if(W_threshold>-1.):
			W = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','W')
			W_selected = (W<=W_threshold)

		p_selected = (gamma_selected & W_selected)

		Px = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Px')
		Px = Px[p_selected]
		Py = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Py')
		Py = Py[p_selected]
		Pz = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Pz')
		Pz = Pz[p_selected]
		X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
		X = X[p_selected]
		Y = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Y')
		Y = Y[p_selected]
		Z = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Z')
		Z = Z[p_selected]

		gamma = np.array(np.sqrt(1. + Px**2 + Py**2 + Pz**2))

		en_spread = np.std(gamma)/np.mean(gamma)
		emittance_y = np.sqrt( np.std(Y)**2*np.std(Py)**2-np.cov(Y,Py)[0][1]**2)/np.mean(gamma);
		emittance_z = np.sqrt( np.std(Z)**2*np.std(Pz)**2-np.cov(X,Pz)[0][1]**2)/np.mean(gamma);
		Charge = 1536.*0.02*1.6e-7*len(X);

		print 'Energy spread: ', en_spread*100 ,'%'
		print 'Normalized Emittance Y: ', ('%3.2e' % emittance_y), 'mm-mrad'
		print 'Normalized Emittance Z: ', ('%3.2e' % emittance_z), 'mm-mrad'
		print 'Charge:',('%3.2e' % Charge), 'pC'

		plt.figure()
		#fig=plt.figure()
		#ax=fig.add_subplot(111,projection='3d')
		#ax.scatter(X,Y,Z,s=.05,edgecolors='None')
		#ax.set_xlabel('X')
		#ax.set_xlabel('Y')
		#ax.set_xlabel('Z')
		plt.subplot(311)
		plt.scatter(X,Px,s=.1,edgecolors='None')
		#plt.scatter(X,Y,Z,s=.1,edgecolors='None')
		#plt.xlabel(r'$X (\mu m) $');plt.ylabel(r'$Z (\mu m)$')
		plt.xlabel(r'$X (\mu m) $',fontsize=24);plt.ylabel(r'$P_x/mc$',fontsize=24)
		#plt.text(90,327, r'$Energy spread = 20,2% $',size=18,ha='left',va='top')
		plt.subplot(312)
		#       plt.scatter(Y,Z,s=.1,edgecolors='None')
		#       plt.xlabel(r'$Y (\mu m) $');plt.ylabel(r'$Z (\mu m)$')
		plt.scatter(Y,Py,s=.1,edgecolors='None')
		plt.xlabel(r'$Y (\mu m)$',fontsize=24);plt.ylabel(r'$P_y/mc$',fontsize=24)
		plt.subplot(313)
		#       plt.scatter(X,Y,s=.1,edgecolors='None')
		#      plt.xlabel(r'$X (\mu m) $');plt.ylabel(r'$Y (\mu m)$')
		#plt.scatter(Pz,Py,s=.1,edgecolors='None')
		#plt.xlabel('Pz/mc');plt.ylabel('Py/mc')
		#plt.subplot(313)
		plt.hist(gamma*0.511,100)
		#       plt.hist(gamma[gamma_selected],50)
		plt.xlabel('$Energy (MeV)$',fontsize=24);plt.ylabel('$dN/dE$',fontsize=24);
		#       plt.axis('tight')
		name_output = 'Phsp_enspect'+'_'+('%2.2i'%i)+'.png'
		#name_output =  '3Dplot' + '_' +('%2.2i'%i)+ '.png'
		plt.savefig( os.path.join(path_write,'data','phasespace',name_output) )
		# plt.close()
		plt.show()
