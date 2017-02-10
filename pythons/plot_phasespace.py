#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_phasespace.py
# Author:       F. Mira, A. Marocchino
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
if(len(sys.argv)<9):
	print 'Usage: python ~/Codes/tools-ALaDyn/ALaDyn_Pythons/ALaDyn_plot_phasespace.py 0 9 30'
	print 'frame begin'
	print 'frame end'
	print 'selection via Weights: -1 for all particles'
	print 'gamma cut: only particle with gamma greater then: -1 for all particles'
	print ' longitudinal threshold'
	print 'transverse threshold '
	print 'select only "one over jump" particles'
	print 'NUMBER OF SLICES'
	exit('---not enough arguments---')
#---
frame_begin = int(sys.argv[1])
frame_end   = int(sys.argv[2])
W_threshold = float(sys.argv[3])
gamma_threshold   = float(sys.argv[4])
X_threshold = float(sys.argv[5])
transverse_threshold = float(sys.argv[6])
jump = int(sys.argv[7])
number_slices     = int(sys.argv[8])

#-path
path = os.getcwd()
#-folder output structure
generate_folder_phasespace(path)

emittance     = np.zeros((frame_end+1-frame_begin,number_slices))
energy_spread = np.zeros((frame_end+1-frame_begin,number_slices))
sigma_trans   = np.zeros((frame_end+1-frame_begin,number_slices))
mean_energy   = np.zeros((frame_end+1-frame_begin,number_slices))
distance      = np.zeros((frame_end+1-frame_begin,number_slices))
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
		if(transverse_threshold > 0.):
			Y = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Y')
			Y_selected = (abs(Y)<=transverse_threshold)
			Z = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Z')
			Z_selected = (abs(Z)<=transverse_threshold)
		if(X_threshold>0.):
			X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
			X_selected = (X>=X_threshold)
		p_selected = (gamma_selected & W_selected & Y_selected & Z_selected & X_selected)

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
		W = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','W')
		W = W[p_selected]

                #Charge  = 1536.*0.01*1.6e-7*len(X);
                Charge  =  1536.*np.sum(W)*1.6e-7
		sigmax  = np.std(X)
		sigmar  = np.std(Y)
		Current = Charge*1e-15*3e8/np.sqrt(2*pi)/sigmax/1e-6
		if (jump > 0):
			Px = Px[0:len(Px):jump]
                        Py = Py[0:len(Py):jump]
                        Pz = Pz[0:len(Pz):jump]
                        X  =  X[0:len(X) :jump]
                        Y  =  Y[0:len(Y) :jump]
                        Z  =  Z[0:len(Z) :jump]
                        W  =  W[0:len(W) :jump]


		gamma = np.array(np.sqrt(1. + Px**2 + Py**2 + Pz**2))

		en_spread = round(np.std(gamma)/np.mean(gamma)*100,1)
		emittance_y = np.sqrt( np.std(Y)**2*np.std(Py)**2-np.cov(Y,Py)[0][1]**2);
		emittance_z = np.sqrt( np.std(Z)**2*np.std(Pz)**2-np.cov(X,Pz)[0][1]**2);

		print 'Whole bunch'
		print 'Energy spread: ', en_spread ,'%'
		print 'Normalized Emittance Y: ', (round(emittance_y,2)), 'mm-mrad'
		print 'Normalized Emittance Z: ', (round(emittance_z,2)), 'mm-mrad'
		print 'Mean Energy: ', round(np.mean(gamma)*0.51,1), 'MeV'
		print 'Charge:',('%3.2e' % Charge), 'pC'
		print 'Peak Current:',('%3.2e' % Current), 'kA'
		print 'Sigmax:',('%3.2e' % sigmax),'mum'
		print 'Sigmar:',('%3.2e' % sigmar),'mum'
		plt.figure(figsize=(20,10))
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
		plt.xlabel(r'$X (\mu m) $',fontsize=24,fontweight='bold');plt.ylabel(r'$P_x/mc$',fontsize=24,fontweight='bold')
		#plt.text(90,327, r'$Energy spread = 20,2% $',size=18,ha='left',va='top')
		plt.subplot(312)
		#       plt.scatter(Y,Z,s=.1,edgecolors='None')
		#       plt.xlabel(r'$Y (\mu m) $');plt.ylabel(r'$Z (\mu m)$')
#		plt.scatter(Y,Py,s=.1,edgecolors='None')
#		plt.xlabel(r'$Y (\mu m)$',fontsize=24,fontweight='bold');plt.ylabel(r'$P_y/mc$',fontsize=24,fontweight='bold')
                plt.scatter(Y,Py,s=.1,edgecolors='None')
                plt.xlabel(r'$Y (\mu m)$',fontsize=24,fontweight='bold');plt.ylabel(r'$P_y/mc$',fontsize=24,fontweight='bold')
		plt.subplot(313)
		#       plt.scatter(X,Y,s=.1,edgecolors='None')
		#      plt.xlabel(r'$X (\mu m) $');plt.ylabel(r'$Y (\mu m)$')
		#plt.scatter(Pz,Py,s=.1,edgecolors='None')
		#plt.xlabel('Pz/mc');plt.ylabel('Py/mc')
		#plt.subplot(313)
		plt.hist(gamma*0.511,300)
		#       plt.hist(gamma[gamma_selected],50)
		plt.xlabel('$Energy (MeV)$',fontsize=24,fontweight='bold');plt.ylabel('$dN/dE$',fontsize=24,fontweight='bold');
		#       plt.axis('tight')
		plt.subplots_adjust(left=0.08,right=0.8,hspace=0.45)
		name_output = 'Phsp_enspect'+'_'+('%2.2i'%i)+'.png'
		#np.savetxt('X_'+('%2.2i'%i)+'.txt',X)
                #np.savetxt('Px_'+('%2.2i'%i)+'.txt',Px)
		#name_output =  '3Dplot' + '_' +('%2.2i'%i)+ '.png'
		plt.savefig( os.path.join(path_write,'data','phasespace',name_output) )
		#plt.close()
		plt.show()

		for j in range (0,number_slices):

			Y_selected = ((Y >= np.mean(Y)-(j+1)*np.std(Y))&(Y <= np.mean(Y)+(j+1)*np.std(Y)))
			Y_slice  = Y[Y_selected]
			Px_slice = Px[Y_selected]
			Py_slice = Py[Y_selected]
			Pz_slice = Pz[Y_selected]
			X_slice = X[Y_selected]
			Z_slice = Z[Y_selected]
			W_slice = W[Y_selected]

			gamma_slice = np.array(np.sqrt(1. + Px_slice**2 + Py_slice**2 + Pz_slice**2))
			en_spread_slice = round(np.std(gamma_slice)/np.mean(gamma_slice)*100,1)
			emittance_y_slice = np.sqrt( np.std(Y_slice)**2*np.std(Py_slice)**2-np.cov(Y_slice,Py_slice)[0][1]**2);
			emittance_z_slice = np.sqrt( np.std(Z_slice)**2*np.std(Pz_slice)**2-np.cov(X_slice,Pz_slice)[0][1]**2);
		        sigmax_slice  = np.std(X_slice)
	                sigmar_slice  = np.std(np.sqrt(Y_slice**2+Z_slice**2))
	                Charge_slice  =  1536.*np.sum(W_slice)*1.6e-7
			Current_slice = Charge_slice*1e-15*3e8/np.sqrt(2*pi)/sigmax_slice/1e-6


	                print 'Particles in',j,' sigma'
        	        print 'Energy spread: ', en_spread_slice ,'%'
        	        print 'Normalized Emittance Y: ', (round(emittance_y_slice,2)), 'mm-mrad'
        	        print 'Normalized Emittance Z: ', (round(emittance_z_slice,2)), 'mm-mrad'
        	        print 'Mean Energy: ', round(np.mean(gamma_slice)*0.51,1), 'MeV'
        	        print 'Charge:',('%3.2e' % Charge_slice), 'pC'
        	        print 'Peak Current:',('%3.2e' % Current_slice), 'kA'
        	        print 'Sigmax:',('%3.2e' % sigmax_slice),'mum'
        	        print 'Sigmar:',('%3.2e' % sigmar_slice),'mum'

 	                plt.figure(figsize=(20,10))
        	        plt.subplot(311)
               		plt.scatter(X_slice,Px_slice,s=.1,edgecolors='None')
      			plt.xlabel(r'$X (\mu m) $',fontsize=24,fontweight='bold');plt.ylabel(r'$P_x/mc$',fontsize=24,fontweight='bold')
                        plt.subplot(312)
                        plt.scatter(Y_slice,Py_slice,s=.1,edgecolors='None')
                        plt.xlabel(r'$Y (\mu m)$',fontsize=24,fontweight='bold');plt.ylabel(r'$P_y/mc$',fontsize=24,fontweight='bold')
                        plt.subplot(313)
                        plt.hist(gamma_slice*0.511,300)
                        plt.xlabel('$Energy (MeV)$',fontsize=24,fontweight='bold');plt.ylabel('$dN/dE$',fontsize=24,fontweight='bold');
                        plt.subplots_adjust(left=0.08,right=0.8,hspace=0.45)

			slice_output = 'Phsp_slice'+'_'+('%2.2i'%j)+'.png'
         	        plt.savefig( os.path.join(path_write,'data','phasespace',slice_output) )

			plt.show()

			energy_spread[i-frame_begin,j] = en_spread
			emittance[i-frame_begin,j]     = emittance_y
			sigma_trans[i-frame_begin,j]   = sigmar
			mean_energy[i-frame_begin,j]   = np.mean(gamma)*0.511
			distance[i-frame_begin,j]      = np.mean(X)


	else:
		print'not working'


#	energy_spread[i-frame_begin] = en_spread
#	emittance[i-frame_begin]     = emittance_y
#	sigma_trans[i-frame_begin]   = sigmar
#	mean_energy[i-frame_begin]   = np.mean(gamma)*0.511
#	distance[i-frame_begin]      = np.mean(X)


np.savetxt(os.path.join(path_write,'data','phasespace',('Energy_spread.txt')),energy_spread)
np.savetxt(os.path.join(path_write,'data','phasespace',('Emittance.txt')),emittance)
np.savetxt(os.path.join(path_write,'data','phasespace',('Transverse_sigma.txt')),sigma_trans)
np.savetxt(os.path.join(path_write,'data','phasespace',('Mean_energy.txt')),mean_energy)
np.savetxt(os.path.join(path_write,'data','phasespace',('Distance.txt')),distance)
