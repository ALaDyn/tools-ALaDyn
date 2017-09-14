#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_phasespace.py
# Author:       A. Marocchino, F. Mira, N. Panzeri
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
if(len(sys.argv)<10):
	print 'Usage :: python ~/Codes/tools-ALaDyn/ALaDyn_Pythons/plot_phasespace.py 0 9 -1 -1 0 1000 10 1 5'
	print 'Not enough arguments >>>'
	print ''
	print 'frame begin'
	print 'frame end'
	print 'selection via Weights: all particles with weight greater then INput :: -1 for all particles'
	print 'selection via Weights: all particles with weight less then INput :: -1 for all particles'
	print 'gamma cut: only particle with gamma greater then INput :: -1 for all particles'
	print 'longitudinal min cut: remove particle with X less than INput :: -1 for no cut'
	print 'longitudinal max cut: remove particle with X greater than INput  :: -1 for no cut'
	print 'transverse cut: remove particle with R greater than INput  :: -1 for no cut'
	print 'particle Jump: one every jump particle'
	print 'slice analysis: slice longitudinal dimension in um'
	print ''
	exit('---not enough arguments---')
#---
frame_begin = int(sys.argv[1])
frame_end   = int(sys.argv[2])
Wmax_threshold = float(sys.argv[3])
Wmin_threshold = float(sys.argv[4])
gamma_threshold   = float(sys.argv[5])
Xmin_threshold = float(sys.argv[6])
Xmax_threshold = float(sys.argv[7])
transverse_threshold = float(sys.argv[8])
jump = int(sys.argv[9])
slice_um  = float(sys.argv[10])

#-path
path = os.getcwd()
#-folder output structure
generate_folder_phasespace(path)

#--- ***  mining and printing ***---#
#--- if frm_begin == frm_end : run and plot ---#
if( frame_begin == frame_end ):
# for i in range(frame_begin, frame_end + 1 ):
	print '-------------------'
	print 'frame_begin == frame_end  :: analysis and plot'
	i=frame_begin
	s='%2.2i'%i
	path_read  = os.path.join(path,'%4.4i'%i)
	path_write = path
	if output_exists(path_read,'phasespace',i) == True:
		print 'Analysing Phase Space --- frame >>> ',i

		X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
		number_of_particles = len(X); del X
		p_selected = np.full((number_of_particles), True, dtype=bool)
		# gamma_selected = np.full((number_of_particles), True, dtype=bool)
		# W_selected = np.full((number_of_particles), True, dtype=bool)
		# X_selected = np.full((number_of_particles), True, dtype=bool)
		# R_selected = np.full((number_of_particles), True, dtype=bool)

		#--- jump ---#
		p_selected[0:jump:number_of_particles]=False

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
			p_selected = (gamma>gamma_threshold)
			del gamma

		if(Wmin_threshold>-1.):
			W = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','W')
			p_selected = (p_selected & (W<=Wmin_threshold) )

		if(Wmax_threshold>-1.):
			W = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','W')
			p_selected = (p_selected & (W>=Wmax_threshold) )

		if(Xmin_threshold> -1):
			X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
			p_selected = (p_selected & (X>=Xmin_threshold) )
			del X
		if(Xmax_threshold> -1):
			X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
			p_selected = (p_selected & (X<=Xmax_threshold) )
			del X

		if(transverse_threshold > -1):
			Y = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Y')
			p_selected = (p_selected & (abs(Y)<=transverse_threshold) )
			del Y
			Z = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Z')
			p_selected = (p_selected & (abs(Z)<=transverse_threshold) )
			del Z

		# --- *** --- #
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

		#--- *** diagnostics *** ---#
		Charge  =  read_electrons_ppm( path_read,'Elpout'+s+'.dat')*np.sum(W)*1.6e-7
		gamma = np.array(np.sqrt(1. + Px**2 + Py**2 + Pz**2))
		mu_x=np.average(X,weights=W)
		mu_y=np.average(Y,weights=W)
		mu_z=np.average(Z,weights=W)
		mu_Px=np.average(Px,weights=W)
		mu_Py=np.average(Py,weights=W)
		mu_Pz=np.average(Pz,weights=W)
		mu_gamma=np.average(gamma,weights=W)
		sigma_x=np.sqrt( np.average((X-mu_x)**2, weights=W) )
		sigma_y=np.sqrt( np.average((Y-mu_y)**2, weights=W) )
		sigma_z=np.sqrt( np.average((Z-mu_z)**2, weights=W) )
		sigma_Px=np.sqrt( np.average((X-mu_Px)**2, weights=W) )
		sigma_Py=np.sqrt( np.average((Y-mu_Py)**2, weights=W) )
		sigma_Pz=np.sqrt( np.average((Z-mu_Pz)**2, weights=W) )
		sigma_gamma=np.sqrt( np.average((gamma-mu_gamma)**2, weights=W) )
		cov_y_Py=np.average((Y-mu_y)*(Py-mu_Py), weights=W)
		cov_z_Pz=np.average((Z-mu_z)*(Pz-mu_Pz), weights=W)

		en_spread = sigma_gamma/mu_gamma*100.
		emittance_y = np.sqrt( sigma_y**2*sigma_Py**2-cov_y_Py**2);
		emittance_z = np.sqrt( sigma_z**2*sigma_Pz**2-cov_z_Pz**2);
		Current = Charge*1e-15*3e8/np.sqrt(2*np.pi)/sigma_x/1e-6
		#--- *** ---#
		print 'Diagnostic for the whole bunch'
		print 'Energy spread: ', ('%3.2e' % en_spread) ,'%'
		print 'mean energy :',mu_gamma
		print 'sigma_energy:',sigma_gamma
		print 'Normalized Emittance Y: ', (round(emittance_y,2)), 'mm-mrad'
		print 'Normalized Emittance Z: ', (round(emittance_z,2)), 'mm-mrad'
		print 'Mean Energy: ',  ('%3.2e' % (mu_gamma*0.511)), 'MeV'
		print 'Charge:',('%3.2e' % Charge), 'pC'
		print 'Peak Current:',('%3.2e' % Current), 'kA'
		print 'longitudinal - Sigmax:',('%3.2e' % sigma_x),'mum'
		print 'transverse   - Sigmay:',('%3.2e' % sigma_y),'mum'

		plt.figure(figsize=(10,5))
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
		plt.hist(gamma*0.511,300,weights=W)
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

		# --- Second Plot :: PS scatter plot---#
		plt.figure(figsize=(10,5))
		plt.subplot(211)
		plt.scatter(X,Y,s=.1,edgecolors='None')
		plt.xlabel(r'$X (\mu m) Longitudinal$',fontsize=8)
		plt.ylabel(r'$Y (\mu m) Transverse$',fontsize=8)

		plt.subplot(212)
        	plt.scatter(Y,Z,s=.1,edgecolors='None')
		plt.xlabel(r'$Y (\mu m) Transverse$',fontsize=8)
		plt.ylabel(r'$Z (\mu m) Transverse$',fontsize=8)
		plt.show()



		#--- *** Slice diagnostics ***---#
		first_particle_x=np.min(X)
		last_particle_x=np.max(X)
		slice_x_min=first_particle_x
		slice_x_max=slice_x_min+slice_um
		slice_counter=0
		print 'slice_counter,Charge, sigma_x,sigma_y,gamma-mean,energy_spread, emittance_y, emittance_z'
		while (last_particle_x>slice_x_max):
			X_selected = (X>slice_x_min & X<slice_x_max)

			Charge  =  read_electrons_ppm(path_read,'Elpout'+s+'.dat')*np.sum(W[X_selected])*1.6e-7
			sigmax  = np.std(X[X_selected])
			sigmay  = np.std(Y[X_selected])
			Current = Charge*1e-15*3e8/np.sqrt(2*np.pi)/slice_um/1e-6
			gamma = np.array(np.sqrt(1. + Px**2 + Py**2 + Pz**2))
			en_spread = round(np.std(gamma[X_selected])/np.mean(gamma[X_selected])*100,1)
			emittance_y = np.sqrt( np.std(Y[X_selected])**2*np.std(Py[X_selected])**2-np.cov(Y[X_selected],Py[X_selected])[0][1]**2);
			emittance_z = np.sqrt( np.std(Z[X_selected])**2*np.std(Pz[X_selected])**2-np.cov(X[X_selected],Pz[X_selected])[0][1]**2);

			print slice_counter,Charge,sigmax,sigmay,np.mean(gamma[X_selected]),en_spread,emittance_y,emittance_z

			slice_x_min+=slice_um
			slice_x_max+=slice_um
			slice_counter+=1

	else:
		print '--- --- ---'

#--- full analysis :: loop from frm_begin to frm_end ---#
if(frame_end > frame_begin):
	print '-------------------'
	print 'full analysis :: loop from frame_begin to frame_end'

	for i in range(frame_begin, frame_end + 1 ):
		s='%2.2i'%i
		path_read  = os.path.join(path,'%4.4i'%i)
		path_write = path
		if output_exists(path_read,'phasespace',i) == True:
			print 'Analysing Phase Space --- frame >>> ',i
			X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
			number_of_particles_before_selection = len(X); del X
			p_selected = np.full((number_of_particles_before_selection), True, dtype=bool)

			#--- jump ---#
			if(number_of_particles_before_selection>1):
				p_selected[0:jump:number_of_particles_before_selection]=False

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
					p_selected = (gamma>gamma_threshold)
					del gamma

				if(Wmin_threshold>-1.):
					W = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','W')
					p_selected = (p_selected & (W<=Wmin_threshold) )

				if(Wmax_threshold>-1.):
					W = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','W')
					p_selected = (p_selected & (W>=Wmax_threshold) )

				if(Xmin_threshold> -1):
					X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
					p_selected = (p_selected & (X>=Xmin_threshold) )
					del X
				if(Xmax_threshold> -1):
					X = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','X')
					p_selected = (p_selected & (X<=Xmax_threshold) )
					del X

				if(transverse_threshold > -1):
					Y = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Y')
					p_selected = (p_selected & (abs(Y)<=transverse_threshold) )
					del Y
					Z = read_particle_phasespace_bycomponent( path_read,'Elpout'+s+'.bin','Z')
					p_selected = (p_selected & (abs(Z)<=transverse_threshold) )
					del Z

				# --- *** --- #
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

				#--- *** diagnostics *** ---#
				number_of_particles_after_selection = len(W)
				if( number_of_particles_after_selection > 1 ):
					Charge  =  read_electrons_ppm( path_read,'Elpout'+s+'.dat')*np.sum(W)*1.6e-7
					gamma = np.array(np.sqrt(1. + Px**2 + Py**2 + Pz**2))
					mu_x=np.average(X,weights=W)
					mu_y=np.average(Y,weights=W)
					mu_z=np.average(Z,weights=W)
					mu_Px=np.average(Px,weights=W)
					mu_Py=np.average(Py,weights=W)
					mu_Pz=np.average(Pz,weights=W)
					mu_gamma=np.average(gamma,weights=W)
					sigma_x=np.sqrt( np.average((X-mu_x)**2, weights=W) )
					sigma_y=np.sqrt( np.average((Y-mu_y)**2, weights=W) )
					sigma_z=np.sqrt( np.average((Z-mu_z)**2, weights=W) )
					sigma_Px=np.sqrt( np.average((Px-mu_Px)**2, weights=W) )
					sigma_Py=np.sqrt( np.average((Py-mu_Py)**2, weights=W) )
					sigma_Pz=np.sqrt( np.average((Pz-mu_Pz)**2, weights=W) )
					sigma_gamma=np.sqrt( np.average((gamma-mu_gamma)**2, weights=W) )
					cov_y_Py=np.average((Y-mu_y)*(Py-mu_Py), weights=W)
					cov_z_Pz=np.average((Z-mu_z)*(Pz-mu_Pz), weights=W)

					en_spread = sigma_gamma/mu_gamma*100.
					emittance_y = np.sqrt( sigma_y**2*sigma_Py**2-cov_y_Py**2);
					emittance_z = np.sqrt( sigma_z**2*sigma_Pz**2-cov_z_Pz**2);
					Current = Charge*1e-15*3e8/np.sqrt(2*np.pi)/sigma_x/1e-6

			if(number_of_particles_before_selection <= 1 or number_of_particles_after_selection <= 1):
				Charge  = 0.0; gamma   = 0.0
				mu_x    = 0.0; mu_y    = 0.0; mu_z    = 0.0
				mu_Px   = 0.0; mu_Py   = 0.0; mu_Pz   = 0.0
				mu_gamma= 0.0; sigma_x = 0.0; sigma_y = 0.0; sigma_z = 0.0
				sigma_Px= 0.0; sigma_Py= 0.0; sigma_Pz= 0.0; sigma_gamma= 0.0
				cov_y_Py = 0.0;	cov_z_Pz = 0.0; en_spread = 0.0
				emittance_y = 0.0; emittance_z = 0.0; Current = 0.0

			#--- *** ---#
			run_distance = read_run_distance( path_read,'Elpout'+s+'.dat')
			f = open('PS_diagnostic.dat','a')
			f.write( '%3i \t %5.4e \t %5.4e \t %5.4e \t %5.4e \t %5.4e \t %5.4e \t %5.4e \t %5.4e \t %5.4e \t %5.4e \t %5.4e\n' % (i,run_distance,sigma_x,sigma_y,sigma_z,emittance_y,emittance_z,mu_gamma*0.511,en_spread,Charge,sigma_Py,cov_y_Py) )
			f.close()
			# 1 :: frame
			# 2 :: run_distance
			# 3,4,5 :: sigma_x,sigma_y,sigma_z
			# 6,7 :: emittance_y,emittance_z
			# 8 :: energy in MeV, mu_gamma*0.511
			# 9 :: en_spread
			# 10 :: Charge in pC
