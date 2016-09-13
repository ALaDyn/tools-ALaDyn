#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:		A. Marocchino
# Date:			2014-11-04
# Purpose:      it nests into 'ALaDyn_read_binary MATRIX' it does it using much less memory than the previous version
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil
###>>>
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/ALaDyn_Pythons'))
###>>>
### --- ###
from read_ALaDyn_bin import *
from ALaDyn_plot_utilities_1 import *
from ALaDyn_plot_utilities_density import *
from ALaDyn_plot_utilities_Efield import *
from ALaDyn_plot_utilities_Bfield import *
from ALaDyn_plot_utilities_ionization import *
from ALaDyn_plot_utilities_energy_density import *
from ALaDyn_plot_utilities_axes import *


### --- ###
# def read_ALaDyn_bin_section(dir_path,file_name,grid_no_grid,axis_to_cut,cell_to_cut):


### --- ### shell inputs
if(len(sys.argv)<2):
	print 'Input  [1]: frame_begin'
	print 'Input  [2]: frame_end'
	print 'Input  [3]: axis to cut'
	print 'Input  [4]: slice position (cell number)'

	exit(0)

if sys.argv[1] == -1:
	frame_begin		  = 0
	frame_end         = last_output(os.getcwd())
	axis_to_cut		  = 'y'
	cell_to_cut		  = 0
else:
	frame_begin 		= int(		sys.argv[1])
	frame_end			= int(		sys.argv[2])
	axis_to_cut			= str(		sys.argv[3])
	cell_to_cut			= int(		sys.argv[4])
### --- ###



#--- *** ---#
if __name__ == '__main__':

	#-path
	path = os.getcwd()

	#-folder output structure
	generate_folder_output_structure(path,'True')


	for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
		print '-------------------'
		s='%2.2i'%i
		path_read  = os.path.join(path,'%4.4i'%i)
		path_write = path

		if output_exists(path_read,'rho',i) == True:
			print 'rho --- frame >>> ',i
			M,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Bdenout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			K,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Edenout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			np.savetxt( os.path.join(path_write,'data','rho',('rho_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),M.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','rho',('rho_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),np.abs(K.transpose()),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','rho',('rho_tot_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),np.abs(M.transpose())+np.abs(K.transpose()),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')

		if output_exists(path_read,'ionization',i) == True:
			print 'ionization rate --- frame >>> ',i
			I,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'H1dnout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			np.savetxt( os.path.join(path_write,'data','ionization',('ionization_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),I.transpose(),fmt='%15.14e')

		if output_exists(path_read,'Energy_density',i) == True:
            		print 'energy density --- frame >>> ',i
            		ED,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Elenout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
            		np.savetxt( os.path.join(path_write,'data','EneDen',('Energy_Density'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),I.transpose(),fmt='%15.14e')

		if output_exists(path_read,'E',i) == True:
			print 'E --- frame >>> ',i
			Exb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Exbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Exf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Exfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			Eyb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Eybout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Eyf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Eyfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			Ezb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Ezbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Ezf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Ezfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			np.savetxt( os.path.join(path_write,'data','E_field',('Ex_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Exb.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','E_field',('Ey_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Eyb.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','E_field',('Ez_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Ezb.transpose(),fmt='%15.14e')

			np.savetxt( os.path.join(path_write,'data','E_field',('Ex_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Exf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','E_field',('Ey_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Eyf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','E_field',('Ez_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Ezf.transpose(),fmt='%15.14e')

			np.savetxt( os.path.join(path_write,'data','E_field',('Ex_tot_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Exb.transpose()+Exf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','E_field',('Ey_tot_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Eyb.transpose()+Eyf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','E_field',('Ez_tot_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Ezb.transpose()+Ezf.transpose(),fmt='%15.14e')

			#np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
			#np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
			#np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')


		if output_exists(path_read,'B',i) == True:
			print 'B --- frame >>> ',i
			#Bxb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Bxbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Bxf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Bxfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			Byb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Bybout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Byf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Byfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			Bzb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Bzbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Bzf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Bzfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			#np.savetxt( os.path.join(path_write,'data','B_field',('Bx_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bxf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','B_field',('By_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Byb.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','B_field',('Bz_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bzb.transpose(),fmt='%15.14e')

			np.savetxt( os.path.join(path_write,'data','B_field',('Bx_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bxf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','B_field',('By_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Byf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','B_field',('Bz_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bzf.transpose(),fmt='%15.14e')

			np.savetxt( os.path.join(path_write,'data','B_field',('Bx_tot_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bxf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','B_field',('By_tot_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Byf.transpose()+Byb.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path_write,'data','B_field',('Bz_tot_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bzf.transpose()+Bzb.transpose(),fmt='%15.14e')

			#np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
			#np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
			#np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')
