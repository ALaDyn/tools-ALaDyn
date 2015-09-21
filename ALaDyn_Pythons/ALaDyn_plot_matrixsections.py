#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:       
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
		
		if output_exists(path,'rho',i) == True:
			print 'rho --- frame >>> ',i
			M,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Bdenout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			K,ax4,ax5,ax6 = read_ALaDyn_bin_section(path,'Edenout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			np.savetxt( os.path.join(path,'data','rho',('rho_bunch_section_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),M.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','rho',('rho_section_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),np.abs(M.transpose())+np.abs(K.transpose()),fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','Moving_window_axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','Moving_window_axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','Moving_window_axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')

		if output_exists(path,'ionization',i) == True:
			print 'ionization rate --- frame >>> ',i
			I,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'H1dnout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			np.savetxt( os.path.join(path,'data','ionization',('ionization_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),I.transpose(),fmt='%15.14e')

		if output_exists(path,'E',i) == True:
			print 'E --- frame >>> ',i
			Exb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Exbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Exf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path,'Exfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			Eyb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Eybout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Eyf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path,'Eyfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			Ezb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Ezbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Ezf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path,'Ezfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			np.savetxt( os.path.join(path,'data','E_field',('Ex_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Exb.transpose()+Exf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','E_field',('Ey_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Eyb.transpose()+Eyf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','E_field',('Ez_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Ezb.transpose()+Ezf.transpose(),fmt='%15.14e')

			np.savetxt( os.path.join(path,'data','Moving_window_axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','Moving_window_axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','Moving_window_axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')


		if output_exists(path,'B',i) == True:
			print 'B --- frame >>> ',i		
			#Bxb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Bxbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Bxf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path,'Bxfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			#Byb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Bybout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Byf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path,'Byfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			#Bzb,ax1,ax2,ax3 = read_ALaDyn_bin_section(path,'Bzbout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
			Bzf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path,'Bzfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

			np.savetxt( os.path.join(path,'data','B_field',('Bx_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bxf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','B_field',('By_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Byf.transpose(),fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','B_field',('Bz_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Bzf.transpose(),fmt='%15.14e')

			np.savetxt( os.path.join(path,'data','Moving_window_axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','Moving_window_axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
			np.savetxt( os.path.join(path,'data','Moving_window_axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')
											###---ionization_output_test---###
											

											###-----------------------------###



# 		
# 		
# 			
# 		if output_exists(path,'Moving_window_axes',i) == True:
# 			print 'Moving window axes --- frame >>> ',i
# 			if savedata == 'True':
#  				print 'Moving Window Coordinates data --- frame >>> ',i		
#  				save_moving_window_coordinates(path,i)

			
























