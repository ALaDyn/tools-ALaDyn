#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:	A. Marocchino
# Date:		2014-11-04
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
	frame_begin	  = 0
	frame_end         = last_output(os.getcwd())
	axis_to_cut	  = 'y'
	cell_to_cut	  = 0
else:
	frame_begin 	= int(		sys.argv[1])
	frame_end	= int(		sys.argv[2])
	axis_to_cut	= str(		sys.argv[3])
	cell_to_cut	= int(		sys.argv[4])
### --- ###



#--- *** ---#
if __name__ == '__main__':

	#-path
	path = os.getcwd()

	#-folder output structure
	generate_folder_output_structure(path,'True')


	for i in range(frame_begin, frame_end+1 ):
		print '*** frm::',i
		s='%2.2i'%i
		path_read  = os.path.join(path,'%4.4i'%i)
		path_write = path

		K,ax1,ax2,ax3 = read_ALaDyn_bin_section(path_read,'Edenout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
		print K.shape
		K = K[:,720:-720]
		np.savetxt( os.path.join(path_write,'data','rho',('rho_section_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),np.abs(K.transpose()),fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','Moving_window_axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')
		#clear variable:
		del K

		Exf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Exfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
		Eyf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Eyfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)
		Exf = Exf[:,720:-720]
		Exf = Exf[:,720:-720]
		#Ezf,ax4,ax5,ax6 = read_ALaDyn_bin_section(path_read,'Ezfout'+s+'.bin','grid',axis_to_cut,cell_to_cut)

		np.savetxt( os.path.join(path_write,'data','E_field',('Ex_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Exf.transpose(),fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','E_field',('Ey_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Eyf.transpose(),fmt='%15.14e')
		#np.savetxt( os.path.join(path_write,'data','E_field',('Ez_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt')),Ezf.transpose(),fmt='%15.14e')
