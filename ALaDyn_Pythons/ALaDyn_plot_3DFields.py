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

						# read_ALaDyn_bin(dir_path,file_name,grid_no_grid)
						# read_ALaDyn_bin_section(dir_path,file_name,grid_no_grid,axis_to_cut,cell_to_cut)
		
		Exf,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Exfout'+s+'.bin','grid')
		Exb,ax4,ax5,ax6 = read_ALaDyn_bin(path_read,'Exbout'+s+'.bin','grid')
		
		Eyf,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Eyfout'+s+'.bin','grid')
		Eyb,ax4,ax5,ax6 = read_ALaDyn_bin(path_read,'Eybout'+s+'.bin','grid')

		Ezf,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Ezfout'+s+'.bin','grid')
		Ezb,ax4,ax5,ax6 = read_ALaDyn_bin(path_read,'Ezbout'+s+'.bin','grid')
		
		
		np.savetxt( os.path.join(path_write,'data','E_field',('Ex_tot_3D_'+('%2.2i'%i)+'.txt')),Exf.transpose()+Exb.transpose(),fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','E_field',('Ey_tot_3D_'+('%2.2i'%i)+'.txt')),Eyf.transpose()+Eyb.transpose(),fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','E_field',('Ez_tot_3D_'+('%2.2i'%i)+'.txt')),Ezf.transpose()+Ezb.transpose(),fmt='%15.14e')

		
		Bxf,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Bxfout'+s+'.bin','grid')
        
		Byf,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Byfout'+s+'.bin','grid')
		Byb,ax4,ax5,ax6 = read_ALaDyn_bin(path_read,'Bybout'+s+'.bin','grid')
       
		Bzf,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Bzfout'+s+'.bin','grid')
		Bzb,ax4,ax5,ax6 = read_ALaDyn_bin(path_read,'Bzbout'+s+'.bin','grid')
        
        
        np.savetxt( os.path.join(path_write,'data','B_field',('Bx_tot_3D_'+('%2.2i'%i)+'.txt')),Exf.transpose(),fmt='%15.14e')
        np.savetxt( os.path.join(path_write,'data','B_field',('By_tot_3D_'+('%2.2i'%i)+'.txt')),Byf.transpose()+Byb.transpose(),fmt='%15.14e')
        np.savetxt( os.path.join(path_write,'data','B_field',('Bz_tot_3D_'+('%2.2i'%i)+'.txt')),Bzf.transpose()+Bzb.transpose(),fmt='%15.14e')