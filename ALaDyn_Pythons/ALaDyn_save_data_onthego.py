#!/usr/bin/python
######################################################################
# Name:         ALaDyn_save_data.py
# Author:       
# Date:			2014-05-05
# Purpose:      it saves ALaDyn outputs in ASCII format 
# Source:       python
#####################################################################









### ---------------------------------------------------- ###
### loading shell commands
import os, os.path, glob, sys, shutil
from matplotlib import colors, ticker, cm
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
### --- ###

#--- *** ---#
path = os.getcwd()

# --- #
# FRAME NUMBER
frame = 2


# --- #
magnification_fig = 3.0

#--- *** ---#
file_name = 'Bdenout'+('%2.2i'%frame)+'.bin'
rho_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Edenout'+('%2.2i'%frame)+'.bin'
rho_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

rho_b=np.abs( rho_b )
rho_w=np.abs( rho_w )

rho = rho_b+rho_w

p = rho.shape
x2=p[0]/2; y2=p[1]/2; z2=p[2]/2;

#--- *** ---#
file_name = 'Exbout'+('%2.2i'%frame)+'.bin'
Ex_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Eybout'+('%2.2i'%frame)+'.bin'
Ey_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Ezbout'+('%2.2i'%frame)+'.bin'
Ez_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

file_name = 'Exfout'+('%2.2i'%frame)+'.bin'
Ex_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Eyfout'+('%2.2i'%frame)+'.bin'
Ey_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
file_name = 'Ezfout'+('%2.2i'%frame)+'.bin'
Ez_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')

Ex = Ex_b+Ex_w
Ey = Ey_b+Ey_w
Ez = Ez_b+Ez_w
E  = np.power( np.power(Ex,2)+np.power(Ey,2)+np.power(Ez,2) , 0.5 )

#--- *** ---#
# file_name = 'Bxbout'+('%2.2i'%frame)+'.bin'
# Bx_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# file_name = 'Bybout'+('%2.2i'%frame)+'.bin'
# By_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# file_name = 'Bzbout'+('%2.2i'%frame)+'.bin'
# Bz_b,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 
# file_name = 'Bxfout'+('%2.2i'%frame)+'.bin'
# Bx_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# file_name = 'Byfout'+('%2.2i'%frame)+'.bin'
# By_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# file_name = 'Bzfout'+('%2.2i'%frame)+'.bin'
# Bz_w,  x,y,z = read_ALaDyn_bin(path,file_name,'grid')
# 
# Bx = Bx_w#+Bx_b
# By = By_w#+By_b
# Bz = Bz_w#+Bz_b
# B  = np.power( np.power(Bx,2)+np.power(By,2)+np.power(Bz,2) , 0.5 )


size1, size2 = figure_dimension_inch(x,y,z,magnification_fig)

#--- rho-sections ---#
np.savetxt( 'rho_section_'+('%2.2i'%frame)+'.dat' ,rho[:,:,z2].T,fmt='%15.14e')
np.savetxt( 'rho_b_section_'+('%2.2i'%frame)+'.dat' ,rho_b[:,:,z2].T,fmt='%15.14e')

#--- E-sections ---#
np.savetxt( 'Ex_z0_section_'+('%2.2i'%frame)+'.dat' ,Ex[:,:,z2].T,fmt='%15.14e')
np.savetxt( 'Ey_z0_section_'+('%2.2i'%frame)+'.dat' ,Ey[:,:,z2].T,fmt='%15.14e')
np.savetxt( 'Ez_z0_section_'+('%2.2i'%frame)+'.dat' ,Ez[:,:,z2].T,fmt='%15.14e')
np.savetxt( 'E_z0_section_'+('%2.2i'%frame)+'.dat'  ,E[:,:,z2].T,fmt='%15.14e')

#--- B-sections ---#
# np.savetxt( 'Bx_z0_section_'+('%2.2i'%frame)+'.dat' ,Bx[:,:,z2].T,fmt='%15.14e')
# np.savetxt( 'By_z0_section_'+('%2.2i'%frame)+'.dat' ,By[:,:,z2].T,fmt='%15.14e')
# np.savetxt( 'Bz_z0_section_'+('%2.2i'%frame)+'.dat' ,Bz[:,:,z2].T,fmt='%15.14e')
# np.savetxt( 'B_z0_section_'+('%2.2i'%frame)+'.dat'  ,B[:,:,z2].T,fmt='%15.14e')

#--- axis ---#
np.savetxt( 'x_'+('%2.2i'%frame)+'.dat' ,x,fmt='%15.14e')
np.savetxt( 'y_'+('%2.2i'%frame)+'.dat' ,y,fmt='%15.14e')
np.savetxt( 'z_'+('%2.2i'%frame)+'.dat' ,z,fmt='%15.14e')