#!/usr/bin/python
######################################################################
# Name:         ALaDyn_plot_section.py
# Author:
# Date:			2014-02-18
# Purpose:      it nests into 'ALaDyn_read_binary' to plot sections
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, re
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
### --- ###


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

#- get last output number
def last_output(path):
	n_last_output = -1
	for root_dir, sub_dirs, files in os.walk(path):
		for file in files:
			if os.path.splitext(file)[1] == '.bin':
				m = re.match(r'\D+(\d+)',os.path.basename(file))
				n_last_output = max( n_last_output, int(m.group(1)) )
	return n_last_output

#- folder structure for outputs -#
def	generate_folder_output_structure(path,savedata):
	directory = os.path.join(path,'plots')
	if not os.path.exists( directory ):
		os.makedirs(directory)

	directory_rho = os.path.join(directory,'rho')
	directory_ionization   = os.path.join(directory,'ionization')
	directory_Energy_Density = os.path.join(directory, 'Energy_density')
	directory_E   = os.path.join(directory,'E_field')
	directory_B   = os.path.join(directory,'B_field')

	if not os.path.exists( directory_rho ):
		os.makedirs(directory_rho)
	if not os.path.exists( directory_E ):
		os.makedirs(directory_E)
	if not os.path.exists( directory_B ):
		os.makedirs(directory_B)

					###---###
	if not os.path.exists( directory_ionization ):
		os.makedirs(directory_ionization)
					###---###
	if not os.path.exists( directory_Energy_Density ):
		os.makedirs(directory_Energy_Density)

	if (savedata == 'True'):
		directory = os.path.join(path,'data')
		if not os.path.exists( directory ):
			os.makedirs(directory)

		directory_rho 		 = os.path.join(directory,'rho')
		directory_ionization     = os.path.join(directory,'ionization')
		directory_Energy_Density = os.path.join(directory, 'Energy_density')
		directory_E   		 = os.path.join(directory,'E_field')
		directory_B   		 = os.path.join(directory,'B_field')
		directory_moving_window  = os.path.join(directory,'axes')


		if not os.path.exists( directory_rho ):
			os.makedirs(directory_rho)
		if not os.path.exists( directory_ionization ):
			os.makedirs( directory_ionization )
		if not os.path.exists( directory_Energy_Density ):
			os.makedirs( directory_Energy_Density )
		if not os.path.exists( directory_E ):
			os.makedirs(directory_E)
		if not os.path.exists( directory_B ):
			os.makedirs(directory_B)
		if not os.path.exists( directory_moving_window ):
			os.makedirs( directory_moving_window )



#- Image inches dimesions -#
def figure_dimension_inch(x,y,z,scale_factor):
	dx = max(x) - min(x)
	dy = max(y) - min(y)
	dz = max(z) - min(z)

	standard_size_x = 6.5
	standard_size_z = 3.0

	size_x = scale_factor*standard_size_x
	size_z = size_x * dz/dx + 1.5

	return size_x, size_z


#- check whether all the files are in place to be plotted -#
def output_exists(path,quantity,frame):
	if quantity == 'rho':
		if os.path.isfile(os.path.join(path,'Bdenout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Edenout'+('%2.2i'%frame)+'.bin')) == True:
			return True
		else:
			return False

			###---###

	if quantity == 'ionization':
		if os.path.isfile(os.path.join(path,'H1dnout'+('%2.2i'%frame)+'.bin')) == True:
			return True
		else:
			return False

			###---###

	if quantity == 'Energy_density':
		if os.path.isfile(os.path.join(path,'Elenout'+('%2.2i'%frame)+'.bin')) == True:
			return True
		else:
			return False

			###---###
	if quantity == 'phasespace':
		if os.path.isfile(os.path.join(path,'Elpout'+('%2.2i'%frame)+'.bin')) == True:
			return True
		else:
			return False


	if quantity == 'E':
		if os.path.isfile(os.path.join(path,'Exfout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Exbout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Eyfout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Eybout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Ezfout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Ezbout'+('%2.2i'%frame)+'.bin')) == True:
			return True
		else:
			return False

	if quantity == 'B':
		if os.path.isfile(os.path.join(path,'Bxfout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Byfout'+('%2.2i'%frame)+'.bin')) == True \
		and \
		os.path.isfile(os.path.join(path,'Bzfout'+('%2.2i'%frame)+'.bin')) == True:
#		and \
#		os.path.isfile(os.path.join(path,'Bxbout'+('%2.2i'%frame)+'.bin')) == True \
#		and \
#		os.path.isfile(os.path.join(path,'Bybout'+('%2.2i'%frame)+'.bin')) == True \
#		and \
#		os.path.isfile(os.path.join(path,'Bzbout'+('%2.2i'%frame)+'.bin')) == True:
			return True
		else:
			return False


	if quantity == 'Moving_window_axes':
# 		print os.path.join(path,'data','Moving_window_axes',('moving_window_axes_'+('%2.2i'%frame)+'.dat')
# 		window_data_file = os.path.join(path,'data','Moving_window_axes',('moving_window_axes_'+('%2.2i'%frame)+'.dat')
# 		print window_data_file
		#
		if os.path.isfile(os.path.join(path,'Bdenout'+('%2.2i'%frame)+'.bin')) == True:
			return True
		else:
			return False



#- folder structure for VTS outputs -#
def	generate_folder_vts(path):
	directory = os.path.join(path,'VTS_files')
	if not os.path.exists( directory ):
		os.makedirs(directory)
