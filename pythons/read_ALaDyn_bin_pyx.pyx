#!/usr/bin/python
######################################################################
# Name:         read_ALaDyn_dat.py
# Author:		A. Marocchino
# Date:			2014-02-18
# Purpose:      reads dued binary frm output
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time
import struct
#from scipy import *
import numpy as np
from cpython cimport array
import array
### --- ###



# - #
def read_ALaDyn_bin(dir_path,file_name,grid_no_grid):

	cdef int counter_z = 0
	cdef int counter_y = 0
	cdef int k = 0
	cdef int j = 0
	cdef int i = 0
	cdef int nx= 0
	cdef int ny= 0
	cdef int nz= 0
	cdef int nproc_y=0
	cdef int nproc_z=0
	cdef int npx =0
	cdef int npy =0
	cdef int npz =0

	# - #
	path     = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	#- vector length -#
	struct.unpack('i', f.read(4))
	N_param = struct.unpack('i', f.read(4))[0]
	print(N_param)
	struct.unpack('i', f.read(4))


	struct.unpack('i', f.read(4))
	int_param=[]
	for i in range(0,N_param):
		int_param.append( struct.unpack('i', f.read(4)) )
	struct.unpack('i', f.read(4))
	nx= int_param[3][0]
	ny= int_param[4][0]
	nz= int_param[6][0]
	nproc_y = int_param[0][0]
	nproc_z = int_param[1][0]

	struct.unpack('i', f.read(4))
	for i in range(0,N_param):
		struct.unpack('f', f.read(4))
	struct.unpack('i', f.read(4))


	#---***---#
	# r = np.zeros((nx,ny,nz)) #original non-cython
	cdef array.array r = array.array('i', [nx, ny, nz])
	print('total grid size: n=(',nx,ny,nz,')')
	rr=[]
	print('number of Np processors: Np_y=',nproc_y,'Np_z=', nproc_z)

	offsetz = 0
	for counter_z in range(0,nproc_z):
		offsety = 0;

		for counter_y in range(0,nproc_y):
			struct.unpack('i', f.read(4))
			npx= struct.unpack('i', f.read(4))[0]
			npy= struct.unpack('i', f.read(4))[0]
			npz= struct.unpack('i', f.read(4))[0]
			struct.unpack('i', f.read(4))

			struct.unpack('i', f.read(4))
			for k in range(0,npz):
				for j in range(0,npy):
					for i in range(0,npx):
						r[i,j+offsety,k+offsetz] = struct.unpack('f', f.read(4))[0]
			struct.unpack('i', f.read(4))
			offsety += npy
		offsetz += npz;

	if grid_no_grid == 'nogrid':
		return r

	#--- * --- * --- * --- * --- * ---#
	#- reading grid -#
	struct.unpack('i', f.read(4))
	X=[]; [X.append(struct.unpack('f', f.read(4))[0]) for i in range(0,nx)]
	struct.unpack('i', f.read(4))

	struct.unpack('i', f.read(4))
	Y=[];  [Y.append(struct.unpack('f', f.read(4))[0]) for i in range(0,ny)]
	struct.unpack('i', f.read(4))

	struct.unpack('i', f.read(4))
	Z=[];  [Z.append(struct.unpack('f', f.read(4))[0]) for i in range(0,nz)]
	struct.unpack('i', f.read(4))

# 	x=np.zeros((nx,ny,nz)); y=z=x;
# 	for k in range(0,nz):
# 		for j in range(0,ny):
# 			for i in range(0,nx):
# 				x[i,j,k] = X[i]
# 				y[i,j,k] = Y[j]
# 				z[i,j,k] = Z[k]
	x=X; y=Y; z=Z;

	return (r,x,y,z)





# --- --- --- --- --- --- --- --- --- --- --- #
# --- --- --- --- --- --- --- --- --- --- --- #
# --- --- --- --- --- --- --- --- --- --- --- #
def read_ALaDyn_bin_section(dir_path,file_name,grid_no_grid,axis_to_cut,cell_to_cut):

	# - #
	path     = os.path.join(os.path.join(dir_path,file_name))
	f        = open(path,'rb')

	#- vector length -#
	struct.unpack('i', f.read(4))
	N_param = struct.unpack('i', f.read(4))[0]
	print(N_param)
	struct.unpack('i', f.read(4))


	struct.unpack('i', f.read(4))
	int_param=[]
	for i in range(0,N_param):
		int_param.append( struct.unpack('i', f.read(4)) )
	struct.unpack('i', f.read(4))
	nx= int_param[3][0]
	ny= int_param[4][0]
	nz= int_param[6][0]
	nproc_y = int_param[0][0]
	nproc_z = int_param[1][0]

	struct.unpack('i', f.read(4))
	for i in range(0,N_param):
		struct.unpack('f', f.read(4))
	struct.unpack('i', f.read(4))


	#---***---#
	if axis_to_cut == 'x':
		r = np.zeros((ny,nz))
	elif axis_to_cut == 'y':
		r = np.zeros((nx,nz))
	elif axis_to_cut == 'z':
		r = np.zeros((nx,ny))
	rr=[]
# 	print 'number of Np processors: Np_y=',nproc_y,'Np_z=', nproc_z

	offsetz = 0
	for counter_z in range(0,nproc_z):
		offsety = 0;

		for counter_y in range(0,nproc_y):
			struct.unpack('i', f.read(4))
			npx= struct.unpack('i', f.read(4))[0]
			npy= struct.unpack('i', f.read(4))[0]
			npz= struct.unpack('i', f.read(4))[0]
			struct.unpack('i', f.read(4))

			struct.unpack('i', f.read(4))
			for k in range(0,npz):
				for j in range(0,npy):
					for i in range(0,npx):
						temp = struct.unpack('f', f.read(4))[0]
						if axis_to_cut == 'x' and i == nx/2+cell_to_cut:
							r[j+offsety,k+offsetz] = temp
						elif axis_to_cut == 'y' and (j+offsety == ny/2+cell_to_cut or ny==1):
							r[i,k+offsetz] = temp
						elif axis_to_cut == 'z' and (k+offsetz == nz/2+cell_to_cut or nz==1):
							r[i,j+offsety] = temp
			struct.unpack('i', f.read(4))
			offsety += npy
		offsetz += npz;

	if grid_no_grid == 'nogrid':
		return r

	#--- * --- * --- * --- * --- * ---#
	#- reading grid -#
	struct.unpack('i', f.read(4))
	X=[]; [X.append(struct.unpack('f', f.read(4))[0]) for i in range(0,nx)]
	struct.unpack('i', f.read(4))

	struct.unpack('i', f.read(4))
	Y=[];  [Y.append(struct.unpack('f', f.read(4))[0]) for i in range(0,ny)]
	struct.unpack('i', f.read(4))

	struct.unpack('i', f.read(4))
	Z=[];  [Z.append(struct.unpack('f', f.read(4))[0]) for i in range(0,nz)]
	struct.unpack('i', f.read(4))

	x=X; y=Y; z=Z;

	return (r,x,y,z)

# 	if axis_to_cut == 'x':
# 		return (r,y,z)
# 	elif axis_to_cut == 'y':
# 		return (r,x,z)
# 	elif axis_to_cut == 'z':
# 		return (r,x,y)
