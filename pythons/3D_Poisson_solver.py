import numpy as np
import math 
#import matplotlib.pyplot as plt
import scipy as sc
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve
import os, os.path, sys


home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/ALaDyn_Code/tools-ALaDyn/pythons'))
###>>>
### --- ###
from read_ALaDyn_bin import *
from utilities_1 import *
from utility_density import *
from utility_Efield import *
from utility_Bfield import *
from utility_ionization import *
from utility_energy_density import *
from utility_axes import *


#PWFA parameters for nondimensionalization of equations
m   	= 9.1e-31   
e   	= 1.6e-19   
c   	= 3e8       
eps 	= 8.85e-12  
IA  	= 17e3      
pi  	= 3.14159   
Ib      = 10e3      

if(len(sys.argv)<8):
	print 'Usage: python ~/Codes/tools-ALaDyn/ALaDyn_Pythons/PWFA_values.py 0 9 30'
	print 'density in cm^-3'
	print 'sigmar in mum'
	print 'sigmaz in mum'
	print 'gamma'
	print 'frame begin'
	print 'frame end'
	print 'want to show ?'
	exit('---not enough arguments---')
#---
n0      	= float(sys.argv[1])*1e6
sr      	= float(sys.argv[2])*1e-6
sz      	= float(sys.argv[3])*1e-6
gamma   	= float(sys.argv[4])
frame_begin = int(	sys.argv[5])
frame_end	= int(	sys.argv[6])
check   	= str(	sys.argv[7])
   
omegap  = np.sqrt(e**2*n0/m/eps)
kp      = omegap/c*1e-6

lambdap = 2*pi/kp
E0      = kp*m*c**2/e
Qb      = Ib*np.sqrt(2*pi)*sz/c
nb      = Qb/np.sqrt((2*pi)**3)/sz/sr**2/e
aplha   = nb/n0
Qtilde  = Qb*kp**3/n0/e

if check=='yes':	
	print 'Plasma frequency in fs^-1: ', omegap*1e-15     
	print 'Plasma wavenumber in mum^-1: ', kp     
	print 'Plasma wavelength in mum: ', lambdap
	print '1D cold wave-breaking limit in GV/m', E0    
	print 'Bunch charge in pC', Qb*1e12     
	print 'Bunch density in cm^-3', nb*1e-6     
	print 'Alpha: ',aplha  
	print 'Qtilde: ', Qtilde*1.e18 


path = os.getcwd()

#-folder output structure
generate_folder_output_structure(path,'True')


for i in range(frame_begin, min(frame_end,last_output(os.getcwd())) + 1 ):
	print '-------------------'
	s='%2.2i'%i
	path_read  = os.path.join(path,'%4.4i'%i)
	path_write = path

	if output_exists(path_read,'rho',i) == True:
		print 'n --- frame >>> ',i
		n_bunch,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Bdenout'+s+'.bin','grid')
		n_bck,ax4,ax5,ax6   = read_ALaDyn_bin(path_read,'Edenout'+s+'.bin','grid')
	
		n                   = n_bunch + n_bck 
	
		np.savetxt( os.path.join(path_write,'data','axes',('x_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','axes',('y_'+('%2.2i'%i)+'.txt')),ax2,fmt='%15.14e')
		np.savetxt( os.path.join(path_write,'data','axes',('z_'+('%2.2i'%i)+'.txt')),ax3,fmt='%15.14e')
	
	if output_exists(path_read,'Energy_density',i) == True:
	        print 'energy density --- frame >>> ',i
	        n_gamma_bck,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Elenout'+s+'.bin','grid')

	#variable declaration

	[nx,ny,nz]=n.shape
	
	xmin = min(ax1)
	xmax = max(ax1)
	ymin = min(ax2)
	ymax = max(ax2)
	zmin = min(ax3)
	zmax = max(ax3)
	
	dx = kp*(xmax - xmin)/(ny-1)
	dy = kp*(ymax - ymin)/(nx-1)
	dz = kp*(zmax - zmin)/(nz-1)
	
	#Initializing all vectors and matrices
	print'box dimensions in mum:' 
	print'Lx: ',(xmax-xmin)
        print'Ly: ',(ymax-ymin)
        print'Lz: ',(zmax-zmin)

	n_gamma_bunch  = n_bunch*gamma
	n_gamma        = n_gamma_bunch+n_gamma_bck
	# rho_gamma = rho_gamma_bck
	
	n_column       = np.reshape(n      ,nx*ny*nz)
	n_gamma_column = np.reshape(n_gamma,nx*ny*nz)
	
	jx    		   = np.zeros((nx,ny,nz))
	# phi  		   = np.zeros((nx,ny,nz))
	psi  		   = np.zeros((nx,ny,nz))
	gamma_matrix   = np.zeros((nx,ny,nz))
	
	for i in range(0,nx):
		for j in range(0,ny):
			for k in range(0,nz):
				if (n[i][j][k]<1.e-3) or (n_gamma[i][j][k]<1.e-3):
					gamma_matrix[i][j][k] = 1.
				else:
					gamma_matrix[i][j][k] = n_gamma[i][j][k]/n[i][j][k]
				jx[i][j][k] = n[i][j][k]*np.sqrt(abs(1.- 1./gamma_matrix[i][j][k]**2))
	
	del(gamma_matrix);
	del(n_gamma);del(n_gamma_bunch);del(n_gamma_bck)
	del(n_bunch);del(n_bck);

	jx_column  	   = np.reshape(jx ,nx*ny*nz)
#	phi_column 	   = np.reshape(phi,nx*ny*nz)
	psi_column    	   = np.reshape(psi,nx*ny*nz)
	print 'matrix initializated'
	
	Delta = -2.*(1./dz**2+1./dy**2)
	z_down_diag = np.ones(nx*nz*ny)*1./dz**2
	z_up_diag   = np.ones(nx*nz*ny)*1./dz**2   
	y_diag      = np.ones(nx*nz*ny)*1./dy**2
	
	print 'diag built'
	
	for i in range(1,nx*ny):
		z_down_diag[i*nz-1]=0.
		z_up_diag[i*nz]=0. 
	
	main_diag = np.ones(nx*nz*ny)*Delta
	print 'diag inserted'
	
	data = np.array([main_diag,z_down_diag,z_up_diag,y_diag,y_diag])
	print 'data matrix built'
	diagonals = np.array([0,-1,1,-ny,ny])
	
	M = spdiags(data,diagonals,nx*nz*ny,nx*nz*ny).tocsr()
	print 'sparse matrix built'
	del(main_diag);del(z_down_diag);del(z_up_diag);del(y_diag);

	# phi_column = spsolve(M,4.*pi*n_column)
	psi_column = spsolve(M,4*pi*(n_column-jx_column))
	print 'matrix solved'
	# phi = np.reshape(phi_column,(nx,ny,nz))
	psi = np.reshape(psi_column,(nx,ny,nz))
	
	psi2D = psi[:][:][nz/2]
	
	# np.savetxt('phi.txt',phi)
	np.savetxt( os.path.join(path_write,'data',('psi2D_'+('%2.2i'%i)+'.txt')),ax1,fmt='%15.14e')

	# np.savetxt('jx.txt',jx)
	# np.savetxt('gamma_matrix.txt',gamma_matrix.transpose())
	print 'finish'
	
