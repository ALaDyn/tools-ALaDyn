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
# from utility_Efield import *
# from utility_Bfield import *
# from utility_ionization import *
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
frame_begin 	= int(	sys.argv[5])
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


for i in range(frame_begin,frame_end+1):
	print '-------------------'
	s='%2.2i'%i
	path_read  = os.path.join(path,'%4.4i'%i)
	path_write = path
	print path_read
	if output_exists(path_read,'rho',i) == True:
		print 'n_bunch --- frame >>> ',i
		n_bunch,ax1,ax2,ax3 = read_ALaDyn_bin(path_read,'Bdenout'+s+'.bin','grid')
                print 'n_bck --- frame >>> ',i
		n_bck,ax4,ax5,ax6   = read_ALaDyn_bin(path_read,'Edenout'+s+'.bin','grid')
                print 'jx_bunch --- frame >>> ',i
		jx_bunch,ax7,ax8,ax9= read_ALaDyn_bin(path_read,'Jxbout'+s+'.bin','grid')     	
                print 'jx_bck --- frame >>> ',i
                jx_bck,ax7,ax8,ax9  = read_ALaDyn_bin(path_read,'Jxfout'+s+'.bin','grid')

		n                   = n_bunch + n_bck 
		jx		    = jx_bunch + jx_bck
	

	#variable declaration

	[nx,ny,nz]=n_bck.shape
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

	psi  = np.zeros((nx,ny,nz))
	
	del(n_bunch);del(n_bck);
	del(jx_bunch);del(jx_bck);
	
	Delta 	    = -2.*(1./dz**2+1./dy**2)
	y_diag 	    = np.ones( ny*nz)*1./dy**2
	z_diag      = np.ones( ny*nz)*1./dz**2
	
	print 'diag built'
	
	for i in range(1,nz):
		y_diag[ny]=0.
	
	main_diag = np.ones(ny*nz)*Delta
	print 'diag inserted'
	
	data = np.array([main_diag,y_diag,y_diag,z_diag,z_diag])
	print 'data matrix built'
	diagonals = np.array([0,-1,1,-ny,ny])

	
	for cell in range(0,nx):	
		
		n_column	= np.reshape(	   n[cell,:,:]      ,ny*nz)
		jx_column	= np.reshape(	  jx[cell,:,:]      ,ny*nz)
		psi_column	= np.reshape(	 psi[cell,:,:]      ,ny*nz)
		print 'matrix initializated'
		
		
		M = spdiags(data,diagonals,ny*nz,ny*nz).tocsr()
		print 'sparse matrix built'
		
		psi_column = spsolve(M,-(n_column-jx_column))
		print 'matrix solved'
		psi[cell,:,:] = np.reshape(psi_column,(ny,nz))
		
		print 'slice ', cell
		
	np.savetxt( os.path.join(path_write,'data',('psi2D_'+('%2.2s'%s)+'.txt')),np.transpose(psi[:,:,nz/2]),fmt='%15.14e')
	print 'finish'
