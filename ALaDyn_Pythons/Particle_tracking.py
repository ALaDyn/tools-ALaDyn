import numpy as np
import math 
import matplotlib.pyplot as plt
import scipy as sc
import os, os.path, sys

#from read_ALaDyn_bin import *
#from ALaDyn_plot_matrixsections import *
#from ALaDyn_plot_utilities_Efield import *
#from ALaDyn_plot_utilities_Bfield import *
#from ALaDyn_plot_utilities_axes import *


def YEE_Ex_deposition(weights, px_index, dx):
		
	if (px_index[i]<dx/2.):
		px_effective = px_index[i]+dx/2.
	else:
		px_effective = px_index[i]-dx/2.
	
	weights[0][0] = (dx - px_effective)*(dy - py_index[i]   )/(dx*dy)	
	weights[1][0] = (px_effective     )*(dy - py_index[i]   )/(dx*dy)	
	weights[0][1] =	(dx - px_effective)*(py_index[i]        )/(dx*dy)
	weights[1][1] =	(px_effective     )*(py_index[i]        )/(dx*dy)	

def YEE_cell_correction_Ex(px_index, px_effective,dx):
		
	px_effective = px_index[i]
	py_effective = py_index[i]
	if(px_index[i]<dx/2.):    
		px_effective = px_index[i] - 1
	
def YEE_Ey_deposition(weights, py_index, dy):
		
	if (py_index[i]<dy/2.):
		py_effective = py_index[i]+dy/2.
	else:
		py_effective = py_index[i]-dy/2.
	
	weights[0][0] = (dx - px_index[i]	 )*(dy - py_effective)/(dx*dy)
	weights[1][0] = (px_index[i]     	 )*(dy - py_effective)/(dx*dy)
	weights[0][1] =	(dx - px_index[i]	 )*(py_effective	 )/(dx*dy)
	weights[1][1] =	(px_index[i]     	 )*(py_effective	 )/(dx*dy)

def YEE_cell_correction_Ey(py_index, py_effective, dy):
		
	px_effective = px_index[i]
	py_effective = py_index[i]
	if(py_index[i]<dy/2.):
		py_effective = py_index[i] - 1
		
def YEE_Bz_deposition(weights, px_index, py_index, dx, dy):

	if (px_index[i]<dx/2.):
		px_effective = px_index[i]+dx/2.
	else:
		px_effective = px_index[i]-dx/2.
	
	if (py_index[i]<dx/2.):
		py_effective = py_index[i]+dy/2.
	else:
		py_effective = py_index[i]-dy/2.
	
	weights[0][0] = (dx - px_effective	  )*(dy - py_effective)/(dx*dy)
	weights[1][0] = (px_effective         )*(dy - py_effective)/(dx*dy)
	weights[0][1] = (dx - px_effective	  )*(py_effective     )/(dx*dy)
	weights[1][1] = (px_effective		  )*(py_effective	  )/(dx*dy)

def YEE_cell_correction_Bz(px_index,py_index, px_effective, py_effective, dx, dy):

	px_effective = px_index[i]
	py_effective = py_index[i]
	if(px_index[i]<dx/2.):    
		px_effective = px_index[i] - 1
	if(py_index[i]<dy/2.):
		py_effective = py_index[i] - 1

def find_left_right(point,mesh_index, axes):
	
	if(axes == 'X'):
		if (cell_index_x[mesh_index] > point):
			return -1
		elif (cell_index_x[mesh_index] == point):
			return 0
		else:
			return 1
	
	if(axes == 'Y'):
		if (cell_index_y[mesh_index] > point):
			return -1
		elif (cell_index_y[mesh_index] == point):
			return 0
		else:
			return 1
	
			
	return -1000000000000		
	
def findMeshIndex(loc_r, axes):
	
	condition = 0
	extl = 0

	
	if (axes == 'X'):
	mesh_index = (NGX-1)/2
	extl = 0
	extr = NGX
	
	if (axes == 'X'):
	mesh_index = (NGY-1)/2
	extl = 0
	extr = NGY

	while (condition == 0 and extl != extr):
	result = find_left_right(loc_r, mesh_index, axes)
	if(result == 1):
		extl = mesh_index
		mesh_index = extl + (extr-extl)/2
	elif(result==-1):
		extr = mesh_index
		mesh_index=extr-(extr-extl)/2
      
    else:            
		return mesh_index;
      
    if (extl==extr-1):
		condition=1;
		return extl;
      
	return extl	

		
def Reassign_particle_position(px_index, cell_index_x, py_index, cell_index_y,dx, dy):

	if (px_index[i] > dx and px_index[i] < 2.*dx and cell_index_x[i] < NGX-1):
		px_index[i] -= dx
		cell_index_x[i] += 1
	elif (px_index[i] < 0. and abs(px_index[i]) < dx and cell_index_x[i]	> 0):
		px_index[i] = px_index[i] + dx
		cell_index_x[i] -= 1
	else:	
		loc_x 		 = px_index[i] + cell_index_x[i]
		cell_index_x[i] = findMeshIndex(loc_x,'X')
		px_index[i] 	 = loc_x - cell_index_x[i] 
		
	if (py_index[i] > dy and py_index[i] < 2.*dy and cell_index_y[i] < NGY-1):
		py_index[i] -= dy
		cell_index_y[i] += 1
	elif (py_index[i] < 0. and abs(py_index[i]) < dy and cell_index_y[i] > 0):
		py_index[i] = py_index[i] + dy
		cell_index_y[i] -= 1
	else:
		loc_y = py_index[i] + cell_index_y[i]
		cell_index_y[i] = findMeshIndex(loc_y,'Y')
		py_index[i] = loc_y - cell_index_y[i]
		
def RK4_Pusher(Exp, Eyp, Bzp, px_index, py_index, ux_p, uy_p, Dt):

	gamma_old = np.sqrt(1+ux_p*ux_p+uy_p*uy_p)

	ux_p_old = ux_p
	uy_p_old = uy_p
	
	px_index_old = px_index
	py_index_old = py_index
	
	k1_ux = Exp + (uy_p*Bzp			  )/gamma_old
	k1_uy = Eyp + (		    - ux_p*Bzp)/gamma_old
	
	k2_ux = Exp + (uy_p*Bzp			  )/gamma_old + Dt/2.*(k1_uy*Bzp		    )
	k2_uy = Eyp + (			- ux_p*Bzp)/gamma_old + Dt/2.*(			 - k1_ux*Bzp)

	k3_ux = Exp + (uy_p*Bzp			  )/gamma_old + Dt/2.*(k2_uy*Bzp		    )
	k3_uy = Eyp + (			- ux_p*Bzp)/gamma_old + Dt/2.*(			 - k2_ux*Bzp)
	
	k4_ux = Exp + (uy_p*Bzp			  )/gamma_old + Dt   *(k3_uy*Bzp		    )
	k4_uy = Eyp + (			- ux_p*Bzp)/gamma_old + Dt   *(			 - k3_ux*Bzp)
	
	
	ux_p 	 += qm * Dt *(k1_ux+2*k2_ux+2*k3_ux+k4_ux)/6.
	uy_p  	 += qm * Dt *(k1_uy+2*k2_uy+2*k3_uy+k4_uy)/6.	
	
	gamma = np.sqrt(1+ux_p*ux_p+uy_p*uy_p)
	
	px_index += Dt * (ux_p/gamma + ux_p_old/gamma_old)/2.
	px_index += Dt * (uy_p/gamma + uy_p_old/gamma_old)/2.
	
	
	Reassign_particle_position(px_index,cell_index_x, py_index,cell_index_y, dx, dy)
	
	
	
###########-------Main-------###########	

#data=np.array(np.genfromtxt('bubble_RK4.dat',dtype=float))


#Ex_b   = np.array(np.genfromtxt('Ex_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt',dtype=float))
#Ey_b   = np.array(np.genfromtxt('Ey_bunch_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt',dtype=float))

#Ex_bkg = np.array(np.genfromtxt('Ex_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt',dtype=float))
#Ey_bkg = np.array(np.genfromtxt('Ey_bck_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt',dtype=float))
 
#Bzp    = np.array(np.genfromtxt('Bz_'+axis_to_cut+'_'+('%2.2i'%i)+'.txt',dtype=float)) 
	
Ex    = np.array(np.genfromtxt(os.path.join('/home/HPC/mira/SCRATCH/2nd_bubble/data/E_field','Ex_tot_y_13.txt')))*1e9
Ey    = np.array(np.genfromtxt(os.path.join('/home/HPC/mira/SCRATCH/2nd_bubble/data/E_field','Ey_tot_y_13.txt')))*1e9
	
Bz    = np.array(np.genfromtxt(os.path.join('/home/HPC/mira/SCRATCH/2nd_bubble/data/B_field','Bz_tot_y_13.txt')))*1e9 
	
#print Exp.shape, Eyp.shape, Bzp.shape	
[NGX,NGY]=Exp.shape
# NGX,NGY
# NGX = 8
# NGY = 8
dx = 1e-6
dy = 1e-6
Np = 1
Dt = 1e-5
qm = 175824175824.1758
px_effective = 0
py_effective = 0
count = 0
iter = 1000




xp = dx*(np.ones(Np)*(NGX-1)-0.5)
yp = np.linspace(2,5,Np)

ux_p = 3.e8*np.ones(Np)
uy_p = np.zeros(Np)


px_index 		= np.zeros(Np)
py_index 		= np.zeros(Np)

px_history = np.zeros(iter)
py_history = np.zeros(iter)

cell_index_x 	= np.zeros(Np)
cell_index_y 	= np.zeros(Np)

weights 		= np.zeros((2,2))



for count in range (0,iter):	
	
	
	
	
	for i in range (0,Np):
		
		px_history[count] = px_index[i]
		py_history[count] = py_index[i]
		
		
		px_index[i] = xp[i] - np.floor(xp[i]) 
		py_index[i] = yp[i] - np.floor(yp[i])
		cell_index_x[i] = np.floor(xp[i])
		cell_index_y[i] = np.floor(yp[i])
	
		YEE_Ex_deposition(weights, px_index, dx)	
		YEE_cell_correction_Ex(px_index, px_effective,dx)
		Ex00 = weights[0][0]*Ex[px_effective  ][py_effective  ]
		Ex10 = weights[1][0]*Ex[px_effective+1][py_effective  ]
		Ex01 = weights[0][1]*Ex[px_effective  ][py_effective+1]
		Ex11 = weights[1][1]*Ex[px_effective+1][py_effective+1]
	
		Exp  = Ex00+Ex01+Ex10+Ex11 
	
		YEE_Ey_deposition(weights, py_index, dy)
		YEE_cell_correction_Ey(py_index, py_effective, dy)
		Ey00 = weights[0][0]*Ey[px_effective  ][py_effective  ]
		Ey10 = weights[1][0]*Ey[px_effective+1][py_effective  ]
		Ey01 = weights[0][1]*Ey[px_effective  ][py_effective+1]
		Ey11 = weights[1][1]*Ey[px_effective+1][py_effective+1]
		
		Eyp  = Ey00+Ey01+Ey10+Ey11
			
		YEE_Bz_deposition(weights, px_index, py_index, dx, dy)
		YEE_cell_correction_Bz(px_index,py_index, px_effective, py_effective, dx, dy)
		Bz00 = weights[0][0]*Bz[px_effective  ][py_effective  ]
		Bz10 = weights[1][0]*Bz[px_effective+1][py_effective  ]
		Bz01 = weights[0][1]*Bz[px_effective  ][py_effective+1]
		Bz11 = weights[1][1]*Bz[px_effective+1][py_effective+1]
	
		Bzp = Bz00+Bz01+Bz10+Bz11
	
		RK4_Pusher(Exp, Eyp, Bzp, px_index, py_index, ux_p, uy_p, Dt)
		
	count +=1	
	
	
	

