#!/usr/bin/python
######################################################################
# Name:         PWFA_bunch_values
# Author:       A. Marocchino
# Date:			24-11-2015
# Purpose:      matching condition for PWFA
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, glob, sys, shutil, time, datetime, scipy, pylab
from scipy import *
import numpy as np
import matplotlib.pyplot as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.constants import codata
#from mpl_toolkits.mplot3d import Axes3D
# - #
home_path = os.path.expanduser('~')
sys.path.append(os.path.join(home_path,'Codes/PlasmaParameters'))
from plasma_basic_parameters import *
# --- #

### --- constants --- ###
physical_const = codata.physical_constants #SI
eps0   = physical_const['electric constant'][0]
mu0    = physical_const['magnetic constant'][0]
c      = physical_const['speed of light in vacuum'][0]
me     = physical_const['electron mass'][0]
qe     = physical_const['elementary charge'][0]


### --- ###
def density_ramp(z,Lramp):
	if z<=Lramp:
		rho = (Lramp-z)/Lramp
	else:
		rho = 0.
	return rho

def fext(z,Lramp,np):
	f_ext = np * qe**2 / (me*eps0*c**2) * density_ramp(z,Lramp)
	return f_ext

def Kext(z,Lramp,np,gamma):
	K_ext = fext(z,Lramp,np) / (2.*gamma)
	return K_ext

def f_rk(z,u1,u2,Lramp,np,gamma,eps):
	return -Kext(z,Lramp,np,gamma)*u1+eps**2/gamma**2/u1**3

def g_rk(z,u1,u2,Lramp,np,gamma,eps):
	return u2


### --- Main --- ###
def matching_with_ramp(sigma_z,Lramp,Lvacuum,np,gamma,emittance):
	h=.5e-6 #integration step
	z=0.
	u1=sigma_z
	u2=0. #D(sigma_z)
	while(z<Lramp+Lvacuum):
		k1=f_rk(z     ,u1        ,u2        ,Lramp,np,gamma,eps)
		g1=g_rk(z     ,u1        ,u2        ,Lramp,np,gamma,eps)

		k2=f_rk(z+h/2.,u1        ,u2+h/2.*k1,Lramp,np,gamma,eps)
		g2=g_rk(z+h/2.,u1+h/2.*g1,u2        ,Lramp,np,gamma,eps)

		k3=f_rk(z+h/2.,u1        ,u2+h/2.*k2,Lramp,np,gamma,eps)
		g3=g_rk(z+h/2.,u1+h/2.*g2,u2        ,Lramp,np,gamma,eps)

		k4=f_rk(z+h   ,u1        ,u2+h*k3   ,Lramp,np,gamma,eps)
		g4=g_rk(z+h   ,u1+h*g3   ,u2        ,Lramp,np,gamma,eps)

		z=z+h
		u2=u2+h*(k1+2.*k2+2.*k3+k4)/6.
		u1=u1+h*(g1+2.*g2+2.*g3+g4)/6.

	sigma=u1
	alphaTwiss=u1*u2/eps*gamma
	betaTwiss=u1**2/eps/1e-6

	return sigma,alphaTwiss,betaTwiss
