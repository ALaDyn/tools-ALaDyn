#Author:           Alberto Albz Marocchino
# Date:                 10-08-2016
# Purpose:   Active Plasma Lens plot with alpha
# Source:     python
#####################################################################


### loading shell commands
import os, os.path, glob, sys, shutil
import time, datetime
import scipy
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
### --- ###
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes/Code_ALaDyn/tools-ALaDyn/pythons'))
#from read_architect_bin import *
from transp import *
from read_particle_phasespace import *
from utilities_1 import *

### --- ###

# --- #
colormap = pyl.cm.gist_ncar
color_matrix = [colormap(i) for i in np.linspace(0.0, 0.9, 10)]
#ne_cmap = plt.get_cmap('Blues_r')
ne_cmap = plt.get_cmap('autumn')
#ne_cmap = plt.get_cmap('cool')
#ne_cmap = plt.get_cmap('Greys')
#ne_cmap=truncate_colormap(plt.get_cmap('Greys'), minval=0.10, maxval=0.75, n=100)
#bunch_cmap = plt.get_cmap('YlOrRd_r')
#bunch_cmap = plt.get_cmap('PuBu_r')
bunch_cmap = plt.get_cmap('hot')
# --- #

#--- read architect files ---#
#path = '/Volumes/mira-ui-hpc/SCRATCH/test_Ar_ionization_layer_weights/data'
path = '/storage/gpfs_maestro/hpc/user/mira/sims/sims_albe/data'
rho_b = np.loadtxt(os.path.join(path,'rho','rho_bunch_y_09.txt'))
rho_bck = np.loadtxt(os.path.join(path,'rho','rho_bck_y_09.txt'))
x = np.loadtxt(os.path.join(path,'Moving_window_axes','x_09.txt'))
y = np.loadtxt(os.path.join(path,'Moving_window_axes','y_09.txt'))
#Etot = np.loadtxt(os.path.join(path,'E_field','E_tot_z_01.txt'))
#---------------- Plots ------------------#

rho_min_bck=0.
rho_max_bck=30.
rho_min_b=1.
rho_max_b=10.


fig = pyl.figure(1) #, figsize=(3.25,3.0))
fig.set_size_inches(3.25, 3.0, forward=True)
#--- sub1
ax  = pyl.subplot(111)
extent = (np.min(x),np.max(x),np.min(y),np.max(y))
im = ax.imshow(rho_bck,interpolation = 'bicubic', cmap = ne_cmap, aspect = 'auto' , extent=extent, clim=[rho_min_bck, rho_max_bck], alpha=0.999)
im2 = Imshow(rho_b,interpolation = 'bicubic', cmap = bunch_cmap, aspect = 'auto' , extent=extent, clim=[rho_min_b, rho_max_b], alpha=1.0, tvmin=0., tvmax=rho_min_b)
#im2 = Imshow(Etot,interpolation = 'bicubic', cmap = bunch_cmap, aspect = 'auto' , extent=extent, clim=[200, 900], alpha=1.0, tvmin=0., tvmax=200)
ax.set_ylim([-25,25])
ax.set_xlim([450,630])
# ax.set_yticks([-20,0,20])
# ax.set_xticks([-100,0,100,200])
# pyl.xticks(fontsize=8)
# pyl.yticks(fontsize=8)

plt.xlabel(r'z $[\mu m]$',fontsize=14);plt.ylabel(r'R $[\mu m]$',fontsize=14);

#--- common label ---#
# fig.text(0.5, 0.04, r'$\xi$=Z-ct ($\mu$m)', ha='center', va='center', fontsize = 9.0)
# fig.text(0.06, 0.5, 'X ($\mu$m)', ha='center', va='center', rotation='vertical', fontsize = 9.0)



cbaxes_bck = fig.add_axes([0.85, 0.60, 0.03, 0.30])
cbar_bck = plt.colorbar(im,ax = ax,cmap=ne_cmap,cax = cbaxes_bck,format = '%1i' ,ticks = [rho_min_bck,rho_max_bck])#,aspect = 5,fraction = 0.07)
cbar_bck.ax.tick_params(labelsize=9)
cbar_bck.ax.set_title(r'$n_{bck}/10^{16}$',rotation=0,color='k', fontsize=14)

cbaxes_b = fig.add_axes([0.85, 0.13, 0.03, 0.30])
cbar_b = plt.colorbar(im2,ax = ax,cmap=bunch_cmap,cax = cbaxes_b,format = '%1i' ,ticks = [rho_min_b,rho_max_b]) #aspect = 5,fraction = 0.07)
cbar_b.ax.tick_params(labelsize=9)
cbar_b.ax.set_title(r'$n_{b}/10^{16}$',rotation=0,color='k', fontsize=14)


pyl.xticks(fontsize=8)
pyl.yticks(fontsize=8)

# plt.rc('axes', linewidth=0.75)
fig.subplots_adjust(bottom=0.13,left=0.13,right=0.80)


# pyl.savefig('/Users/albz/I_images/im_PWFA/im_PWFA_256.pdf', format='pdf')
pyl.savefig('/storage/gpfs_maestro/hpc/user/mira/sims/sims_albe/shock_injection_density.eps', format='eps')
pyl.show()

