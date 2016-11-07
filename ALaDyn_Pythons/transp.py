#########################################
# Author : Remi Lehe
# transp.py : a wrapper for imshow,     #
# that allows to plot matrices with     #
# pixel-dependent transparency          #
#########################################

import matplotlib.pyplot as plt
import numpy as np

#def Imshow( data, gam=1., cmap=plt.cm.Blues, tvmin=None, tvmax=None,tmax=None,**kwargs ) :
def Imshow( data, gam=1., tvmin=None, tvmax=None,tmax=None,**kwargs ) :
   # "Shows the data with pixel-dependent transparency"
   # "Key word arguments : vmin and vmax set the range of the colormap, "
   # "while tvmin and tvmax set the range of the transparency"

    # Determine the min and max value of the plot
    if 'vmax' in kwargs :
        vmax = kwargs['vmax']
    else :
        vmax = data.max()
    if 'vmin' in kwargs :
        vmin = kwargs['vmin']
    else :
        vmin = data.min()
    if tvmax is None :
        tvmax = vmax
    if tvmin is None :
        tvmin = vmin    
    if tmax is None :
        tmax = 1.

    cmap=kwargs['cmap']
    #print 'Hello'
    #print vmin,vmax,data.min(),data.max()
    # Rescale the data to get the transparency and color
    color = (data-vmin)/(vmax-vmin)
    color[color >1.] = 1.
    color[color <0.] = 0.
    transparency = tmax*(data-tvmin)/(tvmax-tvmin)
    transparency[transparency > 1.] = 1
    transparency[transparency < 0.] = 0.
    # Application of a gamma
    transparency = tmax*transparency**gam

    # Create an rgba stack of the data, using the colormap 
    rgba_data = cmap( color )
    # Modify the transparency
    rgba_data[:,:,3] = transparency

    im=plt.imshow( rgba_data, **kwargs )
    


    return im
