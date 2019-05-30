#!python
#cython: language_level=3

import os
import struct
import numpy as np

from cpython cimport array
from ctypes import *
import array

def total_phase_space_read(dir_path,file_name,n_dimensions):
    
    # - #

    path     = os.path.join(os.path.join(dir_path,file_name+'.bin'))
    f        = open(str(path),'rb')
    path     = c_char_p(os.path.join(os.path.join(dir_path,file_name+'.bin')))
    cpath=os.path.dirname(sys.modules[__name__].__file__)
    cpath=os.path.join(cpath,'fastlib','clibs')
    read_bin= CDLL(os.path.join(cpath,'lib_read_phase_space.so'))

    #---***---#
    #totlen=(c_int*1)
    #part_number=totlen(0)
    part_number=read_bin.count_particles(n_dimensions,path)
    jump=2*n_dimensions+1
    #part_number=np.int(part_number[0])
    partarray=(c_float*jump*part_number)    

    ps=partarray()    
    
    read_bin.read_phasespace(byref(ps), n_dimensions, path)
    ps=np.ndarray(buffer=ps,dtype=np.float32,shape=(jump,part_number),order='F')


    return ps
