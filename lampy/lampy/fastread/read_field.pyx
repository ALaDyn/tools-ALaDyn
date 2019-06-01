#!python
#cython: language_level=3

import os
import struct
import numpy as np
from cpython cimport array
from ctypes import *
from ..utilities.Utility import module_path

def read_ALaDyn_bin(file_path,params):
    
    # - #
    
    path     = file_path+'.bin'
    f        = open(str(path), 'rb')
    path     = c_char_p(path.encode('utf-8'))
    cpath = module_path
    cpath = os.path.join(cpath, 'fastread', 'clibs')
    read_bin = CDLL(os.path.join(cpath, 'lib_read_binary.so'))

    struct.unpack('i', f.read(4))
    N_param = struct.unpack('i', f.read(4))[0]

    struct.unpack('i', f.read(4))
    struct.unpack('i', f.read(4))
    int_param = list()
    for i in range(0, N_param):
        int_param.append(struct.unpack('i', f.read(4)))
    struct.unpack('i', f.read(4))
    nproc_y = int(int_param[0][0])
    nproc_z = int(int_param[1][0])
    ndimension = int(int_param[14][0])
    struct.unpack('i', f.read(4))
    for i in range(0, N_param):
        struct.unpack('f', f.read(4))
    struct.unpack('i', f.read(4))
    nx = params['nx']
    ny = params['ny']
    nz = params['nz']
    #---***---#
    totlen = nx*ny*nz
    fieldarr = (c_float*totlen)    
    gridx = (c_float*nx)
    gridy = (c_float*ny)
    gridz = (c_float*nz)
    r = fieldarr()    
    x = gridx()
    y = gridy()
    z = gridz()

    read_bin.read_binary(byref(r), byref(x), byref(y), byref(z), path)
    r = np.ndarray(buffer=r, dtype=np.float32, shape=(nx,ny,nz), order='F')
    x = np.ndarray(buffer=x, dtype=np.float32, shape=(nx))
    y = np.ndarray(buffer=y, dtype=np.float32, shape=(ny))
    z = np.ndarray(buffer=z, dtype=np.float32, shape=(nz))

    return (r, x, y, z)
