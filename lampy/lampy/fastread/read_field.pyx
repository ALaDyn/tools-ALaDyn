#!python
#cython: language_level=3
#setuptools: sources = clibs/lib_read_binary.c

cdef extern from "clibs/lib_read_binary.c":
    pass
cdef extern from 'clibs/lib_read_binary.h':
    void read_binary(float* field, float* x, float* y, float* z, char* file_name)

from cython.view cimport array as cvarray
import struct
import numpy as np
cimport numpy as np

def read_ALaDyn_bin(file_path,params):
    
    # - #
    cdef int nx
    cdef int ny
    cdef int nz
    
    path = file_path+'.bin'
    f = open(path, 'rb')
    uni_path = path.encode('UTF-8')
    cdef char* c_path = uni_path

    # Here it reads the header via a standard struct unpack procedure
    f.seek(4)
    Nparam = struct.unpack('=i', f.read(4))[0]
    f.seek(16)
    integerdata_temp = struct.unpack('='+Nparam*'i', f.read(Nparam*4))
    f.seek(104)
    realdata_temp = struct.unpack('='+Nparam*'f', f.read(Nparam*4))
    f.seek(4)
    nproc_y = integerdata_temp[0]
    nproc_z = integerdata_temp[1]
    ndimension = integerdata_temp[14]
    nx = params['nx']
    ny = params['ny']
    nz = params['nz']
    # End of the header read

    totlen = nx*ny*nz
    f = cvarray(shape=(totlen,), itemsize=sizeof(float), format="f")  
    x = cvarray(shape=(nx,), itemsize=sizeof(float), format="f")    
    y = cvarray(shape=(ny,), itemsize=sizeof(float), format="f")   
    z = cvarray(shape=(nz,), itemsize=sizeof(float), format="f") 

    read_file(f, x, y, z, c_path)

    f = np.reshape(f, (nx,ny,nz), order='F')
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    return (f, x, y, z)

def read_file(f, x, y, z, f_path):

    cdef float [::1] f_view = f
    cdef float [::1] x_view = x
    cdef float [::1] y_view = y
    cdef float [::1] z_view = z

    read_binary(&f_view[0], &x_view[0], &y_view[0], &z_view[0], f_path)

    return (f, x, y, z)
