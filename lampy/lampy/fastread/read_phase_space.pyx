#!python
#cython: language_level=3
#setuptools: sources = clibs/lib_read_phase_space.c

cdef extern from "clibs/lib_read_phase_space.c":
    pass
cdef extern from "clibs/lib_read_phase_space.h":
    void read_phasespace(float* phase_space, int n_dimensions, char* file_path)
    int count_particles(int n_dimensions, char* file_path)
import numpy as np
cimport numpy as np
from cython.view cimport array as cvarray

def total_phase_space_read(file_path, params):
    
    # - #

    path     = file_path+'.bin'
    f        = open(str(path), 'rb')
    uni_path = path.encode('UTF-8')

    cdef char* c_path = uni_path
    cdef int n_dimensions = params['n_dimensions']
    cdef int part_number

    part_number = count_particles(n_dimensions, c_path)
    jump = 2*n_dimensions+1
    tot_dimension = jump*part_number

    ps = cvarray(shape=(tot_dimension,), itemsize=sizeof(float), format="f")    

    read_phase_space(ps, n_dimensions, c_path)

    ps = np.reshape(ps, (jump, part_number), order='F')
    ps = np.asarray(ps)

    return (ps, part_number)

def read_phase_space(ps, n_dimensions, f_path):

    cdef float[::1] ps_view = ps

    read_phasespace(& ps_view[0], n_dimensions, f_path)

    return ps