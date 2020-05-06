#!python
#cython: language_level=3
#setuptools: sources = clibs/lib_read_tracking.c

cdef extern from "clibs/lib_read_tracking.c":
    pass
cdef extern from "clibs/lib_read_tracking.h":
    void read_tracking_phasespace(float* phase_space, int n_dimensions, char* file_path)
    int count_tracked_particles(int n_dimensions, char* file_path)
import numpy as np
cimport numpy as np
from cython.view cimport array as cvarray

def total_tracking_read(file_path, params, species):
    
    # - #

    path     = file_path+'.bin'
    f        = open(str(path), 'rb')
    uni_path = path.encode('UTF-8')

    cdef char* c_path = uni_path
    cdef int n_dimensions = params['n_dimensions']
    cdef int non_phase_space_outputs = 3
    cdef int jump
    cdef int part_number

    if params['a_on_particles'][species - 1]:
        non_phase_space_outputs = non_phase_space_outputs + 1

    jump = 2*n_dimensions + non_phase_space_outputs
    part_number = count_tracked_particles(jump, c_path)
    tot_dimension = jump*part_number

    ps = cvarray(shape=(tot_dimension,), itemsize=sizeof(float), format="f")    

    read_tracking(ps, jump, c_path)

    ps = np.reshape(ps, (jump, part_number), order='F')
    ps = np.asarray(ps)

    return (ps, part_number)

def read_tracking(ps, jump, f_path):

    cdef float[::1] ps_view = ps

    read_tracking_phasespace(& ps_view[0], jump, f_path)

    return ps