import os
import sys
from scipy.constants import pi
from numpy import sqrt

lampy_version = '0.1.dev0'
module_path = os.path.dirname(sys.modules[__name__].__file__)
module_path = os.path.join(module_path, os.pardir)

_translated_filenames = dict()
_total_filenamelist = list()
_translated_filenames['Ex'] = 'Exfout'
_translated_filenames['Ey'] = 'Eyfout'
_translated_filenames['Ez'] = 'Ezfout'
_translated_filenames['Bx'] = 'Bxfout'
_translated_filenames['By'] = 'Byfout'
_translated_filenames['Bz'] = 'Bzfout'
_translated_filenames['Jx'] = 'Jxfout'
_translated_filenames['Jy'] = 'Jyfout'
_translated_filenames['Jz'] = 'Jzfout'
_translated_filenames['A'] = 'Aenvout'
_translated_filenames['a'] = 'aenvout'
_translated_filenames['rho_fluid'] = 'Fdenout'
_translated_filenames['rho'] = 'Edenout'
for n in range(1, 6):
    _translated_filenames['rho_ion'+str(n)] = 'H'+str(n)+'dnout'
    _translated_filenames['energy_ion'+str(n)] = 'H'+str(n)+'enout'

_translated_filenames['energy'] = 'Elenout'
_translated_filenames['phase_space'] = 'Elpout'
_translated_filenames['phase_space_ionization'] = 'Eionzout'
_translated_filenames['phase_space_high_energy'] = 'E_hg_out'

_total_filenamelist += ['Exfout']
_total_filenamelist += ['Eyfout']
_total_filenamelist += ['Ezfout']
_total_filenamelist += ['Bxfout']
_total_filenamelist += ['Byfout']
_total_filenamelist += ['Bzfout']
_total_filenamelist += ['Jxfout']
_total_filenamelist += ['Jyfout']
_total_filenamelist += ['Jzfout']
_total_filenamelist += ['Renvout']
_total_filenamelist += ['Ienvout']
_total_filenamelist += ['Aenvout']
_total_filenamelist += ['aenvout']
_total_filenamelist += ['Fdenout']
_total_filenamelist += ['Edenout']
_total_filenamelist += ['Elenout']
_total_filenamelist += ['Elpout']
_total_filenamelist += ['Eionzout']
_total_filenamelist += ['E_hg_out']
_total_filenamelist += ['H'+str(n)+'dnout' for n in range(1, 6)]
_total_filenamelist += ['H'+str(n)+'enout' for n in range(1, 6)]


def _read_simulation(path):
    """
    Utility that reads the 'input.nml' file
    in the given folder to assign the simulation
    parameters
    """
    param_dic = dict()
    elements = '&!/'
    parameter_file = 'input.nml'
    with open(os.path.join(path, parameter_file), 'r') as f:
        lines = f.readlines()
        for line in lines:
            if(len(line) > 1):
                line = line.replace('=', ' = ')
                line = line.replace(',', ' , ')
                splitline = line.split()
                check = True
                for element in elements:
                    check = check and element not in splitline[0]
                if check:
                    param_dic[splitline[0]] = splitline[2]
    intparamlist = ['nx', 'ny', 'nz', 'ny_targ', 'nprocx', 'nprocy', 'nprocz']
    for key, item in param_dic.items():
        if 'true' not in item and 'false' not in item:
            param_dic[key] = float(item)
            if key in intparamlist:
                param_dic[key] = int(item)
        elif 'true' in item:
            param_dic[key] = True
        elif 'false' in item:
            param_dic[key] = False

    return param_dic


def _compute_physical_parameters(dictionary):
    """
    Utility that computes the relevant physical parameters
    starting from the input file.
    """
    speed_of_light = 0.3

    if dictionary['ny'] < 2 and dictionary['nz'] < 2:
        dictionary['n_dimensions'] = 1
    elif dictionary['ny'] > 1 and dictionary['nz'] < 2:
        dictionary['n_dimensions'] = 2
    else:
        dictionary['n_dimensions'] = 3

    if 'lam0' in dictionary.keys():
        dictionary['omega_0'] = 2*pi/dictionary['lam0']
    if 'lam1' in dictionary.keys():
        dictionary['omega_1'] = 2*pi/dictionary['lam1']
    if 'n_over_nc' in dictionary.keys():
        dictionary['omega_p'] =\
            dictionary['omega_0']*sqrt(dictionary['n_over_nc'])


def _compute_simulation_parameters(dictionary):
    """
    Utility that computes the relevant physical parameters
    starting from the input file.
    """
    dictionary['dx'] = 1./dictionary['k0']
    if dictionary['n_dimensions'] == 1:
        return
    if 'yx_rat' in dictionary.keys():
        dictionary['dy'] = dictionary['yx_rat']*dictionary['dx']
    else:
        dictionary['dy'] = dictionary['dx']

    if dictionary['n_dimensions'] == 2:
        return
    if 'zx_rat' in dictionary.keys():
        dictionary['dz'] = dictionary['zx_rat']*dictionary['dx']
    else:
        dictionary['dz'] = dictionary['dy']


def _translate_filename(fname):
    """
    Utility that connects the fields and particles informations with the
    corresponding output file
    """

    return(_translated_filenames[fname])


def _grid_convert(timestep_dic, timestep, Directories, path, params, **kwargs):
    from ..datas.Field import _Electromagnetic_fields
    from ..fastread.parameter_read import _read_box_limits

    for key, value in timestep_dic.items():
        if timestep == value:
            folder = key
    for elem in Directories._filelist([folder]):
        for em_field in _Electromagnetic_fields:
            if em_field in elem:
                file_path = os.path.join(path, folder, elem)
                break
    box_limits = _read_box_limits(file_path)
    grid_point = list()
    if 'x' in kwargs:
        x = kwargs['x']
        dx = params['dx']
        grid_point += [int((x-box_limits['x_min'])/(dx*params['jump']))]
    if 'y' in kwargs:
        y = kwargs['y']
        dy = params['dy']
        grid_point += [int((y-box_limits['y_min'])/(dy*params['jump']))]
    if params['n_dimensions'] == 3:
        if 'z' in kwargs:
            z = kwargs['z']
            dz = params['dz']
            grid_point += [int((z-box_limits['z_min'])/(dz*params['jump']))]
    return grid_point
