import os
import sys
from scipy.constants import pi
from numpy import sqrt
import numpy as np

module_path = os.path.dirname(sys.modules[__name__].__file__)
module_path = os.path.join(module_path, os.pardir)

speed_of_light = 0.3
r_e = 2.81794033

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
_translated_filenames['rho_electrons'] = 'Edenout'
_translated_filenames['rho_protons'] = 'Pdenout'
for n in range(1, 6):
    _translated_filenames['rho_ion'+str(n)] = 'H'+str(n)+'dnout'
    _translated_filenames['energy_ion'+str(n)] = 'H'+str(n)+'enout'

_translated_filenames['electron_energy'] = 'Elenout'
_translated_filenames['proton_energy'] = 'Prenout'
_translated_filenames['phase_space_electrons'] = 'Elpout'
_translated_filenames['phase_space_protons'] = 'Prpout'
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


def _compute_physical_parameters(dictionary):
    """
    Utility that computes the relevant physical parameters
    starting from the input file.
    """
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
        if 'lam0' in dictionary.keys():
            dictionary['n_crit'] = pi/(r_e*dictionary['lam0']**2)
            dictionary['n0'] = dictionary['n_crit']*dictionary['n_over_nc']
            dictionary['n0'] = dictionary['n0']*1.E3*1.E18
            dictionary['n_crit'] = dictionary['n_crit']*1.E21


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


def _grid_convert(box_limits, params, **kwargs):

    grid_point = dict()
    if params['str_flag'] == 0:
        stretched = False
    elif (params['str_flag'] == 1) or (params['str_flag'] == 2):
        stretched = True
    if 'x' in kwargs:
        x = kwargs['x']
        dx = params['dx']
        comp = 'x'
        grid_point['x'] = int((x-box_limits['x_min'])/(dx*params['jump']))
    if 'y' in kwargs:
        y = kwargs['y']
        if not stretched:
            dy = params['dy']
            grid_point['y'] = int((y-box_limits['y_min'])/(dy*params['jump']))
        else:
            comp = 'y'
            grid_point['y'] = int(_transverse_stretch(y, params, comp))

    if params['n_dimensions'] == 3:
        if 'z' in kwargs:
            z = kwargs['z']
            if not stretched:
                dz = params['dz']
                grid_point['z'] = int((z-box_limits['z_min']) /
                                      (dz*params['jump']))
            else:
                comp = 'z'
                grid_point['z'] = int(_transverse_stretch(z, params, comp))
    return grid_point


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


def _translate_filename(fname):
    """
    Utility that connects the fields and particles informations with the
    corresponding output file
    """

    return(_translated_filenames[fname])


def _translate_timestep(timestep, timestep_dic):
    for key, value in timestep_dic.items():
        if timestep == value:
            folder = key
    return folder


def _transverse_stretch(x, params, comp):
    """
    Function that wraps the inverse stretch functions.
    It takes as argument a coordinate and computes the grid index.
    """
    nl = params['n'+comp]
    dl = params['d'+comp]

    if params['str_flag'] == 1:
        ns = nl/6
    elif params['str_flag'] == 2:
        ns = nl/4

    return _inverse_stretch_function(x, nl, ns, dl)


def _stretch_function(x, nx, ns, dx):
    """
    Stretch function.

        / (5 dx ns)/(2 pi) tan(2 pi/5 (n-ns)/ns) + ns dx + C,   if 0 < n < ns

f(x)= <  dx n + C,  if ns < n < nt-ns

        \ dx (nt-ns) + (5 dx ns)/(2 pi) tan(2 pi/5 (n-(nt-ns))/ns) + C , else

    nt : total number of cells in the transverse direction.
    ns : total number of stretched cells
        ns = nt/6 if str_flag = 1
        ns = nt/4 if str_flag = 2
    """
    max_angle = 2*pi/5  # Arbitrary choice, can be changed
    const = -nx*dx/2    # Assumes that grid is centered around zero
    if type(x) is np.ndarray:
        x = x.astype(float)
    elif type(x) is int:
        x = float(x)

    def f(x):

        stretch = dx*ns/max_angle*np.tan(max_angle*(x-ns)/ns)+ns*dx+const
        return stretch

    def g(x):

        stretch = dx*x+const
        return stretch

    symm_center = g(nx/2)
    str_func =\
        np.piecewise(x, [x < ns, (x >= ns) & (x <= nx-ns)],
                     [lambda x: f(x), lambda x: g(x),
                      lambda x: 2*symm_center-f(nx-x)])

    return str_func


def _inverse_stretch_function(x, nx, ns, dx):
    """
    Inverse of the stretch function.
    """
    max_angle = 2*pi/5  # Arbitrary choice, can be changed
    const = nx*dx/2    # Assumes that grid is centered around zero
    xs = ns*dx

    def f(x):

        stretch = ns/max_angle*np.arctan(max_angle*((x+const-xs)/xs))+ns
        if type(stretch) is np.ndarray:
            stretch.astype(int)
        elif type(stretch) is float:
            stretch = int(stretch)
        return stretch

    def g(x):

        stretch = (x+const)/dx
        if type(stretch) is np.ndarray:
            stretch.astype(int)
        elif type(stretch) is float:
            stretch = int(stretch)
        return stretch

    symm_center = nx*dx/2-const
    alpha = (nx-2*ns)*dx-const/2

    inv_str_func =\
        np.piecewise(x, [x < -alpha, (x >= -alpha) & (x <= alpha), x > alpha],
                     [lambda x: f(x), lambda x: g(x),
                      lambda x: nx-f(2*symm_center-x)])

    return inv_str_func
