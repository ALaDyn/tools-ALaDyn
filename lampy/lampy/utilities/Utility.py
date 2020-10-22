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
_translated_filenames['Jx'] = 'denvxout'
_translated_filenames['Jy'] = 'denvyout'
_translated_filenames['Jz'] = 'denvzout'
_translated_filenames['A'] = 'Aenvout'
_translated_filenames['a'] = 'aenvout'
_translated_filenames['ReA'] = 'Renvout'
_translated_filenames['ImA'] = 'Ienvout'
_translated_filenames['rho_fluid'] = 'Fdenout'
_translated_filenames['Px_fluid'] = 'Flpxout'
_translated_filenames['Py_fluid'] = 'Flpyout'
_translated_filenames['Pz_fluid'] = 'Flpzout'
_translated_filenames['rho_electrons'] = 'Edenout'
_translated_filenames['rho_protons'] = 'Pdenout'
for n in range(1, 6):
    _translated_filenames['rho_ion'+str(n)] = 'H'+str(n)+'dnout'
    _translated_filenames['energy_ion'+str(n)] = 'H'+str(n)+'enout'

_translated_filenames['electrons_energy'] = 'Elenout'
_translated_filenames['protons_energy'] = 'Prenout'
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
_total_filenamelist += ['denvxout']
_total_filenamelist += ['denvyout']
_total_filenamelist += ['denvzout']
_total_filenamelist += ['Fdenout']
_total_filenamelist += ['Flpxout']
_total_filenamelist += ['Flpyout']
_total_filenamelist += ['Flpzout']
_total_filenamelist += ['Edenout']
_total_filenamelist += ['Elenout']
_total_filenamelist += ['Elpout']
_total_filenamelist += ['Eionzout']
_total_filenamelist += ['E_hg_out']
_total_filenamelist += ['H'+str(n)+'dnout' for n in range(1, 6)]
_total_filenamelist += ['H'+str(n)+'enout' for n in range(1, 6)]
_total_filenamelist += ['Track_']

_tracking_directory = 'tracking'
_tracking_dictionary = dict()
_tracking_basename = dict()
for n in range(1, 6):
    _tracking_dictionary[n] = 'tracking_dictionary_'+str(n)+'.dat'

# Warning, to be extended for more tracked species
for n in range(1, 6):
    _tracking_basename[n] = 'Track_'+str(n)+'_'

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

    dictionary['particle_dimensions'] = dictionary['n_dimensions']
    if dictionary['model_id'] == 2 or dictionary['model_id'] == 3:
        dictionary['particle_dimensions'] = 3
    if 'lam0' in dictionary.keys():
        dictionary['omega_0'] = 2*pi/dictionary['lam0']
    if 'lam1' in dictionary.keys() and dictionary['lam1'] > 0:
        dictionary['omega_1'] = 2*pi/dictionary['lam1']
    if 'n_over_nc' in dictionary.keys():
        dictionary['omega_p'] =\
            dictionary['omega_0']*sqrt(dictionary['n_over_nc'])
        if 'lam0' in dictionary.keys():
            dictionary['n_crit'] = pi/(r_e*dictionary['lam0']**2)
            dictionary['n0'] = dictionary['n_crit']*dictionary['n_over_nc']
            dictionary['n0'] = dictionary['n0']*1.E3*1.E18
            dictionary['n_crit'] = dictionary['n_crit']*1.E21
            dictionary['n_reference'] = dictionary['n0']
    if 'n0_ref' in dictionary.keys():
        try:
            lambda_p = 33*np.sqrt(1/dictionary['n0_ref'])
            dictionary['omega_p'] = 2*pi/lambda_p
        except ZeroDivisionError:
            dictionary['omega_p'] = 0.
        dictionary['n_reference'] = 1.e18
        if 'lam0' in dictionary.keys():
            dictionary['n_crit'] = pi/(r_e*dictionary['lam0']**2)
            dictionary['n0'] = dictionary['n0_ref']*dictionary['n_reference']
    dictionary['dt'] = (1./dictionary['k0'])*dictionary['cfl']


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

    if dictionary['str_flag'] > 0:
        dictionary['stretched'] = True
    if 'a_on_particles' not in dictionary.keys():
        dictionary['a_on_particles'] = [False]*dictionary['nsp']



def _convert_component_to_index(params, component):

    ndim = params['particle_dimensions']
    if ndim == 3:
        switch = {
            'x': 0,
            'y': 1,
            'z': 2,
            'px': 3,
            'py': 4,
            'pz': 5,
            'weight': 6,
            'gamma': 7,
            'index': 8,
            'a': 9
        }
    elif ndim == 2:
        switch = {
            'x': 0,
            'y': 1,
            'px': 2,
            'py': 3,
            'weight': 4,
            'gamma': 5,
            'index': 6,
            'a': 7
        }
    return switch[component]


def _find_inputs_json(path):
    """
    Utility that finds and classifies the input files found
    in the folder.
    """
    inputs = []
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            if 'input_' in f and '.json' in f:
                inputs.append(os.path.join(dirpath, f))
        break

    if len(inputs) == 0:
        raise FileNotFoundError('No input JSON have been found in {}'
                                .format(path))

    return inputs

def _find_inputs_nml(path):
    """
    Utility that finds and classifies the input files found
    in the folder.
    """
    inputs = []
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            if 'input_' in f and '.nml' in f:
                inputs.append(os.path.join(dirpath, f))
        break

    if len(inputs) == 0:
        raise FileNotFoundError('No input namelist have been found in {}'
                                .format(path))

    return inputs


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
            grid_point['y'] = int(_transverse_stretch(y, params, comp)/params['jump'])

    if params['n_dimensions'] == 3:
        if 'z' in kwargs:
            z = kwargs['z']
            if not stretched:
                dz = params['dz']
                grid_point['z'] = int((z-box_limits['z_min']) /
                                      (dz*params['jump']))
            else:
                comp = 'z'
                grid_point['z'] = int(_transverse_stretch(z, params, comp)/params['jump'])
    return grid_point

def _nearest_particle( phase_space, component_dict):

    import numpy as np

    nparts = len(phase_space[1])
    dist = np.zeros(nparts)    
    index_comp = component_dict.pop('index')
    for comp, coord in component_dict.items():
        dist += (phase_space[comp] - coord)**2
    dist = np.sqrt(dist)

    minloc = np.argmin(dist)

    return phase_space[index_comp][minloc]

def _read_simulation_json(path):
    """
    Utility that reads the 'input_??.json' file
    in the given folder to assign the simulation
    parameters.
    It requires the package json to be installed.
    """
    import json

    param_json = dict()
    param_dic = dict()
    inputs = _find_inputs_json(path)
    if len(inputs) > 1:
        print("LAMPy not yet prepared to read more outputs.")
        for name in inputs[1:]:
            inputs.remove(name)

    for f in inputs:
        name = f[-6:-4]
        with open(f, 'r') as fp:
            param_json[name] = json.load(fp)

    names_view = param_json.keys()
    names_iterator = iter(names_view)
    first_name = next(names_iterator)
    for key in param_json[first_name].keys():
        for (key2, value2) in param_json[first_name][key].items():
            param_dic[key2] = value2

    return param_dic

def _read_simulation_nml(path):
    """
    Utility that reads the 'input_??.nml' file
    in the given folder to assign the simulation
    parameters.
    It requires the package f90nml to be installed.
    """
    import f90nml

    param_dic = dict()
    inputs = _find_inputs_nml(path)
    if len(inputs) > 1:
        print("LAMPy not yet prepared to read more outputs.")
        for name in inputs[1:]:
            inputs.remove(name)

    n = dict()
    for f in inputs:
        name = f[-6:-4]
        n[name] = f90nml.read(f)

    for input_name in n.keys():
        for namelist_name in n[input_name]:
            for key, item in n[input_name][namelist_name].items():
                param_dic[key] = item

    return param_dic


def _read_simulation_without_nml(path):
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


def _sort_particles(phase_space, component):

    sorted_index = np.argsort(phase_space[component])
    ps_sorted = dict()
    for key in phase_space.keys():
        ps_sorted[key] = phase_space[key][sorted_index]

    return ps_sorted

def _sort_tracked_particles(phase_space, index_in):

    sorted_index = np.argsort(phase_space[index_in])
    return phase_space[:,sorted_index]

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

        ( (5 dx ns)/(2 pi) tan(2 pi/5 (n-ns)/ns) + ns dx + C,   if 0 < n < ns
       (
f(x)= <  dx n + C,  if ns < n < nt-ns
       (
        ( dx (nt-ns) + (5 dx ns)/(2 pi) tan(2 pi/5 (n-(nt-ns))/ns) + C , else

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

def moving_average(times, values, simulation):

    dt = simulation.params['dt']
    every = simulation.params['every_track'][0]
    dt = every*dt
    lam0 = simulation.params['lam0']
    npoints_in_lam = int(lam0/dt + 0.5)
    vals = values.copy()
    N = npoints_in_lam
    offset_left = N//2
    offset_right = N-N//2
    vals = np.insert(vals, 0, offset_left*[vals[0]])
    vals = np.append(vals, offset_right*[vals[-1]])
    averaged = np.zeros_like(vals)
    cumsum = np.zeros_like(vals)
    cumsum = np.append(cumsum, [0])
    for i, x in enumerate(vals, 1):
        cumsum[i] = cumsum[i-1] + x
        if i >= N:
            averaged[i-offset_right] = (cumsum[i] - cumsum[i-N])/N
    averaged = averaged[offset_left:-offset_right]

    return averaged