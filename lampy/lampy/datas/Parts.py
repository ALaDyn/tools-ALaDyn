from ..compiled_cython.read_phase_space import total_phase_space_read
import matplotlib.pylab as plt
import numpy as np

_phase_spaces = ['phase_space', 'phase_space_ionization',
                 'phase_space_high_energy']
m_e = 0.511
e = 1.6E-7


class Particles(object):

    def __init__(self, Simulation):

        self._Simulation = Simulation
        self._params = Simulation.params
        self._directories = Simulation.directories
        self._timesteps = Simulation._timesteps
        self._dimensions = Simulation._dimensions
        self._path = Simulation.path
        self.normalized = False
        self.comoving = False
        self._Directories = Simulation._Directories
        self._stored_phase_space = dict()

    def _search_ps_by_timestep(self, timestep):

        phase_space_list = list()
        for key in self._stored_phase_space.keys():
            if timestep in key and key[0] not in phase_space_list:
                phase_space_list += [key[0]]

        return phase_space_list

    def _search_ps_by_ps_name(self, phase_space_name):

        phase_space_list = list()
        for key in self._stored_fields.keys():
            if phase_space_name in key and key[1] not in phase_space_list:
                phase_space_list += [key[1]]

        return phase_space_list

    def _return_phase_space(self, phase_space_name, timestep):

        phase_space_list = self._search_ps_by_timestep(timestep)
        if phase_space_name in phase_space_list:
            pass
        else:
            self._phase_space_read(phase_space_name, timestep)

    def _phase_space_read(self, phase_space_name, timestep):

        sim = self._Simulation
        ndim = self._params['n_dimensions']
        file_path = sim._derive_file_path(phase_space_name, timestep)
        phase_space = dict()
        ps = total_phase_space_read(file_path, self._params)

        if ndim == 3:
            phase_space['x'] = ps[0]
            phase_space['y'] = ps[1]
            phase_space['z'] = ps[2]
            phase_space['px'] = ps[3]
            phase_space['py'] = ps[4]
            phase_space['pz'] = ps[5]
            phase_space['weight'] = ps[6]
            gamma = np.sqrt(1+ps[3]**2+ps[4]**2+ps[5]**2)
            phase_space['gamma'] = gamma
        elif ndim == 2:
            phase_space['x'] = ps[0]
            phase_space['y'] = ps[1]
            phase_space['px'] = ps[2]
            phase_space['py'] = ps[3]
            phase_space['weight'] = ps[4]
            gamma = np.sqrt(1+ps[3]**2+ps[2]**2)
            phase_space['gamma'] = gamma

        self._stored_phase_space[(phase_space_name, timestep)] = phase_space

    def scatter(self, phase_space, timestep=None, component1='x',
                component2='y', comoving=False, **kwargs):

        accepted_types = [str, dict]

        if type(phase_space) not in accepted_types:
            print("""
        Input phase_space must be either a string with the phase_space name
        or a dictionary
                  """)
            return

        if type(phase_space) is str:
            if timestep is None:
                print("""
        Time not known, impossible to plot datas.
        Please specify a time variable.
                      """)
                return
            file_path = self._Simulation._derive_file_path(phase_space,
                                                           timestep)
            file_path = file_path+'.bin'

        if type(phase_space) is str:
            self._return_phase_space(phase_space, timestep)
            ps = self._stored_phase_space[(phase_space, timestep)]
        elif type(phase_space) is dict:
            ps = phase_space

        if comoving:
            if timestep is None:
                print("Time not known, impossible to shift datas.")
            else:
                v = self._params['w_speed']
                ps['x'] = ps['x']-v*timestep

        plt.scatter(ps[component1], ps[component2], s=1, **kwargs)

    def get_data(self, phase_space_name, timestep):

        ps = dict()
        self._return_phase_space(phase_space_name, timestep)
        ps['data'] = self._stored_fields[(phase_space_name, timestep)]
        ps['time'] = timestep

        return ps

    def bunch_analysis(self, phase_space, timestep, **kwargs):

        n_dimensions = self._params['n_dimensions']
        dx = self._params['dx']
        dy = self._params['dy']
        w0 = self._params['w0_y']
        n0 = self._params['n0']*1.E-6

        pi = np.pi

        if n_dimensions == 3:
            dz = self._params['dz']

        all_particles = True

        accepted_types = [str, dict]

        if type(phase_space) not in accepted_types:
            print("""
        Input phase_space must be either
        a string with the phase_space name or a dictionary.
                  """)
            return

        if type(phase_space) is str:
            file_path = self._Simulation._derive_file_path(phase_space,
                                                           timestep)
            file_path = file_path+'.bin'

        if type(phase_space) is str:
            self._return_phase_space(phase_space, timestep)
            ps = self._stored_phase_space[(phase_space, timestep)]
        elif type(phase_space) is dict:
            ps = phase_space

        possible_kwargs = ['gamma_min', 'gamma_max', 'x_min', 'x_max',
                           'y_min', 'y_max', 'z_min', 'z_max',
                           'weight_min', 'weight_max']

        for kw in possible_kwargs:
            if kw in kwargs:
                all_particles = False
                break

        if(not all_particles):
            ps = _select_particles(ps, self._params, **kwargs)

        n_parts = len(ps['x'])

        weight_sum = np.sum(ps['weight'])
        x_ave = np.sum(ps['weight']*ps['x'])/weight_sum
        y_ave = np.sum(ps['weight']*ps['y'])/weight_sum
        px_ave = np.sum(ps['weight']*ps['px'])/weight_sum
        py_ave = np.sum(ps['weight']*ps['py'])/weight_sum
        y_diff = ps['y']-y_ave
        py_diff = ps['py']-py_ave
        x_square = np.sum(ps['weight']*ps['x']**2)/weight_sum
        y_square = np.sum(ps['weight']*ps['y']**2)/weight_sum
        px_square = np.sum(ps['weight']*ps['px']**2)/weight_sum
        py_square = np.sum(ps['weight']*ps['py']**2)/weight_sum
        gamma_ave = np.sum(ps['weight']*ps['gamma'])/weight_sum
        gamma_square = np.sum(ps['weight']*ps['gamma']**2)/weight_sum

        if n_dimensions == 3:
            z_ave = np.sum(ps['weight']*ps['z'])/weight_sum
            z_square = np.sum(ps['weight']*ps['z']**2)/weight_sum
            z_diff = ps['z']-z_ave
            pz_ave = np.sum(ps['weight']*ps['pz'])/weight_sum
            pz_diff = ps['pz']-pz_ave
            pz_square = np.sum(ps['weight']*ps['pz']**2)/weight_sum

        sigma_x = x_square-x_ave**2
        sigma_y = y_square-y_ave**2
        sigma_px = px_square-px_ave**2
        sigma_py = py_square-py_ave**2
        sigma_gamma = gamma_square-gamma_ave**2
        if n_dimensions == 3:
            sigma_z = z_square-z_ave**2
            sigma_pz = pz_square-pz_ave**2
        y_py_corr = np.sum(ps['weight']*y_diff*py_diff)/weight_sum
        if n_dimensions == 3:
            z_pz_corr = np.sum(ps['weight']*z_diff*pz_diff)/weight_sum

        emittance_y = np.sqrt(sigma_y*sigma_py-y_py_corr**2)

        energy_spread = np.sqrt(sigma_gamma)/gamma_ave

        if n_dimensions == 2:
            charge = e*n0*dx*dy*w0*weight_sum*pi/2
        elif n_dimensions == 3:
            charge = e*n0*dx*dy*dz*weight_sum
            emittance_z = np.sqrt(sigma_z*sigma_pz-z_pz_corr**2)
        print("""
        The total number of particles is {}.
        The sum on the weights of all the selected particles is {}.
        The total charge measured is Q={:4.1e}pC.
        The normalized emittance along the first perpendicular axis
        is eps={} mm mrad
              """.format(n_parts, weight_sum, charge, emittance_y))
        if n_dimensions == 3:
            print("""
        The normalized emittance along the second perpendicular
        axis is eps={}mm mrad
                  """.format(emittance_z))
        print("""
        The mean energy is E={}MeV, while the bunch energy spread is
        rms(E)/mean(E)={}%
              """.format(m_e*gamma_ave, energy_spread*100))
        if(n_dimensions == 3):
            print("""
        <x>={} <y>={}, <z>={}
        <(x-<x>)^2>={}, <(y-<y>)^2>={}, <(z-<z>)^2>={}
        <px>={}, <py>={}, <pz>={},
        <(px-<px>)^2>={}, <(py-<py>)^2>={}, <(pz-<pz>)^2>={}
                  """.format(x_ave, y_ave, z_ave, sigma_x, sigma_y,
                             sigma_z, px_ave, py_ave, pz_ave,
                             sigma_px, sigma_py, sigma_pz))

        if(n_dimensions == 2):
            print("""
        <x>={} <y>={},
        <(x-<x>)^2>={}, <(y-<y>)^2>={},
        <px>={}, <py>={},
        <(px-<px>)^2>={}, <(py-<py>)^2>={},
                  """.format(x_ave, y_ave, sigma_x, sigma_y, px_ave, py_ave,
                             sigma_px, sigma_py))

    def get_data(self, phase_space_name, timestep):

        ps = dict()
        temp = self._return_phase_space(phase_space_name, timestep)
        ps['data'] = temp
        ps['time'] = timestep


def _select_particles(ps, params, **kwargs):

    n_dimensions = params['n_dimensions']

    possible_kwargs = ['gamma_min', 'gamma_max', 'x_min', 'x_max',
                       'y_min', 'y_max', 'z_min', 'z_max',
                       'weight_min', 'weight_max']

    var_dic_min = dict()
    var_dic_max = dict()
    var_dic_min['gamma_min'] = 'gamma'
    var_dic_max['gamma_max'] = 'gamma'
    var_dic_min['x_min'] = 'x'
    var_dic_min['y_min'] = 'y'
    var_dic_min['z_min'] = 'z'
    var_dic_max['x_max'] = 'x'
    var_dic_max['y_max'] = 'y'
    var_dic_max['z_max'] = 'z'
    var_dic_min['weight_min'] = 'weight'
    var_dic_max['weight_max'] = 'weight'

    if 'z_min' in kwargs and n_dimensions == 2:
        del kwargs['z_min']
    if 'z_max' in kwargs and n_dimensions == 2:
        del kwargs['z_max']

    index = list()
    kw_list = list(set(possible_kwargs).intersection(kwargs.keys()))

    if not kw_list:
        return ps

    for kw in kw_list:
        if kw in var_dic_min.keys():
            index.append(np.where(ps[var_dic_min[kw]] > kwargs[kw])[0])
        elif kw in var_dic_max.keys():
            index.append(np.where(ps[var_dic_max[kw]] < kwargs[kw])[0])

    tot_index = index[0]
    for i in range(len(index)):
        tot_index = set(tot_index).intersection(index[i])
    ps_selected = dict()
    tot_index = list(tot_index)
    for key in ps.keys():
        ps_selected[key] = ps[key][tot_index]
    return ps_selected
