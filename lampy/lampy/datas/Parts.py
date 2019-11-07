from ..compiled_cython.read_phase_space import total_phase_space_read
import matplotlib.pylab as plt
import numpy as np
from ..utilities.Utility import speed_of_light
import os

_phase_spaces = ['phase_space_electrons', 'phase_space_ionization',
                 'phase_space_high_energy', 'phase_space_protons']
m_e = 0.511
e = 1.6E-7


class Particles(object):

    def __init__(self, Simulation):
        """
        Class that contains all the methods and datas related to the
        particles phase space.

        With the available methods it is possible to access, manipulate and
        plot all the produced phase spaces.

        Any phase space ps is a dictionary of 8 (3D) or 6 (2D) components.
        ps['x'] particles coordinates on the x axis
        ps['y'] particles coordinates on the y axis
        ps['z'] particles coordinates on the z axis (only if 3D simulation)
        ps['px'] particles momenta along the x direction
        ps['py'] particles momenta along the y direction
        ps['pz'] particles momenta along the z direction (only in 3D)
        ps['weight'] particles computational weights
        ps['gamma'] particles Lorentz factor.

        Possible phase spaces are:
        - phase_space_electrons: phase space of all the electrons
        - phase_space_ionization: phase space of all the particles produced
          via the ionization process
        - phase_space_high_energy: phase space of all the particles with
          gamma > gamma_min.
        - phase_space_protons: phase space of all the protons
        """
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
        self._dx = self._params['dx']
        self._selected_index = dict()
        self._selected_percentage = 1

    def _search_ps_by_timestep(self, timestep):

        phase_space_list = list()
        for key in self._stored_phase_space.keys():
            if timestep in key and key[0] not in phase_space_list:
                phase_space_list += [key[0]]

        return phase_space_list

    def _search_ps_by_ps_name(self, phase_space_name):

        phase_space_list = list()
        for key in self._stored_phase_space.keys():
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

        if phase_space_name not in sim.output:
            print("""
        {} is not available.
        Available output are {}.
        """.format(phase_space_name, sim.output))
            return

        file_path = sim._derive_file_path(phase_space_name, timestep)
        try:
            phase_space = dict()
            ps, part_number = total_phase_space_read(file_path, self._params)
        except FileNotFoundError:
            print("""
        Phase space {} not available, impossible to read.
            """.format(phase_space_name))
            raise

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

        self._selected_index[(phase_space_name, timestep)] =\
            [part_number, np.arange(part_number)]
        self._stored_phase_space[(phase_space_name, timestep)] = phase_space

    def scatter(self, phase_space, time=None, component1='x',
                component2='y', comoving=False,
                selected_percentage=None, **kwargs):
        """
        Method that produces a scatter plot of the given phase space.

        It takes as input the phase space type (e.g. all electrons, energetic
        electrons, etc.) and the time of interest.

        Parameters
        --------
        phase_space : str or dict
            If it is a string, is the phase space name.
            Otherwise, a phase space dictionary (i.e. a phase space collected
            via a get_data()) can be given.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        time : float, optional (default None)
            Time at the phase space is plotted. It is necessary if the phase
            space name is passed, it is optional if the phase space dictionary
            is passed. However, it may be necessary if the longitudinal
            axis is set on 'comoving'.
        component1 : str, optional
            First component of the scatter plot. Default is taken as 'x'.
        component2 : str, optional
            Second component of the scatter plot. Default is taken as 'y'
        comoving : bool, optional
            If True, the longitudinal axis is transformed as xi=x-v_w t.
            Remember that, to obtain the comoving axis, the time is needed
            even when the phase space dictionary is passed.
        selected_percentage : float, optional
            A number between 0 and 1. It determines the percentage of particles
            that have to be selected from the phase space before to be plotted.
            It can be important for heavily populated phase spaces, that are
            difficult to load. Once the particles have been selected, they are
            always kept between various plots, even if the component is
            changed until the percentage is changed and a new selection happens

        Kwargs
        --------
        All the kwargs for the pyplot.scatter function.
        """
        accepted_types = [str, dict]

        if type(phase_space) not in accepted_types:
            print("""
        Input phase_space must be either a string with the phase_space name
        or a dictionary
                  """)
            return

        if type(phase_space) is str:
            if time is None:
                print("""
        Time not known, impossible to plot datas.
        Please specify a time variable.
                      """)
                return
            time = self._Simulation._nearest_time(time)

        if type(phase_space) is str:
            self._return_phase_space(phase_space, time)
            ps = self._stored_phase_space[(phase_space, time)].copy()
        elif type(phase_space) is dict:
            ps = phase_space.copy()

        if comoving:
            if time is None:
                print("Time not known, impossible to shift datas.")
                return
            else:
                v = self._params['w_speed']
                ps['x'] = ps['x']-v*time

        if (selected_percentage != self._selected_percentage)\
                and (selected_percentage is not None):

            self._selected_percentage = selected_percentage
            part_number = self._selected_index[(phase_space, time)][0]
            sel_part_numb = int(selected_percentage * part_number)
            sel_index =\
                np.random.choice(np.arange(part_number), size=sel_part_numb,
                                 replace=False)
            sel_index.sort()
            self._selected_index[(phase_space, time)][1] = sel_index

        inds = self._selected_index[(phase_space, time)][1]
        plt.scatter(ps[component1][inds], ps[component2][inds], s=1, **kwargs)

    def get_data(self, phase_space_name, timestep):
        """
        Method that retrieves any phase space data
        and returns it as a dictionary.

        Parameters
        --------
        phase_space_name : str
            Name of the phase space that is to be collected
        time : float
            Instant at which dictionary is retrieved

        Returns
        --------
        phase space : dict
            Data are returned as a dictionary.
            phase space ['data'] is the dictionary containing the data
            phase space['time'] is the corresponding time .

            phase space ['data'] is again a dictionary, with possible keys:
            'x', 'y', 'z', 'px', 'py', 'pz', 'weight', 'gamma'
        """
        ps = dict()
        timestep = self._Simulation._nearest_time(timestep)
        self._return_phase_space(phase_space_name, timestep)
        ps['data'] = self._stored_phase_space[(phase_space_name, timestep)]
        ps['time'] = timestep

        return ps

    def bunch_analysis(self, phase_space, time, **kwargs):
        """
        Method that analyzed the given phase space, computing all it's
        statistical properties.

        Parameters
        --------
        phase_space : str or dict
            If it is a string, is the phase space name.
            Otherwise, a phase space dictionary (i.e. a phase space collected
            via a get_data()) can be given.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        time : float
            Time at which the phase space is analyzed.

        Kwargs
        --------
        List of possible kwargs:

            'gamma_min', 'gamma_max', 'x_min', 'x_max', 'y_min', 'y_max',
            'z_min', 'z_max', 'weight_min', 'weight_max'

            It is possible to select parts of phase space to analyze via
            the input kwargs.
        """
        n_dimensions = self._params['n_dimensions']
        dx = self._params['dx']
        dy = self._params['dy']
        w0 = self._params['w0_y']
        n0 = self._params['n0']*1.E-12

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
            self._return_phase_space(phase_space, time)
            ps = self._stored_phase_space[(phase_space, time)]
        elif type(phase_space) is dict:
            ps = phase_space

        possible_kwargs = ['gamma_min', 'gamma_max', 'x_min', 'x_max',
                           'y_min', 'y_max', 'z_min', 'z_max',
                           'weight_min', 'weight_max']

        for kw in possible_kwargs:
            if kw in kwargs:
                all_particles = False
                break

        if not all_particles:
            ps = self.select_particles(ps, **kwargs)

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
        The total charge measured is Q={:4.2e}pC.
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

    def slice_analysis(self, phase_space, time, number_of_slices=100,
                       filename='slice_analysis.dat', **kwargs):
        """
        Method that computes the slice analysis of a given particle bunch.

        Parameters
        --------
        phase_space : str or dict
            If it is a string, is the phase space name.
            Otherwise, a phase space dictionary (i.e. a phase space collected
            via a get_data()) can be given.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        time : float
            Time at which the slice analysis is performed.
        number_of_slices : int, optional
            Number of slices in which the bunch is divided.
            Defatul value is 100.
        filename : str, optional
            Name of the file on which results are written.
            Default is 'slice_analysis.dat'.

        Kwargs
        --------
        List of possible kwargs:

            'x_min', 'x_max', 'slice_length'

        x_min : float
            x value from which the analysis starts.
        x_max : float
            x value in which the analysis ends.
        If neither x_min nor x_max are specified, the analysis is performed on
        the whole bunch.

        slice_length : float
            lenght of every slice
        """
        n_dimensions = self._params['n_dimensions']
        dx = self._params['dx']
        dy = self._params['dy']
        w0 = self._params['w0_y']
        n0 = self._params['n0']*1.E-12
        if n_dimensions == 3:
            dz = self._params['dz']
        accepted_types = [str, dict]

        path = os.path.join(self._path, filename)
        f = open(path, 'w')

        if type(phase_space) not in accepted_types:
            print("""
        Input phase_space must be either a string with the phase_space name
        or a dictionary
                  """)
            return

        if type(phase_space) is str:
            if time is None:
                print("""
        Time not known, impossible to plot datas.
        Please specify a time variable.
                      """)
                return

        if type(phase_space) is str:
            self._return_phase_space(phase_space, time)
            ps = self._stored_phase_space[(phase_space, time)]
        elif type(phase_space) is dict:
            ps = phase_space

        if 'x_min' in kwargs:
            x_min = kwargs['x_min']
        else:
            x_min = np.amin(ps['x'])
        if 'x_max' in kwargs:
            x_max = kwargs['x_max']
        else:
            x_max = np.amax(ps['x'])
        if 'slice_length' in kwargs:
            slice_length = kwargs['slice_length']
        else:
            slice_length = 4*self._params['dx']

        delta = float(x_max-x_min)/number_of_slices
        slice_length_fs = slice_length/speed_of_light

        sorted_index = np.argsort(ps['x'])
        ps_sorted = dict()
        for key in ps.keys():
            ps_sorted[key] = ps[key][sorted_index]

        for i in range(number_of_slices):

            x_center = x_min+delta*i
            x0 = x_center-slice_length/2.
            x1 = x_center+slice_length/2.
            min_index = np.where(ps_sorted['x'] >= x0)[0]
            max_index = np.where(ps_sorted['x'] <= x1)[0]

            if len(min_index) > 0 and len(max_index) > 0:
                min_index = min_index[0]
                max_index = max_index[-1]

                yp = ps_sorted['y'][min_index:max_index]
                pyp = ps_sorted['py'][min_index:max_index]
                wgh = ps_sorted['weight'][min_index:max_index]
                gamma_single_particle = ps_sorted['gamma'][min_index:max_index]

                if n_dimensions == 3:

                    zp = ps_sorted['z'][min_index:max_index]
                    pzp = ps_sorted['pz'][min_index:max_index]

                weight_sum = np.sum(wgh)
                y_ave = np.sum(wgh*yp)/weight_sum
                py_ave = np.sum(wgh*pyp)/weight_sum
                y_diff = yp-y_ave
                py_diff = pyp-py_ave
                y_square = np.sum(wgh*yp**2)/weight_sum
                py_square = np.sum(wgh*pyp**2)/weight_sum

                gamma_ave = np.sum(wgh*gamma_single_particle)/weight_sum
                gamma_square = np.sum(wgh*gamma_single_particle**2)/weight_sum
                if n_dimensions == 3:

                    z_ave = np.sum(wgh*zp)/weight_sum
                    z_square = np.sum(wgh*zp**2)/weight_sum
                    z_diff = zp-z_ave
                    pz_ave = np.sum(wgh*pzp)/weight_sum
                    pz_diff = pzp-pz_ave
                    pz_square = np.sum(wgh*pzp**2)/weight_sum

                sigma_y = y_square-y_ave**2
                sigma_py = py_square-py_ave**2
                sigma_gamma = gamma_square-gamma_ave**2
                if n_dimensions == 3:
                    sigma_z = z_square-z_ave**2
                    sigma_pz = pz_square-pz_ave**2
                y_py_corr = np.sum(wgh*y_diff*py_diff)/weight_sum
                if n_dimensions == 3:
                    z_pz_corr = np.sum(wgh*z_diff*pz_diff)/weight_sum

                emittance_y = np.sqrt(sigma_y*sigma_py-y_py_corr**2)

                energy_spread = np.sqrt(sigma_gamma)/gamma_ave
                if n_dimensions == 2:
                    charge = e*n0*dx*dy*w0*weight_sum*np.pi/2.
                if n_dimensions == 3:
                    emittance_z = np.sqrt(sigma_z*sigma_pz-z_pz_corr**2)
                    charge = e*n0*dx*dy*dz*weight_sum

                current = charge/slice_length_fs

                if n_dimensions == 3:
                    f.writelines('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\
                                 \n'.format(x_center, emittance_y, emittance_z,
                                            sigma_y, sigma_z, current,
                                            gamma_ave, energy_spread))
                if n_dimensions == 2:
                    f.writelines('{}\t{}\t{}\t{}\t{}\t{}\
                                 \n'.format(x_center, emittance_y, sigma_y,
                                            current, gamma_ave, energy_spread))

        f.close()
        slice_analysis = self.read_slice_analysis(path)

        return slice_analysis

    def read_slice_analysis(self, file_name):
        """
        Method that reads the file containing a slice analysis.

        Parameters
        --------
        file_name : str
            Name of the file
        """
        f = open(file_name, 'r')
        lines = f.readlines()

        n_dimensions = self._params['n_dimensions']

        slice_analysis = dict()
        x_center = list()
        emittance_y = list()
        sigma_y = list()
        current = list()
        gamma_ave = list()
        energy_spread = list()
        emittance_z = list()
        sigma_z = list()

        if n_dimensions == 3:
            for line in lines:
                x_center.append(float(line.split()[0]))
                emittance_y.append(float(line.split()[1]))
                emittance_z.append(float(line.split()[2]))
                sigma_y.append(float(line.split()[3]))
                sigma_z.append(float(line.split()[4]))
                current.append(float(line.split()[5]))
                gamma_ave.append(float(line.split()[6]))
                energy_spread.append(float(line.split()[7]))

        elif n_dimensions == 2:
            for line in lines:
                x_center.append(float(line.split()[0]))
                emittance_y.append(float(line.split()[1]))
                sigma_y.append(float(line.split()[2]))
                current.append(float(line.split()[3]))
                gamma_ave.append(float(line.split()[4]))
                energy_spread.append(float(line.split()[5]))

        slice_analysis['x_center'] = x_center
        slice_analysis['emittance_y'] = emittance_y
        slice_analysis['sigma_y'] = sigma_y
        slice_analysis['J'] = current
        slice_analysis['gamma'] = gamma_ave
        slice_analysis['e_spread'] = energy_spread

        if n_dimensions == 3:
            slice_analysis['emittance_z'] = emittance_z
            slice_analysis['sigma_z'] = sigma_z

        return slice_analysis

    def select_particles(self, phase_space, time=None, **kwargs):
        """
        Function that takes in input a phase space dictionary, and the array of
        parameters and selects particles according to the given conditions.

        Parameters
        --------
        phase_space : dict or str
            If it is a string, is the phase space name.
            Otherwise, a phase space dictionary (i.e. a phase space collected
            via a get_data()) can be given.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        time : float, optional
            Time variable is needed when phase_space is a string

        Kwargs
        --------
        List of possible kwargs:

                'gamma_min', 'gamma_max', 'x_min', 'x_max', 'y_min', 'y_max',
                'z_min', 'z_max', 'weight_min', 'weight_max'

                It is possible to select parts of phase space to analyze via
                the input kwargs.

        Returns
        --------
        ps_selected : dict
            Phase space dictionary with the selected particles

        """
        n_dimensions = self._params['n_dimensions']

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

        accepted_types = [str, dict]

        if type(phase_space) not in accepted_types:
            print("""
        Input phase_space must be either a string with the phase_space name
        or a dictionary
                  """)
            return

        if type(phase_space) is str:
            if time is None:
                print("""
        Time not known, impossible to plot datas.
        Please specify a time variable.
                      """)
                return
            time = self._Simulation._nearest_time(time)

        if type(phase_space) is str:
            self._return_phase_space(phase_space, time)
            ps = self._stored_phase_space[(phase_space, time)].copy()
        elif type(phase_space) is dict:
            ps = phase_space.copy()

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

    def histogram(self, phase_space, time=None, component1='x',
                  component2=None, **kwargs):
        """
        Method that generates either a 1D or a 2D histogram of a given set and
        for given components of the phase space.
        It is based on numpy histogram and, in the 2D case, represents results
        via a pyplot.pcolormesh plot, masking away all points with a zero
        contribution to the histogram.

        Parameters
        --------
        phase_space : dict or str
            If it is a string, is the phase space name.
            Otherwise, a phase space dictionary (i.e. a phase space collected
            via a get_data()) can be given.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        time : float, optional
            Time variable is needed when phase_space is a string
        component1 : str, optional
            First component of the scatter plot. Default is taken as 'x'.
        component2 : str, optional
            Second component of the scatter plot. Default is taken as None

        Kwargs
        --------
        List of possible kwargs:

                'gamma_min', 'gamma_max', 'x_min', 'x_max', 'y_min', 'y_max',
                'z_min', 'z_max', 'weight_min', 'weight_max'

                It is possible to select parts of phase space to analyze via
                the input kwargs.

        For the histogram and the plot, possible kwargs are:

            alpha, bins, cmap, density

        """
        cmap = 'Reds'
        bins = 1000
        density = True
        alpha = 1

        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
            del kwargs['cmap']

        if 'bins' in kwargs:
            bins = kwargs['bins']
            del kwargs['bins']

        if 'density' in kwargs:
            density = kwargs['density']
            del kwargs['density']

        if 'alpha' in kwargs:
            alpha = kwargs['alpha']
            del kwargs['alpha']

        ps = self.select_particles(phase_space, time=time, **kwargs)

        if component2 is None:
            _ = plt.hist(ps[component1], weights=ps['weight'],
                         bins=bins, density=density, alpha=alpha)

        else:
            H, xedge, yedge = np.histogram2d(ps[component1], ps[component2],
                                             bins=bins, weights=ps['weight'],
                                             density=density)
            H = H.T
            X, Y = np.meshgrid(xedge, yedge)
            H = np.ma.masked_where(H == 0, H)
            plt.pcolormesh(X, Y, H, cmap=cmap, alpha=alpha)
