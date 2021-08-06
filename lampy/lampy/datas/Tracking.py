from compiled_cython.read_tracking import total_tracking_read
import matplotlib.pylab as plt
import numpy as np
from ..utilities.Utility import _convert_component_to_index,\
    _tracking_directory, _tracking_dictionary, _nearest_particle
import os

m_e = 0.511
e = 1.6E-7
_can_dict = dict()
_can_dict['y'] = 'canonical momentum y'


class Tracking(object):

    def __new__(cls, Simulation):

        if not Simulation._tracking:
            if Simulation._verbose_warning:
                print("Tracking not available")
            return None
        else:
            return super(Tracking, cls).__new__(cls)

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
        self._stored_tracking = dict()
        self.iter_dictionary = dict()
        self._available_iterations = list()
        self._available_times = list()
        self._available_species = list()
        self._tracking = self._Simulation._tracking
        self._read_iter_dic()
        self._tracking_dictionary = _tracking_dictionary
        if self._params['particle_dimensions'] == 2:
            self._array_dimensions = 7
        elif self._params['particle_dimensions'] == 3:
            self._array_dimensions = 9

    def _check_index_in_ps(self, phase_space, index):

        comp = _convert_component_to_index(self._params, 'index')
        if index == 'all':
            return True
        mask = np.isin(phase_space[comp], index, assume_unique=True)
        return mask

    def _get_stored_tracking(self, timestep, species):

        if self._Simulation._save_data:
            ps = self._stored_tracking[(timestep, species)]
        else:
            ps = self._stored_tracking.pop((timestep, species))

        return ps

    def _read_iter_dic(self):

        if not self._tracking:
            return

        for key, value in _tracking_dictionary.items():
            track_dic_path =\
                os.path.join(self._path, _tracking_directory, value)

            if os.path.isfile(track_dic_path):
                self._available_species += [key]
                self.iter_dictionary[key] = dict()
            elif not os.path.isfile(track_dic_path):
                continue

            with open(track_dic_path) as fp:
                lines = fp.readlines()

            # Removing header
            lines.pop(0)
            for line in lines:
                self._available_iterations += [int(line.split()[0])]
                self._available_times += [float(line.split()[1])]
                self.iter_dictionary[key][float(line.split()[1])] =\
                    int(line.split()[0])

    def _return_tracking_phase_space(self, timestep, species):

        track_list = self._search_track_by_timestep(timestep, species)
        if track_list is not None:
            pass
        else:
            self._track_read(timestep, species)

    def _search_track_by_timestep(self, timestep, species):

        track_list = None
        for key in self._stored_tracking.keys():
            if (timestep, species) == key:
                track_list = key

        return track_list

    def _track_read(self, timestep, species):

        from ..utilities.Utility import _sort_tracked_particles
        sim = self._Simulation

        file_path = sim._derive_tracking_file_path(timestep, species)
        try:
            ps, part_number = total_tracking_read(file_path, self._params,
                                                  species)
        except FileNotFoundError:
            if self._Simulation._verbose_error:
                print("""
        Tracking at time {} not available, impossible to read.
                """.format(timestep))
            raise

        index_in = _convert_component_to_index(self._params, 'index')

        ps_sorted = _sort_tracked_particles(ps, index_in)

        self._stored_tracking[(timestep, species)] = ps_sorted

    def comp_to_index(self, index):
        return _convert_component_to_index(self._params, index)

    def get_data(self, species=1, time=0):
        """
        Method that retrieves any tracked particle data and returns
        phase space at a given time.

        Parameters
        --------
        field_name : str
            Name of the array that is to be collected
        time : float
            Instant at which array is retrieved

        Returns
        --------
        field : dict
            Data are returned as a dictionary.
            field['data'] is the array containing the data
            field['time'] is the corresponding time
        """
        f = dict()
        time = self._Simulation._nearest_tracking_time(time, species)
        self._return_tracking_phase_space(time, species)
        f['data'] = self._get_stored_tracking(time, species)
        f['time'] = time

        return f

    def get_trajectory(self, species=1, index=1):
        """
        Method that retrieves any field data and returns the 2D or 3D array.

        Parameters
        --------
        field_name : str
            Name of the array that is to be collected
        time : float
            Instant at which array is retrieved

        Returns
        --------
        field : dict
            Data are returned as a dictionary.
            field['data'] is the array containing the data
            field['index'] is the corresponding index
        """
        try:
            val = int(index)
        except ValueError:
            if self._Simulation._verbose_error:
                print("Index must be a number")
            return
        del(index)
        f = dict()
        timedict = self._Simulation._iter_dictionary[species]
        lasttime = list(timedict.keys())[-1]
        time_index = self._available_times.index(lasttime)
        time_array = self._available_times[:time_index + 1]

        pscut = list()

        # Warning, only working for single particle for now
        for instant in time_array:
            self._return_tracking_phase_space(instant, species)
            ps0 = self._get_stored_tracking(instant, species)
            ps = ps0.copy()
            mask = self._check_index_in_ps(ps, val)
            if not mask.any():
                continue
            pstemp = ps[:, mask]
            pstemp = np.append(pstemp, instant)
            pstemp = pstemp.reshape(-1, 1)
            pscut += [pstemp]

        f['data'] = np.column_stack(pscut)
        f['index'] = val

        return f

    def select_index(self, time=None, species=1, **kwargs):
        """
        Function that takes in input a phase space and
        selects particles according to the given conditions.

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
                'z_min', 'z_max', 'weight_min', 'weight_max', 'nearest',
                'x', 'y', 'z'

                It is possible to select parts of phase space to analyze via
                the input kwargs.

        Returns
        --------
        ps_selected : np.array
            Phase space with the selected particles

        """
        n_dimensions = self._params['particle_dimensions']

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

        nearest_part = False
        if 'nearest' in kwargs:
            if kwargs['nearest']:
                nearest_part = True
            if 'x' not in kwargs and 'y' not in kwargs and 'z' not in kwargs:
                if self._Simulation._verbose_error:
                    print("""
        Error, when you ask 'nearest' you should specify the component
        'x', 'y' or 'z'
                        """)
                return
        comp_dict = dict()
        if 'x' in kwargs:
            comp = _convert_component_to_index(self._params, 'x')
            comp_dict[comp] = kwargs['x']
        if 'y' in kwargs:
            comp = _convert_component_to_index(self._params, 'y')
            comp_dict[comp] = kwargs['y']
        if 'z' in kwargs:
            comp = _convert_component_to_index(self._params, 'z')
            comp_dict[comp] = kwargs['z']
        comp_dict['index'] = _convert_component_to_index(self._params, 'index')

        index = list()
        kw_list = list(set(possible_kwargs).intersection(kwargs.keys()))

        time = self._Simulation._nearest_tracking_time(time, species)
        self._return_tracking_phase_space(time, species)

        ps0 = self._get_stored_tracking(time, species)
        ps = ps0.copy()

        if nearest_part:
            ind = _nearest_particle(ps, comp_dict)
            return ind

        for kw in kw_list:
            if kw in var_dic_min.keys():
                comp = _convert_component_to_index(self._params,
                                                   var_dic_min[kw])
                index.append(np.where(ps[comp] > kwargs[kw])[0])
            elif kw in var_dic_max.keys():
                comp = _convert_component_to_index(self._params,
                                                   var_dic_max[kw])
                index.append(np.where(ps[comp] < kwargs[kw])[0])

        tot_index = index[0]
        for i in range(len(index)):
            tot_index = set(tot_index).intersection(index[i])
        tot_index = list(tot_index)

        index_comp = _convert_component_to_index(self._params, 'index')

        return np.array(ps[index_comp][tot_index])

    def scatterplot(self, time=None, component1='x', component2='y', species=1,
                    index='all', comoving=False, s=1, **kwargs):
        """
        Method that produces a scatter plot of the given tracked phase space.

        Parameters
        --------
        time : float, optional
            Time of the given plot, default is None
        component1: str, optional
            Phase space component that goes on the x axis.
            Default value is 'x'.
        component2: str, optional
            Phase space component that goes on the y axis.
            Default value is 'y'.
        species: int, optional
            Tracked species number, default is 1
        index: str, optional
            Particle index to be plotted, default is 'all'.
        comoving: bool, optional
            Flag that shifts the longitudinal axis in the comoving
            reference frame. Default is False.
        s: int, optional
            Scatterplot element size. Default is 1.

        Kwargs
        ---------
        Accepts all the kwargs for matplotlib.scatterplot
        """
        if time is None and self._Simulation._verbose_warning:
            print("""
        Time not known, plotting last available dataset.
                """)
            time = self._available_times[-1]
        else:
            time = self._Simulation._nearest_tracking_time(time, species)

        self._return_tracking_phase_space(time, species)
        ps0 = self._get_stored_tracking(time, species)
        ps = ps0.copy()

        if comoving:
            v = self._params['w_speed']
            t0 = self._params['wi_time']
            xcomp = _convert_component_to_index(self._params, 'x')
            ps[xcomp] = ps[xcomp]-v*(time - t0)
        possible_components = ['x', 'y', 'z', 'px', 'py', 'pz', 'weight',
                               'gamma', 'index', 'a']

        mask = self._check_index_in_ps(ps, index)
        if component1 in possible_components:
            component1 = _convert_component_to_index(self._params, component1)
            psx = ps[component1][mask]
        if component2 in possible_components:
            component2 = _convert_component_to_index(self._params, component2)
            psy = ps[component2][mask]

        if component1 == _can_dict['y'] or component2 == _can_dict['y']:
            cmppy = _convert_component_to_index(self._params, 'py')
            pspy = ps[cmppy][mask]
            cmpa = _convert_component_to_index(self._params, 'a')
            psa = ps[cmpa][mask]
            if component1 == _can_dict['y']:
                psx = pspy - psa
            if component2 == _can_dict['y']:
                psy = pspy - psa

        plt.scatter(psx, psy, s=s, **kwargs)

    def trajectory_plot(self, time=None, component1='x', component2='y',
                        species=1, index='all', **kwargs):
        """
        Method that plots the trajectory of a given particle over time.

        Parameters
        ---------
        time : float, optional
            Final time of the given plot, default is None
        component1: str, optional
            Phase space component that goes on the x axis.
            component1 can also be 'time', which generates a
            component2 plot vs time.
            Default value is 'x'.
        component2: str, optional
            Phase space component that goes on the y axis.
            Default value is 'y'.
        species: int, optional
            Tracked species number, default is 1
        index: str, optional
            Particle index to be plotted, default is 'all'.
        comoving: bool, optional
            Flag that shifts the longitudinal axis in the comoving
            reference frame. Default is False.
        s: int, optional
            Scatterplot element size. Default is 1.

        Kwargs
        ---------
        Accepts all the kwargs for matplotlib.plot
        """

        time = self._Simulation._nearest_tracking_time(time, species)
        time_index = self._available_times.index(time)
        time_array = self._available_times[:time_index + 1]
        plot_list_abscissa = dict()
        plot_list_ordinate = dict()

        if component1 != 'time' and component1 != _can_dict['y']:
            component1 = _convert_component_to_index(self._params, component1)
        if component2 != 'time' and component2 != _can_dict['y']:
            component2 = _convert_component_to_index(self._params, component2)
        ind_comp = _convert_component_to_index(self._params, 'index')

        # cmap_name = 'copper'
        # norm = plt.Normalize(min_time, max_time)
        # Performing the first iteration explicitly to generate the numpy
        # arrays to be plotted
        instant = time_array.pop(0)
        self._return_tracking_phase_space(instant, species)
        ps0 = self._get_stored_tracking(instant, species)
        ps = ps0.copy()
        mask = self._check_index_in_ps(ps, index)
        if component1 == _can_dict['y'] or component2 == _can_dict['y']:
            pycomp = _convert_component_to_index(self._params, 'py')
            acomp = _convert_component_to_index(self._params, 'a')
        if component1 != 'time' and component1 != _can_dict['y']:
            temp_comp1 = ps[component1][mask]
        if component2 != _can_dict['y']:
            temp_comp2 = ps[component2][mask]
        if component1 == _can_dict['y']:
            temp_comp1 = ps[pycomp][mask] - ps[acomp][mask]
        if component2 == _can_dict['y']:
            temp_comp2 = ps[pycomp][mask] - ps[acomp][mask]
        temp_index = ps[ind_comp][mask]
        particle_number = np.count_nonzero(mask)
        if component1 == 'time':
            for i in range(particle_number):
                plot_list_abscissa[temp_index[i]] = np.array([instant])
                plot_list_ordinate[temp_index[i]] = np.array(temp_comp2[i])
        else:
            for i in range(particle_number):
                plot_list_abscissa[temp_index[i]] = np.array(temp_comp1[i])
                plot_list_ordinate[temp_index[i]] = np.array(temp_comp2[i])

        for instant in time_array:
            self._return_tracking_phase_space(instant, species)
            ps0 = self._get_stored_tracking(instant, species)
            ps = ps0.copy()
            mask = self._check_index_in_ps(ps, index)
            if component1 != 'time' and component1 != _can_dict['y']:
                temp_comp1 = ps[component1][mask]
            if component2 != _can_dict['y']:
                temp_comp2 = ps[component2][mask]
            if component1 == _can_dict['y']:
                temp_comp1 = ps[pycomp][mask] - ps[acomp][mask]
            if component2 == _can_dict['y']:
                temp_comp2 = ps[pycomp][mask] - ps[acomp][mask]
            temp_index = ps[ind_comp][mask]
            particle_number = np.count_nonzero(mask)
            for i in range(particle_number):
                if component1 == 'time':
                    plot_list_abscissa[temp_index[i]] = \
                        np.append(plot_list_abscissa[temp_index[i]], [instant])
                else:
                    plot_list_abscissa[temp_index[i]] = \
                        np.append(plot_list_abscissa[temp_index[i]],
                                  temp_comp1[i])
                plot_list_ordinate[temp_index[i]] = \
                    np.append(plot_list_ordinate[temp_index[i]], temp_comp2[i])

        for key in plot_list_abscissa.keys():
            plt.plot(plot_list_abscissa[key],
                     plot_list_ordinate[key], **kwargs)
