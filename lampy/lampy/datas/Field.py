from compiled_cython.read_field import read_ALaDyn_bin
from ..utilities.Utility import _grid_convert, _translate_timestep
from ..fastread.parameter_read import _read_box_limits
import matplotlib.pyplot as plt
import os
import numpy as np

_axis_names = ['x', 'y', 'z']
_Electromagnetic_fields = ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']
_Currents = ['Jx', 'Jy', 'Jz', 'Px_fluid', 'Py_fluid', 'Pz_fluid']
_Envelope = ['A1', 'A2', 'Re1A', 'Im1A', 'Re2A', 'Im2A']
_Densities = ['rho_electrons', 'rho_fluid', 'rho_protons']
_Lorentz_force = ['Fy', 'Fz']
_Energies = ['electrons_energy', 'protons_energy']
_Laser1_from_envelope = ['E1_laser', 'E1_envelope']
_Laser2_from_envelope = ['E2_laser', 'E2_envelope']
# Lorentz force is assumed for a relativistic particle travelling along the x
# direction

for i in range(6):
    _Densities += ['rho_ion'+str(i)]

_Field_list = _Electromagnetic_fields + _Envelope + _Densities + _Energies +\
    _Currents


class Field(object):
    """
    Class that contains every physical dimension which
    is defined on a grid.

    The methods inside 'Field' allow to manipulate and plot
    datas related to every electromagnetic field or fluid
    quantity generated.
    """

    def __init__(self, Simulation, field_list):

        self._Simulation = Simulation
        self._params = Simulation.params
        self._directories = Simulation.directories
        self._timesteps = Simulation._timesteps
        self._dimensions = Simulation._dimensions
        self._path = Simulation.path
        self._Field_list = field_list
        self._stored_fields = dict()
        self._stored_axis = dict()
        self.normalized = False
        self.comoving = False
        self._E0 = 1
        self._A = 1
        self._n0 = 1
        self._omegap = self._params['omega_p']
        self._Directories = Simulation._Directories
        self._box_limits = self._read_all_box_limits()
        self._field_dictionary = self._initialize_field_dic()

    def _field_read(self, field, timestep):

        f = 0
        field_list = self._field_dictionary[field]
        for field_name in field_list:
            file_path = self._Simulation._derive_file_path(field_name,
                                                           timestep)
            try:
                f_temp, x, y, z = read_ALaDyn_bin(file_path, self._params)
            except FileNotFoundError:
                print("""
            Field {} not available, impossible to read.
                """.format(field_name))
                raise

            self._stored_fields[(field_name, timestep)] = f_temp
            if self._params['n_dimensions'] == 2:
                self._stored_fields[(field_name, timestep)] =\
                    self._stored_fields[(field_name, timestep)][:, :, 0]
            self._stored_axis[('x', timestep)] = x
            if self._params['n_dimensions'] >= 2:
                self._stored_axis[('y', timestep)] = y
            if self._params['n_dimensions'] == 3:
                self._stored_axis[('z', timestep)] = z

            if len(field_list) == 1:
                return
            f += f_temp

        self._stored_fields[(field, timestep)] = f
        if self._params['n_dimensions'] == 2:
            self._stored_fields[(field, timestep)] =\
                self._stored_fields[(field, timestep)][:, :, 0]

    def _get_stored_field(self, field_name, timestep):

        if self._Simulation._save_data:
            f = self._stored_fields[(field_name, timestep)]
        else:
            f = self._stored_fields.pop((field_name, timestep))

        return f

    def _initialize_field_dic(self):

        field_dic = dict()
        for field in self._Field_list:
            field_dic[field] = [field]

        if 'Re1A' in self._Simulation.outputs and \
                'Im1A' in self._Simulation.outputs:
            del field_dic['A1']
        if 'Re2A' in self._Simulation.outputs and \
                'Im2A' in self._Simulation.outputs:
            del field_dic['A2']

        return field_dic

    def _read_all_box_limits(self):

        box_limits = dict()
        file_path = None
        for key in self._timesteps.keys():
            folder = key
            for elem in self._Directories._filelist([folder]):
                for em_field in _Electromagnetic_fields:
                    if em_field in elem:
                        file_path = os.path.join(self._path, folder, elem)
                        break
                if file_path is not None:
                    break
            box_limits[folder] = _read_box_limits(file_path)
            file_path = None
        if not box_limits:
            print(""" No EM field output have been found in this
                  folder""")
            return

        return box_limits

    def _return_field(self, field_name, timestep):

        Lorentz_force = False
        Complex_envelope = False
        Las_from_env = False
        field_list = self._search_field_by_timestep(timestep)

        if field_name in _Lorentz_force:
            Lorentz_force = True
            component = field_name[1]
            if component == 'y':
                fields = ['Ey', 'Bz']
            elif component == 'z':
                fields = ['Ez', 'By']

        if self._Simulation._a_from_imaginary and field_name == 'A1':
            fields = ['Re1A', 'Im1A']
            Complex_envelope = True

        if self._Simulation._a_from_imaginary and field_name == 'A2':
            fields = ['Re2A', 'Im2A']
            Complex_envelope = True

        if field_name in _Laser1_from_envelope:
            fields = ['Re1A', 'Im1A']
            Las_from_env = True

        if field_name in _Laser2_from_envelope:
            fields = ['Re2A', 'Im2A']
            Las_from_env = True

        if field_name in field_list:
            pass

        elif field_name in self._field_dictionary.keys():
            self._field_read(field_name, timestep)

        elif Lorentz_force:
            stor_fie = list()
            for field in fields:
                if field in field_list:
                    pass
                else:
                    self._field_read(field, timestep)
                stor_fie += [self._stored_fields[(field, timestep)]]
            temp = np.zeros_like(stor_fie[0])
            if component == 'y':
                temp[:-1, ...] = \
                    0.5*(stor_fie[0][1:, ...] + stor_fie[0][:-1, ...])
                self._stored_fields[(field_name, timestep)] = \
                    temp - stor_fie[1]
            elif component == 'z':
                temp[:, :-1, ...] = \
                    0.5*(stor_fie[0][:, 1:, ...] + stor_fie[0][:, :-1, ...])
                self._stored_fields[(field_name, timestep)] = \
                    temp + stor_fie[1]

        elif Complex_envelope:
            stor_fie = list()
            for field in fields:
                if field in field_list:
                    pass
                else:
                    self._field_read(field, timestep)
                stor_fie += [self._stored_fields[(field, timestep)]]
            self._stored_fields[(field_name, timestep)] = \
                np.sqrt(stor_fie[0]**2+stor_fie[1]**2)

        elif Las_from_env:

            from ..utilities.FieldUtil import convert_a_in_e, \
                convert_a_in_e_envelope

            stor_fie = list()
            for field in fields:
                if field in field_list:
                    pass
                else:
                    self._field_read(field, timestep)
                stor_fie += [self._stored_fields[(field, timestep)]]

            if field_name == 'E1_laser':
                x = self._stored_axis[('x', timestep)]
                self._stored_fields[(field_name, timestep)] = \
                    convert_a_in_e(stor_fie[0], stor_fie[1], x, timestep,
                                   self._params, 1)
            elif field_name == 'E1_envelope':
                x = self._stored_axis[('x', timestep)]
                self._stored_fields[(field_name, timestep)] = \
                    convert_a_in_e_envelope(stor_fie[0], stor_fie[1],
                                            x, timestep,
                                            self._params, 1)

            if field_name == 'E2_laser':
                x = self._stored_axis[('x', timestep)]
                self._stored_fields[(field_name, timestep)] = \
                    convert_a_in_e(stor_fie[0], stor_fie[1], x, timestep,
                                   self._params, 2)
            elif field_name == 'E2_envelope':
                x = self._stored_axis[('x', timestep)]
                self._stored_fields[(field_name, timestep)] = \
                    convert_a_in_e_envelope(stor_fie[0], stor_fie[1],
                                            x, timestep,
                                            self._params, 2)

    def _search_field_by_field(self, field_name):

        time_list = list()
        for key in self._stored_fields.keys():
            if field_name in key and key[1] not in time_list:
                time_list += [key[1]]

        return time_list

    def _search_field_by_timestep(self, timestep):

        field_list = list()
        for key in self._stored_fields.keys():
            if timestep in key and key[0] not in field_list:
                field_list += [key[0]]

        return field_list

    def map_2d(self, field, timestep,
               plane='xy', normalized=False, comoving=False,
               mask=None, mask_argument=None, **kwargs):
        """
        Method that generates a 2D map.

        This method produces a 2D map of the given field
        at the given timestep.
        It is based on pyplot.pcolormesh, so it accepts all its **kwargs.

        Parameters
        --------
        field : str or numpy array type
            If a string is given, it is the name of the plotted field.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        timestep : float
            Variable that defines the instant at which
            the field should be plotted.
        plane : str, optional
            Possible values are 'xy', 'xz', 'yz'. Default value is taken
            plane = 'xy'.
            Plane where to execute the cut. Every plane is always defined
            at the center of the orthogonal axis.
        normalized : bool, optional
            Set normalized = True to plot the field normalized to the
            'unit_field' quantity.
            By default, fields are displayed in code units.
        comoving : bool, optional
            Set comoving = True to plot the fields respect to the longitudinal
            variable xi=x-v_w t, where v_w is the moving window velocity.
        mask : function, optional
            If a mask function is given, only the non-masked values of the
            output field will be plotted. Tipically, a mask function can be
            passed as an anonymus lambda function, e.g.

                mask = lambda x, y: (x-130)**2 + y**2 <= 200

            that will be interpreted using the numpy method
            numpy.ma.masked_where.
        mask_argument : str, optional
            Possible values are None (default), 'field' or 'axes'.
            If None, an explicit mask is used for the field array via

                f = mask.

            If 'field' or 'axes', the variables of the anonymus function will
            be interpreted as the plotted field itself or the axes
            respectively.

        kwargs
        --------
        unit_field : float
            Fields are normalized to the given unit_field.
            To normalize electromagnetic fields to the cold wavebreaking limit,
            call unit_field = s.params['omega_p']
        x : float
            Moves the cutting plane along the x axis to a given x value
        y : float
            Moves the cutting plane along the y axis to a given y value
        z : float
            Moves the cutting plane along the z axis to a given z value
        """

        accepted_types = [str, np.ndarray]
        if type(field) not in accepted_types:
            print("""
        Input field must be either a string with the field name or a numpy
        array """)
            return

        if 'unit_field' in kwargs:
            norm = kwargs['unit_field']
            del kwargs['unit_field']
        else:
            norm = None

        if field not in self._Simulation.outputs:
            print("""
        {} is not available.
        Available output are {}.
        """.format(field, self._Simulation.outputs))
            return

        timestep = self._Simulation._nearest_time(timestep)

        folder = _translate_timestep(timestep, self._timesteps)
        box_limits = self._box_limits[folder]

        if 'x' in kwargs:
            x_plane = kwargs['x']
            del kwargs['x']
        else:
            x_plane = (box_limits['x_max']+box_limits['x_min'])/2
        if 'y' in kwargs:
            y_plane = kwargs['y']
            del kwargs['y']
        else:
            y_plane = 0
        if 'z' in kwargs:
            z_plane = kwargs['z']
            del kwargs['z']
        else:
            z_plane = 0

        error = ''
        if x_plane > box_limits['x_max'] or x_plane < box_limits['x_min']:
            axis_warning = 'x'
            error += 'WARNING: '
            error += axis_warning
            error += """value chosen for the map
                is not in the computational box\n"""
        if y_plane > box_limits['y_max'] or y_plane < box_limits['y_min']:
            axis_warning = 'y'
            error += 'WARNING: '
            error += axis_warning
            error += """value chosen for the map
                is not in the computational box\n"""
        if z_plane > box_limits['z_max'] or z_plane < box_limits['z_min']:
            axis_warning = 'z'
            error += 'WARNING: '
            error += axis_warning
            error += """value chosen for the map
                is not in the computational box\n"""

        if error != '':
            print(error)

        if self._params['n_dimensions'] == 3:
            map_plane = _grid_convert(box_limits, self._params, x=x_plane,
                                      y=y_plane, z=z_plane)

        if type(field) is str:
            self._return_field(field, timestep)
            f = self._get_stored_field(field, timestep)
        elif type(field) is np.ndarray:
            f = field

        if self.normalized or normalized:
            if field in _Electromagnetic_fields or field in _Lorentz_force:
                if norm is None:
                    norm = self._E0

            if field in _Densities:
                if norm is None:
                    norm = self._n0

            if field in _Envelope:
                if norm is None:
                    norm = self._A

            f = f/norm

        if mask is not None and mask_argument == 'field':
            f = np.ma.masked_where(mask(f), f)
        elif mask is not None and mask_argument is None:
            f = mask

        shading_type = 'auto'
        if 'shading' in kwargs:
            shading_type = kwargs['shading']
            del kwargs['shading']

        if self._params['n_dimensions'] == 2:
            if plane != 'xy':
                print("""WARNING: output data is in two dimensions:
                        map will be plotted in the x-y plane""")
            x = self._stored_axis[('x', timestep)]
            if self.comoving or comoving:
                x = x-x[0]
            y = self._stored_axis[('y', timestep)]
            if mask is not None and mask_argument == 'axes':
                X, Y = np.meshgrid(x, y)
                f = np.ma.masked_where(mask(X.transpose(), Y.transpose()), f)
            plt.pcolormesh(x, y, f.transpose(), shading=shading_type, **kwargs)

        elif self._params['n_dimensions'] == 3:

            if plane == 'xy' or plane == 'yx':
                x = self._stored_axis[('x', timestep)]
                if self.comoving or comoving:
                    x = x-x[0]
                y = self._stored_axis[('y', timestep)]
                nz_map = map_plane['z']
                if mask is not None and mask_argument == 'axes':
                    X, Y = np.meshgrid(x, y)
                    f = np.ma.masked_where(mask(X.transpose(),
                                           Y.transpose()), f)
                plt.pcolormesh(x, y, f[..., nz_map].transpose(),
                               shading=shading_type, **kwargs)

            elif plane == 'xz' or plane == 'zx':
                x = self._stored_axis[('x', timestep)]
                if self.comoving or comoving:
                    x = x-x[0]
                z = self._stored_axis[('z', timestep)]
                ny_map = map_plane['y']
                if mask is not None and mask_argument == 'axes':
                    X, Z = np.meshgrid(x, z)
                    f = np.ma.masked_where(mask(X.transpose(),
                                           Z.transpose()), f)
                plt.pcolormesh(x, z, f[:, ny_map, :].transpose(),
                               shading=shading_type, **kwargs)

            elif plane == 'zy' or plane == 'yz':
                y = self._stored_axis[('y', timestep)]
                z = self._stored_axis[('z', timestep)]
                nx_map = map_plane['x']
                if mask is not None and mask_argument == 'axes':
                    Y, Z = np.meshgrid(y, z)
                    f = np.ma.masked_where(mask(Y.transpose(),
                                           Z.transpose()), f)
                plt.pcolormesh(y, z, f[nx_map, ...].transpose(),
                               shading=shading_type, **kwargs)

    def lineout(self, field, timestep, axis='x',
                normalized=False, comoving=False, **kwargs):
        """
        Method that generates a lineout (plot along a single axis).

        This method produces a lineout of the given field
        at the given timestep along one axis.
        It is based on pylab.plot, so it accepts all its **kwargs.

        Parameters
        --------
        field : str or numpy array type
            If a string is given, it is the name of the plotted field.
            To know the available field in the simulation,
            check the s.show_outputs() variable.
        timestep : float
            Variable that defines the instant at which
            the field should be plotted.
        axis : str, optional
            Possible values are 'x', 'y', 'z'. Default value is taken
            axis = 'x'.
            Axis along which the lineout is plotted. The plotting axis
            always intersects the orthogonal plane in its center.
        normalized : bool, optional
            Set normalized = True to plot the field normalized to the
            'unit_field' quantity.
            By default, fields are displayed in code units.
        comoving : bool, optional
            Set comoving = True to plot the fields respect to the longitudinal
            variable xi=x-v_w t, where v_w is the moving window velocity.

        kwargs
        --------
        unit_field : float
            Fields are normalized to the given unit_field.
            To normalize electromagnetic fields to the cold wavebreaking limit,
            call unit_field = s.params['omega_p']
        x : float
            Moves the lineout axis along the x axis to a given x value
        y : float
            Moves the lineout axis along the y axis to a given y value
        z : float
            Moves the lineout axis along the z axis to a given z value
        """

        accepted_types = [str, np.ndarray]
        if type(field) not in accepted_types:
            if self._Simulation._verbose_error:
                print("""Input field must be either a string with the field name
                or a numpy array """)
            return

        if 'unit_field' in kwargs:
            norm = kwargs['unit_field']
            del kwargs['unit_field']
        else:
            norm = None

        if field not in self._Simulation.outputs:
            if self._Simulation._verbose_error:
                print("""
        {} is not available.
        Available output are {}.
            """.format(field, self._Simulation.outputs))
            return

        timestep = self._Simulation._nearest_time(timestep)

        folder = _translate_timestep(timestep, self._timesteps)
        box_limits = self._box_limits[folder]

        if 'x' in kwargs:
            x_line = kwargs['x']
            del kwargs['x']
        else:
            x_line = (box_limits['x_max']+box_limits['x_min'])/2
        if 'y' in kwargs:
            y_line = kwargs['y']
            del kwargs['y']
        else:
            y_line = 0
        if 'z' in kwargs:
            z_line = kwargs['z']
            del kwargs['z']
        else:
            z_line = 0

        error = ''
        if x_line > box_limits['x_max'] or x_line < box_limits['x_min']:
            axis_warning = 'x'
            error += 'WARNING: '
            error += axis_warning
            error += """value chosen for the lineout
                is not in the computational box\n"""
        if y_line > box_limits['y_max'] or y_line < box_limits['y_min']:
            axis_warning = 'y'
            error += 'WARNING: '
            error += axis_warning
            error += """value chosen for the lineout
                is not in the computational box\n"""
        if z_line > box_limits['z_max'] or z_line < box_limits['z_min']:
            axis_warning = 'z'
            error += 'WARNING: '
            error += axis_warning
            error += """value chosen for the lineout
                is not in the computational box\n"""
        if error != '':
            if self._Simulation._verbose_error:
                print(error)

        if type(field) is str:
            self._return_field(field, timestep)
            f = self._get_stored_field(field, timestep)
        elif type(field) is np.ndarray:
            f = field

        if self.normalized or normalized:
            if field in _Electromagnetic_fields or field in _Lorentz_force:
                if norm is None:
                    norm = self._E0

            if field in _Densities:
                if norm is None:
                    norm = self._n0

            if field in _Envelope:
                if norm is None:
                    norm = self._A

            f = f/norm

        line = _grid_convert(box_limits, self._params, x=x_line, y=y_line,
                             z=z_line)
        if axis == 'z' and self._params['n_dimensions'] == 2:
            if self._Simulation._verbose_warning:
                print("""No lineout along the z axis is possible
                     in 2 dimensions.
                     Lineout will be performed along the x axis""")
            axis = 'x'
        if axis == 'x':
            x = self._stored_axis[('x', timestep)]

            if self.comoving or comoving:
                x = x-x[0]
            if self._params['n_dimensions'] == 2:
                ny_lineout = line['y']
                plt.plot(x, f[:, ny_lineout], **kwargs)

            if self._params['n_dimensions'] == 3:
                ny_lineout = line['y']
                nz_lineout = line['z']
                plt.plot(x, f[:, ny_lineout, nz_lineout], **kwargs)
        elif axis == 'y':
            y = self._stored_axis[('y', timestep)]
            if self._params['n_dimensions'] == 2:
                nx_lineout = line['x']
                plt.plot(y, f[nx_lineout, :], **kwargs)
            if self._params['n_dimensions'] == 3:
                nx_lineout = line['x']
                nz_lineout = line['z']
                plt.plot(y, f[nx_lineout, :, nz_lineout], **kwargs)
        elif axis == 'z':
            z = self._stored_axis[('z', timestep)]
            nx_lineout = line['x']
            ny_lineout = line['y']
            plt.plot(z, f[nx_lineout, ny_lineout, :], **kwargs)

    def normalize(self, field_type=None, **kwargs):
        """
        Method that normalizes the electromagnetic fields.

        When the s.normalize() method is called, all the electromagnetic
        fields will be normalized to a given norm.
        By defalut, they will be normalized to the cold wavebreaking limit.

        Parameters
        --------
        field_type : str, optional
            Type of field to be normalized. It is needed to know which is
            the fallback value if nothing is specified. Possibilities are
            'Electromagnetic', 'Envelope', 'Density', 'All'
        kwargs
        --------
        unit_field : float
            Defines the normalizing value
        """
        norm = None
        if 'unit_field' in kwargs:
            norm = kwargs['unit_field']

        if field_type == 'Electromagnetic' or field_type == 'All':
            if norm is not None:
                self._E0 = norm
            else:
                self._E0 = self._omegap

        if field_type == 'Envelope' or field_type == 'All':
            if norm is not None:
                self._A = norm

        if field_type == 'Density' or field_type == 'All':
            if norm is not None:
                self._n0 = norm
            else:
                self._n0 = self._params['n0']

        self.normalized = True

    def code_units(self):
        """
        Method that normalized the fields to the code units.
        """
        self.normalized = False

    def set_comoving(self):
        """
        Method used to change the reference frame
        from the lab fixed to the comoving longitudinal axis.

        The transformation is xi=x-v_w t, where
        v_w is the moving window velocity.
        """
        self.comoving = True

    def lab_frame(self):
        """
        Method used to change the reference frame
        from the comoving to the lab fixed longitudinal axis.

        The transformation is x=xi+v_w t, where
        v_w is the moving window velocity.
        """
        self.comoving = False

    def get_data(self, field_name, time):
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
            field['time'] is the corresponding time
        """
        f = dict()
        time = self._Simulation._nearest_time(time)
        self._return_field(field_name, time)
        f['data'] = self._get_stored_field(field_name, time)
        f['time'] = time

        return f

    def get_axis(self, axis, time):
        """
        Method that retrieves any axis data and returns it as a 1D array.

        Parameters
        --------
        axis : str
            Name of the axis that is to be collected
        time : float
            Instant at which array is retrieved

        Returns
        --------
        axis : dict
            Data are returned as a dictionary.
            axis['data'] is the array containing the data
            axis['time'] is the corresponding time
        """
        time = self._Simulation._nearest_time(time)
        field_list = self._search_field_by_timestep(time)
        ax = dict()
        axis_is_not_empty = bool(self._stored_axis)
        if axis_is_not_empty or \
            (axis, time) not in self._stored_axis.keys() \
                and self._Simulation._verbose_error:
            print("""No axis have been read yet at timestep {}.
            Call first
            >>> f = s.Field.get_data(any_available_field_name,{})
            to load the axis data""".format(time, time))
            return
        ax['data'] = self._stored_axis[(axis, time)]
        ax['time'] = time

        return ax

    def sum(self, field_list, new_name):
        """
        Method that defines a new custom field as the sum of
        other two available fields.
        WARNING: sum is performed on the available grid, so only fields
        sharing the same grid should be summed.

        Parameters
        --------
        field_list: list
            List of strings containing the fields that should be summed
            (e.g. ['rho_electrons', 'rho_fluid'])
        new_name: str
            Name of the newly defined field
        """
        self._field_dictionary[new_name] = field_list
        self._Simulation.outputs += [new_name]
        self._Field_list += [new_name]
