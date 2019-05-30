from ..compiled_cython.read_field import read_ALaDyn_bin
from ..utilities.Utility import _grid_convert
import matplotlib.pylab as plt

_axis_names = ['x', 'y', 'z']
_Electromagnetic_fields = ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']


class Field(object):
    """
    Class that contains every physical dimension which
    is defined on a grid.

    The methods inside 'Field' allow to manipulate and plot
    datas related to every electromagnetic field or fluid
    quantity generated.
    """

    def __init__(self, Simulation):

        self._Simulation = Simulation
        self._params = Simulation.params
        self._directories = Simulation.directories
        self._timesteps = Simulation._timesteps
        self._dimensions = Simulation._dimensions
        self._path = Simulation.path
        self._stored_fields = dict()
        self._stored_axis = dict()
        self.normalized = False
        self.comoving = False
        self._E0 = 1
        self._omegap = self._params['omega_p']
        self._Directories = Simulation._Directories

    def _field_read(self, field_name, timestep):

        file_path = self._Simulation._derive_file_path(field_name, timestep)
        f, x, y, z = read_ALaDyn_bin(file_path, self._params)

        self._stored_fields[(field_name, timestep)] = f
        self._stored_axis[('x', timestep)] = x
        if self._params['n_dimensions'] >= 2:
            self._stored_axis[('y', timestep)] = y
        if self._params['n_dimensions'] == 3:
            self._stored_axis[('z', timestep)] = z
        if self._params['n_dimensions'] == 2:
            self._stored_fields[(field_name, timestep)] =\
                self._stored_fields[(field_name, timestep)][:, :, 0]

    def _search_field_by_timestep(self, timestep):

        field_list = list()
        for key in self._stored_fields.keys():
            if timestep in key and key[0] not in field_list:
                field_list += [key[0]]

        return field_list

    def _search_field_by_field(self, field_name):

        time_list = list()
        for key in self._stored_fields.keys():
            if field_name in key and key[1] not in time_list:
                time_list += [key[1]]

        return time_list

    def _return_field(self, field_name, timestep):

        field_list = self._search_field_by_timestep(timestep)
        if field_name in field_list:
            pass

        else:
            self._field_read(field_name, timestep)

    def map_2d(self, field_name, timestep,
               plane='xy', normalized=False, comoving=False, **kwargs):
        """
        Method that generates a 2D map.

        This method produces a 2D map of the given field
        at the given timestep.
        It is based on pylab.imshow, so it accepts all its **kwargs.

        Parameters
        --------
        field_name : str
            Name of the plotted field.
            To know the available field in the simulation,
            check the s.outputs variable.
        timestep : float
            Variable that defines the instant at which
            the field should be plotted.
        plane : str, optional
            Possible values are 'xy', 'xz', 'yz'. Default value is taken
            plane = 'xy'.
            Plane where to execute the cut. Every plane is always defined
            at the center of the orthogonal axis.
        normalized : bool, optional
            Set normalized = True to plot the field normalized to the 'norm'
            quantity. By default, fields are displayed in code units.
        comoving : bool, optional
            Set comoving = True to plot the fields respect to the longitudinal
            variable xi=x-v_w t, where v_w is the moving window velocity.

        kwargs
        --------
        norm : float
            Fields are normalized to the given norm.
            To normalize electromagnetic fields to the cold wavebreaking limit,
            call norm = s.params['omega_p']
        x : float
            Moves the cutting plane along the x axis to a given x value
        y : float
            Moves the cutting plane along the y axis to a given y value
        z : float
            Moves the cutting plane along the z axis to a given z value
        """
        from ..fastread.parameter_read import _read_box_limits

        if 'norm' in kwargs:
            norm = kwargs['norm']
            del kwargs['norm']
        else:
            norm = None

        file_path = self._Simulation._derive_file_path(field_name, timestep)
        file_path = file_path+'.bin'
        box_limits = _read_box_limits(file_path)

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
            plane = _grid_convert(self._timesteps, timestep, self._Directories,
                                  self._path, self._params, x=x_plane,
                                  y=y_plane, z=z_plane)
        self._return_field(field_name, timestep)
        f = self._stored_fields[(field_name, timestep)]

        if (self.normalized or normalized) and\
           field_name in _Electromagnetic_fields:
            if norm is None:
                E0 = self._E0
            else:
                E0 = norm
            f = f/E0

        if self._params['n_dimensions'] == 2:
            if plane != 'xy':
                print("""WARNING: output data is in two dimensions:
                        map will be plotted in the x-y plane""")
            x = self._stored_axis[('x', timestep)]
            if self.comoving or comoving:
                x = x-x[0]
            y = self._stored_axis[('y', timestep)]
            plt.imshow(f.transpose(), origin='low',
                       extent=(x[0], x[-1], y[0], y[-1]), **kwargs)

        elif self._params['n_dimensions'] == 3:

            if plane == 'xy' or plane == 'yx':
                x = self._stored_axis[('x', timestep)]
                if self.comoving or comoving:
                    x = x-x[0]
                y = self._stored_axis[('y', timestep)]
                nz_map = plane[2]
                plt.imshow(f[..., nz_map].transpose(), origin='low',
                           extent=(x[0], x[-1], y[0], y[-1]), **kwargs)

            elif plane == 'xz' or plane == 'zx':
                x = self._stored_axis[('x', timestep)]
                if self.comoving or comoving:
                    x = x-x[0]
                z = self._stored_axis[('z', timestep)]
                ny_map = plane[1]
                plt.imshow(f[:, ny_map, :].transpose(), origin='low',
                           extent=(x[0], x[-1], z[0], z[-1]), **kwargs)

            elif plane == 'zy' or plane == 'yz':
                y = self._stored_axis[('y', timestep)]
                z = self._stored_axis[('z', timestep)]
                nx_map = plane[0]
                plt.imshow(f[nx_map, ...].transpose(), origin='low',
                           extent=(y[0], y[-1], z[0], z[-1]), **kwargs)

    def lineout(self, field_name, timestep, axis='x',
                normalized=False, comoving=False, **kwargs):
        """
        Method that generates a lineout (plot along a single axis).

        This method produces a lineout of the given field
        at the given timestep along one axis.
        It is based on pylab.plot, so it accepts all its **kwargs.

        Parameters
        --------
        field_name : str
            Name of the plotted field.
            To know the available field in the simulation,
            check the s.outputs variable.
        timestep : float
            Variable that defines the instant at which
            the field should be plotted.
        axis : str, optional
            Possible values are 'x', 'y', 'z'. Default value is taken
            axis = 'x'.
            Axis along which the lineout is plotted. The plotting axis
            always intersects the orthogonal plane in its center.
        normalized : bool, optional
            Set normalized = True to plot the field normalized to the 'norm'
            quantity. By default, fields are displayed in code units.
        comoving : bool, optional
            Set comoving = True to plot the fields respect to the longitudinal
            variable xi=x-v_w t, where v_w is the moving window velocity.

        kwargs
        --------
        norm : float
            Fields are normalized to the given norm.
            To normalize electromagnetic fields to the cold wavebreaking limit,
            call norm = s.params['omega_p']
        x : float
            Moves the lineout axis along the x axis to a given x value
        y : float
            Moves the lineout axis along the y axis to a given y value
        z : float
            Moves the lineout axis along the z axis to a given z value
        """

        from ..fastread.parameter_read import _read_box_limits

        if 'norm' in kwargs:
            norm = kwargs['norm']
            del kwargs['norm']
        else:
            norm = None

        file_path = self._Simulation._derive_file_path(field_name, timestep)
        file_path = file_path+'.bin'
        box_limits = _read_box_limits(file_path)

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
            print(error)

        self._return_field(field_name, timestep)
        f = self._stored_fields[(field_name, timestep)]

        if (self.normalized or normalized)\
           and field_name in _Electromagnetic_fields:
            if norm is None:
                E0 = self._E0
            else:
                E0 = norm
            f = f/E0
        line = _grid_convert(self._timesteps, timestep, self._Directories,
                             self._path, self._params, x=x_line, y=y_line)
        if axis == 'z' and self._params['n_dimensions'] == 2:
            print("""WARNING: No lineout along the z axis is possible
                     in 2 dimensions.
                     Lineout will be performed along the x axis""")
            axis = 'x'
        if axis == 'x':
            x = self._stored_axis[('x', timestep)]

            if self.comoving or comoving:
                x = x-x[0]
            if self._params['n_dimensions'] == 2:
                ny_lineout = line[1]
                plt.plot(x, f[:, ny_lineout], **kwargs)
            if self._params['n_dimensions'] == 3:
                ny_lineout = line[1]
                nz_lineout = line[2]
                plt.plot(x, f[:, ny_lineout, nz_lineout], **kwargs)
        elif axis == 'y':
            y = self._stored_axis[('y', timestep)]
            if self._params['n_dimensions'] == 2:
                nx_lineout = line[0]
                plt.plot(y, f[nx_lineout, :], **kwargs)
            if self._params['n_dimensions'] == 3:
                nx_lineout = line[0]
                nz_lineout = line[2]
                plt.plot(y, f[nx_lineout, :, nz_lineout], **kwargs)
        elif axis == 'z':
            z = self._stored_axis[('z', timestep)]
            ny_lineout = line[1]
            nz_lineout = line[2]
            plt.plot(z, f[nx_lineout, ny_lineout, :], **kwargs)

    def normalize(self, **kwargs):
        """
        Method that normalizes the electromagnetic fields.

        When the s.normalize() method is called, all the electromagnetic
        fields will be normalized to a given norm.
        By defalut, they will be normalized to the cold wavebreaking limit.

        kwargs
        --------
        norm : float
            Defines the normalizing value
        """
        self._E0 = self._omegap
        self.normalized = True
        if 'norm' in kwargs:
            self._E0 = kwargs['norm']

    def code_units(self):
        """
        Method that normalized the fields to the code units.
        """
        self.normalized = False

    def comoving(self):
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
