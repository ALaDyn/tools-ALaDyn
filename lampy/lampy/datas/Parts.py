from ..compiled_cython.read_phase_space import total_phase_space_read

_phase_spaces = ['phase_space', 'phase_space_ionization',
                 'phase_space_high_energy']


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
        elif ndim == 2:
            phase_space['x'] = ps[0]
            phase_space['y'] = ps[1]
            phase_space['px'] = ps[2]
            phase_space['py'] = ps[3]
            phase_space['weight'] = ps[4]

        self._stored_phase_space[(phase_space_name, timestep)] = phase_space

#    def scatter(self, phase_space, timestep, component1='x',
#                component2='y', **kwargs):