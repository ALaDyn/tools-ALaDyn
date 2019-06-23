import os
import numpy as np


class Diagnostics(object):

    def __init__(self, Simulation):

        self._Simulation = Simulation
        self._params = Simulation.params
        self._directories = Simulation.directories
        self._timesteps = Simulation._timesteps
        self._dimensions = Simulation._dimensions
        self._path = Simulation.path
        self._Directories = Simulation._Directories

    def _derive_diag_path(self, diag_type, time):

        time = self._Simulation._nearest_time(time)

        for key, value in self._timesteps.items():
            if time == value:
                folder = key
        diag_name = diag_type+folder[-2:]
        path = os.path.join(self._path, 'diagnostics', diag_name)
        file_path = path+'.dat'
        return file_path

    def read_ionz_diagnostic(self, time):

        path = self._derive_diag_path('ionz_emitt', time)

        file_read = open(path, 'r')

        lines = file_read.readlines()

        number_of_outputs = int(lines[17].split()[1])
        sample = np.zeros(number_of_outputs)
        time = sample
        charge = sample
        xmean = sample
        ymean = sample
        zmean = sample
        pxmean = sample
        pymean = sample
        pzmean = sample
        xrms = sample
        yrms = sample
        zrms = sample
        pxrms = sample
        pyrms = sample
        pzrms = sample
        emitx = sample
        emity = sample
        dgamma = sample
        gamma = sample
        offset = lines.index('time')+1
        offset += (number_of_outputs-1)/5+1+1+2*((number_of_outputs-1)/6+1)+1+1
        for i in range(number_of_outputs):
            time[i] = float(lines[20+i/5].split()[i % 5])
            charge[i] = float(lines[20+(number_of_outputs-1)/5+1+1 +
                              (number_of_outputs-1)/6+1+1+i/6].split()[i % 6])

        ind = 0
        for line in lines[offset:offset+number_of_outputs]:
            xmean[ind] = float(line.split()[0])
            ymean[ind] = float(line.split()[1])
            zmean[ind] = float(line.split()[2])
            pxmean[ind] = float(line.split()[3])
            pymean[ind] = float(line.split()[4])
            pzmean[ind] = float(line.split()[5])
            ind += 1

        offset += number_of_outputs+1
        ind = 0
        for line in lines[offset:offset+number_of_outputs]:
            xrms[ind] = float(line.split()[0])
            yrms[ind] = float(line.split()[1])
            zrms[ind] = float(line.split()[2])
            pxrms[ind] = float(line.split()[3])
            pyrms[ind] = float(line.split()[4])
            pzrms[ind] = float(line.split()[5])
            ind += 1
        offset += number_of_outputs+1
        ind = 0
        for line in lines[offset:offset+number_of_outputs]:
            emitx[ind] = np.sqrt(float(line.split()[0]))
            emity[ind] = np.sqrt(float(line.split()[1]))
            gamma[ind] = float(line.split()[2])
            dgamma[ind] = float(line.split()[3])
            ind += 1

        totdata = dict()
        totdata['time'] = time
        totdata['charge'] = charge
        totdata['zmean'] = zmean
        totdata['xmean'] = xmean
        totdata['ymean'] = ymean
        totdata['pxmean'] = pxmean
        totdata['pymean'] = pymean
        totdata['pzmean'] = pzmean
        totdata['xrms'] = xrms
        totdata['yrms'] = yrms
        totdata['zrms'] = zrms
        totdata['pxrms'] = pxrms
        totdata['pyrms'] = pyrms
        totdata['pzrms'] = pzrms
        totdata['emittance_y'] = emitx
        totdata['emittance_z'] = emity
        totdata['gamma'] = gamma
        totdata['e_spread'] = dgamma

        return totdata
