import numpy as np


def convert_a_in_e(real_part, imaginary_part, x_axis, time,
                   parameters, **kwargs):

    dx = parameters['dx']

    omega_0 = parameters['omega_0']

    phi = omega_0*(x_axis - time)

    e_field = (real_part.transpose()*np.sin(phi) +
               imaginary_part.transpose()*np.cos(phi)).transpose()

    e_field = np.diff(e_field, axis=0)/dx

    e_field = np.append(e_field, [e_field[-1, ...]], axis=0)

    return e_field


def convert_a_in_e_envelope(real_part, imaginary_part, x_axis, time,
                            parameters, **kwargs):

    from scipy.signal import hilbert

    dx = parameters['dx']

    omega_0 = parameters['omega_0']

    phi = omega_0*(x_axis - time)

    e_field = (real_part.transpose()*np.cos(phi) -
               imaginary_part.transpose()*np.sin(phi)).transpose()

    e_field = np.diff(e_field, axis=0)/dx

    e_field = np.append(e_field, [e_field[-1, ...]], axis=0)

    e_field = np.abs(hilbert(e_field, axis=0))

    return e_field
