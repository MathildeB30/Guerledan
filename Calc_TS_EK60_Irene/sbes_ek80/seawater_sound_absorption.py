"""

.. module:: water_sound_absorption
.. moduleauthor:: MLG

Computation of water sound absorption

"""

from numpy import *

def Francois_Garrison( ff, pH, S, T, z ):
    """
    Water sound absorption in dB/m

    :param f: frequency (Hz)
    :param pH: sea water pH
    :param S: salinity (psu)
    :param T: temperature (°C)
    :param z: depth (m)
    :return: water sound absorption in dB/m

    """

    f = ff / 1000. # Original formula is in kHz
    c = 1412. + 3.21 * T + 1.19 * S + 0.0167 * z

    # Borique acid
    A1 = 8.86 / c* (10.**(0.78 * pH - 5))
    P1 = 1
    f1 = 2.8 * sqrt(S / 35.) * (10.**(4. - 1245. / (T + 273.)))

    # Magnesium sulfate
    A2 = (21.44 * S / c) * ( 1 + 0.025 * T)
    P2 = 1 - (1.37e-4) * z + (6.2e-9) * (z**2)
    f2 = (8.17 * 10.**(8. - 1990./(T + 273.)))/(1. + 0.0018 * (S - 35.))

    # Water viscosity
    P3=1-3.83e-5 * z + 4.9e-10 * (z**2)
    A3= (T < 20) * \
        (4.937e-4 - 2.59e-5 * T + 9.11e-7 * (T**2) - 1.5e-8 * (T**3)) +\
        (T >= 20) * \
        (3.964e-4 - 1.146e-5 * T + 1.45e-7 * (T**2) - 6.5e-10 * (T**3))

    alpha = A1 * P1 * f1 * f * f / (f * f + f1 * f1) + A2 * P2 * f2 * f * f/( f * f + f2 * f2) + A3 * P3 * f * f
    return alpha / 1000.
