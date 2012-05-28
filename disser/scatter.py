from collections import namedtuple
import numpy as np
import scattering
from disser.units import to_dB, to_dBz, angle, exp_to_dB
import quantities as pq

ScatterResults = namedtuple('ScatterResults',
    'z zdr atten diff_atten kdp delta')
def bulk_scatter(wavelength, temp, dist, diameters):
    csec_units = 'mm**2'

    scatt = scattering.scatterer(wavelength, temp, 'water',
            diameters=diameters, shape='oblate')
    scatt.set_scattering_model('tmatrix')

    z = to_dBz(scatt.get_reflectivity_factor(dist, polar='h'))
    zdr = (z - to_dBz(scatt.get_reflectivity_factor(dist,
        polar='v'))).rescale('dB')
    atten = (scatt.get_attenuation(dist, polar='h').rescale('1/km') * exp_to_dB)
    diff_atten = (atten -
            scatt.get_attenuation(dist, polar='v').rescale('1/km') * exp_to_dB)
    kdp = ((scatt.get_propagation_wavenumber(dist, polar='h') -
            scatt.get_propagation_wavenumber(dist, polar='v')) * pq.rad).rescale(
                    'deg/km')
    delta = (scatt.get_backscatter_differential_phase(dist) * pq.rad).rescale('deg')

    return ScatterResults(z, zdr, atten, diff_atten, kdp, delta)
