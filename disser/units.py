import numpy as np
import quantities as pq

def log10(x, out=None):
    """
    Raises a ValueError if input cannot be rescaled to radians.

    Returns a dimensionless quantity.
    """
    if not isinstance(x, pq.Quantity):
        return np.log10(x, out)

    return pq.Quantity(np.log10(x.rescale(pq.dimensionless).magnitude, out),
                    copy=False)

c = pq.constants.natural_unit_of_velocity.simplified

zUnit = pq.CompoundUnit('mm**6/m**3')

dB = pq.UnitQuantity('decibel', pq.dimensionless, 'dB')
dBz = pq.UnitQuantity('dB_relative_to_1_mm^6_m^-3',
    pq.dimensionless, 'dBz', aliases=['dBZ'])
dBm = pq.UnitQuantity('dB_relative_to_1_mW', pq.dimensionless, 'dBm')
dBW = pq.UnitQuantity('dB_relative_to_1_watt', pq.dimensionless, 'dBW')

def to_linear(dB):
    '''Convert a linear dimensionless value to decibels.'''
    return 10. ** (dB.rescale(pq.dimensionless) / 10.)

def to_dB(lin):
    '''Convert a value in (dimensionless) decibels to linear units.'''
    return 10. * log10(lin)

def make_dB(log):
    '''Take a log difference and attach dB units.'''
    return log.magnitude * dB

def to_dBz(z):
    '''Convert reflectivity factor in linear units to dBz.'''
    return to_dB(z.rescale(zUnit).magnitude) * dBz

def to_dBm(w):
    '''Convert power in linear units to dBm.'''
    return to_dB(w.rescale(pq.watt).magnitude) * dBm

def dBW_to_dBm(dbw):
    return (dbw.magnitude + 30.) * dBm
