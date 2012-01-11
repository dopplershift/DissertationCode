import numpy as np
import quantities as pq

c = pq.constants.natural_unit_of_velocity.simplified

zUnit = pq.CompoundUnit('mm**6/m**3')

dB = pq.UnitQuantity('decibel', pq.dimensionless, 'dB')
dBz = pq.UnitQuantity('decibels relative to 1 mm^6 m^-3',
    pq.dimensionless, 'dBz')
dBm = pq.UnitQuantity('decibels relative to 1 mW', pq.dimensionless, 'dBm')
dBW = pq.UnitQuantity('decibels relative to 1 W', pq.dimensionless, 'dBW')

def to_linear(dB):
    '''Convert a linear dimensionless value to decibels.'''
    return 10. ** (dB.rescale(pq.dimensionless) / 10.)

def to_dB(lin):
    '''Convert a value in (dimensionless) decibels to linear units.'''
    return 10. * np.log10(lin.rescale(pq.dimensionless))

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
