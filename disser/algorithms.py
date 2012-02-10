from .plugintools import PluginRegistry
from .datatypes import DataType
AlgorithmRegistry = PluginRegistry()
import numpy as np


class Algorithm(object):

    def __init__(self, name, func, products, sources):
        self.name = name
        self._calc_func = func
        self.products = products
        self.sources = sources

    def __call__(self, *args, **kwargs):
        return self._calc_func(*args, **kwargs)


def algorithm(makes, uses, name=None):

    def dec(func):
        alg = Algorithm(name if name is None else func.func_name, func, makes,
            uses)
        for m in makes:
            m.add_source(alg)
        AlgorithmRegistry.register(alg) # Is this necessary?
        return func
    return dec

#def snr(ref, rng, radConst):
#    return snr

#def tempK(tempC):
#    return tempK

WindUComp = DataType(name='U Wind Component', units='m / s', abbr='U')
WindVComp = DataType(name='V Wind Component', units='m / s', abbr='V')
WindDir = DataType(name='Wind Direction', units='deg', abbr='TH')
WindSpd = DataType(name='Wind Speed', units='m / s', abbr='SP')


@algorithm(makes=(WindUComp, WindVComp), uses=(WindSpd, WindDir))
def windcomps(spd, direc):
    return spd * np.sin(direc), spd * np.cos(direc)

WindCompAlg = Algorithm("WindComponents", windcomps,
    makes=(WindUComp, WindVComp), uses=(WindSpd, WindDir))
