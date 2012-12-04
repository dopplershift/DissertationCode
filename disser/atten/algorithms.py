import numpy as np
import scipy.integrate as si
import quantities as pq

from .. import datatypes
from ..plugintools import PluginRegistry
from ..units import dB
from fit_coeffs import za_coeffs, ka_coeffs, sc_coeffs

class AttenuationRegistry(PluginRegistry):
    def runAll(self, data, var):
        for alg in self:
            alg(data, var)
    def register(self, name, dt, kinds, coeffs):
        def wrapper(func):
            alg = AttenuationAlgorithm(name, func, dt, kinds, coeffs)
            PluginRegistry.register(self, alg)
            return func
        return wrapper
attenAlgs = AttenuationRegistry()

class AttenuationAlgorithm(object):
    typeInfo = {'H':(datatypes.Attenuation, 'H'),
                'V':(datatypes.Attenuation, 'V'),
                'diff':(datatypes.DiffAtten, None)}

    def __init__(self, name, alg, dt, kinds, coeffs):
        self.name = name
        self.alg = alg
        self.dt = dt
        self.kinds = kinds
        self.coeffs = coeffs

    def __call__(self, data, var='H'):
        # TODO: Need to handle automatically calculation of diffAtten from
        # H and V, if necessary
        coeffs = self.coeffs[data.waveBand, var]
        args = [data.fields.grabData(f, pol=var) for f in self.dt]
        atten = self.alg(*(args + [coeffs]))
        dt, pol = self.typeInfo[var]
        data.addField(atten, dt, pol=pol, source=self.name)
        return atten


# This is the main calculation, as outlined in equation (24) in the Testud
# et al. (2000) paper
def zphi(z, dr, delta_phi, b=0.7644, gamma=3.2e-1):
    phi_factor = np.expm1(np.log(10) * b * gamma * delta_phi / 10.)
    return phi_factor * 10**((b / 10.) * z) / (I0(z, dr, b)
            + phi_factor * I(z, dr, b))

# This is the implementation of the I function. We have two versions, since
# one needs to be implemented as a cumulative function in order to maintain
# a dependence on r.

# This is what the 0.46 really is. Using just 0.46 introduces a 1% error
db_conv = (2. / 10.) * np.log(10)

# Z needs to be in dB here. This is why the formula differs from the paper,
# which is linear. Keeping Z in dB let's us only need a single pow (**) call.
def I(z, dr, b):
    # Cumulatively integrate away from the last value. Return these in the
    # reverse order since it's supposed to be in order of increasing r.
    return db_conv * b * si.cumtrapz(10 ** ((b / 10.)* z[::-1]), dx=dr,
        initial=0)[::-1]

# And here too
def I0(z, dr, b):
    return db_conv * b * np.trapz(10 ** ((b / 10.) * z), dx=dr)


bringi_coeffs = {('C', 'H') : 0.08, ('C', 'diff') : 0.02, ('X', 'H') : 0.233,
        ('X', 'diff') : 0.033}

@attenAlgs.register('LinearBringi', [datatypes.PhiDP], ('H', 'diff'), bringi_coeffs)
@attenAlgs.register('Linear', [datatypes.PhiDP], ('H', 'diff'), ka_coeffs)
def linear(phi, coeff=0.08):
    try:
        coeff.magnitude
    except AttributeError:
        coeff = coeff * dB / pq.degrees

    return coeff * phi
