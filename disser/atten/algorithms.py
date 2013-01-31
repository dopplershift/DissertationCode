import numpy as np
import scipy.integrate as si
import scipy.optimize as so
import quantities as pq

from .. import datatypes
from ..plugintools import PluginRegistry
from ..units import dB, dBz
from ..cbook import is_string_like
from fit_coeffs import za_coeffs, ka_coeffs, sc_coeffs

class AttenuationRegistry(PluginRegistry):
    def __init__(self):
        super(AttenuationRegistry, self).__init__()
        self.nameLookup = dict()
    def run(self, algs, data, var):
        if is_string_like(algs):
            self.nameLookup[algs](data, var)
        else:
            for alg in algs:
                self.nameLookup[alg](data, var)
    def runAll(self, data, var):
        for alg in self:
            alg(data, var)
    def register(self, name, dt, kinds, coeffs, **kwargs):
        def wrapper(func):
            alg = AttenuationAlgorithm(name, func, dt, kinds, coeffs, **kwargs)
            PluginRegistry.register(self, alg)
            self.nameLookup[name] = alg
            return func
        return wrapper
attenAlgs = AttenuationRegistry()

class AttenuationAlgorithm(object):
    typeInfo = {'H':(datatypes.Attenuation, 'H'),
                'V':(datatypes.Attenuation, 'V'),
                'diff':(datatypes.DiffAtten, None)}

    def __init__(self, name, alg, dt, kinds, coeffs, **kwargs):
        self.name = name
        self.alg = alg
        self.dt = dt
        self.kinds = kinds
        self.coeffs = coeffs
        self.data_args = kwargs

    def __call__(self, data, var='H'):
        # TODO: Need to handle automatically calculation of diffAtten from
        # H and V, if necessary
        kwargs = {k:f(data) for k,f in self.data_args.items()}

        if var == 'diff' and var not in self.kinds:
            atten_h = self.__call__(data, var='H')
            atten_v = self.__call__(data, var='V')
            atten = atten_h - atten_v
        else:
            coeffs = self.coeffs[data.waveBand, var]
            try:
                coeffs = list(coeffs)
            except TypeError:
                coeffs = [coeffs]
            args = [data.fields.grabData(f, pol=var) for f in self.dt]
            atten = self.alg(*(args + coeffs), **kwargs)
        dt, pol = self.typeInfo[var]
        data.addField(atten, dt, pol=pol, source=self.name)
        return atten

zphi_coeffs = {k:(v[1], ka_coeffs[k]) for k,v in za_coeffs.items()}

@attenAlgs.register('SC',
        [datatypes.SNR, datatypes.Reflectivity, datatypes.PhiDP],
        ('H', 'V'), zphi_coeffs, dr=lambda d: d.gate_length,
        phi0=lambda d: d.phi_offset)
def self_consistent(snr, z, phi, b, gamma, dr, phi0):
    #b = pq.Quantity(b, units=1./dBz)
    #gamma = pq.Quantity(gamma, units=dB / pq.degrees)
    z = z.rescale(dBz).magnitude
    snr = snr.rescale(dB).magnitude
    phi = phi.rescale(pq.degree).magnitude
    dr = dr.rescale(pq.kilometer).magnitude
    phi0 = phi0.rescale(pq.degree).magnitude
    atten = np.zeros_like(z)
    for ray in range(atten.shape[0]):
        #good_snr = np.argwhere(snr[ray] > 1.0)
        #if not good_snr.size > 0:
            #continue
        #begin = good_snr[0]
        #end = good_snr[-1]
        #if end < phi.shape[-1] - 1:
            #end += 1
        mask = (~np.isnan(phi[ray])) & (snr[ray] > 0.)
        if not np.any(mask):
            continue

        # Use a 1.5 km half-width for smoothing. Use the gatewidth to
        # convert this to a number of bins.
        smooth_phi = running_mean(phi[ray, mask], int(1.5 / dr))
        delta_phi = smooth_phi[-1] - phi0
        try:
            atten[ray, mask] = tuned_zphi(z[ray, mask],
                    running_mean(phi[ray, mask], int(1.0 / dr)), dr,
                    delta_phi, b, gamma)
            atten[ray] = (2 * dr * atten[ray]).cumsum()
        except RuntimeError:
            print "Minimizer failed on ray: %d" % ray
            atten[ray] = np.nan
    atten[snr < 0.0] = np.nan
    atten[np.isnan(snr)] = np.nan
    return atten * dB

def running_mean(data, halfWidth):
    numPoints = 2 * halfWidth + 1
    # Pad the boundaries with nans so we don't need special handling
    tempData = np.concatenate([[np.nan] * halfWidth, data,
            [np.nan] * halfWidth])

    # Play with strides so that we can easily get the moving window
    tempData = np.lib.stride_tricks.as_strided(tempData,
            (data.size, numPoints), (tempData.itemsize, tempData.itemsize))

    # Mean accounting for nans is the nansum() divided by number of
    # non-nan items
    return np.nansum(tempData, axis=1) / np.sum(~np.isnan(tempData), axis = 1)

def zphi_error(gamma, z, phi, dr, delta_phi, b):
    A = zphi_atten(z, dr, delta_phi, b, gamma)
    phi_calc = (2. / gamma) * si.cumtrapz(A, dx=dr, initial=0)
    return np.abs(phi - phi_calc).mean()

def tuned_zphi(z, phi, dr, delta_phi, b=0.7644, gamma_default=0.1):
    gamma_min = 0.01 * gamma_default
    gamma_max = 10 * gamma_default
    ret = so.minimize_scalar(zphi_error, bounds=(gamma_min, gamma_max),
                             args=(z, phi, dr, delta_phi, b), method='Bounded')
    if ret.success:
        gamma = ret.x
        #print ret.message
    else:
        raise RuntimeError(ret.message)
        gamma = gamma_default
    #print 'Running zphi with gamma: %f' % gamma
    return zphi_atten(z, dr, delta_phi, b, gamma)


@attenAlgs.register('ZPHI',
        [datatypes.SNR, datatypes.Reflectivity, datatypes.PhiDP],
        ('H', 'V'), zphi_coeffs, dr=lambda d: d.gate_length,
        phi0=lambda d: d.phi_offset)
def zphi(snr, z, phi, b, gamma, dr, phi0):
    #b = pq.Quantity(b, units=1./dBz)
    #gamma = pq.Quantity(gamma, units=dB / pq.degrees)
    z = z.rescale(dBz).magnitude
    snr = snr.rescale(dB).magnitude
    phi = phi.rescale(pq.degree).magnitude
    dr = dr.rescale(pq.kilometer).magnitude
    phi0 = phi0.rescale(pq.degrees).magnitude
    atten = np.zeros_like(z)
    for ray in range(atten.shape[0]):
        #good_snr = np.argwhere(snr[ray] > 1.0)
        #if not good_snr.size > 0:
            #continue
        #begin = good_snr[0]
        #end = good_snr[-1]
        #if end < phi.shape[-1] - 1:
            #end += 1
        mask = (~np.isnan(phi[ray])) & (snr[ray] > 0.0)
        if not np.any(mask):
            continue

        smooth_phi = running_mean(phi[ray, mask], int(1.5 / dr))
        delta_phi = smooth_phi[-1] - phi0
        atten[ray, mask] = zphi_atten(z[ray, mask], dr, delta_phi,
                b, gamma)
        atten[ray] = (2 * dr * atten[ray]).cumsum()
    atten[snr < 0.0] = np.nan
    atten[np.isnan(snr)] = np.nan
    return atten * dB

# This is the main calculation, as outlined in equation (24) in the Testud
# et al. (2000) paper
def zphi_atten(z, dr, delta_phi, b=0.7644, gamma=3.2e-1):
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


# These come from the Bring et al. (1990) paper that originally introduced
# the linear relation between attenuation and phi.
bringi_coeffs = {('S', 'H') : 0.016, ('S', 'diff') : 0.00367,
        ('C', 'H') : 0.054, ('C', 'diff') : 0.0157,
        ('X', 'H') : 0.25, ('X', 'diff') : 0.05}

@attenAlgs.register('Linear-Bringi', [datatypes.PhiDP], ('H', 'diff'), bringi_coeffs)
@attenAlgs.register('Linear', [datatypes.PhiDP], ('H', 'diff'), ka_coeffs)
def linear(phi, coeff=0.08):
    try:
        coeff.magnitude
    except AttributeError:
        coeff = coeff * dB / pq.degrees

    return coeff * phi
