import numpy as np
import scipy.integrate as si

from fit_coeffs import za_coeffs, ka_coeffs, sc_coeffs

# This is the main calculation, as outlined in equation (24) in the Testud et al. (2000) paper
def zphi(z, dr, delta_phi, b=0.7644, gamma=3.2e-1):
    phi_factor = np.expm1(np.log(10) * b * gamma * delta_phi / 10.)
    return phi_factor * 10**((b / 10.) * z) / (I0(z, dr, b) + phi_factor * I(z, dr, b))

# This is the implementation of the I function. We have two versions, since one needs to be
# implemented as a cumulative function in order to maintain a dependence on r.

# This is what the 0.46 really is. Using just 0.46 introduces a 1% error needlessly.
db_conv = (2. / 10.) * np.log(10)

# Z needs to be in dB here. This is why the formula differs from the paper, which is linear.
# Keeping Z in dB let's us only need a single pow (**) call.
def I(z, dr, b):
    # Cumulatively integrate away from the last value. Return these in the reverse order since
    # it's supposed to be in order of increasing r.
    return db_conv * b * si.cumtrapz(10 ** ((b / 10.)* z[::-1]), dx=dr, initial=0)[::-1]

# And here too
def I0(z, dr, b):
    return db_conv * b * np.trapz(10 ** ((b / 10.) * z), dx=dr)


def linear_phi(phi, coeff=None, lam='C', var='H'):
    if coeff is None:
        coeff = ka_coeffs[lam, var]
    return coeff * phi
