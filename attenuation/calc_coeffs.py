import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import kilo
from scattering import scatterer
import dsd

if __name__ == '__main__':
    T = 15
    diam = np.linspace(0.00001, 0.01, 250).reshape(-1, 1)
    x = scatterer(0.054, T, shape='oblate', diameters=diam)
    x.set_scattering_model('tmatrix')
#    c = scatterer(0.054, T, shape='oblate', diameters=diam)
#    c.set_scattering_model('tmatrix')
#    s = scatterer(0.10, T, shape='oblate', diameters=diam)
#    s.set_scattering_model('tmatrix')

    d0 = np.linspace(0.000001, 0.0025, 100)
    shape_param = np.linspace(-1, 4, 20)

    log_factor = np.log10(np.e)
    N0_size = 100
    final_shape = (d0.size, shape_param.size, N0_size)
    Kdp = np.zeros(final_shape)
    Ah = np.zeros(final_shape)
    for d_ind, D0 in enumerate(d0):
        for n_ind, n in enumerate(shape_param):
            lowN0 = (2.8 * n * log_factor + 4.2) / 2
            upperN0 = (3.57 * n * log_factor + 5.5) * 2
            N0 = np.logspace(lowN0, upperN0, N0_size)
            dsd_weights = N0 * diam**n * np.exp(-(3.67 + n) * diam / D0)
            Kdp[d_ind, n_ind, :] = (x.get_propagation_wavenumber(dsd_weights,
                polar='h') - x.get_propagation_wavenumber(dsd_weights,
                polar='v'))
            Ah[d_ind, n_ind, :] = x.get_attenuation(dsd_weights)
            
    plt.scatter(np.rad2deg(Kdp.flatten()) * kilo,
        10 * log_factor * kilo * Ah.flatten(), c=d0)
    plt.xlabel('$K_{DP} (deg km^{-1})$')
    plt.ylabel('$A_h (dB km^{-1})')
    plt.xlim(xmin=0, xmax=12)
    plt.ylim(ymin=0, ymax=0.6)
    plt.show()
