import scattering
import scipy.constants as consts
import numpy as np
import matplotlib.pyplot as plt

def liebe(mat, lam, temp=20.):
    f = consts.c / (lam * consts.giga)
    th = 1 - 300. / consts.C2K(temp)
    eps0 = 77.66 - 103.3 * th
    eps1 = 0.0671 * eps0
    gamma1 = 20.20 + 146.4 * th + 316 * th * th
    eps2 = 3.52 + 7.52 * th
    gamma2 = 39.8 * gamma1
    epsM = (eps0 - eps1) / (1 - 1.0j * (f / gamma1)) + (eps1 - eps2) / (
            1 - 1.0j * (f / gamma2)) + eps2
    return np.sqrt(epsM)

temps = np.linspace(-10, 40)
fig, (ax1, ax2) = plt.subplots(2, 1)

for lam,color in zip([0.03, 0.05, 0.1], ['r', 'g', 'b']):
    eps_ray = scattering.refractive_index('water', lam, temps)
    eps_liebe = liebe('water', lam, temps)

    ax1.plot(temps, eps_ray.real, color + '-', label='%.2f' % lam)
    ax1.plot(temps, eps_liebe.real, color + '--', label='%.2f' % lam)
    ax1.set_title('Real part of index of refraction')
    ax1.set_xlabel(u'Temperature (\N{DEGREE SIGN}C)')
    ax1.grid(True)

    ax2.plot(temps, eps_ray.imag, color + '-', label='%.2f' % lam)
    ax2.plot(temps, eps_liebe.imag, color + '--', label='%.2f' % lam)
    ax2.set_title('Imaginary part of index of refraction')
    ax2.set_xlabel(u'Temperature (\N{DEGREE SIGN}C)')
    ax2.grid(True)

plt.tight_layout()
plt.show()
