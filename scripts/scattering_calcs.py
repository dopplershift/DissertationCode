import matplotlib.pyplot as plt
import numpy as np
import scattering
import dsd
import scipy.constants as consts

d = (np.linspace(0.01, 0.6, 300).reshape(-1, 1, 1).astype(np.float32) 
        * consts.centi)
d0_lut = (np.linspace(0.01, .5, 400).reshape(1, -1, 1).astype(np.float32)
        * consts.centi)
nr_lut = np.logspace(-2, 6, 800).reshape(1, 1,-1).astype(np.float32)
dist = dsd.gamma(d, d0_lut, nr_lut, nu=-0.8)

wavelength = 0.1
temps = np.arange(0, 30.0, 5.0)

db_factor = 10.0 * np.log10(np.e)
ref_adjust = 180

sigma_bh = np.empty((temps.size, d.size), dtype=np.float32)
sigma_eh = np.empty_like(sigma_bh)
sigma_ev = np.empty_like(sigma_bh)
phase = np.empty_like(sigma_bh)

for t in temps:
    scatt = scattering.scatterer(wavelength, t, 'water', diameters=d)
    scatt.set_scattering_model('tmatrix')
    ref = 10.0 * np.log10(scatt.get_reflectivity_factor(dist)) + ref_adjust
    atten = scatt.get_attenuation(dist) * consts.kilo * db_factor

d = d.squeeze() / consts.centi

plt.subplot(2, 2, 1)
plt.semilogy(d, s_ray.sigma_b, 'b--', label = 'Rayleigh (10cm)')
plt.semilogy(d, s_mie.sigma_b, 'b-', label = 'Mie (10cm)')
plt.semilogy(d[::5], s_tmat.sigma_b[::5], 'bx', label = 'T-matrix (10cm)')
plt.semilogy(d, x_ray.sigma_b, 'r--', label = 'Rayleigh (3.21cm)')
plt.semilogy(d, x_mie.sigma_b, 'r-', label = 'Mie (3.21cm)')
plt.semilogy(d[::5], x_tmat.sigma_b[::5], 'rx', label = 'T-matrix (3.21cm)')
plt.legend(loc = 'lower right')
plt.xlabel('Diameter (cm)')
plt.ylabel(r'$\sigma_b \rm{(m^2)}$')

plt.subplot(2, 2, 2)
plt.plot(l, s_ray_ref, 'b--', label = 'Rayleigh (10cm)')
plt.plot(l, s_mie_ref, 'b-', label = 'Mie (10cm)')
plt.plot(l[::5], s_tmat_ref[::5], 'bx', label = 'T-matrix (10cm)')
plt.plot(l, x_ray_ref, 'r--', label = 'Rayleigh (3.21cm)')
plt.plot(l, x_mie_ref, 'r-', label = 'Mie (3.21cm)')
plt.plot(l[::5], x_tmat_ref[::5], 'rx', label = 'T-matrix (3.21cm)')
plt.xlabel('Rain Content (g/m^3)')
plt.ylabel(r'Z$_{e}$ (dBZ)')

plt.subplot(2, 2, 3)
plt.semilogy(d, s_ray.sigma_e, 'b--', label = 'Rayleigh (10cm)')
plt.semilogy(d, s_mie.sigma_e, 'b-', label = 'Mie (10cm)')
plt.semilogy(d[::5], s_tmat.sigma_e[::5], 'bx', label = 'T-matrix (10cm)')
plt.semilogy(d, x_ray.sigma_e, 'r--', label = 'Rayleigh (3.21cm)')
plt.semilogy(d, x_mie.sigma_e, 'r-', label = 'Mie (3.21cm)')
plt.semilogy(d[::5], x_tmat.sigma_e[::5], 'rx', label = 'T-matrix (3.21cm)')
plt.xlabel('Diameter (cm)')
plt.ylabel(r'$\sigma_e \rm{(m^2)}$')

plt.subplot(2, 2, 4)
plt.plot(l, s_ray_atten, 'b--', label = 'Rayleigh (10cm)')
plt.plot(l, s_mie_atten, 'b-', label = 'Mie (10cm)')
plt.plot(l[::5], s_tmat_atten[::5], 'bx', label = 'T-matrix (10cm)')
plt.plot(l, x_ray_atten, 'r--', label = 'Rayleigh (3.21cm)')
plt.plot(l, x_mie_atten, 'r-', label = 'Mie (3.21cm)')
plt.plot(l[::5], x_tmat_atten[::5], 'rx', label = 'T-matrix (3.21cm)')
plt.xlabel('Rain Content (g/m^3)')
plt.ylabel('1-way Attenuation (db/km)')
plt.gcf().text(0.5,0.95,'Comparison of Rayleigh and Mie Scattering models',
  horizontalalignment='center',fontsize=16)
plt.show()
