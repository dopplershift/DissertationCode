import matplotlib.pyplot as plt
import numpy as np
import scattering
import dsd
import scipy.constants as consts
from disser.units import to_dB

d = (np.linspace(0.01, 0.8, 100).reshape(-1, 1, 1).astype(np.float32) 
        * consts.centi)
d0_lut = (np.linspace(0.01, .5, 40).reshape(1, -1, 1).astype(np.float32)
        * consts.centi)
nr_lut = np.logspace(-2, 6, 110).reshape(1, 1,-1).astype(np.float32)
dist = dsd.gamma(d, d0_lut, nr_lut, nu=-0.8)

wavelength = 0.053
temps = np.arange(0, 35.0, 10.0)

db_factor = 10.0 * np.log10(np.e)
ref_adjust = 180
d_plot = d.squeeze() / consts.centi

z = np.empty((temps.size,) + dist.shape[1:], dtype=np.float32)
zdr = np.empty_like(z)
atten = np.empty_like(z)
diff_atten = np.empty_like(z)
kdp = np.empty_like(z)
delta = np.empty_like(z)

fig, ax = plt.subplots(2, 3, sharex=True)

qr = dsd.lwc(d, dist)
fig2, ax2 = plt.subplots(2, 3, sharex=True)

fig3, ax3 = plt.subplots(1, 2, sharex=True)

for ind, t in enumerate(temps):
    scatt = scattering.scatterer(wavelength, t, 'water', diameters=d,
            shape='oblate')
    scatt.set_scattering_model('tmatrix')
    plot_label = u'%d\u00B0C' % t
    ax[0,0].plot(d_plot, scatt.sigma_bh.squeeze(), label=plot_label)
    ax[0,1].plot(d_plot, scatt.sigma_eh.squeeze(), label=plot_label)
    ax[0,2].plot(d_plot, -np.angle(-scatt.S_bkwd[0,0].conj() *
        scatt.S_bkwd[1,1]).squeeze(), label=plot_label)
    ax[1,0].plot(d_plot, to_dB(scatt.sigma_bh / scatt.sigma_bv).squeeze(),
            label=plot_label)
    ax[1,1].plot(d_plot, to_dB(scatt.sigma_eh / scatt.sigma_ev).squeeze(),
            label=plot_label)
    diff_phase = (scatt.S_frwd[0, 0] - scatt.S_frwd[1, 1]).real.squeeze()
    ax[1,2].plot(d_plot, diff_phase, label=plot_label)

    z[ind] = to_dB(scatt.get_reflectivity_factor(dist, polar='h'))
    zdr[ind] = z[ind] - scatt.get_reflectivity_factor(dist, polar='v')
    atten[ind] = scatt.get_attenuation(dist, polar='h') * db_factor
    diff_atten[ind] = (atten[ind] -
            scatt.get_attenuation(dist, polar='v') * db_factor)
    kdp[ind] = (scatt.get_propagation_wavenumber(dist, polar='h') -
            scatt.get_propagation_wavenumber(dist, polar='v'))
    delta[ind] = scatt.get_backscatter_differential_phase(dist)
    temp_color = plt.get_cmap('RdBu')(plt.Normalize(0,30)(t))
    scatter_kwargs = dict(c=temp_color, edgecolor='none', alpha=0.4)
    ax2[0,0].scatter(qr, z[ind] + ref_adjust, **scatter_kwargs)
    ax2[0,1].scatter(qr, atten[ind], **scatter_kwargs)
    ax2[0,2].scatter(qr, delta[ind], **scatter_kwargs)
    ax2[1,0].scatter(qr, zdr[ind], **scatter_kwargs)
    ax2[1,1].scatter(qr, diff_atten[ind], **scatter_kwargs)
    ax2[1,2].scatter(qr, kdp[ind], **scatter_kwargs)
#    ref = to_dB(scatt.get_reflectivity_factor(dist)) + ref_adjust
#    atten = scatt.get_attenuation(dist) * consts.kilo * db_factor
    ax3[0].scatter(kdp[ind], atten[ind], **scatter_kwargs)
    ax3[1].scatter(kdp[ind], diff_atten[ind], **scatter_kwargs)


ax[0,0].set_xlabel('Diameter')
ax[0,0].set_ylabel('Backscatter Cross-Section')
ax[0,0].legend(loc='upper left')
ax[0,0].grid()
ax[0,0].set_xlim(0, d_plot.max())

ax[0,1].set_xlabel('Diameter')
ax[0,1].set_ylabel('Horizontal Extinction Cross-Section')
ax[0,1].grid()

ax[0,2].set_xlabel('Diameter')
ax[0,2].set_ylabel('Differential Backscatter Phase')
ax[0,2].grid()

ax[1,0].set_xlabel('Diameter')
ax[1,0].set_ylabel('Differential Backscatter Cross-Section')
ax[1,0].grid()

ax[1,1].set_xlabel('Diameter')
ax[1,1].set_ylabel('Differential Extinction Cross-Section')
ax[1,1].grid()

ax[1,2].set_xlabel('Diameter')
ax[1,2].set_ylabel('Differential Wavenumber')
ax[1,2].grid()

ax2[0,0].set_xlabel('Rain Water')
ax2[0,0].set_ylabel('Reflectivity Factor')
ax2[0,0].grid()
ax2[0,0].set_xlim(0, 35)

ax2[0,1].set_xlabel('Rain Water')
ax2[0,1].set_ylabel('Attenuation')
ax2[0,1].grid()

ax2[0,2].set_xlabel('Rain Water')
ax2[0,2].set_ylabel('Differential Backscatter Phase')
ax2[0,2].grid()

ax2[1,0].set_xlabel('Rain Water')
ax2[1,0].set_ylabel('Differential Reflectivity')
#ax2[1,0].set_ylim(-6,6)
ax2[1,0].grid()

ax2[1,1].set_xlabel('Rain Water')
ax2[1,1].set_ylabel('Differential Attenuation')
ax2[1,1].grid()

ax2[1,2].set_xlabel('Rain Water')
ax2[1,2].set_ylabel('KDP')
ax2[1,2].grid()

ax3[0].set_xlabel('KDP')
ax3[0].set_ylabel('Attenuation')
ax3[0].grid()
ax3[0].set_xlim(0,3)

ax3[1].set_xlabel('KDP')
ax3[1].set_ylabel('Differential Attenuation')
ax3[1].grid()

plt.show()
