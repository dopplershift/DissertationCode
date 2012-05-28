import matplotlib.pyplot as plt
import numpy as np
import scattering
import dsd
from disser.units import to_dB, to_dBz, angle, exp_to_dB
from disser.scatter import bulk_scatter
import quantities as pq

d = (np.linspace(0.01, 8, 100).reshape(-1, 1, 1).astype(np.float32) 
        * pq.mm)
qr = (np.linspace(0, 15.0, 40).reshape(1, -1, 1).astype(np.float32)
        * pq.g / pq.m**3)
n = (np.logspace(-2, 6, 110).reshape(1, 1,-1).astype(np.float32)
        / pq.m**3)
dist = dsd.constrained_gamma_from_moments(n, qr, d)

qr_calc = dsd.lwc(d, dist).rescale('g/m**3')
print np.abs(qr - qr_calc).max()
rr = dsd.rainrate(d, dist, dsd.rain_fallspeed(d)).rescale('mm/hr')
print rr.min(), rr.max()

wavelength = 0.053 * pq.m
temp_colors = {0:'b', 10:'g', 20:'y', 30:'r'}
temps = np.array(temp_colors.keys())

d_plot = d.squeeze()

z = np.empty((temps.size,) + dist.shape[1:], dtype=np.float32)
zdr = np.empty_like(z)
atten = np.empty_like(z)
diff_atten = np.empty_like(z)
kdp = np.empty_like(z)
delta = np.empty_like(z)

max_kdp = 100
max_atten = 40
max_diff_atten = 10

fig2, ax2 = plt.subplots(2, 3, sharex=True)
fig3, ax3 = plt.subplots(1, 2, sharex=True)

for ind, t in enumerate(temps):
    scatter_kwargs = dict(c=temp_colors[t], edgecolor='none', alpha=0.4)
    scatt = bulk_scatter(wavelength, t, dist, d)
    ax2[0,0].scatter(qr, scatt.z, **scatter_kwargs)
    ax2[0,1].scatter(qr, scatt.atten, **scatter_kwargs)
    ax2[0,2].scatter(qr, scatt.delta, **scatter_kwargs)
    ax2[1,0].scatter(qr, scatt.zdr, **scatter_kwargs)
    ax2[1,1].scatter(qr, scatt.diff_atten, **scatter_kwargs)
    ax2[1,2].scatter(qr, scatt.kdp, **scatter_kwargs)
    ax3[0].scatter(scatt.kdp, scatt.atten, **scatter_kwargs)
    ax3[1].scatter(scatt.kdp, scatt.diff_atten, **scatter_kwargs)

ax2[0,0].set_xlabel('Rain Water (%s)' % qr.dimensionality.latex)
ax2[0,0].set_ylabel('Reflectivity Factor (%s)' % z.dimensionality.latex)
ax2[0,0].grid()
ax2[0,0].set_xlim(0, 35)
ax2[0,0].set_ylim(0, 90)

ax2[0,1].set_xlabel('Rain Water (%s)' % qr.dimensionality.latex)
ax2[0,1].set_ylabel('Attenuation (%s)' % atten.dimensionality.latex)
ax2[0,1].grid()
ax2[0,1].set_ylim(0, max_atten)

ax2[0,2].set_xlabel('Rain Water (%s)' % qr.dimensionality.latex)
ax2[0,2].set_ylabel('Differential Backscatter Phase (%s)'
        % delta.dimensionality.latex)
ax2[0,2].grid()

ax2[1,0].set_xlabel('Rain Water (%s)' % qr.dimensionality.latex)
ax2[1,0].set_ylabel('Differential Reflectivity (%s)' % zdr.dimensionality.latex)
#ax2[1,0].set_ylim(-6,6)
ax2[1,0].grid()

ax2[1,1].set_xlabel('Rain Water (%s)' % qr.dimensionality.latex)
ax2[1,1].set_ylabel('Differential Attenuation (%s)'
        % diff_atten.dimensionality.latex)
ax2[1,1].grid()
ax2[1,1].set_ylim(0, max_diff_atten)

ax2[1,2].set_xlabel('Rain Water (%s)' % qr.dimensionality.latex)
ax2[1,2].set_ylabel('KDP (%s)' % kdp.dimensionality.latex)
ax2[1,2].grid()
ax2[1,2].set_ylim(0, max_kdp)

ax3[0].set_xlabel('KDP (%s)' % kdp.dimensionality.latex)
ax3[0].set_ylabel('Attenuation (%s)' % atten.dimensionality.latex)
ax3[0].grid()
ax3[0].set_xlim(0, max_kdp)
ax3[0].set_ylim(0, max_atten)

ax3[1].set_xlabel('KDP (%s)' % kdp.dimensionality.latex)
ax3[1].set_ylabel('Differential Attenuation (%s)'
        % diff_atten.dimensionality.latex)
ax3[1].grid()
ax3[1].set_xlim(0, max_kdp)
ax3[1].set_ylim(0, max_diff_atten)

plt.show()
