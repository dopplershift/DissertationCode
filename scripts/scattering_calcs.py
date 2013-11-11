import matplotlib.pyplot as plt
import numpy as np
import scattering
import dsd
import scipy.constants as consts
from disser.units import to_dB, to_dBz, angle, exp_to_dB
import quantities as pq

d = (np.linspace(0.01, 8, 100).reshape(-1, 1, 1).astype(np.float32) 
        * pq.mm)
lam = (np.linspace(0.01, 15, 40).reshape(1, -1, 1).astype(np.float32)
        / pq.mm)
nr_lut = (np.logspace(-2, 6, 110).reshape(1, 1,-1).astype(np.float32)
        / (pq.m**3 * pq.cm))
d0 = (np.linspace(0.01, 4, 40).reshape(1, -1, 1).astype(np.float32)
        * pq.mm)
#dist = dsd.volume_gamma(d, d0, nr_lut * pq.cm, nu=-0.8)
#dist = dsd.modified_gamma(d, lam, nr_lut, dsd.constrained_gamma_shape(lam))
#dist = dsd.modified_gamma(d, lam, nr_lut, 0.8)
#dist = dsd.mp_from_lwc(d, pq.Quantity(10.0, 'g/m^3'))
qr = np.linspace(0.01, 15, 100).reshape(1, -1) * pq.g / pq.m**3
dist = dsd.mp_from_lwc(d, qr)
fallspeed = dsd.rain_fallspeed(d)
print fallspeed.min(), fallspeed.max()
rr = dsd.rainrate(d, dist, fallspeed).rescale('mm/hr')
print rr.min(), rr.max()

for title,wavelength in [('S-band', 10.8*pq.cm),('C-band', 5.3*pq.cm),
        ('X-band', 3.3*pq.cm)]:
    temps = np.arange(0, 35, 10)
    temp_colors = {0:'b', 10:'g', 20:'y', 30:'r'}

    d_plot = d.squeeze()

    z = np.empty((temps.size,) + dist.shape[1:], dtype=np.float32)
    zdr = np.empty_like(z)
    atten = np.empty_like(z)
    diff_atten = np.empty_like(z)
    kdp = np.empty_like(z)
    delta = np.empty_like(z)

    max_kdp, max_atten, max_diff_atten = {'C-band':(30,5,2),
            'S-band':(15, 0.5, 0.25), 'X-band':(50,15,5)}[title]

    csec_units = 'mm**2'

    fig, ax = plt.subplots(2, 3, sharex=True, figsize=(10,8))
    fig.suptitle(title)
    fig.subplots_adjust(wspace=0.4)
    fig2, ax2 = plt.subplots(2, 3, sharex=True, figsize=(10,8))
    fig2.suptitle(title)
    fig2.subplots_adjust(wspace=0.4)
    fig3, ax3 = plt.subplots(1, 2, sharex=True, figsize=(10,8))
    fig3.suptitle(title)
    fig4, ax4 = plt.subplots(2, 3, sharex=True, figsize=(10,8))
    fig4.suptitle(title)
    fig4.subplots_adjust(wspace=0.4)

    for ind, t in enumerate(temps):
        scatt = scattering.scatterer(wavelength, t, 'water', diameters=d,
                shape='oblate')
        scatt.set_scattering_model('tmatrix')
        plot_label = u'%d\u00B0C' % t
        temp_color = temp_colors[int(t)]

        sig_b = scatt.sigma_bh.squeeze().rescale(csec_units)
        ax[0,0].plot(d_plot, sig_b, temp_color, label=plot_label)

        sig_e = scatt.sigma_eh.squeeze().rescale(csec_units)
        ax[0,1].plot(d_plot, sig_e, temp_color, label=plot_label)

        delta_d = -angle(-scatt.S_bkwd[0,0].conj() * scatt.S_bkwd[1,1],
                deg=True).squeeze()
        ax[0,2].plot(d_plot, delta_d, temp_color, label=plot_label)

        diff_sig_b = to_dB(scatt.sigma_bh / scatt.sigma_bv).squeeze()
        ax[1,0].plot(d_plot, diff_sig_b, temp_color, label=plot_label)

        diff_sig_e = to_dB(scatt.sigma_eh / scatt.sigma_ev).squeeze()
        ax[1,1].plot(d_plot, diff_sig_e, temp_color, label=plot_label)

        diff_phase = (scatt.S_frwd[0, 0] -
                scatt.S_frwd[1, 1]).real.squeeze().rescale('mm')
        ax[1,2].plot(d_plot, diff_phase, temp_color, label=plot_label)

        z = to_dBz(scatt.get_reflectivity_factor(dist, polar='h'))
        zdr = (z - to_dBz(scatt.get_reflectivity_factor(dist,
            polar='v'))).rescale('dB')
        atten = (scatt.get_attenuation(dist, polar='h').rescale('1/km') * exp_to_dB)
        diff_atten = (atten -
                scatt.get_attenuation(dist, polar='v').rescale('1/km') * exp_to_dB)
        kdp = ((scatt.get_propagation_wavenumber(dist, polar='h') -
                scatt.get_propagation_wavenumber(dist, polar='v')) * pq.rad).rescale(
                        'deg/km')
        delta = (scatt.get_backscatter_differential_phase(dist) * pq.rad).rescale('deg')
    #    temp_color = plt.get_cmap('RdBu')(plt.Normalize(0,30)(t))
        scatter_kwargs = dict(c=temp_color, edgecolor='none', alpha=0.4)
        ax2[0,0].scatter(qr, z, **scatter_kwargs)
        ax2[0,1].scatter(qr, atten, **scatter_kwargs)
        ax2[0,2].scatter(qr, delta, **scatter_kwargs)
        ax2[1,0].scatter(qr, zdr, **scatter_kwargs)
        ax2[1,1].scatter(qr, diff_atten, **scatter_kwargs)
        ax2[1,2].scatter(qr, kdp, **scatter_kwargs)
        ax3[0].scatter(kdp, atten, **scatter_kwargs)
        ax3[1].scatter(kdp, diff_atten, **scatter_kwargs)
        ax4[0,0].scatter(rr, z, **scatter_kwargs)
        ax4[0,1].scatter(rr, atten, **scatter_kwargs)
        ax4[0,2].scatter(rr, delta, **scatter_kwargs)
        ax4[1,0].scatter(rr, zdr, **scatter_kwargs)
        ax4[1,1].scatter(rr, diff_atten, **scatter_kwargs)
        ax4[1,2].scatter(rr, kdp, **scatter_kwargs)

    ax[0,0].set_xlabel('Diameter (%s)' % d_plot.dimensionality.latex)
    ax[0,0].set_ylabel('Backscatter Cross-Section (%s)' % sig_b.dimensionality.latex)
    ax[0,0].legend(loc='upper left')
    ax[0,0].grid()
    ax[0,0].set_xlim(0, d_plot.max().magnitude)

    ax[0,1].set_xlabel('Diameter (%s)' % d_plot.dimensionality.latex)
    ax[0,1].set_ylabel('Horizontal Extinction Cross-Section (%s)'
            % sig_e.dimensionality.latex)
    ax[0,1].grid()

    ax[0,2].set_xlabel('Diameter (%s)' % d_plot.dimensionality.latex)
    ax[0,2].set_ylabel('Differential Backscatter Phase (%s)'
            % delta_d.dimensionality.latex)
    ax[0,2].grid()

    ax[1,0].set_xlabel('Diameter (%s)' % d_plot.dimensionality.latex)
    ax[1,0].set_ylabel('Differential Backscatter Cross-Section (%s)'
            % diff_sig_b.dimensionality.latex)
    ax[1,0].grid()

    ax[1,1].set_xlabel('Diameter (%s)' % d_plot.dimensionality.latex)
    ax[1,1].set_ylabel('Differential Extinction Cross-Section (%s)'
            % diff_sig_e.dimensionality.latex)
    ax[1,1].grid()

    ax[1,2].set_xlabel('Diameter (%s)' % d_plot.dimensionality.latex)
    ax[1,2].set_ylabel('Differential Wavenumber (%s)' % diff_phase.dimensionality.latex)
    ax[1,2].grid()

    ax2[0,0].set_xlabel('Rain Water (%s)' % qr.dimensionality.latex)
    ax2[0,0].set_ylabel('Reflectivity Factor (%s)' % z.dimensionality.latex)
    ax2[0,0].grid()
    ax2[0,0].set_xlim(0, 15)
    ax2[0,0].set_ylim(0, 70)

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

    ax4[0,0].set_xlabel('Rain Rate (%s)' % rr.dimensionality.latex)
    ax4[0,0].set_ylabel('Reflectivity Factor (%s)' % z.dimensionality.latex)
    ax4[0,0].grid()
    ax4[0,0].xaxis.set_ticks(np.arange(0, 500, 100))
    ax4[0,0].set_xlim(0, 400)
    ax4[0,0].set_ylim(0, 70)

    ax4[0,1].set_xlabel('Rain Rate (%s)' % rr.dimensionality.latex)
    ax4[0,1].set_ylabel('Attenuation (%s)' % atten.dimensionality.latex)
    ax4[0,1].grid()
    ax4[0,1].set_ylim(0, max_atten)

    ax4[0,2].set_xlabel('Rain Rate (%s)' % rr.dimensionality.latex)
    ax4[0,2].set_ylabel('Differential Backscatter Phase (%s)'
            % delta.dimensionality.latex)
    ax4[0,2].grid()

    ax4[1,0].set_xlabel('Rain Rate (%s)' % rr.dimensionality.latex)
    ax4[1,0].set_ylabel('Differential Reflectivity (%s)' % zdr.dimensionality.latex)
    #ax4[1,0].set_ylim(-6,6)
    ax4[1,0].grid()

    ax4[1,1].set_xlabel('Rain Rate (%s)' % rr.dimensionality.latex)
    ax4[1,1].set_ylabel('Differential Attenuation (%s)'
            % diff_atten.dimensionality.latex)
    ax4[1,1].grid()
    ax4[1,1].set_ylim(0, max_diff_atten)

    ax4[1,2].set_xlabel('Rain Rate (%s)' % rr.dimensionality.latex)
    ax4[1,2].set_ylabel('KDP (%s)' % kdp.dimensionality.latex)
    ax4[1,2].grid()
    ax4[1,2].set_ylim(0, max_kdp)

    fig.savefig('%s-cross-sections.png' % title, dpi=150)
    fig2.savefig('%s-vars-rainwater.png' % title, dpi=150)
    fig4.savefig('%s-vars-rainrate.png' % title, dpi=150)
    fig3.savefig('%s-atten-kdp.png' % title, dpi=150)
