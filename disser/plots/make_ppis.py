#!/usr/bin/env python
import numpy as np
from numpy.ma import masked_array
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.constants import degree, kilo, centi
from mpl_toolkits.axes_grid import make_axes_locatable, Divider
from mpl_toolkits.axes_grid.axes_size import Scaled, Fixed
#from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.offsetbox import AnchoredText


from metpy.vis import ctables

def remainder(x, div):
    n = np.round(x / div)
    return x - n * div


#
# TODO: Need to make these time series estimators/correlation routine
# transpose arrays or do calculation over different axes so that the lag
# is the first axis. This will remove a *ton* of ellipsis.
#

def correlate(x, y=None, maxlag=None, scale=None):
    if y is None:
        y = x

    x = x.conj()

    M = max(x.shape[-1], y.shape[-1])
    if maxlag is None:
        maxlag = M - 1

    R = np.empty(x.shape[:-1] + (maxlag + 1,), dtype=y.dtype)
    R[..., 0] = (x*y).sum(axis=-1)
    if maxlag > 0:
        for l in xrange(1, maxlag+1):
            R[..., l] = (x[..., :-l] * y[..., l:]).sum(axis=-1)

    if scale == 'biased':
        R /= M
    elif scale == 'unbiased':
        R /= np.arange(M, M - maxlag - 1, -1)
    elif scale is not None:
        raise ValueError("Scale must be 'biased', 'unbiased', or None")

    return R

def auto_moments(ts, nyq):
    nyq = nyq[:, np.newaxis]
    ac = correlate(ts, maxlag=2, scale='biased')

    pwr = np.abs(ac[..., 0])
    vel = np.angle(ac[..., 1]) * -nyq / np.pi

    log_rat = np.log(np.abs(ac[..., 1] / ac[..., 2]))
    sig = 2 * nyq / (np.sqrt(6) * np.pi) * np.sqrt(np.abs(log_rat))

    return pwr, vel, sig

def auto_dual_pol(H_series, V_series, noise_h=0, noise_v=0, lag_base=0):
    maxlag = lag_base
    Rhh = correlate(H_series, maxlag=maxlag, scale='biased')
    Rvv = correlate(V_series, maxlag=maxlag, scale='biased')
    Rhv1 = correlate(H_series, V_series, maxlag=maxlag, scale='biased')
    if lag_base == 0:
        return dual_pol_covar(Rhh, Rvv, Rhv1, noise_h, noise_v)
    else:
        Rhv2 = correlate(V_series, H_series, maxlag=maxlag, scale='biased')
        return dual_pol_covar1(Rhh, Rvv, Rhv1, Rhv2)

def dual_pol_covar(Rhh, Rvv, Rhv, noise_h, noise_v):
    Zdr = 10.0 * np.log10(np.abs(Rhh[..., 0] - noise_h)
        / np.abs(Rvv[..., 0] - noise_v))

    rho_hv = Rhv[..., 0] / np.sqrt(np.abs(Rhh[..., 0] * Rvv[..., 0]))

    phi_dp = np.angle(Rhv[..., 0])

    return Zdr, rho_hv, phi_dp

def dual_pol_covar1(Rhh, Rvv, Rhv1, Rhv2):
    Zdr = np.abs(Rhh[..., 1]) / np.abs(Rvv[..., 1])

    rho_hv = (np.abs(Rhv1[..., 1]) + np.abs(Rhv2[..., 1])) / (2 *
        np.sqrt(np.abs(Rhh[..., 1] * Rvv[..., 1])))

    phi_dp = np.angle(Rhv1[..., 0])

    return Zdr, rho_hv, phi_dp

def ppi_plot(x, y, data, cmap=None, norm=None, ax=None, rings=None):
    if ax is None:
        ax = plt.gca()

    mesh = ax.pcolormesh(x, y, data, cmap=cmap, norm=norm)

    if rings is None:
        rings = range(10, 70, 10)

    ring_patches = []
    for rng in rings:
        c = Circle(xy=(0,0), radius=rng, fill=False)
        ax.add_patch(c)
        ring_patches.append(c)

    return mesh, ring_patches

def create_multipanel_plot(size, dpi, shape, layout, var_info, cmap, lims):
        fig = plt.figure(figsize=size, dpi=dpi)
        rings = []

        # the rect parameter will be ignore as we will set axes_locator
        rect = (0.08, 0.08, 0.9, 0.9)
        nrow,ncol = shape

        # divide the axes rectangle into grid whose size is specified
        # by horiz * vert
        horiz = [Scaled(1.)]
        for i in range(ncol - 1):
            horiz.extend([Fixed(.2), Scaled(1.)])

        vert = [Scaled(.1), Fixed(.35), Scaled(1.)]
        for i in range(nrow - 1):
            vert.extend([Fixed(.1), Scaled(1.)])

        divider = Divider(fig, rect, horiz, vert, aspect=False)

#        ax0 = fig.add_axes(rect, label="0")
#        ax0.set_aspect('equal', 'datalim')
#        ax = [ax0] + [fig.add_axes(rect, label="%d"%i, sharex=ax0, sharey=ax0)
#            for i in range(1,6)]
        ax = [fig.add_axes(rect, label="%d"%i) for i in range(len(layout))]
        cax = [fig.add_axes(rect, label='cb%d'%i) for i in range(ncol)]

        for i,a in enumerate(ax):
#            a.set_axes_locator(divider.new_locator(nx=(i // nrow) * 2,
#                ny=((i%nrow) + 1) * 2))
            a.set_axes_locator(divider.new_locator(nx=(i % ncol) * 2,
                ny=(nrow - (i // ncol)) * 2))
            a.set_aspect('equal', 'datalim')

        for i,a in enumerate(cax):
            a.set_axes_locator(divider.new_locator(nx=2 * i, ny=0))

        for num,(a,(data, label, var)) in enumerate(zip(ax, layout)):
            norm,ticks,units = var_info[var]
            ppi_plot(init_data.xlocs, init_data.ylocs, data, norm=norm,
                cmap=cmap, ax=a, rings=rings)
#            a.set_title('%s (%s)' % (moment, units))

            if num >= ncol:
                a.set_xlabel('X Distance (km)')
                cbar = ColorbarBase(ax=cax[num%ncol], norm=norm, cmap=cmap,
                    orientation='horizontal')
                cbar.set_label('%s (%s)' % (label, units))
                cbar.set_ticks(ticks)
            else:
                a.xaxis.set_major_formatter(plt.NullFormatter())

            if num % ncol == 0:
                a.set_ylabel('Y Distance (km)')
            else:
                a.yaxis.set_major_formatter(plt.NullFormatter())

            if lims:
                a.xaxis.set_major_locator(plt.MultipleLocator(lims[0]))
                a.yaxis.set_major_locator(plt.MultipleLocator(lims[0]))
                a.set_xlim(*lims[1:3])
                a.set_ylim(*lims[3:])

            # loc = 2 is upper left. TODO: Should patch matplotlib to use
            # same strings as legend
            at = AnchoredText("%s)" % chr(97 + num), loc=2, prop=dict(size='large'),
                frameon=True)
#            at.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            a.add_artist(at)

        return fig

def data_to_fname(data):
    return '%s_%s_%s' % (data.radar, data.scattering, data.drop_model)

if __name__ == '__main__':
    import sys
    import os.path
    import glob
    from matplotlib.colors import Normalize
    from matplotlib.colorbar import ColorbarBase
    from matplotlib import rcParams

    cmap = ctables.Carbone42
    pwr_norm = Normalize(-115, -25)
    pwr_ticks = np.arange(-115, -20, 15)
    ref_norm = Normalize(-20, 70)
    ref_ticks = np.arange(-20, 80, 10)
    vel_norm = Normalize(-30, 30)
    vel_ticks = np.arange(-25, 30, 10)
    spw_norm = Normalize(0, 15)
    spw_ticks = np.arange(0, 20, 5)
    zdr_norm = Normalize(-5, 5)
    zdr_ticks = np.arange(-4, 12, 2)
    kdp_norm = Normalize(-5, 25)
    kdp_ticks = np.arange(-25, 30, 5)
    rhohv_norm = Normalize(0.98, 1.0)
    diff_atten_norm = Normalize(0, 5)
    diff_atten_ticks = np.arange(0, 6, 1)

    phidp_norms = dict(COMMAS=Normalize(-180, 180), ARPS=Normalize(-180, 180))
    phidp_ticks = dict(COMMAS=np.arange(-180, 180, 30), ARPS=np.arange(-180, 200, 60))
    atten_norms = dict(COMMAS=Normalize(0, 10), ARPS=Normalize(0, 30))
    atten_ticks = dict(COMMAS=np.arange(0, 12, 2), ARPS=np.arange(0, 35, 5))

    limits = dict(COMMAS=(5, -15, 15, 0, 30), ARPS=(10, -30, 0, -20, 10))
    limits = dict(COMMAS=(5, -24, -9, -14, 1), ARPS=(10, -30, 0, -20, 10))

    rcParams['font.size'] = 8

    np.seterr(all='ignore')
#    for fname in glob.glob('data/*.nc'):
#    for fname in glob.glob('*.nc'):
    for fname in glob.glob('../../svn_docs/prospectus/figures/data/Sband_COMMAS_20100416_153716.nc'):
        print 'Processing:', fname
        init_data = process_data(fname)

        moment_info = {'power':(pwr_norm, pwr_ticks, 'dBm$$'),
                       'ref':(ref_norm, ref_ticks, 'dBZ$$'),
                       'vel':(vel_norm, vel_ticks, r'$m\,s^{-1}$'),
                       'spw':(spw_norm, spw_ticks, r'$m\,s^{-1}$'),
                       'zdr':(zdr_norm, zdr_ticks, 'dB$$'),
                       'phidp':(phidp_norms[init_data.model], phidp_ticks[init_data.model], 'deg$$'),
                       'kdp':(kdp_norm, kdp_ticks, '$deg\,km^{-1}$'),
                       'atten':(atten_norms[init_data.model], atten_ticks[init_data.model], 'dB$$'),
                       'diff_atten':(diff_atten_norm, diff_atten_ticks, 'dB$$')}

        size = (6.3, 4.)
        screen_dpi = 107
        save_dpi = 400

        info = [(init_data.pwr_h, 'Power', 'power'),
                (init_data.vel, 'Velocity', 'vel'),
                (init_data.spw, 'Spectrum Width', 'spw'),
                (init_data.pwr_h_ts_dbm, 'Power', 'power'),
                (init_data.vel_h_ts, 'Velocity', 'vel'),
                (init_data.spw_h_ts, 'Spectrum Width', 'spw')]

        fig = create_multipanel_plot(size, screen_dpi, (2, 3), info,
            moment_info, cmap, limits.get(init_data.model, None))

        save_name = data_to_fname(init_data) + '_basic.png'
        fig.savefig(save_name, dpi=save_dpi, bbox_inches='tight')
        fig.clf()

        info = [(init_data.ref_h, '$Z_H$', 'ref'),
                (init_data.zdr, '$Z_{DR}$', 'zdr'),
                (init_data.atten, '$A_H\,or\,A_D$', 'atten'),
                (init_data.ref_h_ts, '$Z_H $', 'ref'),
                (init_data.zdr_ts, '$Z_{DR}$', 'zdr'),
                (init_data.diff_atten, '$A_H\,or\,A_D$', 'atten')]

        fig = create_multipanel_plot(size, screen_dpi, (2, 3), info,
            moment_info, cmap, limits.get(init_data.model, None))

        save_name = data_to_fname(init_data) + '_atten.png'
        fig.savefig(save_name, dpi=save_dpi, bbox_inches='tight')
        fig.clf()

        info = [(init_data.phi_dp, '$\phi_{DP}$', 'phidp'),
                (init_data.kdp, '$K_{DP}$', 'kdp'),
                (init_data.phi_dp_ts, '$\phi_{DP}$', 'phidp'),
                (init_data.kdp_ts, '$K_{DP}$', 'kdp')]

        fig = create_multipanel_plot(size, screen_dpi, (2, 2), info,
            moment_info, cmap, limits.get(init_data.model, None))

        save_name = data_to_fname(init_data) + '_phase.png'
        fig.savefig(save_name, dpi=save_dpi, bbox_inches='tight')
        fig.clf()
