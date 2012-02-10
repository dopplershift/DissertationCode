#!/usr/bin/env python
import numpy as np
from numpy.ma import masked_array
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm
from matplotlib.patches import Circle
from matplotlib.colorbar import ColorbarBase
from scipy.constants import degree, kilo

from netCDF4 import Dataset

from metpy.vis import ctables

class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

def remainder(x, div):
    n = np.round(x / div)
    return x - n * div

def nc_get_masked(ncdf_file, var):
    fill_val = ncdf_file.variables[var]._FillValue
    var = ncdf_file.variables[var][:]
    return masked_array(var, mask=(var==fill_val))

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

def auto_dual_pol(H_series, V_series):
    maxlag = 1
    Rhh = correlate(H_series, maxlag=maxlag, scale='biased')
    Rvv = correlate(V_series, maxlag=maxlag, scale='biased')
    Rhv1 = correlate(H_series, V_series, maxlag=maxlag, scale='biased')
    Rhv2 = correlate(V_series, H_series, maxlag=maxlag, scale='biased')
    return dual_pol_covar(Rhh, Rvv, Rhv1, Rhv2)

def dual_pol_covar(Rhh, Rvv, Rhv1, Rhv2):
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

def process_data(fname):
    data = Bunch()

    nc = Dataset(fname, 'r')
    data.nc = nc

    data.az = nc.variables['Azimuth'][:]
    data.gate_length = nc.variables['GateLength'].getValue()
    data.pulse_length = nc.variables['PulseDuration'][0] * 3e8 / 2.
    data.rng = (np.arange(len(nc.dimensions['gates']))
        * data.gate_length + nc.variables['RangeToFirstGate'].getValue())

    data.nyquist = nc.variables['NyquistVelocity'][:]
    data.noise_pwr_db = nc.variables['MinimumDetectableSignal'].getValue()
    data.noise_pwr = 10**(data.noise_pwr_db / 10.)

    data.xlocs = data.rng * np.sin(data.az[:, np.newaxis] * degree) / kilo
    data.ylocs = data.rng * np.cos(data.az[:, np.newaxis] * degree) / kilo

    data.pwr_h = nc_get_masked(nc, 'Power_H')
    data.vel = nc_get_masked(nc, 'Velocity')
    data.spw = nc_get_masked(nc, 'SpectrumWidth')

    if 'I_H' in nc.variables:
        data.time_series = True
        data.iq_h = nc.variables['I_H'][:] + 1.0j * nc.variables['Q_H'][:]
        data.pwr_h_ts, data.vel_h_ts, data.spw_h_ts = auto_moments(data.iq_h,
            data.nyquist)
        data.pwr_h_ts_dbm = 10.0 * np.log10(data.pwr_h_ts) + 30.
        data.rad_const_h = nc.variables['RefCalibration_H'].getValue()
        data.ref_h_ts = data.pwr_h_ts_dbm + 20 * np.log10(data.rng) + data.rad_const_h - 30.
        data.ref_h = nc_get_masked(nc, 'Reflectivity_H')

        data.iq_v = nc.variables['I_V'][:] + 1.0j * nc.variables['Q_V'][:]
        data.pwr_v_ts, data.vel_v_ts, data.spw_v_ts = auto_moments(data.iq_v,
            data.nyquist)
        data.pwr_v_ts_dbm = 10.0 * np.log10(data.pwr_v_ts) + 30.

        data.rad_const_v = nc.variables['RefCalibration_V'].getValue()
        data.ref_v_ts = data.pwr_v_ts_dbm + 20 * np.log10(data.rng) + data.rad_const_v - 30.
        data.zdr_ts = 10.0 * np.log10(np.abs(data.pwr_h_ts - data.noise_pwr)
            / np.abs(data.pwr_v_ts - data.noise_pwr))
        data.zdr_ts2, data.rho_hv, data.phi_dp = auto_dual_pol(data.iq_h,
            data.iq_v)
        data.rho_hv = np.ma.array(data.rho_hv, mask=data.rho_hv>1.)
        data.phi_dp /= degree
    else:
        data.time_series = False

    if 'Power_V' in nc.variables:
        data.dual_pol = True

        data.pwr_v = nc_get_masked(nc, 'Power_V')
        data.ref_v = nc_get_masked(nc, 'Reflectivity_V')

        data.zdr = data.ref_h - data.ref_v
        data.kdp = nc_get_masked(nc, 'KDP')
#        kdp[np.isnan(kdp)] = np.ma.masked
#        kdp[np.abs(kdp)>1e10] = np.ma.masked
        data.kdp = data.kdp * (kilo / degree)
#        print data.kdp
#        data.phi_dp2 = 2 * data.kdp.cumsum(axis=1) * (data.pulse_length / kilo)
        data.phi_dp2 = nc.variables['PhiDP'][:] / degree
        data.phi_dp2 = remainder(data.phi_dp2, 360.)
        thgrad, rgrad = np.gradient(data.phi_dp2, 1, data.gate_length)
        data.kdp2 = rgrad * kilo / 2.
    else:
        data.dual_pol = False

    if 'MeanQR' in nc.variables:
        data.means = True
        
        data.mean_atten = -np.log(nc.variables['MeanAtten'][:].mean(axis=-1)) * 10.0 * np.log10(np.e)
        data.mean_diff_atten = -np.log(nc.variables['MeanDiffAtten'][:].mean(axis=-1)) * 10.0 * np.log10(np.e)
#        print data.mean_diff_atten.min(), data.mean_diff_atten.max()
        data.mean_phidp = nc.variables['MeanPhiDP'][:].mean(axis=-1) / degree
        data.mean_qr = nc.variables['MeanQR'][:].mean(axis=-1) * kilo
        data.std_qr = nc.variables['StdQR'][:].mean(axis=-1)
        data.mean_nr = nc.variables['MeanNR'][:].mean(axis=-1)
        data.mean_temp = nc.variables['MeanTemp'][:].mean(axis=-1)
        data.std_temp = nc.variables['StdTemp'][:].mean(axis=-1)
        print data.std_qr.min(), data.std_qr.max()
#        data.mean_nr[data.mean_nr<0] = 0.
    else:
        data.means = False

    return data

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    opt_parser = OptionParser(usage='usage: %prog [options] datafile')

    opts,args = opt_parser.parse_args()
    if not args:
        print "Must specify a data file to use."
        opt_parser.print_help()
        sys.exit(-1)

    fname = args[0]

    init_data = process_data(fname)
    nrow = 2
    ncol = 3

    pwr_norm = Normalize(-115, -25)
    ref_norm = Normalize(-20, 70)
    vel_norm = Normalize(-init_data.nyquist.max(), init_data.nyquist.max())
    spw_norm = Normalize(0, init_data.nyquist.max())
    cmap = ctables.Carbone42
    rings = []

    fig = plt.figure()
    fig.suptitle(fname, fontsize=14)
    fig.canvas.manager.set_window_title('Basic Moments')

    # Plot the powers
    ax = fig.add_subplot(nrow, ncol, 1)
    ax.set_aspect('equal', 'datalim')
    ppi_plot(init_data.xlocs, init_data.ylocs, init_data.pwr_v_ts_dbm, norm=pwr_norm, cmap=cmap, ax=ax, rings=rings)
    ax.set_title('Power (dBm)')

    # Plot the vels
    ax3 = fig.add_subplot(nrow, ncol, 2, sharex=ax, sharey=ax)
    ppi_plot(init_data.xlocs, init_data.ylocs, init_data.vel_v_ts, norm=vel_norm, cmap=cmap, ax=ax3, rings=rings)
    ax3.set_title('Velocity (m/s)')

    # Plot the spectrum widths
    ax5 = fig.add_subplot(nrow, ncol, 3, sharex=ax, sharey=ax)
    ppi_plot(init_data.xlocs, init_data.ylocs, init_data.spw_v_ts, norm=spw_norm, cmap=cmap, ax=ax5, rings=rings)
    ax5.set_title('Spectrum Width (m/s)')

    # Plot time_series data if present
    if init_data.time_series:
        ax2 = fig.add_subplot(nrow, ncol, 4, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.pwr_h_ts_dbm, norm=pwr_norm, cmap=cmap, ax=ax2, rings=rings)
        ax2.set_title('TS Power (dBm)')

        ax4 = fig.add_subplot(nrow, ncol, 5, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.vel_h_ts, norm=vel_norm, cmap=cmap, ax=ax4, rings=rings)
        ax4.set_title('TS Velocity (m/s)')

        ax6 = fig.add_subplot(nrow, ncol, 6, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.spw_h_ts, norm=spw_norm, cmap=cmap, ax=ax6, rings=rings)
        ax6.set_title('TS Spectrum Width (m/s)')

    fig.subplots_adjust(bottom=0.15)
    lpos = list(ax.get_position().bounds)
    lpos[1] = 0.04
    lpos[3] = 0.04
    cax_left = fig.add_axes(lpos)

    cbar_left = ColorbarBase(ax=cax_left,
        norm=pwr_norm, cmap=cmap, orientation='horizontal')
    cbar_left.set_label('Power (dBm)')

    mpos = list(ax3.get_position().bounds)
    mpos[1] = 0.04
    mpos[3] = 0.04
    cax_mid = fig.add_axes(mpos)
    cbar_mid = ColorbarBase(ax=cax_mid,
        norm=vel_norm, cmap=cmap, orientation='horizontal')
    cbar_mid.set_label('Velocity (m/s)')

    rpos = list(ax5.get_position().bounds)
    rpos[1] = 0.04
    rpos[3] = 0.04
    cax_right = fig.add_axes(rpos)
    cbar_right = ColorbarBase(ax=cax_right,
        norm=spw_norm, cmap=cmap, orientation='horizontal')
    cbar_right.set_label('Spectrum Width (m/s)')

    # Plot dual pol data if present
    if init_data.dual_pol and init_data.time_series:
        zdr_norm = Normalize(-5, 5)
        phidp_norm = Normalize(0, 60)
        rhohv_norm = Normalize(0.9, 1.0)
        kdp_norm = Normalize(0, 10)

        fig = plt.figure()
        fig.suptitle(fname, fontsize=14)
        fig.canvas.manager.set_window_title('Polarimetric Moments')

        ax = fig.add_subplot(2, 3, 1, sharex=ax, sharey=ax)
        ax.set_aspect('equal', 'datalim')
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.ref_h_ts, norm=ref_norm, cmap=cmap, ax=ax, rings=rings)
        ax.set_title('$Z_H$ (dBZ)')

#        ax = fig.add_subplot(2, 3, 1, sharex=ax, sharey=ax)
#        ax.set_aspect('equal', 'datalim')
#        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.zdr, norm=zdr_norm, cmap=cmap, ax=ax, rings=rings)
#        ax.set_title('$Z_{DR}$ (from ref)')

        ax4 = fig.add_subplot(2, 3, 2, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.phi_dp, norm=phidp_norm, cmap=cmap, ax=ax4, rings=rings)
        ax4.set_title('$\phi_{DP}$ (deg)')

#        ax5 = fig.add_subplot(2, 3, 3, sharex=ax, sharey=ax)
#        ppi_plot(data.xlocs, data.ylocs, rho_hv, norm=rhohv_norm, cmap=cmap, ax=ax5, rings=rings)
#        ax5.set_title(r'$\rho_{HV}$')

#        ax5 = fig.add_subplot(2, 3, 3, sharex=ax, sharey=ax)
#        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.phi_dp2, norm=phidp_norm, cmap=cmap, ax=ax5, rings=rings)
#        ax5.set_title('$\phi_{DP}$ Calculated (deg)')

        ax5 = fig.add_subplot(2, 3, 3, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.rho_hv, norm=rhohv_norm, cmap=cmap, ax=ax5, rings=rings)
        ax5.set_title(r'$\rho_{HV}$')

        ax3 = fig.add_subplot(2, 3, 4, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.zdr_ts, norm=zdr_norm, cmap=cmap, ax=ax3, rings=rings)
        ax3.set_title('$Z_{DR}$ (from ts)')

        ax2 = fig.add_subplot(2, 3, 5, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.kdp, norm=kdp_norm, cmap=cmap, ax=ax2, rings=rings)
        ax2.set_title('$K_{DP}$ (deg/km)')

        ax6 = fig.add_subplot(2, 3, 6, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.kdp2, cmap=cmap, norm=kdp_norm, ax=ax6, rings=rings)
        ax6.set_title('$K_{DP}$ 2')

#        ax5 = fig.add_subplot(2, 3, 6, sharex=ax, sharey=ax)
#        ppi_plot(data.xlocs, data.ylocs, back_phase, cmap=cmap, ax=ax5, rings=rings)
#        ax5.set_title('$\delta$')

        fig.subplots_adjust(bottom=0.15)
        lpos = list(ax.get_position().bounds)
        lpos[1] = 0.04
        lpos[3] = 0.04
        cax_left = fig.add_axes(lpos)

        cbar_left = ColorbarBase(ax=cax_left,
            norm=ref_norm, cmap=cmap, orientation='horizontal')
        cbar_left.set_label('Reflectivity (dBZ)')

        mpos = list(ax4.get_position().bounds)
        mpos[1] = 0.04
        mpos[3] = 0.04
        cax_mid = fig.add_axes(mpos)
        cbar_mid = ColorbarBase(ax=cax_mid,
            norm=phidp_norm, cmap=cmap, orientation='horizontal')
        cbar_mid.set_label('$\phi_{DP}$ (deg)')

        rpos = list(ax5.get_position().bounds)
        rpos[1] = 0.04
        rpos[3] = 0.04
        cax_right = fig.add_axes(rpos)
        cbar_right = ColorbarBase(ax=cax_right,
            norm=kdp_norm, cmap=cmap, orientation='horizontal')
#            norm=Normalize(0,3), cmap=cmap, orientation='horizontal')
        cbar_right.set_label('$K_{DP}$ (deg/km)')


    # Plot mean analysis if present
    if init_data.means:
        qr_norm = Normalize(0, 10)
        nr_norm = LogNorm(1e-4, 1e6)
        atten_norm = Normalize(0, 1)
        diff_atten_norm = Normalize(0, .1)

        fig = plt.figure()
        fig.suptitle(fname, fontsize=14)
        fig.canvas.manager.set_window_title('Means')

        ax = fig.add_subplot(2, 3, 1, sharex=ax, sharey=ax)
        ax.set_aspect('equal', 'datalim')
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.mean_qr, norm=qr_norm, cmap=cmap, ax=ax, rings=rings)
        ax.set_title('$q_r$ (g m^-3)')

        ax4 = fig.add_subplot(2, 3, 2, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.mean_atten, norm=atten_norm, cmap=cmap, ax=ax4, rings=rings)
        ax4.set_title('atten')

        ax5 = fig.add_subplot(2, 3, 3, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.mean_diff_atten, norm=diff_atten_norm, cmap=cmap, ax=ax5, rings=rings)
        ax5.set_title('diff atten')

#        ax3 = fig.add_subplot(2, 3, 4, sharex=ax, sharey=ax)
#        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.mean_nr, norm=nr_norm, cmap=cmap, ax=ax3, rings=rings)
#        ax3.set_title('Nr')

        ax3 = fig.add_subplot(2, 3, 4, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.mean_temp, norm=Normalize(273, 303), cmap=cmap, ax=ax3, rings=rings)
        ax3.set_title('Mean T')

        ax2 = fig.add_subplot(2, 3, 5, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.mean_phidp, norm=phidp_norm, cmap=cmap, ax=ax2, rings=rings)
        ax2.set_title('phidp')

        ax6 = fig.add_subplot(2, 3, 6, sharex=ax, sharey=ax)
        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.std_qr, norm=Normalize(0, 0.1), cmap=cmap, ax=ax6, rings=rings)
        ax6.set_title('std. dev. qr')

#        ax6 = fig.add_subplot(2, 3, 6, sharex=ax, sharey=ax)
#        ppi_plot(init_data.xlocs, init_data.ylocs, init_data.std_temp, norm=Normalize(0, 25), cmap=cmap, ax=ax6, rings=rings)
#        ax6.set_title('std. dev. T')

        fig.subplots_adjust(bottom=0.15)
        lpos = list(ax.get_position().bounds)
        lpos[1] = 0.04
        lpos[3] = 0.04
        cax_left = fig.add_axes(lpos)

        cbar_left = ColorbarBase(ax=cax_left,
            norm=qr_norm, cmap=cmap, orientation='horizontal')
        cbar_left.set_label('Qr')

        mpos = list(ax4.get_position().bounds)
        mpos[1] = 0.04
        mpos[3] = 0.04
        cax_mid = fig.add_axes(mpos)
        cbar_mid = ColorbarBase(ax=cax_mid,
            norm=phidp_norm, cmap=cmap, orientation='horizontal')
        cbar_mid.set_label('$\phi_{DP}$ (deg)')

        rpos = list(ax5.get_position().bounds)
        rpos[1] = 0.04
        rpos[3] = 0.04
        cax_right = fig.add_axes(rpos)
        cbar_right = ColorbarBase(ax=cax_right,
            norm=kdp_norm, cmap=cmap, orientation='horizontal')
#            norm=Normalize(0,3), cmap=cmap, orientation='horizontal')
        cbar_right.set_label('$K_{DP}$ (deg/km)')


    ax.set_xlim(-40, 10)
    ax.set_ylim(-20, 30)

    plt.show()
