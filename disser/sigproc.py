import numpy as np
from quantities import degrees

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

    return Zdr, np.abs(rho_hv), phi_dp


def dual_pol_covar1(Rhh, Rvv, Rhv1, Rhv2):
    Zdr = np.abs(Rhh[..., 1]) / np.abs(Rvv[..., 1])

    rho_hv = (np.abs(Rhv1[..., 1]) + np.abs(Rhv2[..., 1])) / (2 *
        np.sqrt(np.abs(Rhh[..., 1] * Rvv[..., 1])))

    phi_dp = np.angle(Rhv1[..., 0])

    return Zdr, rho_hv, phi_dp


def shift_phi(phi):
    from quantities import degrees
    phi[phi < -15. * degrees] += 360 * degrees
