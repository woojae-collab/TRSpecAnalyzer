"""Kinetic basis functions: exponentials and stretched exponentials with a
Gaussian instrument response (IRF).

Port of expIRFConv, stretchedIRFConv, gaussLegendre, buildGaBasis and
buildGaBasisExt of TRSpecAnalyzer.m.
"""
import numpy as np
from scipy.special import erfc, erfcx

_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
_GL = None


def irf_mode_norm(mode):
    """'none' stays 'none'; every other name (also names of older versions)
    means the Gaussian IRF, 'numerical'."""
    if isinstance(mode, str) and mode.strip().lower() == 'none':
        return 'none'
    return 'numerical'


def matlab_round(x):
    """Round half away from zero (MATLAB round)."""
    return float(np.sign(x) * np.floor(np.abs(x) + 0.5))


def ca_options(ca=None):
    """Settings of the coherent-artifact terms with defaults filled in.
    Port of caOptions. ca: None or a dict with any of
      order : number of terms (0 = none, 1..3: Gaussian and its first
              and second derivative)
      mode  : 'irf' (width of the IRF), 'free' (own width, fitted),
              'fixed' (own width, kept)
      fwhm  : the own width (NaN: start from the IRF width)
    """
    d = {'order': 0, 'mode': 'irf', 'fwhm': np.nan}
    if not isinstance(ca, dict) or not ca:
        return d
    out = dict(d)
    for k in d:
        v = ca.get(k)
        if v is not None and not (hasattr(v, '__len__') and len(v) == 0):
            out[k] = v
    try:
        out['order'] = int(max(0, min(3, matlab_round(float(out['order'])))))
    except (TypeError, ValueError):
        out['order'] = 0
    if out['mode'] not in ('irf', 'free', 'fixed'):
        out['mode'] = 'irf'
    try:
        out['fwhm'] = float(out['fwhm'])
    except (TypeError, ValueError):
        out['fwhm'] = np.nan
    return out


def ca_basis(t, t0, fwhm, order):
    """Coherent-artifact terms at time zero: a Gaussian of the given
    FWHM and its first two derivatives (columns 1..order). Port of
    caBasis."""
    t = np.asarray(t, dtype=float).ravel()
    sg = fwhm * _FWHM_TO_SIGMA
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        x = (t - t0) / sg
        G = np.exp(-x ** 2 / 2)
        B = np.column_stack([G, -x * G, (x ** 2 - 1) * G])
    return B[:, :int(order)]


def irf_delay_step(t, t0):
    """Delay step at time zero: median of the (up to four) intervals
    next to t0. NaN if it cannot be given. Port of irfDelayStep."""
    t = np.asarray(t, dtype=float).ravel()
    d = np.unique(t[np.isfinite(t)])
    if d.size < 2 or not np.isfinite(t0):
        return np.nan
    j = int(np.argmin(np.abs(d - t0)))
    k = d[max(0, j - 2):min(d.size, j + 3)]
    return float(np.median(np.diff(k)))


def irf_unresolved(fw, fw_res):
    """True when an IRF width is at or below what the delays resolve
    (fw_res = half the delay step at time zero)."""
    return bool(fw_res > 0 and fw <= 1.2 * fw_res)


def irf_width_lost(fw, fw_std, fw_res):
    """True when a fitted IRF width carries no information: unresolved,
    or its uncertainty is not finite or exceeds the value."""
    return bool(irf_unresolved(fw, fw_res) or not np.isfinite(fw_std)
                or fw_std > fw)


def exp_irf_conv(t, tau, t0, fwhm):
    """exp(-(t - t0)/tau) for t >= t0, convolved with a normalized Gaussian
    of the given FWHM. tau = inf gives the step response; fwhm <= 0 gives
    the bare decay."""
    sigma = fwhm * _FWHM_TO_SIGMA
    dt = np.asarray(t, dtype=float).ravel() - t0
    c = np.zeros(dt.shape)
    if not np.isfinite(sigma) or sigma <= 0:
        if np.isinf(tau):
            c = (dt >= 0).astype(float)
        else:
            pos = dt >= 0
            with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
                c[pos] = np.exp(-dt[pos] / tau)
        c[~np.isfinite(c)] = 0
        return c
    if np.isinf(tau):
        c = 0.5 * erfc(-dt / (sigma * np.sqrt(2)))
        c[~np.isfinite(c)] = 0
        return c
    if not np.isfinite(tau) or tau <= 0:
        return c
    thresh = 5 * sigma + sigma ** 2 / tau
    near = np.abs(dt) <= thresh
    far_pos = dt > thresh
    with np.errstate(over='ignore', invalid='ignore', under='ignore'):
        if far_pos.any():
            dtf = dt[far_pos]
            c[far_pos] = np.exp(sigma ** 2 / (2 * tau ** 2) - dtf / tau)
        if near.any():
            dtn = dt[near]
            b = (sigma / tau - dtn / sigma) / np.sqrt(2)
            c[near] = 0.5 * np.exp(-dtn ** 2 / (2 * sigma ** 2)) * erfcx(b)
    c[~np.isfinite(c)] = 0
    return c


def gauss_legendre():
    """Nodes x (in (-1, 1)) and weights w of the 48-point Gauss-Legendre rule."""
    global _GL
    if _GL is None:
        x, w = np.polynomial.legendre.leggauss(48)
        _GL = (x.copy(), w.copy())
    return _GL


def stretched_irf_conv(t, tau, beta, t0, fwhm):
    """exp(-((t - t0)/tau)^beta) for t >= t0, convolved with a normalized
    Gaussian of the given FWHM (two-panel Gauss-Legendre quadrature)."""
    dt = np.asarray(t, dtype=float).ravel() - t0
    c = np.zeros(dt.shape)
    if not np.isfinite(tau) or tau <= 0 or not np.isfinite(beta) or beta <= 0:
        return c
    sigma = fwhm * _FWHM_TO_SIGMA
    if not np.isfinite(sigma) or sigma <= 0:
        return np.exp(-(np.maximum(dt, 0) / tau) ** beta) * (dt >= 0)
    gx, gw = gauss_legendre()
    lo = np.maximum(dt - 6 * sigma, 0)
    hi = dt + 6 * sigma
    act = hi > 0
    if act.any():
        xa = dt[act][:, None]
        lo_a = lo[act][:, None]
        hi_a = hi[act][:, None]
        spl = lo_a + np.minimum(4 * tau, (hi_a - lo_a) / 2)
        acc = np.zeros(xa.shape[0])
        for a, b in ((lo_a, spl), (spl, hi_a)):
            half = (b - a) / 2
            mid = (b + a) / 2
            u = np.maximum(mid + half * gx[None, :], 0)
            f = np.exp(-(u / tau) ** beta)
            g = np.exp(-(xa - u) ** 2 / (2 * sigma ** 2))
            acc = acc + half[:, 0] * ((f * g) @ gw)
        c[act] = acc / (sigma * np.sqrt(2 * np.pi))
    c[~np.isfinite(c)] = 0
    return c


def build_ga_basis(t, tau, t0, fwhm, has_inf=False):
    """Columns exp(-t/tau_j) (x) IRF; a last step column if has_inf."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float)).ravel()
    t = np.asarray(t, dtype=float).ravel()
    c = np.zeros((t.size, tau.size + int(bool(has_inf))))
    for j in range(tau.size):
        c[:, j] = exp_irf_conv(t, tau[j], t0, fwhm)
    if has_inf:
        c[:, -1] = exp_irf_conv(t, np.inf, t0, fwhm)
    return c


def build_ga_basis_ext(t, tau, beta, stretch_on, t0, fwhm, has_inf=False, irf_mode='numerical'):
    """Like build_ga_basis, with optional stretched columns.
    irf_mode 'none' gives the bare kinetics (no IRF at all)."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float)).ravel()
    t = np.asarray(t, dtype=float).ravel()
    n_tau = tau.size
    beta = np.ones(n_tau) if beta is None else np.atleast_1d(np.asarray(beta, dtype=float)).ravel()
    stretch_on = (np.zeros(n_tau, dtype=bool) if stretch_on is None
                  else np.atleast_1d(np.asarray(stretch_on, dtype=bool)).ravel())
    c = np.zeros((t.size, n_tau + int(bool(has_inf))))
    if isinstance(irf_mode, str) and irf_mode.lower() == 'none':
        tt = t - t0
        step = (tt >= 0).astype(float)
        for j in range(n_tau):
            if stretch_on[j]:
                c[:, j] = step * np.exp(-(np.maximum(tt, 0) / tau[j]) ** beta[j])
            else:
                c[:, j] = step * np.exp(-np.maximum(tt, 0) / tau[j])
        if has_inf:
            c[:, -1] = step
        return c
    for j in range(n_tau):
        if stretch_on[j]:
            c[:, j] = stretched_irf_conv(t, tau[j], beta[j], t0, fwhm)
        else:
            c[:, j] = exp_irf_conv(t, tau[j], t0, fwhm)
    if has_inf:
        c[:, -1] = exp_irf_conv(t, np.inf, t0, fwhm)
    return c
