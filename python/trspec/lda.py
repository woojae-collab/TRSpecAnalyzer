"""Lifetime density analysis (Tikhonov-regularized) and its L-curve.

Port of computeLDA and computeLcurve of TRSpecAnalyzer.m.
"""
import numpy as np
from ._util import Result, lstsq_min_norm
from .kinetics import exp_irf_conv


def _t_select(t, t_range):
    if t_range is None or np.size(t_range) == 0 or not np.all(np.isfinite(t_range)):
        return np.ones(t.size, dtype=bool)
    sel = (t >= np.min(t_range)) & (t <= np.max(t_range))
    return sel if sel.any() else np.ones(t.size, dtype=bool)


def compute_lda(D, t, tau_grid, t0=0.0, fwhm=0.15, alpha=0.05, reg_type='l2',
                irf_mode='irf', t_range=None):
    """Decompose D(lambda, t) = sum_k A(lambda, tau_k) exp(-(t - t0)/tau_k) (x) IRF
    on a fixed lifetime grid.

    alpha    : regularization strength
    reg_type : 'l2' (penalty on A) or 'l2deriv' (penalty on its first
               derivative along log tau)
    irf_mode : 'irf' or 'none'

    Returns a Result with amap (n_wavelength x n_tau), drec (reconstruction),
    rms, ratio (rms / rms of the data), tau_grid, alpha, reg_type.
    """
    D = np.asarray(D, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    tau_grid = np.asarray(tau_grid, dtype=float).ravel()
    K, N, M = tau_grid.size, t.size, D.shape[0]
    if irf_mode is None or irf_mode == '':
        irf_mode = 'irf'
    t_sel = _t_select(t, t_range)
    E = np.zeros((N, K))
    if isinstance(irf_mode, str) and irf_mode.lower() == 'none':
        tt = t - t0
        stp = (tt >= 0).astype(float)
        for k in range(K):
            E[:, k] = stp * np.exp(-np.maximum(tt, 0) / tau_grid[k])
    else:
        for k in range(K):
            E[:, k] = exp_irf_conv(t, tau_grid[k], t0, fwhm)
    if str(reg_type).lower() == 'l2deriv' and K >= 2:
        ws = np.mean(np.diff(np.log(tau_grid)))
        L = (-np.eye(K) + np.diag(np.ones(K - 1), 1)) / ws
        L = L[:-1, :]
    else:
        L = np.eye(K)
    e_aug = np.vstack([E[t_sel, :], alpha * L])
    d_clean = D[:, t_sel].copy()
    d_clean[~np.isfinite(d_clean)] = 0
    d_aug = np.vstack([d_clean.T, np.zeros((L.shape[0], M))])
    try:
        a_t = lstsq_min_norm(e_aug, d_aug)
    except Exception:
        a_t = np.linalg.pinv(e_aug) @ d_aug
    amap = a_t.T
    amap[~np.isfinite(amap)] = 0
    drec = amap @ E.T
    R = D[:, t_sel] - drec[:, t_sel]
    fin = np.isfinite(R)
    rms = float(np.sqrt(np.mean(R[fin] ** 2))) if fin.any() else np.nan
    dw = D[:, t_sel]
    fd = np.isfinite(dw)
    denom = float(np.sqrt(np.mean(dw[fd] ** 2))) if fd.any() else np.nan
    ratio = rms / denom if denom > 0 else np.nan
    return Result(amap=amap, drec=drec, rms=rms, ratio=ratio, tau_grid=tau_grid,
                  alpha=alpha, reg_type=reg_type, K=K, M=M, N=N)


def compute_lcurve(D, t, tau_grid, t0=0.0, fwhm=0.15, reg_type='l2', alpha_list=None,
                   irf_mode='irf', t_range=None):
    """L-curve of the lifetime density analysis.

    Returns (alpha_list, res_norm, sol_norm, alpha_corner); the corner is
    the point of largest curvature of log(sol_norm) against log(res_norm).
    """
    D = np.asarray(D, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    if alpha_list is None or np.size(alpha_list) == 0:
        alpha_list = np.logspace(-6, 2, 25)
    alpha_list = np.asarray(alpha_list, dtype=float).ravel()
    t_sel = _t_select(t, t_range)
    n = alpha_list.size
    res_norm = np.zeros(n)
    sol_norm = np.zeros(n)
    for ii in range(n):
        r = compute_lda(D, t, tau_grid, t0, fwhm, alpha_list[ii], reg_type, irf_mode, t_range)
        R = D[:, t_sel] - r.drec[:, t_sel]
        res_norm[ii] = np.sqrt(np.sum(R[np.isfinite(R)] ** 2))
        sol_norm[ii] = np.sqrt(np.sum(r.amap[np.isfinite(r.amap)] ** 2))
    with np.errstate(divide='ignore', invalid='ignore'):
        x = np.log(res_norm)
        y = np.log(sol_norm)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    alpha_ok = alpha_list[ok]
    if x.size < 4:
        return alpha_list, res_norm, sol_norm, alpha_list[int(np.floor(n / 2 + 0.5)) - 1]
    dx, dy = np.gradient(x), np.gradient(y)
    ddx, ddy = np.gradient(dx), np.gradient(dy)
    with np.errstate(divide='ignore', invalid='ignore'):
        curv = np.abs(dx * ddy - dy * ddx) / (dx ** 2 + dy ** 2) ** 1.5
    curv[~np.isfinite(curv)] = 0
    return alpha_list, res_norm, sol_norm, alpha_ok[int(np.argmax(curv))]
