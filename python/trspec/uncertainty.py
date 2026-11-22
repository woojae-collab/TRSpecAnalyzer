"""1-sigma parameter uncertainties from a numerical Jacobian.

Port of fitParamUncertainty of TRSpecAnalyzer.m.
"""
import numpy as np


def fit_param_uncertainty(res_fun, p, rss, n_data):
    """sqrt(diag(s^2 inv(J'J))) with s^2 = rss / (n_data - M).

    res_fun : p -> residual vector (physical parameters, not log-scaled)
    p       : optimum (M values);  rss : residual sum of squares there
    n_data  : number of data points, already reduced by the number of
              linear amplitudes

    Returns an array of M values. inf marks a parameter the data do not
    determine (no effect on the residual, or acting only in combination
    with another one); NaN = could not be evaluated.
    """
    p = np.asarray(p, dtype=float).ravel()
    m = p.size
    sig = np.full(m, np.nan)
    if m == 0 or n_data <= m:
        return sig
    try:
        r0 = np.asarray(res_fun(p), dtype=float).ravel()
        jac = np.zeros((r0.size, m))
        for j in range(m):
            h = max(abs(p[j]), 1e-12) * 1e-4
            pp = p.copy()
            pp[j] = p[j] + h
            pm = p.copy()
            pm[j] = p[j] - h
            jac[:, j] = (np.asarray(res_fun(pp)).ravel() - np.asarray(res_fun(pm)).ravel()) / (2 * h)
        s2 = rss / max(1, n_data - m)
        cn = np.sqrt(np.sum(jac ** 2, axis=0))
        dead = ~(cn > 0) | ~np.isfinite(cn)
        und = dead.copy()
        live = np.flatnonzero(~dead)
        if live.size:
            jn = jac[:, live] / cn[live]
            _, sv, vt = np.linalg.svd(jn, full_matrices=False)
            v = vt.T
            flat = sv < 1e-5 * sv[0]
            if flat.any():
                w_null = np.sum(v[:, flat] ** 2, axis=1)
                und[live[w_null > 1e-6]] = True
            vg = v[:, ~flat]
            dg = np.sum((vg / sv[~flat]) ** 2, axis=1)
            with np.errstate(invalid='ignore', divide='ignore'):
                sg = np.sqrt(s2 * dg) / cn[live]
            ok = np.isfinite(sg)
            sig[live[ok]] = sg[ok]
        sig[und] = np.inf
    except Exception:
        pass
    return sig
