"""Chirp (wavelength-dependent time zero) model and fit.

Port of chirpModel and fitChirpParams of TRSpecAnalyzer.m.
"""
import numpy as np
from ._optim import fminsearch


def chirp_model(p, w):
    """t0(w) = a * sqrt((b*w^2 - 1) / (c*w^2 - 1)) + d, p = [a, b, c, d]."""
    w = np.asarray(w, dtype=float).ravel()
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = (p[1] * w ** 2 - 1) / (p[2] * w ** 2 - 1)
        ratio = np.where((ratio < 0) | ~np.isfinite(ratio), np.nan, ratio)
        return p[0] * np.sqrt(ratio) + p[3]


def fit_chirp_params(pts):
    """Fit the chirp formula to points (wavelength, time zero).

    pts : array (n, 2), n >= 4. Returns (p, rms).
    Six starting points, Nelder-Mead each; wavelengths are scaled by their
    median during the search and b, c are scaled back at the end.
    """
    pts = np.asarray(pts, dtype=float)
    w = pts[:, 0]
    t_obs = pts[:, 1]
    if w.size < 4:
        raise ValueError('At least 4 points are required.')
    w0 = np.median(w)
    ws = w / w0

    def obj(p):
        m = chirp_model(p, ws)
        return np.nansum((t_obs - m) ** 2) + 1e6 * float(np.any(~np.isfinite(m)))

    # MATLAB max/min skip NaN
    t_range = np.nanmax(t_obs) - np.nanmin(t_obs)
    t_mid = (np.nanmax(t_obs) + np.nanmin(t_obs)) / 2
    if t_range == 0:
        t_range = 1.0
    starts = [
        [t_range, 1.2, 1.5, t_mid - t_range],
        [-t_range, 1.2, 1.5, t_mid + t_range],
        [t_range / 2, 2.0, 3.0, t_mid],
        [t_range, 1.1, 2.0, t_mid - t_range / 2],
        [2 * t_range, 0.8, 1.2, t_mid - t_range],
        [t_range, 5.0, 10.0, t_mid - t_range],
    ]
    best_p, best_r = None, np.inf
    for s in starts:
        try:
            pk, rk, _, _ = fminsearch(obj, s, tol_x=1e-9, tol_fun=1e-11,
                                      max_iter=20000, max_fun_evals=20000)
        except Exception:
            continue
        if np.all(np.isfinite(pk)) and rk < best_r:
            best_r, best_p = rk, pk
    if best_p is None:
        raise RuntimeError('Chirp fit did not converge. Try adjusting points.')
    best_p = best_p.copy()
    best_p[1] /= w0 ** 2
    best_p[2] /= w0 ** 2
    return best_p, float(np.sqrt(best_r / w.size))
