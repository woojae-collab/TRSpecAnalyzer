"""Chirp (wavelength-dependent time zero): model, fit, and the
interpolation along the delays that the correction and the time shifts use.

Port of chirpModel, fitChirpParams, delayRuns, splineCurv, interpRuns and
chirpShift of TRSpecAnalyzer.m.
"""
import numpy as np
from scipy.linalg import solve_banded
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


def delay_runs(t):
    """Stretches of the delay axis with equal steps.

    Returns a list of (first, last) indices (both included) of the
    stretches of at least 4 delays whose steps all lie within 10 % of the
    first step of the stretch. Port of delayRuns."""
    t = np.asarray(t, dtype=float).ravel()
    n = t.size
    runs = []
    if n < 4:
        return runs
    d = np.diff(t)
    a = 0
    for i in range(1, n):                # step i is d[i]; i = n - 1 closes
        if i == n - 1 or not (abs(d[i] - d[a]) <= 0.10 * abs(d[a])):
            if i - a >= 3:
                runs.append((a, i))
            a = i
    return runs


def spline_curv(x, Y):
    """Second derivatives at the delays x of the cubic spline through the
    rows of Y (rows x len(x), at least 4 delays): not-a-knot ends, the
    spline of MATLAB's spline() and of scipy's CubicSpline.
    Port of splineCurv."""
    x = np.asarray(x, dtype=float).ravel()
    Y = np.atleast_2d(np.asarray(Y, dtype=float))
    n = x.size
    h = np.diff(x)
    D = np.diff(Y, axis=1) / h                    # slopes of the chords
    R = np.zeros((n, Y.shape[0]))
    R[1:n - 1, :] = 6.0 * (D[:, 1:] - D[:, :-1]).T
    ab = np.zeros((5, n))                         # ab[2 + i - j, j] = A[i, j]
    ab[3, 0:n - 2] = h[0:n - 2]                   # A[i, i - 1], i = 1 ... n - 2
    ab[2, 1:n - 1] = 2.0 * (h[0:n - 2] + h[1:n - 1])
    ab[1, 2:n] = h[1:n - 1]                       # A[i, i + 1]
    # the third derivative is continuous at the second and at the last but
    # one delay
    ab[2, 0] = h[1]
    ab[1, 1] = -(h[0] + h[1])
    ab[0, 2] = h[0]
    ab[4, n - 3] = h[n - 2]
    ab[3, n - 2] = -(h[n - 3] + h[n - 2])
    ab[2, n - 1] = h[n - 3]
    return solve_banded((2, 2), ab, R).T


def interp_runs(t, Y, Tq):
    """The rows of Y (rows x delays, no empty values) at the delays Tq.

    Tq : one row of delays for all rows, or one row per row of Y.
    Between two delays the values are interpolated linearly. Inside every
    stretch of equally spaced delays (delay_runs) they come from a cubic
    spline through the delays of that stretch instead: a spline keeps
    oscillations and sharp features that linear interpolation flattens,
    but it must not run across a jump of the delay step, where it would
    overshoot. Empty (NaN) outside the measured delays.
    Port of interpRuns."""
    t = np.asarray(t, dtype=float).ravel()
    Y = np.atleast_2d(np.asarray(Y, dtype=float))
    n = t.size
    nr = Y.shape[0]
    Tq = np.asarray(Tq, dtype=float)
    if Tq.ndim == 1:
        Tq = np.broadcast_to(Tq, (nr, Tq.size))
    Yq = np.full(Tq.shape, np.nan)
    if n < 2 or nr == 0:
        return Yq
    inside = (Tq >= t[0]) & (Tq <= t[-1])
    if not inside.any():
        return Yq
    ri = np.nonzero(inside)[0]
    q = Tq[inside]
    # interval of every query: t[k] <= q <= t[k + 1]
    k = np.clip(np.searchsorted(t, q, side='right') - 1, 0, n - 2)
    h = np.diff(t)
    Yf = np.ascontiguousarray(Y).ravel()
    i0 = ri * n + k                               # Y[ri, k] is Yf[i0]
    b = (q - t[k]) / h[k]
    y = Yf[i0] + b * (Yf[i0 + 1] - Yf[i0])
    for a, e in delay_runs(t):
        s = np.flatnonzero((q >= t[a]) & (q <= t[e]))
        if s.size == 0:
            continue
        M = spline_curv(t[a:e + 1], Y[:, a:e + 1])
        Mf = np.ascontiguousarray(M).ravel()
        w = e - a + 1
        ks = np.clip(k[s], a, e - 1)
        rs = ri[s]
        hs = h[ks]
        bb = (q[s] - t[ks]) / hs
        aa = 1.0 - bb
        j0 = rs * n + ks
        m0 = rs * w + (ks - a)
        y[s] = (aa * Yf[j0] + bb * Yf[j0 + 1]
                + ((aa ** 3 - aa) * Mf[m0]
                   + (bb ** 3 - bb) * Mf[m0 + 1]) * hs ** 2 / 6.0)
    Yq[inside] = y
    return Yq


def chirp_shift(t, X, t0):
    """The chirp correction of the rows of X (wavelength x delay):
    row i at the delays t + t0[i], by interp_runs. Empty (NaN) where
    t + t0 lies outside the measured delays, where t0 is not a number and,
    for a row with empty values, where a delay next to t + t0 is empty.
    Port of chirpShift."""
    t = np.asarray(t, dtype=float).ravel()
    X = np.asarray(X, dtype=float)
    t0 = np.asarray(t0, dtype=float).ravel()
    Y = np.full(X.shape, np.nan)
    okT = np.isfinite(t0)
    full = okT & np.all(np.isfinite(X), axis=1)
    if full.any():
        Y[full, :] = interp_runs(t, X[full, :], t[None, :] + t0[full, None])
    for i in np.flatnonzero(okT & ~full):
        row = X[i, :]
        m = np.isfinite(row)
        if np.count_nonzero(m) < 2:
            continue
        tq = t + t0[i]
        lin = np.interp(tq, t, row, left=np.nan, right=np.nan)
        sp = interp_runs(t[m], row[m], tq)[0]
        use = np.isfinite(lin)
        lin[use] = sp[use]
        Y[i, :] = lin
    return Y
