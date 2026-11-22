"""Multivariate curve resolution by alternating least squares (MCR-ALS).

Port of computeMCR of TRSpecAnalyzer.m. D = C S' with constraints applied
after every half step.
"""
import time
import numpy as np
from ._util import Result

_EPS = np.finfo(float).eps


def _solve_s(d, c):
    try:
        s = np.linalg.lstsq(c, d, rcond=None)[0].T
    except Exception:
        s = (np.linalg.pinv(c) @ d).T
    s[~np.isfinite(s)] = 0
    return s


def _solve_c(d, s):
    try:
        c = np.linalg.lstsq(s, d.T, rcond=None)[0].T
    except Exception:
        c = d @ np.linalg.pinv(s.T)
    c[~np.isfinite(c)] = 0
    return c


def _init_svd(d0, n):
    u0 = np.linalg.svd(d0, full_matrices=False)[0][:, :n].copy()
    for jj in range(u0.shape[1]):
        idx = int(np.argmax(np.abs(u0[:, jj])))
        if u0[idx, jj] < 0:
            u0[:, jj] = -u0[:, jj]
    c0 = np.abs(u0)
    if np.any(~np.isfinite(c0)):
        c0 = np.abs(np.random.default_rng().standard_normal((d0.shape[0], n)))
    return c0


def _start_spectra(d0, s_in, n):
    """Starting spectra from given ones (linearly dependent or empty columns
    are dropped); missing components are filled from the SVD of what the
    given spectra do not describe."""
    if s_in is None or np.size(s_in) == 0:
        return None, 0
    s_in = np.array(s_in, dtype=float)
    if s_in.ndim == 1:
        s_in = s_in[:, None]
    if s_in.shape[0] != d0.shape[1]:
        return None, 0
    s_in[~np.isfinite(s_in)] = 0
    s_in = s_in[:, :min(n, s_in.shape[1])]
    keep = np.zeros(s_in.shape[1], dtype=bool)
    Q = np.zeros((s_in.shape[0], 0))
    for jj in range(s_in.shape[1]):
        v = s_in[:, jj]
        nv = np.linalg.norm(v)
        if nv == 0:
            continue
        r = v - Q @ (Q.T @ v)
        r = r - Q @ (Q.T @ r)
        if np.linalg.norm(r) > 1e-6 * nv:
            keep[jj] = True
            Q = np.column_stack([Q, r / np.linalg.norm(r)])
    n_used = int(np.count_nonzero(keep))
    if n_used == 0:
        return None, 0
    s0 = s_in[:, keep]
    if n_used < n:
        dr = d0 - (d0 @ Q) @ Q.T
        v0 = np.linalg.svd(dr, full_matrices=False)[2][:n - n_used, :].T.copy()
        if v0.shape[1] < n - n_used or np.any(~np.isfinite(v0)):
            return None, 0
        for jj in range(v0.shape[1]):
            p = dr @ v0[:, jj]
            ip = int(np.argmax(np.abs(p)))
            if p[ip] < 0:
                v0[:, jj] = -v0[:, jj]
        s0 = np.column_stack([s0, v0 * np.mean(np.sqrt(np.sum(s0 ** 2, axis=0)))])
    return s0, n_used


def _unimodal(c):
    for jj in range(c.shape[1]):
        col = c[:, jj].copy()
        ipk = int(np.argmax(col))
        for kk in range(ipk - 1, -1, -1):
            if col[kk] > col[kk + 1]:
                col[kk] = col[kk + 1]
        for kk in range(ipk + 1, col.size):
            if col[kk] > col[kk - 1]:
                col[kk] = col[kk - 1]
        c[:, jj] = col
    return c


def compute_mcr(D, n, nn_c=True, uni_c=False, nn_s=False, tol=1e-4, max_iter=200,
                init='svd', c_init=None, s_init=None, seed=None):
    """MCR-ALS of D (n_wavelength x n_delay) with n components.

    nn_c, uni_c : non-negative / unimodal concentration profiles
    nn_s        : non-negative spectra
    init        : 'svd', 'random' or 'custom' (c_init: n_delay x n, or
                  s_init: n_wavelength x up to n starting spectra, e.g. the
                  EADS of a global fit)
    seed        : seed of the random start (numpy generator; the numbers
                  differ from the random start of the MATLAB program)

    Returns a Result with C (n_delay x n), S (n_wavelength x n), lof (lack
    of fit in %), iters, lof_history, converged, seed, n_s_init.
    """
    D = np.asarray(D, dtype=float).T          # n_delay x n_wavelength, as in the program
    n = int(n)
    n_t, _ = D.shape
    valid = np.isfinite(D)
    d_clean = np.where(valid, D, 0.0)
    S = None
    seed_used = np.nan
    n_s_init = np.nan
    mode = str(init).lower()
    if mode == 'custom':
        if c_init is not None and np.ndim(c_init) == 2 and np.shape(c_init)[1] == n:
            C = np.array(c_init, dtype=float)
        else:
            S, n_used = _start_spectra(d_clean, s_init, n)
            if S is None:
                C = _init_svd(d_clean, n)
            else:
                C = _solve_c(d_clean, S)
                n_s_init = n_used
    elif mode == 'random':
        seed_used = seed
        if (seed_used is None or not np.isscalar(seed_used) or not np.isfinite(seed_used)
                or seed_used < 1):
            seed_used = 1 + int(time.time() * 1000) % 999999
        seed_used = int(min(np.floor(float(seed_used)), 2 ** 32 - 1))
        C = np.abs(np.random.default_rng(seed_used).standard_normal((n_t, n)))
    else:
        C = _init_svd(d_clean, n)
    if S is None:
        S = _solve_s(d_clean, C)

    lof_history = []
    prev_lof = np.inf
    converged = False
    lof = np.nan
    sum_d = np.sum(D[valid] ** 2)
    it = 0
    for it in range(1, int(max_iter) + 1):
        S = _solve_s(d_clean, C)
        if nn_s:
            S[S < 0] = 0
        C = _solve_c(d_clean, S)
        if nn_c:
            C[C < 0] = 0
        if uni_c:
            C = _unimodal(C)
        R = D - C @ S.T
        with np.errstate(divide='ignore', invalid='ignore'):
            lof = 100 * np.sqrt(np.sum(R[valid] ** 2) / sum_d)
        lof_history.append(lof)
        with np.errstate(invalid='ignore'):       # (inf / inf in the first round)
            rel = abs(prev_lof - lof) / max(prev_lof, _EPS)
        if rel < tol:
            converged = True
            break
        prev_lof = lof
    return Result(C=C, S=S, lof=float(lof), iters=it, lof_history=np.asarray(lof_history),
                  converged=converged, seed=seed_used, n_s_init=n_s_init)
