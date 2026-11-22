"""Sequential model: EADS from DADS, species populations and SAS.

Port of computeEADSFromDADS, seqClockConc, computeSASsequential and
suggestNComponents of TRSpecAnalyzer.m.
"""
import numpy as np
from scipy.special import erf
from ._util import lstsq_min_norm, colon as _colon
from .kinetics import exp_irf_conv, stretched_irf_conv

_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def compute_eads_from_dads(dads, tau, has_inf=False):
    """Convert DADS (parallel decays) to EADS of the sequential model
    1 -> 2 -> ... -> n (-> non-decaying species if has_inf).

    Returns (eads, tau_sorted, b_matrix); columns ordered by increasing tau.
    """
    dads = np.asarray(dads, dtype=float)
    tau = np.atleast_1d(np.asarray(tau, dtype=float)).ravel()
    n_decay = tau.size
    sort_idx = np.argsort(tau, kind='stable')
    tau_sorted = tau[sort_idx]
    if has_inf:
        k_vec = np.concatenate([1.0 / tau_sorted, [0.0]])
        perm = np.concatenate([sort_idx, [n_decay]])
    else:
        k_vec = 1.0 / tau_sorted
        perm = sort_idx
    n = k_vec.size
    b = np.zeros((n, n))
    for ii in range(n):
        prod_k = 1.0 if ii == 0 else np.prod(k_vec[:ii])
        for jj in range(ii + 1):
            denom = 1.0
            for mm in range(ii + 1):
                if mm != jj:
                    denom *= (k_vec[mm] - k_vec[jj])
            with np.errstate(divide='ignore', invalid='ignore'):
                b[ii, jj] = prod_k / denom
    dads_perm = dads[:, perm]
    try:
        eads = np.linalg.solve(b.T, dads_perm.T).T
    except Exception:
        eads = np.full(dads_perm.shape, np.nan)
    if np.any(~np.isfinite(eads)):
        eads = dads_perm @ np.linalg.pinv(np.nan_to_num(b, nan=0.0, posinf=0.0, neginf=0.0))
    return eads, tau_sorted, b




def seq_clock_conc(t, tau_s, beta_s, str_s, has_inf, t0, fwhm, use_irf):
    """Populations of the sequential model in which the clock of a
    stretched step starts at t0 ('t0' model). tau_s ascending."""
    t = np.asarray(t, dtype=float).ravel()
    tau_s = np.asarray(tau_s, dtype=float).ravel()
    beta_s = np.asarray(beta_s, dtype=float).ravel()
    str_s = np.asarray(str_s, dtype=bool).ravel()
    n_t = t.size
    n_decay = tau_s.size
    n_sp = n_decay + int(bool(has_inf))
    conc = np.zeros((n_t, n_sp))
    xd = t - t0
    if use_irf:
        sigma = fwhm * _FWHM_TO_SIGMA
        x_end = xd.max() + 6 * sigma
    else:
        sigma = 0.0
        x_end = xd.max()
    if not (x_end > 0):
        return conc
    growth = 0.01
    x_lo = 1e-40 * tau_s[0]
    n_geo = int(np.ceil(np.log(x_end / x_lo) / np.log(1 + growth)))
    xg = x_lo * (1 + growth) ** np.arange(n_geo + 1)
    xg = np.unique(np.concatenate([[0.0], xg, xd[xd > 0]]))
    n_g = xg.size
    dx = np.diff(xg)
    xm = xg[:-1] + dx / 2
    H = np.zeros((n_g, n_sp))
    icum = np.ones(n_g)
    with np.errstate(over='ignore', invalid='ignore', divide='ignore', under='ignore'):
        for ii in range(n_decay):
            tau_i = tau_s[ii]
            beta_i = beta_s[ii]
            if ii == 0:
                hi = np.exp(-(xg / tau_i) ** beta_i) if str_s[0] else np.exp(-xg / tau_i)
            else:
                d_i = np.diff(icum)
                if str_s[ii]:
                    g = (xg / tau_i) ** beta_i
                    dec = np.exp(-np.diff(g))
                    arr = d_i * np.exp(-(g[1:] - (xm / tau_i) ** beta_i))
                else:
                    z = dx / tau_i
                    dec = np.exp(-z)
                    phi = np.ones(z.shape)
                    big = z > 1e-8
                    phi[big] = -np.expm1(-z[big]) / z[big]
                    arr = d_i * phi
                hi = np.zeros(n_g)
                acc = 0.0
                for nn in range(1, n_g):
                    acc = acc * dec[nn - 1] + arr[nn - 1]
                    hi[nn] = acc
            hi[~np.isfinite(hi)] = 0
            H[:, ii] = hi
            icum = icum - hi
    if has_inf:
        H[:, n_sp - 1] = icum
    if not use_irf:
        pos = xd >= 0
        loc = np.searchsorted(xg, xd[pos])
        conc[pos, :] = H[loc, :]
        return conc
    if str_s[0]:
        conc[:, 0] = stretched_irf_conv(t, tau_s[0], beta_s[0], t0, fwhm)
    else:
        conc[:, 0] = exp_irf_conv(t, tau_s[0], t0, fwhm)
    a = xg[:-1]
    b = xg[1:]
    Ha = H[:-1, 1:n_decay]
    Hb = H[1:, 1:n_decay]
    g_n = 1.0 / (sigma * np.sqrt(2 * np.pi))
    s2 = sigma * np.sqrt(2)
    wide = dx >= 1e-3 * sigma
    for n in range(n_t):
        X = xd[n]
        if not (X + 6 * sigma > 0):
            continue
        k1 = int(np.searchsorted(b, X - 6 * sigma, side='right'))
        k2 = int(np.searchsorted(a, X + 6 * sigma, side='left')) - 1
        if k1 >= b.size or k2 < 0 or k2 < k1:
            continue
        kk = np.arange(k1, k2 + 1)
        va = X - a[kk]
        vb = X - b[kk]
        wA = np.zeros(kk.size)
        wB = np.zeros(kk.size)
        L = wide[kk]
        P = 0.5 * (erf(va[L] / s2) - erf(vb[L] / s2))
        Q = (va[L] * P + sigma ** 2 * g_n *
             (np.exp(-va[L] ** 2 / (2 * sigma ** 2)) - np.exp(-vb[L] ** 2 / (2 * sigma ** 2)))) / dx[kk[L]]
        wB[L] = Q
        wA[L] = P - Q
        N = ~L
        w = 0.5 * g_n * dx[kk[N]] * np.exp(-(X - xm[kk[N]]) ** 2 / (2 * sigma ** 2))
        wA[N] = w
        wB[N] = w
        conc[n, 1:n_decay] = wA @ Ha[kk, :] + wB @ Hb[kk, :]
    if has_inf:
        phi_inf = exp_irf_conv(t, np.inf, t0, fwhm)
        conc[:, n_sp - 1] = phi_inf - conc[:, :n_decay].sum(axis=1)
    conc[~np.isfinite(conc)] = 0
    return conc


def compute_sas_sequential(D, t, tau, beta=None, stretch_on=None, has_inf=False,
                           t0=0.0, fwhm=0.0, irf_mode='numerical', t_range=None,
                           seq_model='age'):
    """Species-associated spectra and populations of the sequential model
    1 -> 2 -> ... -> n, also with stretched steps.

    seq_model : 'age' (the clock of a stretched step starts when the
                population enters that state) or 't0' (it starts at t0)

    Returns (sas, species_conc, tau_sorted); sas is NaN where it cannot be
    determined.
    """
    D = np.asarray(D, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    tau = np.atleast_1d(np.asarray(tau, dtype=float)).ravel()
    n_t = t.size
    n_decay = tau.size
    beta = np.ones(n_decay) if beta is None else np.atleast_1d(np.asarray(beta, dtype=float)).ravel()
    stretch_on = (np.zeros(n_decay, dtype=bool) if stretch_on is None
                  else np.atleast_1d(np.asarray(stretch_on, dtype=bool)).ravel())
    sort_idx = np.argsort(tau, kind='stable')
    tau_sorted = tau[sort_idx]
    beta_sorted = beta[sort_idx]
    stretch_sorted = stretch_on[sort_idx]
    n_species = n_decay + int(bool(has_inf))
    if seq_model is None or seq_model == '':
        seq_model = 'age'
    sigma = fwhm * _FWHM_TO_SIGMA
    has_sig = np.isfinite(sigma) and sigma > 0
    no_irf = isinstance(irf_mode, str) and irf_mode.lower() == 'none'
    use_irf = bool(has_sig and not no_irf)
    clock_model = (str(seq_model).lower() == 't0' and n_decay > 1 and bool(np.any(stretch_sorted[1:])))

    sas = np.full((D.shape[0], n_species), np.nan)
    species_conc = np.zeros((n_t, n_species))
    if n_decay < 1 or np.any(~np.isfinite(tau_sorted)) or np.any(tau_sorted <= 0):
        return sas, species_conc, tau_sorted

    if clock_model:
        species_conc = seq_clock_conc(t, tau_sorted, beta_sorted, stretch_sorted,
                                      has_inf, t0, fwhm, use_irf)
    else:
        if use_irf:
            t_start = t0 - 6 * sigma
            h0 = min(0.02 * tau_sorted[0], 0.1 * sigma)
            h0 = max(h0, 6 * sigma / 1e4)
        else:
            t_start = t0
            h0 = 0.02 * tau_sorted[0]
        t_end = t.max()
        if not (t_end > t_start):
            return sas, species_conc, tau_sorted
        growth = 0.01
        x_end = t_end - t0
        x_uni_end = h0 / growth
        x_uni = _colon(t_start - t0, h0, min(x_uni_end, x_end))
        if x_end > x_uni_end:
            n_geo = int(np.ceil(np.log(x_end / x_uni_end) / np.log(1 + growth)))
            x_geo = x_uni_end * (1 + growth) ** np.arange(1, n_geo + 1)
        else:
            x_geo = np.zeros(0)
        in_range = t >= t_start
        t_all = np.unique(np.concatenate([t0 + x_uni, t0 + x_geo, t[in_range]]))
        t_all = t_all[t_all >= t_start]
        loc = np.searchsorted(t_all, t[in_range])
        n_g = t_all.size
        x = t_all - t0
        dt_g = np.diff(t_all)
        t_mid = t_all[:-1] + dt_g / 2
        icum = exp_irf_conv(t_all, np.inf, t0, fwhm) if use_irf else np.ones(n_g)
        c = np.zeros((n_g, n_species))
        with np.errstate(over='ignore', invalid='ignore', divide='ignore', under='ignore'):
            for ii in range(n_decay):
                tau_i = tau_sorted[ii]
                beta_i = beta_sorted[ii]
                if ii == 0:
                    if use_irf and stretch_sorted[0]:
                        ci = stretched_irf_conv(t_all, tau_i, beta_i, t0, fwhm)
                    elif use_irf:
                        ci = exp_irf_conv(t_all, tau_i, t0, fwhm)
                    elif stretch_sorted[0]:
                        ci = np.exp(-(np.maximum(x, 0) / tau_i) ** beta_i)
                    else:
                        ci = np.exp(-np.maximum(x, 0) / tau_i)
                    ci = np.array(ci, dtype=float)
                else:
                    d_i = np.diff(icum)
                    ci = np.zeros(n_g)
                    if stretch_sorted[ii]:
                        for nn in range(1, n_g):
                            age = t_all[nn] - t_mid[:nn]
                            ci[nn] = np.sum(d_i[:nn] * np.exp(-(age / tau_i) ** beta_i))
                    else:
                        z = dt_g / tau_i
                        dec = np.exp(-z)
                        phi = np.ones(z.shape)
                        big = z > 1e-8
                        phi[big] = -np.expm1(-z[big]) / z[big]
                        acc = 0.0
                        for nn in range(1, n_g):
                            acc = acc * dec[nn - 1] + d_i[nn - 1] * phi[nn - 1]
                            ci[nn] = acc
                ci[~np.isfinite(ci)] = 0
                c[:, ii] = ci
                icum = icum - ci
        if has_inf:
            c[:, n_species - 1] = icum
        species_conc[in_range, :] = c[loc, :]

    species_conc[~np.isfinite(species_conc)] = 0
    col_norms = np.sqrt(np.sum(species_conc ** 2, axis=0))
    keep = col_norms > 1e-12 * col_norms.max()
    n_k = int(np.count_nonzero(keep))
    t_sel = np.ones(n_t, dtype=bool)
    if t_range is not None and np.size(t_range) and np.all(np.isfinite(t_range)):
        t_sel = (t >= np.min(t_range)) & (t <= np.max(t_range))
        if not t_sel.any():
            t_sel = np.ones(n_t, dtype=bool)
    if n_k == 0 or np.count_nonzero(t_sel) <= n_k:
        return sas, species_conc, tau_sorted

    def solve_ls(cm, y):
        try:
            xs = lstsq_min_norm(cm, y)
        except Exception:
            xs = np.linalg.pinv(cm) @ y
        if np.any(~np.isfinite(xs)):
            xs = np.linalg.pinv(cm) @ y
        return xs

    c_use = species_conc[:, keep]
    d_sel = D[:, t_sel]
    s_kept = np.full((D.shape[0], n_k), np.nan)
    if np.all(np.isfinite(d_sel)):
        s_kept = solve_ls(c_use[t_sel, :], d_sel.T).T
    else:
        for ii in range(D.shape[0]):
            m = t_sel & np.isfinite(D[ii, :])
            if np.count_nonzero(m) <= n_k:
                continue
            s_kept[ii, :] = solve_ls(c_use[m, :], D[ii, m])
    sas[:, keep] = s_kept
    return sas, species_conc, tau_sorted


def suggest_n_components(D):
    """Number of singular values that stand above the noise floor
    (more than 5 times the median of the smaller half), between 1 and 8.

    Returns (n_suggested, singular_values, text).
    """
    D = np.array(D, dtype=float)
    if D.size == 0:
        return 1, np.zeros(0), ''
    D[~np.isfinite(D)] = 0
    sv = np.linalg.svd(D, compute_uv=False)
    n = sv.size
    if n < 4:
        return max(1, min(n, 2)), sv, ''
    noise = np.median(sv[int(np.ceil(n / 2)) - 1:])
    if not (noise > 0):
        noise = np.finfo(float).eps
    n_sug = int(np.count_nonzero(sv > 5 * noise))
    n_sug = max(1, min(n_sug, 8))
    k_show = min(8, n)
    parts = ['%d:%.3g' % (k + 1, sv[k] / sv[0]) for k in range(k_show)]
    txt = 'sigma_k/sigma_1 = %s  (noise ~ %.2g)' % ('  '.join(parts), noise / sv[0])
    return n_sug, sv, txt
