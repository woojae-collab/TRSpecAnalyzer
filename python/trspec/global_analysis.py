"""Global analysis: sum of (stretched) exponentials shared by all wavelengths.

Port of fitGlobalAnalysis of TRSpecAnalyzer.m. Variable projection: the
time constants (and optionally beta, t0, FWHM) are optimized by Nelder-Mead,
the amplitude spectra (DADS) are solved by linear least squares at each step.
"""
import numpy as np
from ._optim import fminsearch
from ._util import Result, lstsq_min_norm, as_row
from .kinetics import (build_ga_basis_ext, ca_basis, ca_options,
                       irf_delay_step, irf_too_wide, irf_unresolved,
                       irf_width_lost)
from .uncertainty import fit_param_uncertainty

_PENALTY = 1e30


def fit_global_analysis(D, t, tau_init, t0_init=0.0, fwhm_init=0.15,
                        tau_fixed=False, t0_fixed=True, fwhm_fixed=True,
                        has_inf=False, beta_init=None, stretch_on=None,
                        irf_mode='numerical', t_range=None, ca=None,
                        fw_mode='auto'):
    """Fit D(lambda, t) = sum_j DADS_j(lambda) * c_j(t).

    D         : data, shape (n_wavelength, n_delay); NaN = not used
    t         : delays (n_delay)
    tau_init  : starting time constants;  tau_fixed : bool per component
    t0_init, fwhm_init, t0_fixed, fwhm_fixed : Gaussian IRF
    has_inf   : add a non-decaying component
    beta_init, stretch_on : stretched-exponential components
    irf_mode  : 'numerical' (Gaussian IRF) or 'none' (no IRF; t0 and FWHM
                are then not fitted)
    t_range   : optional (tmin, tmax) fit window
    ca        : coherent-artifact terms at time zero (see ca_options)
    fw_mode   : how a fitted IRF FWHM is handled:
                'auto'  plain fit; if the width is lost (unresolved, or
                        its uncertainty exceeds it) the fit is repeated in
                        two steps (IRF fixed, then everything free with a
                        lower wall) and the better one is kept
                'plain' plain fit only;  'limit' fit with the wall

    Returns a Result with tau, beta, t0, fwhm, dads (n_wavelength x n_comp),
    fit (same shape as D), rss, rms, iters, initial_loss, initial_rms,
    n_used, tau_std, beta_std, t0_std, fwhm_std, ca_order, ca_spec,
    ca_fwhm, ca_fwhm_std, ca_mode, fwhm_lo (half the delay step at time
    zero, NaN if FWHM was fixed), fwhm_at_limit, irf_retry (0 none,
    1 repeat taken, 2 repeat tried without gain).
    """
    D = np.asarray(D, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    tau_init = np.atleast_1d(np.asarray(tau_init, dtype=float)).ravel()
    n_cmp = tau_init.size
    if beta_init is None or np.size(beta_init) == 0:
        beta_init = np.ones(n_cmp)
    if stretch_on is None or np.size(stretch_on) == 0:
        stretch_on = np.zeros(n_cmp, dtype=bool)
    if irf_mode is None or irf_mode == '':
        irf_mode = 'numerical'
    try:
        stretch_on = as_row(stretch_on, n_cmp, bool)
        beta_init = as_row(beta_init, n_cmp, float)
        tau_fixed = as_row(tau_fixed, n_cmp, bool)
    except ValueError:
        raise ValueError('fit_global_analysis: tau_fixed, beta_init and stretch_on need '
                         'one entry per time constant (%d).' % n_cmp)
    no_irf = isinstance(irf_mode, str) and irf_mode.lower() == 'none'
    if not fw_mode:
        fw_mode = 'auto'
    ca = ca_options(ca)
    if no_irf:
        ca['order'] = 0
    ca_k = ca['order']
    ca_own = ca_k > 0 and ca['mode'] != 'irf'
    ca_free = ca_k > 0 and ca['mode'] == 'free'
    ca_w = ca['fwhm']
    if ca_own and not (ca_w > 0):
        ca_w = float(fwhm_init)
    tau_init_arg, beta_init_arg = tau_init.copy(), beta_init.copy()
    t0_fixed_arg, fwhm_fixed_arg = bool(t0_fixed), bool(fwhm_fixed)

    if t_range is None or np.size(t_range) == 0 or not np.all(np.isfinite(t_range)):
        t_sel = np.ones(t.size, dtype=bool)
    else:
        t_sel = (t >= np.min(t_range)) & (t <= np.max(t_range))
        if not t_sel.any():
            t_sel = np.ones(t.size, dtype=bool)

    state = {'tau': tau_init.copy(), 'beta': beta_init.copy(),
             't0': float(t0_init), 'fwhm': float(fwhm_init), 'n_loss': 0,
             'ca_w': ca_w, 'fw_wall': False}
    if no_irf:
        t0_fixed = True
        fwhm_fixed = True
    idx_free_tau = np.flatnonzero(~tau_fixed)
    idx_free_beta = np.flatnonzero(stretch_on)
    n_free_tau = idx_free_tau.size
    n_free_beta = idx_free_beta.size

    x0 = []
    x0.extend(np.log(tau_init[idx_free_tau]))
    x0.extend(np.log(np.maximum(beta_init[idx_free_beta], 1e-3)))
    if not t0_fixed:
        x0.append(float(t0_init))
    if not fwhm_fixed:
        x0.append(np.log(fwhm_init))
    if ca_free:
        x0.append(np.log(ca_w))
    x0 = np.asarray(x0, dtype=float)

    fit_mask = t_sel
    if not fit_mask.any():
        raise ValueError('No delay points left to fit (check the fit t range and t0).')
    fw_res, fw_lo = 0.0, 0.0
    if not fwhm_fixed:
        dz_irf = irf_delay_step(t[fit_mask], t0_init)
        if np.isfinite(dz_irf) and dz_irf > 0:
            fw_res = 0.5 * dz_irf
            fw_lo = min(fw_res, fwhm_init) * (1 - 1e-9)
            state['fw_wall'] = fw_mode == 'limit'
    if ca_k > 0:
        w0 = ca_w if ca_own else fwhm_init
        n_near = int(np.count_nonzero(fit_mask & (np.abs(t - t0_init) <= 1.5 * w0)))
        if n_near < ca_k + 2:
            raise ValueError('Coherent artifact: %d fitted delay(s) lie within 1.5 FWHM of '
                             'time zero, %d are needed for %d term(s). Use fewer terms or '
                             'switch them off.' % (n_near, ca_k + 2, ca_k))
    m_rows, n_cols = D.shape
    all_finite = bool(np.all(np.isfinite(D)))
    t_span = t.max() - t.min()

    def objective(x):
        x = np.asarray(x, dtype=float).ravel()
        idx = 0
        for j in range(n_free_tau):
            state['tau'][idx_free_tau[j]] = np.exp(x[idx])
            idx += 1
        for j in range(n_free_beta):
            state['beta'][idx_free_beta[j]] = np.exp(x[idx])
            idx += 1
        if not t0_fixed:
            state['t0'] = x[idx]
            idx += 1
        if not fwhm_fixed:
            state['fwhm'] = np.exp(x[idx])
            idx += 1
        if ca_free:
            state['ca_w'] = np.exp(x[idx])
        tau_cur, beta_cur = state['tau'], state['beta']
        t0_cur, fwhm_cur = state['t0'], state['fwhm']
        ca_w_cur = state['ca_w']
        sigma_cur = 0.0 if no_irf else fwhm_cur / (2 * np.sqrt(2 * np.log(2)))
        k = n_cmp + int(bool(has_inf)) + ca_k
        if (np.any(tau_cur < sigma_cur / 4) or np.any(tau_cur > 10 * t_span)
                or np.any(tau_cur <= 0)
                or (not no_irf and (fwhm_cur <= 0 or fwhm_cur > t_span))
                or (state['fw_wall'] and fwhm_cur < fw_lo)
                or (ca_free and (ca_w_cur < fwhm_cur / 5 or ca_w_cur > 5 * fwhm_cur))
                or np.any(beta_cur[stretch_on] <= 0) or np.any(beta_cur[stretch_on] > 2)):
            return _PENALTY, np.zeros((m_rows, k)), np.zeros(D.shape)
        C = build_ga_basis_ext(t, tau_cur, beta_cur, stretch_on, t0_cur, fwhm_cur, has_inf, irf_mode)
        if ca_k > 0:
            C = np.column_stack([C, ca_basis(t, t0_cur, ca_w_cur if ca_own else fwhm_cur, ca_k)])
        col_norms = np.sqrt(np.nansum(C ** 2, axis=0))
        if np.any(col_norms < 1e-12) or np.any(~np.isfinite(C)):
            return _PENALTY, np.zeros((m_rows, C.shape[1])), np.zeros(D.shape)
        k = C.shape[1]
        As = np.zeros((m_rows, k))
        fit_m = np.zeros((m_rows, n_cols))
        if all_finite:
            if fit_mask.all():
                At = lstsq_min_norm(C, D.T)
            else:
                At = lstsq_min_norm(C[fit_mask, :], D[:, fit_mask].T)
            if np.any(~np.isfinite(At)) or np.max(np.abs(At)) > 1e10:
                return _PENALTY, As, fit_m
            As = At.T
            fit_m = As @ C.T
        else:
            for ii in range(m_rows):
                row = D[ii, :]
                mm = np.isfinite(row) & fit_mask
                if np.count_nonzero(mm) <= k:
                    continue
                ai = lstsq_min_norm(C[mm, :], row[mm])
                if np.any(~np.isfinite(ai)) or np.max(np.abs(ai)) > 1e8:
                    continue
                As[ii, :] = ai
                fit_m[ii, :] = ai @ C.T
        R = D - fit_m
        mfin = np.isfinite(R) & fit_mask[None, :]
        state['n_loss'] = int(np.count_nonzero(mfin))
        if not mfin.any():
            return _PENALTY, As, fit_m
        return float(np.sum(R[mfin] ** 2)), As, fit_m

    x_opt = x0
    if x0.size == 0:
        loss, A, fit_mat = objective(x0)
        iters = 0
        init_loss = loss
        n_init = n_opt = state['n_loss']
    else:
        init_loss, _, _ = objective(x0)
        n_init = state['n_loss']
        x_opt, loss, _, iters = fminsearch(lambda x: objective(x)[0], x0, tol_x=1e-8,
                                           tol_fun=1e-10, max_iter=5000, max_fun_evals=20000)
        _, A, fit_mat = objective(x_opt)
        n_opt = state['n_loss']

    tau_out = state['tau'].copy()
    beta_out = state['beta'].copy()
    t0_out = float(state['t0'])
    fwhm_out = float(state['fwhm'])
    state['fw_wall'] = False        # the search is over

    tau_std = np.full(n_cmp, np.nan)
    beta_std = np.full(n_cmp, np.nan)
    t0_std = np.nan
    fwhm_std = np.nan
    ca_std = np.nan
    if x0.size:
        p_opt = []
        p_opt.extend(tau_out[idx_free_tau])
        p_opt.extend(beta_out[idx_free_beta])
        if not t0_fixed:
            p_opt.append(t0_out)
        if not fwhm_fixed:
            p_opt.append(fwhm_out)
        if ca_free:
            p_opt.append(float(state['ca_w']))
        p_opt = np.asarray(p_opt, dtype=float)
        n_lin = A.shape[1] * int(np.count_nonzero(np.any(A != 0, axis=1)))
        n_used = n_opt - n_lin
        res_mask = np.isfinite(D) & fit_mask[None, :]

        def residual_at_phys(p):
            xx = np.zeros(p.size)
            i2 = 0
            for _ in range(n_free_tau + n_free_beta):
                xx[i2] = np.log(max(p[i2], 1e-12))
                i2 += 1
            if not t0_fixed:
                xx[i2] = p[i2]
                i2 += 1
            if not fwhm_fixed:
                xx[i2] = np.log(max(p[i2], 1e-12))
                i2 += 1
            if ca_free:
                xx[i2] = np.log(max(p[i2], 1e-12))
            _, _, fm = objective(xx)
            r = (D - fm)[res_mask]
            r[~np.isfinite(r)] = 0
            return r

        sig = fit_param_uncertainty(residual_at_phys, p_opt, loss, n_used)
        idx = 0
        for j in range(n_free_tau):
            tau_std[idx_free_tau[j]] = sig[idx]
            idx += 1
        for j in range(n_free_beta):
            beta_std[idx_free_beta[j]] = sig[idx]
            idx += 1
        if not t0_fixed:
            t0_std = sig[idx]
            idx += 1
        if not fwhm_fixed:
            fwhm_std = sig[idx]
            idx += 1
        if ca_free:
            ca_std = sig[idx]
        # the Jacobian probing changed the internal state: restore the optimum
        _, A, fit_mat = objective(x_opt)

    k_kin = n_cmp + int(bool(has_inf))
    ca_spec = A[:, k_kin:].copy()
    A = A[:, :k_kin]
    ca_w_out = np.nan
    if ca_k > 0:
        ca_w_out = float(state['ca_w']) if ca_own else fwhm_out

    res = Result(tau=tau_out, beta=beta_out, t0=t0_out, fwhm=fwhm_out,
                 dads=A, fit=fit_mat, rss=loss, iters=iters,
                 rms=float(np.sqrt(loss / max(1, n_opt))),
                 initial_loss=init_loss,
                 initial_rms=float(np.sqrt(init_loss / max(1, n_init))),
                 n_used=n_opt, tau_std=tau_std, beta_std=beta_std,
                 t0_std=t0_std, fwhm_std=fwhm_std, has_inf=bool(has_inf),
                 stretch_on=stretch_on.copy(), irf_mode=irf_mode,
                 ca_order=ca_k, ca_spec=ca_spec, ca_fwhm=ca_w_out,
                 ca_fwhm_std=ca_std, ca_mode=ca['mode'],
                 fwhm_lo=np.nan, fwhm_at_limit=False, irf_retry=0,
                 fwhm_hi=np.nan, fwhm_too_wide=False)
    if not fwhm_fixed and fw_res > 0:
        res.fwhm_lo = fw_res
        res.fwhm_at_limit = irf_unresolved(fwhm_out, fw_res)
    if not fwhm_fixed:
        res.fwhm_hi = 0.5 * t_span
        res.fwhm_too_wide = irf_too_wide(fwhm_out, t_span)
    if (not fwhm_fixed and fw_mode == 'auto' and x0.size
            and irf_width_lost(fwhm_out, fwhm_std, fw_res, t_span)):
        # the direct fit lost the IRF width: IRF fixed first, then all free
        fw_s0 = float(fwhm_init)
        if irf_unresolved(fwhm_init, fw_res):
            fw_s0 = 4 * fw_res
        try:
            s1 = fit_global_analysis(D, t, tau_init_arg, t0_init, fw_s0, tau_fixed,
                                     True, True, has_inf, beta_init_arg, stretch_on,
                                     irf_mode, t_range, ca, 'plain')
            s2 = fit_global_analysis(D, t, s1.tau, t0_init, fw_s0, tau_fixed,
                                     t0_fixed_arg, fwhm_fixed_arg, has_inf, s1.beta,
                                     stretch_on, irf_mode, t_range, ca, 'limit')
            it_all = res.iters + s1.iters + s2.iters
            if (s2.rss < res.rss * (1 - 1e-6)
                    or ((res.fwhm_at_limit or res.fwhm_too_wide)
                        and s2.rss <= res.rss * (1 + 1e-6))):
                s2.initial_loss = res.initial_loss
                s2.initial_rms = res.initial_rms
                res = s2
                res.irf_retry = 1
            else:
                res.irf_retry = 2
            res.iters = it_all
        except Exception:
            pass
    return res
