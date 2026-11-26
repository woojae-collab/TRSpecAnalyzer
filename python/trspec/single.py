"""Single-wavelength kinetic fit.

Port of fitSingleTrace of TRSpecAnalyzer.m.
"""
import numpy as np
from ._optim import fminsearch
from ._util import Result, lstsq_min_norm, as_row
from .kinetics import (exp_irf_conv, irf_delay_step, irf_too_wide,
                       irf_unresolved, irf_width_lost,
                       stretched_irf_conv)
from .uncertainty import fit_param_uncertainty

_EPS = np.finfo(float).eps


def fit_single_trace(t, y, tau_init, tau_fixed=False, beta_init=1.0, beta_fixed=True,
                     stretch_on=False, t0_init=0.0, t0_fixed=True, fwhm_init=0.15,
                     fwhm_fixed=True, has_inf=False, irf_mode='numerical', t_range=None,
                     fw_mode='auto'):
    """Fit one kinetic trace with a sum of (stretched) exponentials.

    Returns a Result with tau, beta, t0, fwhm, amp (amplitudes; the last one
    belongs to the non-decaying term if has_inf), fit, residual (NaN where
    not fitted), rms, fval, exit_flag, iters, mask, tau_std, beta_std,
    t0_std, fwhm_std, fwhm_lo, fwhm_at_limit, irf_retry (fw_mode: see
    fit_global_analysis).
    """
    fw_mode = fw_mode or 'auto'
    t0_fixed_arg, fwhm_fixed_arg = bool(t0_fixed), bool(fwhm_fixed)
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if t.size != y.size:
        raise ValueError('fit_single_trace: t and y must have the same length.')
    tau_init = np.atleast_1d(np.asarray(tau_init, dtype=float)).ravel()
    n_cmp = tau_init.size
    try:
        tau_fixed = as_row(tau_fixed, n_cmp, bool)
        beta_init = as_row(beta_init, n_cmp, float)
        beta_fixed = as_row(beta_fixed, n_cmp, bool)
        stretch_on = as_row(stretch_on, n_cmp, bool)
    except ValueError:
        raise ValueError('fit_single_trace: tau_fixed, beta_init, beta_fixed and stretch_on '
                         'need one entry per time constant (%d).' % n_cmp)
    no_irf = isinstance(irf_mode, str) and irf_mode.lower() == 'none'
    if no_irf:
        t0_fixed = True
        fwhm_fixed = True
    mask = np.isfinite(y)
    if t_range is not None and np.size(t_range) and np.all(np.isfinite(t_range)):
        mask = mask & (t >= np.min(t_range)) & (t <= np.max(t_range))
    n_nl = n_cmp
    n_lin = n_nl + int(bool(has_inf))
    if np.count_nonzero(mask) < n_lin + 1:
        raise ValueError('Not enough data points for the fit (need at least %d, have %d after masking).'
                         % (n_lin + 1, np.count_nonzero(mask)))
    free_tau = np.flatnonzero(~tau_fixed)
    free_beta = np.flatnonzero(stretch_on & ~beta_fixed)
    t0_free = not t0_fixed
    fwhm_free = not fwhm_fixed
    x0 = []
    x0.extend(np.log(np.maximum(tau_init[free_tau], _EPS)))
    x0.extend(np.log(np.maximum(beta_init[free_beta], 1e-3)))
    if t0_free:
        x0.append(float(t0_init))
    if fwhm_free:
        x0.append(np.log(max(fwhm_init, _EPS)))
    x0 = np.asarray(x0, dtype=float)
    t_span = t.max() - t.min()
    ym = y[mask]
    fw_res, fw_lo = 0.0, 0.0
    wall = {'on': False}
    if fwhm_free:
        dz_irf = irf_delay_step(t[mask], t0_init)
        if np.isfinite(dz_irf) and dz_irf > 0:
            fw_res = 0.5 * dz_irf
            fw_lo = min(fw_res, fwhm_init) * (1 - 1e-9)
            wall['on'] = fw_mode == 'limit'

    # (a search that runs off overflows here and there; such values are
    # rejected below, so the warnings of numpy are not shown)
    @np.errstate(over='ignore', invalid='ignore', divide='ignore')
    def unpack(x):
        tau_c = tau_init.copy()
        beta_c = beta_init.copy()
        t0_c = float(t0_init)
        fw_c = float(fwhm_init)
        idx = 0
        if free_tau.size:
            tau_c[free_tau] = np.exp(x[idx:idx + free_tau.size])
            idx += free_tau.size
        if free_beta.size:
            beta_c[free_beta] = np.exp(x[idx:idx + free_beta.size])
            idx += free_beta.size
        if t0_free:
            t0_c = float(x[idx])
            idx += 1
        if fwhm_free:
            fw_c = float(np.exp(x[idx]))
        fw_bad = (not no_irf) and (not np.isfinite(fw_c) or fw_c <= 0 or fw_c > t_span
                                   or (wall['on'] and fw_c < fw_lo))
        if (np.any(~np.isfinite(tau_c)) or np.any(tau_c <= 0) or np.any(~np.isfinite(beta_c))
                or np.any(beta_c <= 0) or np.any(beta_c > 2) or fw_bad):
            return tau_c, beta_c, t0_c, fw_c, np.zeros(n_lin), np.zeros(t.size)
        M = np.zeros((t.size, n_lin))
        if no_irf:
            tt = t - t0_c
            stp = (tt >= 0).astype(float)
            for jj in range(n_nl):
                if stretch_on[jj]:
                    M[:, jj] = stp * np.exp(-(np.maximum(tt, 0) / tau_c[jj]) ** beta_c[jj])
                else:
                    M[:, jj] = stp * np.exp(-np.maximum(tt, 0) / tau_c[jj])
            if has_inf:
                M[:, -1] = stp
        else:
            for jj in range(n_nl):
                if stretch_on[jj]:
                    M[:, jj] = stretched_irf_conv(t, tau_c[jj], beta_c[jj], t0_c, fw_c)
                else:
                    M[:, jj] = exp_irf_conv(t, tau_c[jj], t0_c, fw_c)
            if has_inf:
                M[:, -1] = exp_irf_conv(t, np.inf, t0_c, fw_c)
        Mm = M[mask, :]
        try:
            a_c = lstsq_min_norm(Mm, ym)
        except Exception:
            a_c = np.linalg.pinv(Mm) @ ym
        if np.any(~np.isfinite(a_c)):
            a_c = np.linalg.pinv(Mm) @ ym
        return tau_c, beta_c, t0_c, fw_c, a_c, M @ a_c

    def loss_fn(x):
        fv = unpack(np.asarray(x, dtype=float).ravel())[5]
        val = float(np.sum((y - fv)[mask] ** 2))
        return val if np.isfinite(val) else 1e30

    if x0.size:
        x_hat, fval, exit_flag, iters = fminsearch(loss_fn, x0, tol_x=1e-7, tol_fun=1e-8,
                                                   max_iter=1500, max_fun_evals=5000)
    else:
        x_hat, fval, exit_flag, iters = x0, loss_fn(x0), 1, 0
    wall['on'] = False          # the search is over
    tau_out, beta_out, t0_out, fwhm_out, a_out, fit_vec = unpack(x_hat)
    res_vec = y - fit_vec
    res_vec[~mask] = np.nan
    rms = float(np.sqrt(np.nanmean(res_vec[mask] ** 2)))

    tau_std = np.full(n_cmp, np.nan)
    beta_std = np.full(n_cmp, np.nan)
    t0_std = np.nan
    fwhm_std = np.nan
    if x0.size:
        p_opt = list(tau_out[free_tau]) + list(beta_out[free_beta])
        if t0_free:
            p_opt.append(t0_out)
        if fwhm_free:
            p_opt.append(fwhm_out)
        p_opt = np.asarray(p_opt, dtype=float)

        def residual_at_phys(p):
            xx = np.zeros(p.size)
            kk = 0
            for _ in range(free_tau.size + free_beta.size):
                xx[kk] = np.log(max(p[kk], 1e-12))
                kk += 1
            if t0_free:
                xx[kk] = p[kk]
                kk += 1
            if fwhm_free:
                xx[kk] = np.log(max(p[kk], 1e-12))
            return (y - unpack(xx)[5])[mask]

        sig = fit_param_uncertainty(residual_at_phys, p_opt, fval,
                                    np.count_nonzero(mask) - a_out.size)
        i2 = 0
        for jj in range(free_tau.size):
            tau_std[free_tau[jj]] = sig[i2]
            i2 += 1
        for jj in range(free_beta.size):
            beta_std[free_beta[jj]] = sig[i2]
            i2 += 1
        if t0_free:
            t0_std = sig[i2]
            i2 += 1
        if fwhm_free:
            fwhm_std = sig[i2]

    res = Result(tau=tau_out, beta=beta_out, t0=t0_out, fwhm=fwhm_out, amp=a_out,
                 fit=fit_vec, residual=res_vec, rms=rms, fval=fval,
                 exit_flag=exit_flag, iters=iters, mask=mask, irf_mode=irf_mode,
                 tau_std=tau_std, beta_std=beta_std, t0_std=t0_std, fwhm_std=fwhm_std,
                 fwhm_lo=np.nan, fwhm_at_limit=False, irf_retry=0,
                 fwhm_hi=np.nan, fwhm_too_wide=False)
    if fwhm_free and fw_res > 0:
        res.fwhm_lo = fw_res
        res.fwhm_at_limit = irf_unresolved(fwhm_out, fw_res)
    if fwhm_free:
        res.fwhm_hi = 0.5 * t_span
        res.fwhm_too_wide = irf_too_wide(fwhm_out, t_span)
    if fwhm_free and fw_mode == 'auto' and irf_width_lost(fwhm_out, fwhm_std, fw_res, t_span):
        # the direct fit lost the IRF width: IRF fixed first, then all free
        fw_s0 = float(fwhm_init)
        if irf_unresolved(fwhm_init, fw_res):
            fw_s0 = 4 * fw_res
        try:
            s1 = fit_single_trace(t, y, tau_init, tau_fixed, beta_init, beta_fixed,
                                  stretch_on, t0_init, True, fw_s0, True, has_inf,
                                  irf_mode, t_range, 'plain')
            s2 = fit_single_trace(t, y, s1.tau, tau_fixed, s1.beta, beta_fixed,
                                  stretch_on, t0_init, t0_fixed_arg, fw_s0,
                                  fwhm_fixed_arg, has_inf, irf_mode, t_range, 'limit')
            it_all = res.iters + s1.iters + s2.iters
            if (s2.fval < res.fval * (1 - 1e-6)
                    or ((res.fwhm_at_limit or res.fwhm_too_wide)
                        and s2.fval <= res.fval * (1 + 1e-6))):
                res = s2
                res.irf_retry = 1
            else:
                res.irf_retry = 2
            res.iters = it_all
        except Exception:
            pass
    return res
