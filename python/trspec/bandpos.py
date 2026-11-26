"""Position of a band along the delays (spectral shift, dynamic Stokes
shift), and the fit of its time course.

Port of the band-position methods of TRSpecAnalyzer.m: lmFit,
lognormShape, bandAxis, bandHalfWidth, bandParabola, bandPosition,
bandModelCurve, bandShiftStart, fitBandShift.
"""
import numpy as np

from ._util import Result
from .uncertainty import fit_param_uncertainty

LN2 = float(np.log(2.0))
BAND_METHODS = ('moment', 'parabola', 'lognormal')


def lm_fit(res_fun, p0, max_iter=60):
    """Least squares by the Levenberg-Marquardt method: the parameters p
    that make sum(res_fun(p)^2) smallest, started at p0. The Jacobian is
    taken by central differences (step 1e-6 of the parameter, at least
    1e-9). Returns (p, rss, n_iter).

    A step is taken when it lowers the sum of squares; the damping is
    divided by 10 after a step taken and multiplied by 10 after one
    refused. The search ends when no step up to a damping of 1e12 lowers
    the sum of squares, when a step lowers it by less than 1e-12 of its
    value, or after max_iter steps."""
    p = np.array(p0, dtype=float).ravel()
    m = p.size
    r = np.asarray(res_fun(p), dtype=float).ravel()
    rss = float(r @ r)
    n_iter = 0
    if m == 0 or not np.isfinite(rss):
        return p, rss, n_iter
    lam = 1e-3
    for n_iter in range(1, max_iter + 1):
        jac = np.zeros((r.size, m))
        for j in range(m):
            h = 1e-6 * max(abs(p[j]), 1e-3)
            pp = p.copy()
            pp[j] = p[j] + h
            pm = p.copy()
            pm[j] = p[j] - h
            jac[:, j] = (np.asarray(res_fun(pp), dtype=float).ravel()
                         - np.asarray(res_fun(pm), dtype=float).ravel()) / (2 * h)
        if not np.all(np.isfinite(jac)):
            break
        a = jac.T @ jac
        g = jac.T @ r
        da = np.diag(a).copy()
        ridge = 1e-12 * float(np.sum(da)) / m + 1e-300
        taken = False
        rel = 0.0
        while lam <= 1e12:
            mat = a + lam * np.diag(da) + ridge * np.eye(m)
            try:
                dp = -np.linalg.solve(mat, g)
            except np.linalg.LinAlgError:
                lam = lam * 10
                continue
            pn = p + dp
            rn = np.asarray(res_fun(pn), dtype=float).ravel()
            rssn = float(rn @ rn)
            if np.isfinite(rssn) and rssn < rss:
                rel = (rss - rssn) / rss
                p, r, rss = pn, rn, rssn
                lam = max(lam / 10, 1e-10)
                taken = True
                break
            lam = lam * 10
        if not taken or rel < 1e-12:
            break
    return p, rss, n_iter


def lognorm_shape(v, vp, width, gam):
    """Log-normal band of height 1 with its maximum at vp:
    exp(-ln2 (ln(1 + a) / gam)^2), a = 2 gam (v - vp) / width, and 0
    where a <= -1. gam is the asymmetry (negative: tail to smaller v;
    0: a Gaussian of FWHM width; kept within +-0.95). The full width at
    half maximum is width sinh(gam) / gam."""
    v = np.asarray(v, dtype=float)
    g = min(max(float(gam), -0.95), 0.95)
    u = 2.0 * (v - vp) / width
    if abs(g) < 1e-9:
        return np.exp(-LN2 * u ** 2)
    a = g * u
    out = np.zeros(v.shape)
    ok = a > -1
    out[ok] = np.exp(-LN2 * (np.log1p(a[ok]) / g) ** 2)
    return out


def band_axis(x, spec_unit, pos_unit):
    """The axis along which the position is given, and the factor that
    turns a signal per unit of the spectral axis into one per unit of
    that axis. x: values of the spectral axis ('nm' or 'cm-1' =
    spec_unit). pos_unit the same: (x, 1). The other one: (1e7 / x,
    x^2 / 1e7), since |dx/dv| = x^2 / 1e7 for v = 1e7 / x."""
    x = np.asarray(x, dtype=float).ravel()
    if pos_unit == spec_unit:
        return x.copy(), np.ones(x.shape)
    with np.errstate(divide='ignore', invalid='ignore'):
        return 1e7 / x, x ** 2 / 1e7


def band_half_width(v, y, m, h):
    """Full width at half maximum of the curve y(v) around its point m,
    for the height h: the two places where y falls to h / 2, by linear
    interpolation between the points. NaN when the curve does not fall
    that far on one side, or does not reach h / 2 at the point m."""
    n = v.size
    half = 0.5 * h
    if not y[m] >= half:
        return np.nan
    i = m
    while i >= 0 and y[i] >= half:
        i -= 1
    if i < 0:
        return np.nan
    v_lo = v[i] + (half - y[i]) * (v[i + 1] - v[i]) / (y[i + 1] - y[i])
    k = m
    while k < n and y[k] >= half:
        k += 1
    if k >= n:
        return np.nan
    v_hi = v[k - 1] + (half - y[k - 1]) * (v[k] - v[k - 1]) / (y[k] - y[k - 1])
    return float(v_hi - v_lo)


def band_parabola(v, y, frac):
    """Maximum of the curve y(v) from a parabola through its top: the
    points around the largest one whose value is at least frac times
    that value (at least three points: the largest and its neighbors),
    fitted by least squares. Returns (pos, height, i1, i2, coef): the
    points used are i1 ... i2, the parabola is
    coef[0] u^2 + coef[1] u + coef[2] with u = (v - v[m]) / s (see the
    code). pos is NaN when the parabola has no maximum between the
    first and the last point used (a curvature of less than 1e-9 of the
    largest value counts as none: a flat top)."""
    n = v.size
    bad = (np.nan, np.nan, 0, -1, np.full(5, np.nan))
    if n < 3 or not np.all(np.isfinite(y)):
        return bad
    m = int(np.argmax(y))
    if not y[m] > 0:
        return bad
    lim = frac * y[m]
    i1 = m
    while i1 > 0 and y[i1 - 1] >= lim:
        i1 -= 1
    i2 = m
    while i2 < n - 1 and y[i2 + 1] >= lim:
        i2 += 1
    if i2 - i1 + 1 < 3:
        i1 = max(m - 1, 0)
        i2 = min(m + 1, n - 1)
        if i2 - i1 + 1 < 3:
            return bad
    s = float(v[i2] - v[i1])
    u = (v[i1:i2 + 1] - v[m]) / s
    design = np.column_stack([u ** 2, u, np.ones(u.shape)])
    c = np.linalg.lstsq(design, y[i1:i2 + 1], rcond=None)[0]
    coef = np.array([c[0], c[1], c[2], v[m], s])
    if not c[0] < -1e-9 * y[m]:
        return np.nan, np.nan, i1, i2, coef
    u0 = -c[1] / (2 * c[0])
    pos = float(v[m] + s * u0)
    if pos < v[i1] or pos > v[i2]:
        return np.nan, np.nan, i1, i2, coef
    return pos, float(c[2] - c[1] ** 2 / (4 * c[0])), i1, i2, coef


def band_position(v, Y, method='moment', frac=0.8, min_rel=0.05):
    """Position of a band in every column (delay) of Y.

    v      : values of the position axis, increasing (n)
    Y      : signal, n x number of delays, the band positive. A delay
             with an empty value (NaN) gives no position
    method : 'moment'     first moment, int v y dv / int y dv (trapezoid)
             'parabola'   maximum of a parabola through the top points
                          (band_parabola, with frac)
             'lognormal'  maximum of a log-normal band fitted to the
                          whole range (height, position, width,
                          asymmetry)
    min_rel: a delay counts only when the largest value of its column is
             above 0 and at least min_rel times the largest size
             (absolute value) of the signal in all columns

    Returns a Result with one value per delay (NaN where there is none):
      pos     the position
      height  height of the band: largest value ('moment'), top of the
              parabola, or fitted height
      width   full width at half maximum: 2.3548 x the standard
              deviation of the band about pos ('moment'), from the two
              half-height points of the data ('parabola'), or of the
              fitted band
      area    integral of the signal over the range ('lognormal': of
              the fitted band over the whole axis)
      asym    asymmetry of the fitted band ('lognormal' only)
      rms     RMS difference between data and fitted band ('lognormal')
      ok      the delay counted
      par     numbers that draw the model of a delay (band_model_curve)
    """
    v = np.asarray(v, dtype=float).ravel()
    Y = np.atleast_2d(np.asarray(Y, dtype=float))
    n, nt = Y.shape
    R = Result(pos=np.full(nt, np.nan), height=np.full(nt, np.nan),
               width=np.full(nt, np.nan), area=np.full(nt, np.nan),
               asym=np.full(nt, np.nan), rms=np.full(nt, np.nan),
               ok=np.zeros(nt, dtype=bool), par=np.full((5, nt), np.nan),
               method=method)
    if n < 3 or nt == 0:
        return R
    full = np.all(np.isfinite(Y), axis=0)
    top = np.full(nt, -np.inf)
    top[full] = np.max(Y[:, full], axis=0)
    if not np.any(full) or not np.max(top) > 0:
        return R
    R.ok = (full & (top >= min_rel * float(np.max(np.abs(Y[:, full]))))
            & (top > 0))
    dv = np.diff(v)
    for j in np.flatnonzero(R.ok):
        y = Y[:, j]
        m = int(np.argmax(y))
        area = float(np.sum(dv * (y[:-1] + y[1:])) / 2)
        if method == 'moment':
            if not area > 0:
                continue
            pos = float(np.sum(dv * (v[:-1] * y[:-1] + v[1:] * y[1:])) / 2) / area
            d2 = (v - pos) ** 2 * y
            var = float(np.sum(dv * (d2[:-1] + d2[1:])) / 2) / area
            R.pos[j] = pos
            R.height[j] = y[m]
            R.area[j] = area
            if var > 0:
                R.width[j] = 2 * np.sqrt(2 * LN2) * np.sqrt(var)
        elif method == 'parabola':
            pos, hgt, i1, i2, coef = band_parabola(v, y, frac)
            R.par[:, j] = coef
            if not np.isfinite(pos):
                continue
            R.pos[j] = pos
            R.height[j] = hgt
            R.area[j] = area
            R.width[j] = band_half_width(v, y, m, hgt)
        else:
            # start: top of the parabola, width of the data at half
            # height (else a quarter of the range), no asymmetry
            p0, h0, _, _, _ = band_parabola(v, y, 0.8)
            if not np.isfinite(p0):
                p0, h0 = float(v[m]), float(y[m])
            w0 = band_half_width(v, y, m, h0)
            if not (np.isfinite(w0) and w0 > 0):
                w0 = (v[-1] - v[0]) / 4

            def split(q, y=y):
                shape = lognorm_shape(v, q[0], np.exp(q[1]), q[2])
                den = float(shape @ shape)
                hh = float(shape @ y) / den if den > 0 else 0.0
                return hh, shape

            def res(q, y=y):
                hh, shape = split(q)
                return y - hh * shape

            q, rss, _ = lm_fit(res, [p0, np.log(w0), 0.0])
            hh, _ = split(q)
            wid = float(np.exp(q[1]))
            gam = min(max(float(q[2]), -0.95), 0.95)
            if not (np.isfinite(rss) and hh > 0 and np.isfinite(wid)
                    and v[0] <= q[0] <= v[-1]):
                continue
            R.pos[j] = q[0]
            R.height[j] = hh
            R.asym[j] = gam
            R.width[j] = wid if abs(gam) < 1e-9 else wid * np.sinh(gam) / gam
            R.area[j] = (hh * wid * np.sqrt(np.pi / (4 * LN2))
                         * np.exp(gam ** 2 / (4 * LN2)))
            R.rms[j] = np.sqrt(rss / n)
            R.par[:, j] = [q[0], wid, gam, hh, 0.0]
    R.ok = R.ok & np.isfinite(R.pos)
    return R


def band_model_curve(v, method, par):
    """The model of one delay on the axis values v, from the column of
    Result.par of that delay: the parabola ('parabola'; NaN outside the
    points it was fitted to is left to the caller) or the fitted band
    ('lognormal'). None for 'moment'."""
    v = np.asarray(v, dtype=float)
    par = np.asarray(par, dtype=float).ravel()
    if method == 'parabola':
        u = (v - par[3]) / par[4]
        return par[0] * u ** 2 + par[1] * u + par[2]
    if method == 'lognormal':
        return par[3] * lognorm_shape(v, par[0], par[1], par[2])
    return None


def band_shift_start(t, y, n_exp):
    """Starting time constants for the fit of a position y(t) that runs
    from its first to its last value: the times (counted from the first
    delay) at which k / (n_exp + 1) of the change has been made,
    k = 1 ... n_exp. The last value is the mean of the last three
    points. Each at least half the first delay step, and each at least
    1.5 times the one before."""
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    n = t.size
    y_end = float(np.mean(y[max(n - 3, 0):]))
    tot = y_end - y[0]
    step = float(t[1] - t[0]) if n > 1 else 1.0
    span = float(t[-1] - t[0]) if n > 1 else 1.0
    tau = np.zeros(n_exp)
    for k in range(1, n_exp + 1):
        tq = span * k / (n_exp + 1.0)
        if tot != 0:
            f = (y - y[0]) / tot
            hit = np.flatnonzero(f >= k / (n_exp + 1.0))
            if hit.size:
                tq = float(t[hit[0]] - t[0])
        tau[k - 1] = max(tq, 0.5 * step)
        if k > 1:
            tau[k - 1] = max(tau[k - 1], 1.5 * tau[k - 2])
    return tau


def fit_band_shift(t, pos, t_range, n_exp, tau0=None):
    """Fit of the position of a band along the delays:

        pos(t) = v_inf + sum_i A_i exp(-t / tau_i)

    in the delays t_range[0] ... t_range[1]. The amplitudes are those at
    delay 0, so v0 = v_inf + sum(A_i) is the position extrapolated to
    delay 0. For given time constants the amplitudes and v_inf follow by
    linear least squares; the time constants are found by lm_fit on
    their logarithms.

    tau0: starting time constants (n_exp values). None: three starts are
    tried and the best fit is kept: band_shift_start, and two sets
    spaced evenly on a logarithmic scale between the first delay step
    and 1/5 or 1/50 of the range. A later start replaces an earlier one
    only when its sum of squares is lower by more than 1e-6 of it. The
    time constants are kept between 1/1000 of the first delay step and
    1000 times the range.

    Returns a Result:
      ok, msg          False and a text when no fit was possible
      tau, A           time constants (increasing) and their amplitudes
      v_inf, v0        position at long delays, and extrapolated to 0
      shift            v0 - v_inf, the sum of the amplitudes
      tau_avg          sum(A tau) / sum(A)
      tau_std, A_std, v_inf_std, v0_std   standard deviations (inf: not
                       determined by the data)
      rms, n           RMS of the residual, number of delays fitted
      mask             the delays fitted
      fit              the fitted curve at every delay of t
      tau_start        the starting values that led to the result
      tau_lim          True for a time constant that lies at one of the
                       two limits: not determined
    """
    t = np.asarray(t, dtype=float).ravel()
    pos = np.asarray(pos, dtype=float).ravel()
    n_exp = int(n_exp)
    R = Result(ok=False, msg='', tau=np.zeros(0), A=np.zeros(0),
               v_inf=np.nan, v0=np.nan, shift=np.nan, tau_avg=np.nan,
               tau_std=np.zeros(0), A_std=np.zeros(0), v_inf_std=np.nan,
               v0_std=np.nan, rms=np.nan, n=0,
               mask=np.zeros(t.shape, dtype=bool),
               fit=np.full(t.shape, np.nan), tau_start=np.zeros(0),
               tau_lim=np.zeros(0, dtype=bool))
    tr = np.asarray(t_range, dtype=float).ravel()
    if tr.size < 2 or not np.all(np.isfinite(tr[:2])) or tr[0] == tr[1]:
        R.msg = 'Give the delay range of the fit.'
        return R
    mask = (t >= min(tr[:2])) & (t <= max(tr[:2])) & np.isfinite(pos)
    R.mask = mask
    R.n = int(np.count_nonzero(mask))
    if n_exp < 1 or R.n < 2 * n_exp + 2:
        R.msg = ('%d positions in the range: too few for %d exponential%s '
                 '(at least %d).' % (R.n, n_exp, '' if n_exp == 1 else 's',
                                     2 * n_exp + 2))
        return R
    tf = t[mask]
    yf = pos[mask]
    step = float(tf[1] - tf[0])
    span = float(tf[-1] - tf[0])
    lo = np.log(1e-3 * step)
    hi = np.log(1e3 * span)

    def basis(tau, tt):
        return np.column_stack([np.exp(-tt.reshape(-1, 1) / tau.reshape(1, -1)),
                                np.ones(tt.size)])

    def linear(q):
        tau = np.exp(np.minimum(np.maximum(q, lo), hi))
        B = basis(tau, tf)
        c = np.linalg.lstsq(B, yf, rcond=None)[0]
        return tau, c, yf - B @ c

    def res(q):
        return linear(q)[2]

    if tau0 is not None and np.size(tau0):
        tau0 = np.asarray(tau0, dtype=float).ravel()
        if tau0.size != n_exp or not np.all(np.isfinite(tau0) & (tau0 > 0)):
            R.msg = ('Give %d starting time constant%s above 0, or leave '
                     'the field empty.' % (n_exp, '' if n_exp == 1 else 's'))
            return R
        starts = [tau0]
    else:
        starts = [band_shift_start(tf, yf, n_exp)]
        for div in (5.0, 50.0):
            a, b = np.log(step), np.log(max(span / div, 2 * step))
            if n_exp == 1:
                starts.append(np.array([np.exp(0.5 * (a + b))]))
            else:
                starts.append(np.exp(np.linspace(a, b, n_exp)))
    best = None
    for s0 in starts:
        q, rss, _ = lm_fit(res, np.log(s0))
        if np.isfinite(rss) and (best is None
                                 or rss < best[1] * (1 - 1e-6)):
            best = (q, rss, s0)
    if best is None:
        R.msg = 'The fit of the position failed.'
        return R
    q, rss, s0 = best
    tau, c, _ = linear(q)
    order = np.argsort(tau, kind='stable')
    tau = tau[order]
    amp = c[:n_exp][order]
    v_inf = float(c[n_exp])
    R.tau = tau
    R.tau_lim = ((tau >= np.exp(hi) * (1 - 1e-6))
                 | (tau <= np.exp(lo) * (1 + 1e-6)))
    R.A = amp
    R.v_inf = v_inf
    R.shift = float(np.sum(amp))
    R.v0 = v_inf + R.shift
    if R.shift != 0:
        R.tau_avg = float(np.sum(amp * tau) / R.shift)
    R.rms = float(np.sqrt(rss / R.n))
    R.fit = basis(tau, t) @ np.r_[amp, v_inf]
    R.tau_start = np.asarray(s0, dtype=float)[np.argsort(s0, kind='stable')]

    # standard deviations: all parameters free; for v0 the last
    # amplitude is written as v0 - v_inf - (the other amplitudes)
    def res_all(p):
        return yf - basis(p[:n_exp], tf) @ p[n_exp:]

    def res_v0(p):
        a_last = p[2 * n_exp - 1] - p[2 * n_exp] - np.sum(p[n_exp:2 * n_exp - 1])
        return yf - basis(p[:n_exp], tf) @ np.r_[p[n_exp:2 * n_exp - 1],
                                                  a_last, p[2 * n_exp]]

    p_all = np.r_[tau, amp, v_inf]
    sig = fit_param_uncertainty(res_all, p_all, rss, R.n)
    R.tau_std = sig[:n_exp]
    R.A_std = sig[n_exp:2 * n_exp]
    R.v_inf_std = float(sig[2 * n_exp])
    p_v0 = np.r_[tau, amp[:-1], R.v0, v_inf]
    R.v0_std = float(fit_param_uncertainty(res_v0, p_v0, rss, R.n)[2 * n_exp - 1])
    R.ok = True
    return R
