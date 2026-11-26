"""Vibrational coherence: FFT map of the residual and time-frequency maps.

Port of computeCoherence and computeTimeFrequency of TRSpecAnalyzer.m.
"""
import numpy as np
from ._util import Result, colon as _colon

_REALMIN = np.finfo(float).tiny
_TU2S = {'fs': 1e-15, 'ps': 1e-12, 'ns': 1e-9, 'us': 1e-6}
_C_CM = 2.99792458e10


def _round(x):
    return int(np.floor(x + 0.5))


def _nextpow2(x):
    return int(np.ceil(np.log2(x))) if x > 1 else 0




def _freq_conv(time_unit, freq_unit):
    tu2s = _TU2S.get(time_unit, 1e-12)
    fu = str(freq_unit).lower()
    if fu == 'cm-1':
        return 1.0 / (tu2s * _C_CM)
    if fu == 'thz':
        return 1.0 / (tu2s * 1e12)
    return 1.0 / tu2s


def _window(name, n):
    k = np.arange(n) / (n - 1)
    name = str(name).lower()
    if name == 'hann':
        return 0.5 * (1 - np.cos(2 * np.pi * k))
    if name == 'hamming':
        return 0.54 - 0.46 * np.cos(2 * np.pi * k)
    if name == 'blackman':
        return 0.42 - 0.5 * np.cos(2 * np.pi * k) + 0.08 * np.cos(4 * np.pi * k)
    return np.ones(n)


def compute_coherence(R, t, t_min=None, t_max=None, apod='hann', zero_pad=2,
                      detrend='linear', freq_unit='cm-1', norm_mode='peak', time_unit='ps'):
    """|FFT|^2 map of the residual R (n_wavelength x n_delay).

    The delays inside [t_min, t_max] are resampled to a uniform grid (step =
    median step of the window) before the FFT.

    apod      : 'rect', 'hann', 'hamming', 'blackman'
    zero_pad  : zero-padding factor (>= 1)
    detrend   : 'none', 'mean', 'linear'
    freq_unit : 'cm-1', 'THz' or 'Hz';  time_unit : unit of t (fs, ps, ns, us)
    norm_mode : 'peak' (largest value = 1), 'perWl' (every wavelength), 'none'

    Returns a Result with pmap (n_wavelength x n_freq), freq, phase, dt,
    n_fft, t_uni, n_win and norm_scale (what the raw |FFT|^2 was divided
    by: 1, the largest value, or one value per wavelength as a column).
    """
    R = np.asarray(R, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    if t_min is None:
        t_min = np.nanmin(t)
    if t_max is None:
        t_max = np.nanmax(t)
    mask = (t >= t_min) & (t <= t_max)
    tw = t[mask]
    Rw = R[:, mask]
    if tw.size < 8:
        raise ValueError('Time window is too short (need >= 8 points, got %d).' % tw.size)
    dt_uni = float(np.median(np.diff(tw)))
    n_uni = abs(tw[-1] - tw[0]) / max(abs(dt_uni), _REALMIN) + 1
    if n_uni * max(1, _round(zero_pad)) * Rw.shape[0] > 4e7:
        raise ValueError('The time window resamples to about %.3g uniform points (step %.3g) for %d rows '
                         '-- too large for the FFT map. Narrow the window to the densely sampled part '
                         'of the delay axis.' % (n_uni, dt_uni, Rw.shape[0]))
    t_uni = _colon(tw[0], dt_uni, tw[-1])
    if t_uni.size < 8:
        raise ValueError('Resampled grid too short (n=%d).' % t_uni.size)
    n_t = t_uni.size
    Ru = np.zeros((Rw.shape[0], n_t))
    for ii in range(Rw.shape[0]):
        row = Rw[ii, :]
        m = np.isfinite(row)
        if np.count_nonzero(m) >= 4:
            Ru[ii, :] = np.interp(t_uni, tw[m], row[m], left=0.0, right=0.0)
    dmode = str(detrend).lower()
    if dmode == 'mean':
        Ru = Ru - np.nanmean(Ru, axis=1, keepdims=True)
    elif dmode == 'linear':
        A = np.column_stack([t_uni, np.ones(n_t)])
        for ii in range(Ru.shape[0]):
            y = Ru[ii, :]
            m = np.isfinite(y)
            if np.count_nonzero(m) < 2:
                continue
            c = np.linalg.lstsq(A[m, :], y[m], rcond=None)[0]
            Ru[ii, :] = y - A @ c
    Ru[~np.isfinite(Ru)] = 0
    Ru = Ru * _window(apod, n_t)[None, :]
    zpad = max(1, _round(zero_pad))
    n_fft = max(2 ** _nextpow2(n_t * zpad), n_t)
    X = np.fft.fft(Ru, n=n_fft, axis=1)
    n_half = n_fft // 2 + 1
    Xh = X[:, :n_half]
    P = np.abs(Xh) ** 2
    phi = np.arctan2(Xh.imag, Xh.real)
    f_nyq = 1.0 / (2 * dt_uni)
    freq = np.linspace(0, f_nyq, n_half) * _freq_conv(time_unit, freq_unit)
    nmode = str(norm_mode).lower()
    # norm_scale: what the raw |FFT|^2 was divided by (1, the largest value,
    # or one value per wavelength). A second map is put on the same scale
    # by dividing its raw |FFT|^2 by it.
    norm_scale = 1.0
    if nmode == 'peak':
        pk = np.nanmax(P) if P.size else np.nan
        if np.isfinite(pk) and pk > 0:
            P = P / pk
            norm_scale = float(pk)
    elif nmode == 'perwl':
        pk_w = np.nanmax(P, axis=1, keepdims=True)
        pk_w[pk_w <= 0] = 1
        P = P / pk_w
        norm_scale = np.array(pk_w, float)
    P[~np.isfinite(P)] = 0
    return Result(pmap=P, freq=freq, phase=phi, dt=dt_uni, n_fft=n_fft, t_uni=t_uni,
                  n_win=n_t, apod=apod, zero_pad=zpad, detrend=detrend,
                  freq_unit=freq_unit, norm_mode=norm_mode, norm_scale=norm_scale)


# ----------------------------------------------------------------------
# Solvent coherence: the vibrational coherence of the solvent, measured
# alone, is taken out of the residual before the FFT
#     residual - A * solvent coherence(t - shift)
# A and the shift are fitted in the frequency domain, on the real and
# imaginary parts of the Fourier transforms, in a band that holds a mode
# of the solvent; they are then applied to the residual itself
# (M. Liebel, thesis, Oxford 2014, Fig. 4.6).
# Port of cohSpectrum, solventVC, solventShifted, solventBand and
# fitSolventCoherence of TRSpecAnalyzer.m.
# ----------------------------------------------------------------------
_SHIFT_K = 16      # points to each side in the interpolation of a shift


def coh_spectrum(R, t, t_min=None, t_max=None, apod='hann', zero_pad=2, detrend='linear',
                 freq_unit='cm-1', time_unit='ps'):
    """Complex one-sided Fourier transform of the residual R, made exactly as
    compute_coherence makes the |FFT|^2 map (time window, uniform grid,
    detrending, apodization, zero padding), without normalization.

    Returns (F, freq, o): F is n_wavelength x n_freq (complex), o the Result
    of compute_coherence (dt, n_fft, n_win, ...).
    """
    o = compute_coherence(R, t, t_min, t_max, apod, zero_pad, detrend, freq_unit, 'none', time_unit)
    F = np.sqrt(o.pmap) * np.exp(1j * o.phase)
    return F, o.freq, o


def solvent_vc(Xs, t, t_min, t_max, max_shift=0.0):
    """Vibrational coherence of the solvent data Xs (n_wavelength x n_delay,
    delays t) for the time window [t_min, t_max].

    The data are put on the uniform grid of the FFT, continued on both sides
    by n_pad = ceil(max_shift / step) + 16 points, and their smooth part (a
    polynomial of order 3 fitted in the window) is removed. With these points
    a shift by up to max_shift is made from measured values everywhere in the
    window. Where the delays beyond the window are missing or sparse, the
    data are continued by point reflection at the last dense delay.

    Returns a Result with t_e (grid), dt, R (n_wavelength x len(t_e)), w
    (weights that let a subtraction fade out beyond the window: 1 inside),
    n_pad, n_u (points of the window).
    """
    Xs = np.asarray(Xs, dtype=float)
    t = np.asarray(t, dtype=float).ravel()
    mask = (t >= t_min) & (t <= t_max)
    tw = t[mask]
    if tw.size < 8:
        raise ValueError('Time window is too short (need >= 8 points, got %d).' % tw.size)
    dt = float(np.median(np.diff(tw)))
    n_u = _colon(tw[0], dt, tw[-1]).size
    if not np.isfinite(max_shift):
        max_shift = 0.0
    n_pad = int(np.ceil(abs(max_shift) / dt)) + _SHIFT_K
    t_e = tw[0] + np.arange(-n_pad, n_u + n_pad) * dt
    n = t_e.size
    # densely sampled delays around the window (steps up to 1.5 x the step of the grid)
    idx = np.flatnonzero(mask)
    iw1, iw2 = int(idx[0]), int(idx[-1])
    while iw1 > 0 and (t[iw1] - t[iw1 - 1]) <= 1.5 * dt:
        iw1 -= 1
    while iw2 < t.size - 1 and (t[iw2 + 1] - t[iw2]) <= 1.5 * dt:
        iw2 += 1
    t_lo, t_hi = t[iw1], t[iw2]
    lo = t_e < t_lo
    hi = t_e > t_hi
    tq = t_e.copy()
    tq[lo] = 2 * t_lo - tq[lo]
    tq[hi] = 2 * t_hi - tq[hi]
    tq = np.minimum(np.maximum(tq, t_lo), t_hi)
    in_w = np.zeros(n, bool)
    in_w[n_pad:n_pad + n_u] = True
    t_mid = (t_e[n_pad] + t_e[n_pad + n_u - 1]) / 2
    t_half = max((t_e[n_pad + n_u - 1] - t_e[n_pad]) / 2, _REALMIN)
    x = (t_e - t_mid) / t_half
    Pm = np.column_stack([np.ones(n), x, x ** 2, x ** 3])
    R = np.zeros((Xs.shape[0], n))
    for ii in range(Xs.shape[0]):
        row = Xs[ii, :]
        m = np.isfinite(row)
        if np.count_nonzero(m) < 4:
            continue
        tm = t[m]
        y = np.interp(tq, tm, row[m])
        y_lo = float(np.interp(t_lo, tm, row[m]))
        y_hi = float(np.interp(t_hi, tm, row[m]))
        y[lo] = 2 * y_lo - y[lo]
        y[hi] = 2 * y_hi - y[hi]
        c = np.linalg.lstsq(Pm[in_w, :], y[in_w], rcond=None)[0]
        R[ii, :] = y - Pm @ c
    w = np.ones(n)
    e = 0.5 * (1 - np.cos(np.pi * np.arange(1, n_pad + 1) / (n_pad + 1)))
    w[:n_pad] = e
    w[n - n_pad:] = e[::-1]
    return Result(t_e=t_e, dt=dt, R=R, w=w, n_pad=n_pad, n_u=n_u)


def solvent_shifted(V, t, d):
    """The solvent coherence V (Result of solvent_vc) shifted by d in time,
    V(t - d), on the delays t (0 outside the grid of V). d > 0 moves it to
    later delays.

    Whole steps of the grid are moved as they are; the fraction of a step is
    made by band-limited interpolation (sinc kernel with a Hann window, 16
    points to each side), which keeps the fast oscillations: linear
    interpolation would damp them. The result fades to zero beyond the
    window (weights V.w).
    """
    X = np.asarray(V.R, dtype=float)
    n = X.shape[1]
    K = _SHIFT_K
    ms = int(np.floor(d / V.dt))
    dl = d / V.dt - ms
    kk = np.arange(-K, K + 1)
    if dl == 0:
        h = (kk == 0).astype(float)
    else:
        x = kk - dl
        h = np.sin(np.pi * x) / (np.pi * x) * (0.5 + 0.5 * np.cos(np.pi * x / (K + 1)))
        h = h / np.sum(h)
    Y = np.zeros_like(X)
    for q in range(kk.size):
        if h[q] == 0:
            continue
        sh = ms + int(kk[q])                 # Y[:, j] takes X[:, j - sh]
        j1 = max(0, sh)
        j2 = min(n, n + sh)
        if j1 < j2:
            Y[:, j1:j2] += h[q] * X[:, j1 - sh:j2 - sh]
    Y = Y * np.asarray(V.w, dtype=float)[None, :]
    t = np.asarray(t, dtype=float).ravel()
    Yt = np.zeros((Y.shape[0], t.size))
    inside = np.flatnonzero((t >= V.t_e[0]) & (t <= V.t_e[-1]))
    if inside.size:
        xq = (t[inside] - V.t_e[0]) / V.dt
        i0 = np.minimum(np.maximum(np.floor(xq), 0), n - 2).astype(int)
        wq = xq - i0
        Yt[:, inside] = Y[:, i0] * (1 - wq) + Y[:, i0 + 1] * wq
    return Yt


def solvent_band(G, freq, rows, n_fft, n_win):
    """Frequency band of the strongest mode of the solvent: around the
    largest value of its power spectrum (G: complex spectrum of the solvent
    coherence, summed over the rows given), as far as the power stays above
    10 % of that value. Frequencies below 3 cycles per time window are not
    considered. Returns (lo, hi); (nan, nan) when there is none.
    """
    rows = np.asarray(rows, bool).ravel()
    if not rows.any():
        return np.nan, np.nan
    freq = np.asarray(freq, dtype=float).ravel()
    P = np.sum(np.abs(np.asarray(G)[rows, :]) ** 2, axis=0)
    n_f = freq.size
    i0 = min(n_f, int(np.ceil(3.0 * n_fft / n_win)) + 1) - 1
    if i0 >= n_f:
        return np.nan, np.nan
    ip = int(np.argmax(P[i0:])) + i0
    pk = P[ip]
    if not (pk > 0):
        return np.nan, np.nan
    i_l = ip
    while i_l > i0 and P[i_l - 1] > 0.1 * pk:
        i_l -= 1
    i_r = ip
    while i_r < n_f - 1 and P[i_r + 1] > 0.1 * pk:
        i_r += 1
    return float(freq[i_l]), float(freq[i_r])


def fit_solvent_coherence(R, V, t, band, rows=None, max_shift=0.0, fit_shift=True, d0=0.0,
                          t_min=None, t_max=None, apod='hann', zero_pad=2, detrend='linear',
                          freq_unit='cm-1', time_unit='ps'):
    """Factor A and time shift d of the solvent coherence V that take it out
    of the residual R best: least squares on the real and imaginary parts of
    the Fourier transforms, over the rows given, in the frequency band
    band[0] ... band[1] (unit freq_unit):

        FT[R] ~ A * FT[V(t - d)]

    fit_shift False: d = d0. Otherwise d is searched in [-max_shift,
    max_shift]: on a grid with 16 points per period of the highest frequency
    of the band (at least 41 points), with FT[V(t - d)] = exp(-2 pi i f d)
    FT[V], then by a golden-section search between the neighbors of the best
    point (60 steps). A follows in closed form, with the coherence shifted as
    it is subtracted.

    Returns a Result with A, d, P0, P1 (power of the residual in the band
    before and after the subtraction), n_freq, n_rows.
    """
    kw = dict(t_min=t_min, t_max=t_max, apod=apod, zero_pad=zero_pad, detrend=detrend,
              freq_unit=freq_unit, time_unit=time_unit)
    F, freq, o = coh_spectrum(R, t, **kw)
    freq = np.asarray(freq, dtype=float).ravel()
    fm = (freq >= min(band)) & (freq <= max(band))
    if not fm.any():
        raise ValueError('The band holds no frequency of the FFT. Widen it.')
    if rows is None:
        rows = np.ones(F.shape[0], bool)
    rows = np.asarray(rows, bool).ravel()
    if not rows.any():
        raise ValueError('The spectral range of the fit holds no point.')
    d = float(d0)
    if fit_shift and max_shift > 0:
        G0 = coh_spectrum(solvent_shifted(V, t, 0.0), t, **kw)[0]
        c = np.sum(np.conj(G0[rows][:, fm]) * F[rows][:, fm], axis=0)
        fB = np.flatnonzero(fm) / (o.n_fft * o.dt)          # cycles per time unit
        nG = 2 * int(np.ceil(16 * max_shift * np.max(fB))) + 1
        nG = min(max(nG, 41), 4001)
        dG = np.linspace(-max_shift, max_shift, nG)
        Sg = np.real(np.exp(2j * np.pi * np.outer(dG, fB)) @ c)
        ig = int(np.argmax(Sg))
        lo = float(dG[max(ig - 1, 0)])
        hi = float(dG[min(ig + 1, nG - 1)])

        def sfun(xv):
            return -float(np.real(np.sum(c * np.exp(2j * np.pi * fB * xv))))

        gr = (np.sqrt(5.0) - 1) / 2
        x1 = hi - gr * (hi - lo)
        x2 = lo + gr * (hi - lo)
        f1, f2 = sfun(x1), sfun(x2)
        for _ in range(60):
            if f1 <= f2:
                hi, x2, f2 = x2, x1, f1
                x1 = hi - gr * (hi - lo)
                f1 = sfun(x1)
            else:
                lo, x1, f1 = x1, x2, f2
                x2 = lo + gr * (hi - lo)
                f2 = sfun(x2)
        d = (lo + hi) / 2
    Gd = coh_spectrum(solvent_shifted(V, t, d), t, **kw)[0]
    a = F[rows][:, fm].ravel()
    b = Gd[rows][:, fm].ravel()
    den = float(np.sum(np.abs(b) ** 2))
    A = float(np.real(np.sum(np.conj(b) * a)) / den) if den > 0 else 0.0
    return Result(A=A, d=float(d), P0=float(np.sum(np.abs(a) ** 2)),
                  P1=float(np.sum(np.abs(a - A * b) ** 2)),
                  n_freq=int(np.count_nonzero(fm)), n_rows=int(np.count_nonzero(rows)))


def compute_time_frequency(y, t, t_min=None, t_max=None, detrend='linear', method='stft',
                           freq_unit='cm-1', norm_mode='peak', win_len=np.nan, hop=np.nan,
                           apod='hann', zero_pad=4, omega0=6, n_freq=200, f_lo=0, f_hi=0,
                           time_unit='ps'):
    """Time-frequency map of one residual trace y(t).

    method : 'stft' (sliding-window FFT; win_len and hop in time units) or
             'cwt' (Morlet wavelet; omega0, n_freq, f_lo, f_hi)

    Returns a Result with ptf (n_freq x n_time), freq, tc (time of every
    column), dt, t_uni, yu (resampled, detrended trace), coi_freq, coi_t
    and the settings used.
    """
    t = np.asarray(t, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    if t_min is None:
        t_min = np.nanmin(t)
    if t_max is None:
        t_max = np.nanmax(t)
    mask = (t >= t_min) & (t <= t_max) & np.isfinite(y)
    tw, yw = t[mask], y[mask]
    if tw.size < 16:
        raise ValueError('Time window is too short (need >= 16 points, got %d).' % tw.size)
    dt_uni = float(np.median(np.diff(tw)))
    n_uni = abs(tw[-1] - tw[0]) / max(abs(dt_uni), _REALMIN) + 1
    if n_uni > 2e5:
        raise ValueError('The time window resamples to about %.3g uniform points (step %.3g) -- too many '
                         'for a time-frequency map. Narrow the window to the densely sampled part of '
                         'the delay axis.' % (n_uni, dt_uni))
    t_uni = _colon(tw[0], dt_uni, tw[-1])
    n_t = t_uni.size
    if n_t < 16:
        raise ValueError('Resampled grid too short (n=%d).' % n_t)
    yu = np.interp(t_uni, tw, yw, left=0.0, right=0.0)
    dmode = str(detrend).lower()
    if dmode == 'mean':
        yu = yu - np.nanmean(yu)
    elif dmode == 'linear':
        A = np.column_stack([t_uni, np.ones(n_t)])
        yu = yu - A @ np.linalg.lstsq(A, yu, rcond=None)[0]
    yu[~np.isfinite(yu)] = 0
    f_conv = _freq_conv(time_unit, freq_unit)
    f_nyq_tu = 1.0 / (2 * dt_uni)
    out = Result(dt=dt_uni, t_uni=t_uni, yu=yu, method=str(method).lower(),
                 freq_unit=freq_unit, detrend=detrend)
    meth = str(method).lower()
    if meth == 'stft':
        wl_ = win_len
        if not np.isfinite(wl_) or wl_ <= 0:
            wl_ = (t_uni[-1] - t_uni[0]) / 4
        n_win = _round(wl_ / dt_uni)
        n_win = max(8, min(n_win, n_t))
        if n_win % 2 == 0:
            n_win -= 1
        half_w = (n_win - 1) // 2
        hop_n = 1 if (not np.isfinite(hop) or hop <= 0) else max(1, _round(hop / dt_uni))
        w = _window(apod if str(apod).lower() in ('rect', 'hamming', 'blackman') else 'hann', n_win)
        zpad = max(1, _round(zero_pad))
        n_fft = max(2 ** _nextpow2(n_win * zpad), n_win)
        n_half = n_fft // 2 + 1
        centers = np.arange(0, n_t, hop_n)
        ypad = np.concatenate([np.zeros(half_w), yu, np.zeros(half_w)])
        ptf = np.zeros((n_half, centers.size))
        for k, c in enumerate(centers):
            X = np.fft.fft(ypad[c:c + n_win] * w, n=n_fft)
            ptf[:, k] = np.abs(X[:n_half]) ** 2
        tc = t_uni[centers]
        freq = np.linspace(0, f_nyq_tu, n_half) * f_conv
        edge_t = half_w * dt_uni
        coi_freq = np.zeros(centers.size)
        coi_freq[(tc < t_uni[0] + edge_t) | (tc > t_uni[-1] - edge_t)] = np.inf
        out.update(n_win=n_win, win_len=n_win * dt_uni, hop_n=hop_n, n_fft=n_fft, apod=apod,
                   zero_pad=zpad, coi_freq=coi_freq,
                   coi_t=np.array([t_uni[0] + edge_t, t_uni[-1] - edge_t]))
    elif meth == 'cwt':
        om0 = omega0
        if not np.isfinite(om0) or om0 < 2:
            om0 = 6
        n_f = max(16, _round(n_freq))
        fh = f_hi
        if not np.isfinite(fh) or fh <= 0:
            fh = f_nyq_tu * f_conv
        fh = min(fh, f_nyq_tu * f_conv)
        fl = f_lo
        f_min_useful = 2.0 / (t_uni[-1] - t_uni[0]) * f_conv
        if not np.isfinite(fl) or fl <= 0:
            fl = f_min_useful
        fl = max(fl, f_min_useful / 4)
        if fl >= fh:
            fl = fh / 50
        freq = np.linspace(fl, fh, n_f)
        f_tu = freq / f_conv
        s_vec = om0 / (2 * np.pi * f_tu)
        n_pad = 2 ** _nextpow2(2 * n_t)
        y_hat = np.fft.fft(yu, n=n_pad)
        omega = 2 * np.pi * np.arange(n_pad) / (n_pad * dt_uni)
        hi = omega > np.pi / dt_uni
        omega[hi] = omega[hi] - 2 * np.pi / dt_uni
        pos = omega > 0
        ptf = np.zeros((n_f, n_t))
        for k in range(n_f):
            psi = np.zeros(n_pad)
            psi[pos] = 2 * np.exp(-0.5 * (s_vec[k] * omega[pos] - om0) ** 2)
            W = np.fft.ifft(y_hat * np.conj(psi))
            ptf[k, :] = np.abs(W[:n_t]) ** 2
        tc = t_uni.copy()
        d = np.minimum(tc - t_uni[0], t_uni[-1] - tc)
        d[d <= 0] = dt_uni / 2
        coi_freq = om0 / (2 * np.pi) * np.sqrt(2) / d * f_conv
        out.update(omega0=om0, n_freq=n_f, scales=s_vec, coi_freq=coi_freq, coi_t=np.zeros(0))
    else:
        raise ValueError('Unknown time-frequency method "%s".' % method)
    if str(norm_mode).lower() == 'peak':
        pk = np.nanmax(ptf) if ptf.size else np.nan
        if np.isfinite(pk) and pk > 0:
            ptf = ptf / pk
    ptf[~np.isfinite(ptf)] = 0
    out.update(ptf=ptf, freq=freq, tc=tc, norm_mode=norm_mode)
    return out
