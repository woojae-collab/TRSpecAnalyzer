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
    n_fft, t_uni, n_win.
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
    if nmode == 'peak':
        pk = np.nanmax(P) if P.size else np.nan
        if np.isfinite(pk) and pk > 0:
            P = P / pk
    elif nmode == 'perwl':
        pk_w = np.nanmax(P, axis=1, keepdims=True)
        pk_w[pk_w <= 0] = 1
        P = P / pk_w
    P[~np.isfinite(P)] = 0
    return Result(pmap=P, freq=freq, phase=phi, dt=dt_uni, n_fft=n_fft, t_uni=t_uni,
                  n_win=n_t, apod=apod, zero_pad=zpad, detrend=detrend,
                  freq_unit=freq_unit, norm_mode=norm_mode)


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
