"""Data file reading / writing and the preprocessing pipeline
(background, chirp, time zero, mask, crop).

Port of parseDataFile, writeDataFile, detectAxisSwap, recompute,
computeRidge, applyCropByRange and the related window logic of
TRSpecAnalyzer.m.
"""
import numpy as np
from ._util import interp1_nan
from .chirp import chirp_model, fit_chirp_params


def read_data_file(path):
    """Read a matrix file: first row = delays, first column = wavelengths
    (or wavenumbers), the rest = signal. Comma, tab or blank separated.

    Returns (wl, t, A) with both axes sorted ascending; A has shape
    (n_wavelength, n_delay).
    """
    with open(path, 'r', encoding='utf-8-sig', errors='replace') as fh:
        first = fh.readline()
    delim = ',' if ',' in first else (';' if ';' in first else None)
    M = np.genfromtxt(path, delimiter=delim, dtype=float, encoding='utf-8-sig')
    M = np.atleast_2d(M)
    if M.shape[0] < 3 or M.shape[1] < 3:
        raise ValueError('Data matrix is too small (must be at least 3x3).')
    t = M[0, 1:].copy()
    wl = M[1:, 0].copy()
    A = M[1:, 1:].copy()
    # rows / columns whose axis value is not a number (headers, units) go
    ok_w, ok_t = np.isfinite(wl), np.isfinite(t)
    if np.count_nonzero(ok_w) < 2 or np.count_nonzero(ok_t) < 2:
        raise ValueError('The first column (wavelengths) or the first row (delays) of '
                         'the file holds no numbers.')
    wl, t, A = wl[ok_w], t[ok_t], A[ok_w, :][:, ok_t]
    iw = np.argsort(wl, kind='stable')
    it = np.argsort(t, kind='stable')
    return wl[iw], t[it], A[iw, :][:, it]


def write_data_file(path, wl, t, A, delimiter=None):
    """Write (wl, t, A) in the layout read_data_file reads. The delimiter is
    a comma for .csv and a tab otherwise unless given."""
    wl = np.asarray(wl, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()
    A = np.asarray(A, dtype=float)
    if A.shape != (wl.size, t.size):
        raise ValueError('Data size mismatch: A is %dx%d, expected %dx%d'
                         % (A.shape[0], A.shape[1], wl.size, t.size))
    if delimiter is None:
        delimiter = ',' if str(path).lower().endswith('.csv') else '\t'
    out = np.zeros((wl.size + 1, t.size + 1))
    out[0, 1:] = t
    out[1:, 0] = wl
    out[1:, 1:] = A
    np.savetxt(path, out, delimiter=delimiter, fmt='%.10g')


def detect_axis_swap(wl, t):
    """True if the two axes of a file look interchanged."""
    wl = np.asarray(wl, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()
    wl_like_delay = bool(np.any(wl <= 0) or (np.min(wl) > 0 and np.max(wl) / np.min(wl) > 1000))
    t_like_wl = bool(np.all(t > 150) and np.all(t < 2500) and (np.max(t) - np.min(t) > 50))
    return wl_like_delay or t_like_wl


def default_chirp_tmax(time_unit='ps'):
    return {'fs': 500, 'ns': 0.002, 'us': 2e-6}.get(time_unit, 2)


def movmean(x, n):
    """Centered moving average with a window that shrinks at the ends
    (MATLAB smoothdata(x, 'movmean', n))."""
    x = np.asarray(x, dtype=float).ravel()
    n = int(n)
    if n <= 1:
        return x.copy()
    kb = (n - 1) // 2 if n % 2 else n // 2
    kf = (n - 1) // 2 if n % 2 else n // 2 - 1
    out = np.empty(x.size)
    for i in range(x.size):
        seg = x[max(0, i - kb):min(x.size, i + kf + 1)]
        out[i] = np.nanmean(seg) if np.any(np.isfinite(seg)) else np.nan
    return out


def compute_ridge(base, wl, t, method='max|dA/dt|', t_min=-0.5, t_max=1.5, smooth_n=5):
    """Estimated time zero of every wavelength: the delay of the steepest
    change ('max|dA/dt|') or of the largest signal ('max|deltaA|') inside
    [t_min, t_max], smoothed over smooth_n wavelengths.

    Returns (ridge_wl, ridge_t); empty arrays if fewer than 3 delays are
    inside the window.
    """
    base = np.asarray(base, dtype=float)
    wl = np.asarray(wl, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()
    m = (t >= t_min) & (t <= t_max)
    if np.count_nonzero(m) < 3:
        return np.zeros(0), np.zeros(0)
    t_sub = t[m]
    sub = base[:, m]
    if method == 'max|dA/dt|':
        rate = np.diff(sub, axis=1) / np.diff(t_sub)
        t_axis = (t_sub[:-1] + t_sub[1:]) / 2
        val = np.abs(rate)
    elif method == 'max|deltaA|':
        t_axis = t_sub
        val = np.abs(sub)
    else:
        return np.zeros(0), np.zeros(0)
    idx = np.argmax(np.where(np.isnan(val), -np.inf, val), axis=1)
    ridge_t = t_axis[idx]
    if smooth_n > 1 and ridge_t.size > smooth_n:
        ridge_t = movmean(ridge_t, smooth_n)
    return wl.copy(), ridge_t


def estimate_scan_shifts(datasets, ref=0, ref_wl=None):
    """Time shift of every scan relative to the reference scan (Auto-align
    to Ref): cross-correlation of the kinetics at ref_wl, with parabolic
    interpolation of the correlation peak.

    datasets : list of Dataset (or of (wl, t, A) tuples)
    ref      : index of the reference scan
    ref_wl   : wavelength of the kinetics used; choose one with a sharp rise
               (default: the middle wavelength, as in the program)

    Returns an array of shifts, to be passed to Dataset.from_average.
    A scan whose time zero is late gets a negative shift.
    """
    def parts(d):
        return (d.wavelength, d.delay, d.raw) if isinstance(d, Dataset) else d

    wl_r, t_r, a_r = parts(datasets[ref])
    if ref_wl is None:
        ref_wl = wl_r[int(np.floor(wl_r.size / 2 + 0.5)) - 1]

    def detrend(y):
        y = y - (np.nanmean(y) if np.any(np.isfinite(y)) else 0.0)
        y[~np.isfinite(y)] = 0
        return y

    from ._util import colon
    dts = np.diff(t_r)
    dt_uni = np.median(dts)
    if dt_uni <= 0 or not np.isfinite(dt_uni):
        dt_uni = np.mean(dts)
    t_uni = colon(t_r[0], dt_uni, t_r[-1])
    y_ref = detrend(np.interp(t_uni, t_r, a_r[int(np.argmin(np.abs(wl_r - ref_wl))), :], left=0.0, right=0.0))
    n_u = y_ref.size
    n_fft = 2 ** int(np.ceil(np.log2(2 * n_u - 1)))
    f_ref = np.fft.fft(y_ref, n_fft)
    lags = np.arange(-(n_u - 1), n_u)
    shifts = np.zeros(len(datasets))
    for i, d in enumerate(datasets):
        if i == ref:
            continue
        wl_i, t_i, a_i = parts(d)
        y_i = detrend(np.interp(t_uni, t_i, a_i[int(np.argmin(np.abs(wl_i - ref_wl))), :], left=0.0, right=0.0))
        r = np.real(np.fft.ifft(f_ref * np.conj(np.fft.fft(y_i, n_fft))))
        r = np.concatenate([r[n_fft - n_u + 1:], r[:n_u]])
        i_max = int(np.argmax(r))
        i_max = max(1, min(r.size - 2, i_max))
        r_l, r_0, r_r = r[i_max - 1], r[i_max], r[i_max + 1]
        denom = r_l - 2 * r_0 + r_r
        delta = 0.5 * (r_l - r_r) / denom if abs(denom) > 0 else 0.0
        shifts[i] = (lags[i_max] + delta) * dt_uni
    return shifts


class Dataset:
    """Time-resolved data set with the preprocessing steps of the program.

    The steps are kept as settings and applied to the raw data in a fixed
    order, whatever the order in which they were set:
    background -> chirp -> mask.  crop() makes the current result the new
    raw data (as the Crop window does).

    Attributes: wavelength, delay, data (current working data,
    n_wavelength x n_delay), raw.
    """

    def __init__(self, wavelength, delay, data, source=''):
        wl = np.asarray(wavelength, dtype=float).ravel()
        t = np.asarray(delay, dtype=float).ravel()
        A = np.asarray(data, dtype=float)
        if A.shape != (wl.size, t.size):
            raise ValueError('data must have shape (n_wavelength, n_delay).')
        self.original = (wl.copy(), t.copy(), A.copy())
        self.source = source
        self.time_unit = 'ps'
        self._load(wl, t, A)

    def _load(self, wl, t, A):
        self.wavelength = wl.copy()
        self.delay = t.copy()
        self.raw = A.copy()
        self.data = A.copy()
        self.bg_applied = False
        self.bg_spectrum = None
        self.bg_n = 0
        self.chirp_applied = False
        self.chirp_points = np.zeros((0, 2))
        self.chirp_params = None
        self.chirp_fit_rms = np.nan
        self.ridge_wl = np.zeros(0)
        self.ridge_t = np.zeros(0)
        self.t_zero_shift = 0.0
        self.masked_regions = []
        self.cropped = False

    # ------------------------------------------------------------ loading
    @classmethod
    def from_file(cls, path, transpose=False):
        """Load a data file. transpose=True interchanges the two axes
        (file stored as rows = delays, columns = wavelengths)."""
        wl, t, A = read_data_file(path)
        if transpose:
            wl, t, A = t, wl, A.T
            iw = np.argsort(wl, kind='stable')
            it = np.argsort(t, kind='stable')
            wl, t, A = wl[iw], t[it], A[iw, :][:, it]
        return cls(wl, t, A, source=str(path))

    @classmethod
    def from_average(cls, paths, shifts=None, transpose=False):
        """Average several scans of the same size (Load & Average).

        shifts : time shift of every file (the Shift column); file i is read
                 at its own time  t - shift_i. NaN where a shifted file has
                 no data, so the mean is taken over the files that do.
        """
        sets = [p if isinstance(p, cls) else cls.from_file(p, transpose) for p in paths]
        if not sets:
            raise ValueError('Check at least one file.')
        wl, t = sets[0].wavelength, sets[0].delay
        shifts = np.zeros(len(sets)) if shifts is None else np.asarray(shifts, dtype=float)
        stack = np.full((wl.size, t.size, len(sets)), np.nan)
        for i, ds in enumerate(sets):
            if ds.raw.shape != (wl.size, t.size):
                raise ValueError('All files must have the same size.')
            if shifts[i] == 0 and np.array_equal(ds.delay, t) and np.array_equal(ds.wavelength, wl):
                stack[:, :, i] = ds.raw
                continue
            tq = t - shifts[i]
            for j in range(wl.size):
                stack[j, :, i] = interp1_nan(ds.delay, ds.raw[j, :], tq)
        with np.errstate(invalid='ignore'):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', category=RuntimeWarning)
                avg = np.nanmean(stack, axis=2)
        return cls(wl, t, avg, source='Avg of %d files' % len(sets))

    @property
    def axes_look_swapped(self):
        return detect_axis_swap(self.wavelength, self.delay)

    # ------------------------------------------------------ preprocessing
    def recompute(self):
        """Rebuild the working data from the raw data and the settings."""
        data = self.raw.copy()
        if self.bg_applied and self.bg_spectrum is not None:
            data = data - self.bg_spectrum[:, None]
        if self.chirp_applied and self.chirp_params is not None:
            t = self.delay
            t0 = chirp_model(self.chirp_params, self.wavelength)
            new = np.full(data.shape, np.nan)
            for i in range(self.wavelength.size):
                new[i, :] = interp1_nan(t, data[i, :], t + t0[i])
            data = new
        for lo, hi, mode in self.masked_regions:
            rows = (self.wavelength >= lo) & (self.wavelength <= hi)
            if not rows.any():
                continue
            data[rows, :] = np.nan if mode == 'nan' else 0.0
        self.data = data
        return self

    def chirp_base(self):
        """Raw data with the background removed (what the chirp window shows)."""
        base = self.raw.copy()
        if self.bg_applied and self.bg_spectrum is not None:
            base = base - self.bg_spectrum[:, None]
        return base

    def subtract_background(self, n_points=5):
        """Subtract the mean spectrum of the first n_points delays."""
        n = int(min(max(n_points, 1), self.delay.size))
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=RuntimeWarning)
            self.bg_spectrum = np.nanmean(self.raw[:, :n], axis=1)
        self.bg_n = n
        self.bg_applied = True
        return self.recompute()

    def compute_ridge(self, method='max|dA/dt|', t_min=-0.5, t_max=1.5, smooth_n=5):
        self.ridge_wl, self.ridge_t = compute_ridge(self.chirp_base(), self.wavelength,
                                                    self.delay, method, t_min, t_max, smooth_n)
        return self.ridge_wl, self.ridge_t

    def auto_place_chirp_points(self, n=8, **ridge_kw):
        """n points on the ridge, evenly spaced in wavelength (Auto-place)."""
        if self.ridge_t.size == 0 or ridge_kw:
            self.compute_ridge(**ridge_kw)
        if self.ridge_t.size == 0:
            raise RuntimeError('Could not compute ridge (check search window).')
        wls = np.linspace(self.wavelength[0], self.wavelength[-1], n + 2)[1:-1]
        pts = np.zeros((n, 2))
        for k in range(n):
            j = int(np.argmin(np.abs(self.ridge_wl - wls[k])))
            pts[k, :] = [self.ridge_wl[j], self.ridge_t[j]]
        self.chirp_points = pts
        return pts

    def correct_chirp(self, points=None):
        """Fit the chirp formula to the points (wavelength, time zero) and
        apply the correction (Fit & Apply)."""
        if points is not None:
            self.chirp_points = np.asarray(points, dtype=float)
        self.chirp_params, self.chirp_fit_rms = fit_chirp_params(self.chirp_points)
        self.chirp_applied = True
        return self.recompute()

    def set_time_zero(self, new_zero):
        """Shift the delay axis so that new_zero becomes 0 (Set t=0)."""
        delta = float(new_zero)
        if not np.isfinite(delta) or delta == 0:
            return self
        self.delay = self.delay - delta
        self.t_zero_shift += delta
        return self

    def mask(self, wl_from, wl_to, mode='nan'):
        """Exclude a wavelength range: mode 'nan' (not used by the fits)
        or 'zero'."""
        lo, hi = (wl_from, wl_to) if wl_from <= wl_to else (wl_to, wl_from)
        if mode not in ('nan', 'zero'):
            raise ValueError("mode must be 'nan' or 'zero'.")
        self.masked_regions.append((float(lo), float(hi), mode))
        return self.recompute()

    def clear_masks(self):
        self.masked_regions = []
        return self.recompute()

    def crop(self, wl_range=None, t_range=None):
        """Keep the given wavelength and delay ranges. The current working
        data become the new raw data; the individual corrections can no
        longer be undone separately (reset() still works)."""
        wl_lo, wl_hi = (self.wavelength.min(), self.wavelength.max()) if wl_range is None else sorted(wl_range)
        t_lo, t_hi = (self.delay.min(), self.delay.max()) if t_range is None else sorted(t_range)
        wm = (self.wavelength >= wl_lo) & (self.wavelength <= wl_hi)
        tm = (self.delay >= t_lo) & (self.delay <= t_hi)
        if not wm.any() or not tm.any():
            raise ValueError('Empty crop range.')
        if np.count_nonzero(wm) < 2 or np.count_nonzero(tm) < 2:
            raise ValueError('The crop range must keep at least 2 points on each axis.')
        cropped = self.data[np.ix_(wm, tm)].copy()
        wl, t = self.wavelength[wm], self.delay[tm]
        self._load(wl, t, cropped)
        self.cropped = True
        return self

    def reset(self):
        """Back to the data as loaded (Reset all)."""
        self._load(*self.original)
        return self

    # ------------------------------------------------------------ helpers
    def kinetics(self, wavelength):
        i = int(np.argmin(np.abs(self.wavelength - wavelength)))
        return self.delay, self.data[i, :], self.wavelength[i]

    def spectrum(self, delay):
        j = int(np.argmin(np.abs(self.delay - delay)))
        return self.wavelength, self.data[:, j], self.delay[j]

    def save(self, path, delimiter=None):
        write_data_file(path, self.wavelength, self.delay, self.data, delimiter)

    def __repr__(self):
        return 'Dataset(%d wavelengths x %d delays, %s)' % (
            self.wavelength.size, self.delay.size, self.source or 'in memory')
