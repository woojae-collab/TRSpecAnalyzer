"""Small shared helpers."""
import numpy as np


class Result(dict):
    """Dictionary whose entries are also attributes (res.tau == res['tau'])."""
    __setattr__ = dict.__setitem__

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            # (AttributeError, so that copy, pickle and hasattr work)
            raise AttributeError(name) from None

    def __dir__(self):
        return list(self.keys())


def lstsq_min_norm(a, b):
    """Minimum-norm least-squares solution of a x = b (MATLAB lsqminnorm)."""
    return np.linalg.lstsq(a, b, rcond=None)[0]


def as_row(x, n=None, dtype=float, name='value'):
    """1-D array; a scalar is repeated n times."""
    x = np.atleast_1d(np.asarray(x, dtype=dtype)).ravel()
    if n is not None and x.size == 1 and n != 1:
        x = np.repeat(x, n)
    if n is not None and x.size != n:
        raise ValueError('%s needs one entry per time constant (%d).' % (name, n))
    return x


def interp1_nan(x, y, xq):
    """Linear interpolation, NaN outside the range of x (MATLAB interp1 with
    'linear', NaN). x must be increasing."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xq = np.asarray(xq, dtype=float)
    return np.interp(xq, x, y, left=np.nan, right=np.nan)


def colon(a, step, b):
    """MATLAB a:step:b (the algorithm C. Moler gives for it). The number
    of intervals n is round((b - a)/step), reduced by one if the last
    point would pass b by more than 2*eps*max(|a|, |b|). The right end
    is b itself when a + n*step lies within that tolerance of it, and
    the points are counted from both ends toward the middle."""
    a, step, b = float(a), float(step), float(b)
    if step == 0 or (step > 0 and a > b) or (step < 0 and a < b):
        return np.zeros(0)
    tol = 2.0 * np.finfo(float).eps * max(abs(a), abs(b))
    sig = 1.0 if step > 0 else -1.0
    if a == np.floor(a) and step == 1:
        n = int(np.floor(b) - a)
    elif a == np.floor(a) and step == np.floor(step):
        q = np.floor(a / step)
        r = a - q * step
        n = int(np.floor((b - r) / step) - q)
    else:
        n = int(np.floor((b - a) / step + 0.5))
        if sig * (a + n * step - b) > tol:
            n -= 1
    c = a + n * step
    if sig * (c - b) > -tol:
        c = b
    v = np.zeros(n + 1)
    k = np.arange(n // 2 + 1)
    v[k] = a + k * step
    v[n - k] = c - k * step
    if n % 2 == 0:
        v[n // 2] = (a + c) / 2
    return v
