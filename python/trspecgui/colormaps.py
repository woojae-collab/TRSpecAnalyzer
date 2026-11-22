"""Color maps of the 2D maps (port of getColormap / divergingColormap).

MATLAB's 'parula' is replaced by matplotlib's 'viridis'."""
import numpy as np
from matplotlib import colormaps
from matplotlib.colors import ListedColormap

# names offered in the menus
NAMES = ['turbo', 'viridis', 'jet', 'hot', 'cool', 'gray',
         'RdBu', 'BWR', 'seismic']


def normName(name):
    """Menu name for a stored name ('parula' of a MATLAB session
    becomes 'viridis'); 'turbo' for anything unknown."""
    if not isinstance(name, str):
        return 'turbo'
    if name.lower() == 'parula':
        return 'viridis'
    for n in NAMES:
        if n.lower() == name.lower():
            return n
    return 'turbo'


def divergingColormap(N, c1, cMid, c2):
    half = N // 2
    t1 = np.linspace(0, 1, half)[:, None]
    t2 = np.linspace(0, 1, N - half)[:, None]
    c1, cMid, c2 = (np.asarray(c, float)[None, :] for c in (c1, cMid, c2))
    part1 = (1 - t1) * c1 + t1 * cMid
    part2 = (1 - t2) * cMid + t2 * c2
    return np.vstack([part1, part2])


def getColormap(name, N=256):
    """Colormap object for a menu name; NaN is drawn transparent."""
    key = normName(name).lower()
    if key == 'rdbu':
        cm = ListedColormap(divergingColormap(
            N, [0.02, 0.19, 0.38], [1, 1, 1], [0.40, 0, 0.12]), 'RdBu_trs')
    elif key == 'bwr':
        cm = ListedColormap(divergingColormap(
            N, [0, 0, 1], [1, 1, 1], [1, 0, 0]), 'BWR_trs')
    elif key == 'seismic':
        cm = ListedColormap(divergingColormap(
            N, [0, 0, 0.30], [1, 1, 1], [0.55, 0, 0]), 'seismic_trs')
    else:
        cm = colormaps[key].resampled(N).copy()
    cm.set_bad((0, 0, 0, 0))
    return cm
