"""Zoom and pan of the plots with the mouse.

What the user set is kept per axes in ax._view and put back before
every drawing (apply), so that a zoomed plot stays zoomed when the
program draws it again (another delay, another wavelength, a pinned
curve). The zoom of an axes is dropped when the program gives it
another x-range or another scale (other data, units, display range) --
for maps also another y-range -- and by reset().

The program's own limits are never switched off: user limits are set
with auto=None, so that autoscaling by the program goes on as before
and its result is remembered as the full view ('home').
"""
import numpy as np
from matplotlib.collections import QuadMesh
from matplotlib.ticker import (AutoLocator, FixedLocator, Formatter,
                               FuncFormatter, Locator, LogLocator,
                               MaxNLocator, NullFormatter, ScalarFormatter)

STEP = 1.2            # one notch of the wheel
MIN_FRACTION = 1e-6   # smallest view, as a part of the full one
AXIS_ZONE_PT = (34.0, 58.0)   # room below / left of an axes that
#                               counts as "on the x-axis / y-axis"


# ----------------------------------------------------------------------
# which axes
# ----------------------------------------------------------------------
def zoomable(ax):
    """Axes with data and visible axis lines; not a color bar."""
    return (ax.get_visible() and ax.axison
            and not hasattr(ax, '_colorbar')
            and not ax.__dict__.get('_noZoom')
            and bool(ax.lines or ax.collections or ax.images or ax.patches))


def isMap(ax):
    return bool(ax.images) or any(isinstance(c, QuadMesh)
                                  for c in ax.collections)


def isZoomed(ax):
    z = ax.__dict__.get('_view')
    return z is not None and (z['x'] is not None or z['y'] is not None)


def anyZoomed(cv):
    return any(isZoomed(ax) for ax in cv.figure.axes)


def axesAt(cv, x, y):
    """(axes, mode) under the pixel (x, y) of the canvas: mode 'xy'
    inside an axes, 'x' on the tick labels below it, 'y' on those to
    its left; (None, '') elsewhere."""
    axs = [a for a in cv.figure.axes if zoomable(a)]
    for ax in axs:
        if ax.bbox.contains(x, y):
            return ax, 'xy'
    k = cv.figure.dpi / 72.0
    zx, zy = AXIS_ZONE_PT[0] * k, AXIS_ZONE_PT[1] * k
    best = (None, '', np.inf)
    for ax in axs:
        b = ax.bbox
        if b.x0 <= x <= b.x1 and b.y0 - zx <= y < b.y0:
            d = b.y0 - y
            if d < best[2]:
                best = (ax, 'x', d)
        if b.y0 <= y <= b.y1 and b.x0 - zy <= x < b.x0:
            d = b.x0 - x
            if d < best[2]:
                best = (ax, 'y', d)
    return best[0], best[1]


def _same(a, b):
    w = max(abs(a[1] - a[0]), abs(b[1] - b[0]), 1e-300)
    return abs(a[0] - b[0]) <= 1e-9 * w and abs(a[1] - b[1]) <= 1e-9 * w


def _home(ax):
    z = ax.__dict__.get('_view')
    if z is not None:
        return z['home']
    return [tuple(ax.get_xlim()), tuple(ax.get_ylim())]


def partners(cv, ax, which):
    """Axes that share a range with ax: which = 'y' -> side by side
    with the same y-range (the linear and the logarithmic part of a
    kinetics plot); 'x' -> above each other with the same x-range (the
    two panels of a split map). Axes named together in cv.zoomLinks
    (the data, fit and residual maps of an analysis window) share
    both. Followed through: the partner of a partner is one too."""
    seen, todo = [ax], [ax]
    while todo:
        a = todo.pop()
        for o in _direct(cv, a, which):
            if not any(o is q for q in seen):
                seen.append(o)
                todo.append(o)
    return seen[1:]


def _direct(cv, ax, which):
    out = []

    def add(o):
        if o is not ax and zoomable(o) and not any(o is q for q in out):
            out.append(o)

    for g in (cv.xGroups if which == 'y' else cv.yGroups):
        if any(a is ax for a in g):
            for a in g:
                add(a)
    for g in getattr(cv, 'zoomLinks', []):
        if any(a is ax for a in g):
            for a in g:
                add(a)
    # (axes that matplotlib itself keeps on one range)
    sh = ax.get_shared_y_axes() if which == 'y' else ax.get_shared_x_axes()
    for a in sh.get_siblings(ax):
        add(a)
    if out:
        return out
    # by their places: neighbors with the same range
    W, H = cv.sizePt()
    if W <= 0 or H <= 0:
        return []
    b = ax.get_position()
    i = 1 if which == 'y' else 0
    scale = ax.get_yscale() if which == 'y' else ax.get_xscale()
    for o in cv.figure.axes:
        if o is ax or not zoomable(o):
            continue
        p = o.get_position()
        if which == 'y':
            band = abs(p.y0 - b.y0) * H < 1.5 and abs(p.y1 - b.y1) * H < 1.5
            gap = max(p.x0 - b.x1, b.x0 - p.x1) * W
            ok = band and -1.0 <= gap <= 16.0 and o.get_yscale() == scale
        else:
            band = abs(p.x0 - b.x0) * W < 1.5 and abs(p.x1 - b.x1) * W < 1.5
            gap = max(p.y0 - b.y1, b.y0 - p.y1) * H
            ok = band and -1.0 <= gap <= 40.0 and o.get_xscale() == scale
        if ok and _same(_home(o)[i], _home(ax)[i]):
            out.append(o)
    return out


# ----------------------------------------------------------------------
# ticks of a zoomed axis
# ----------------------------------------------------------------------
class ZoomLogLocator(Locator):
    """Ticks of a logarithmic axis the user zoomed into: the densest of
    1-2-3-5-7, 1-2-5, 1-3, decades whose labels do not touch, for a
    view of a decade and more; round numbers as on a linear axis for a
    narrower one."""

    def __call__(self):
        return self.tick_values(*self.axis.get_view_interval())

    def room(self):
        """(length of the axis, room one label of n characters needs),
        in pixels."""
        try:
            ax = self.axis.axes
            px = self.axis.majorTicks[0].label1.get_size() * \
                ax.figure.dpi / 72.0
            if self.axis.axis_name == 'y':
                return ax.bbox.height, lambda n: 1.5 * px
            return ax.bbox.width, lambda n: (0.62 * n + 1.0) * px
        except Exception:
            return 300.0, lambda n: 40.0

    def tick_values(self, vmin, vmax):
        lo, hi = sorted((vmin, vmax))
        if not (lo > 0 and hi > lo):
            return np.array([])
        length, need = self.room()
        dec = np.log10(hi / lo)
        k0, k1 = int(np.floor(np.log10(lo))), int(np.ceil(np.log10(hi)))

        def inView(t):
            t = np.asarray(t, float)
            return t[(t >= lo * (1 - 1e-12)) & (t <= hi * (1 + 1e-12))]

        def fits(t):
            if t.size < 2:
                return False
            gap = length * float(np.min(np.diff(np.log10(t)))) / dec
            return gap >= need(max(len(_plain(v)) for v in t))

        if dec >= 0.7:
            sets = [[q * 10.0 ** k for k in range(k0, k1 + 1) for q in m]
                    for m in ((1, 2, 3, 5, 7), (1, 2, 5), (1, 3), (1,))]
            sets += [list(10.0 ** np.arange(k0, k1 + 1, st))
                     for st in (2, 3, 4, 6)]
            for t in sets:
                t = inView(t)
                if fits(t):
                    return t
        t = np.array([])
        for nb in (5, 4, 3, 2):
            t = inView(MaxNLocator(nbins=nb, steps=[1, 2, 2.5, 5, 10])
                       .tick_values(lo, hi))
            t = t[t > 0]
            if fits(t):
                return t
        while t.size > 1 and not fits(t):   # (a very short axis)
            t = t[::2]
        return t


class ZoomAutoLocator(AutoLocator):
    """Automatic ticks in place of fixed ones that the zoomed view no
    longer holds."""


def _plain(v, _pos=None):
    return '%g' % float('%.6g' % v)


class ZoomLogFormatter(Formatter):
    """Plain numbers on a zoomed logarithmic y-axis. Where one of them
    would need more than five characters (300000, 0.0002), all are
    written with a common factor (x 10^n) above the axis, as on the
    linear axes: wide numbers push the label of the axis off the
    picture."""

    def __init__(self):
        self._exp = 0

    def set_locs(self, locs):
        self.locs = locs
        self._exp = 0
        v = np.asarray(locs, float)
        v = v[np.isfinite(v) & (v > 0)]
        if self.axis is not None and v.size:
            lo, hi = sorted(self.axis.get_view_interval())
            v = v[(v >= lo * (1 - 1e-9)) & (v <= hi * (1 + 1e-9))]
        if v.size and max(len(_plain(x)) for x in v) > 5:
            self._exp = int(np.floor(np.log10(float(np.min(v))) + 1e-9))

    def __call__(self, v, pos=None):
        return _plain(v / 10.0 ** self._exp)

    def get_offset(self):
        if self._exp == 0:
            return ''
        return r'$\times\mathdefault{10^{%d}}$' % self._exp


def plainLogY(ax):
    """Plain numbers (1.37, 1.38) at the ticks of a logarithmic y-axis,
    as on a zoomed one. For values that stay within a decade: there
    matplotlib writes every tick as 1.375 x 10^0, which is wider than
    the room left of the plot, so that the numbers were cut off and
    the label of the axis fell off the picture."""
    axis = ax.yaxis
    axis.set_major_locator(ZoomLogLocator())
    axis.set_major_formatter(ZoomLogFormatter())
    axis.get_offset_text().set_size(max(axis.label.get_size() - 1.0, 6.0))
    axis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)))
    axis.set_minor_formatter(NullFormatter())
    z = ax.__dict__.get('_view')
    if z is not None:           # (ticks kept aside from an earlier drawing:
        z['ticks'].pop('y', None)    # a reset of the zoom leaves these)


def _mine(axis):
    return isinstance(axis.get_major_locator(),
                      (ZoomLogLocator, ZoomAutoLocator))


def _ticks(ax, z):
    """Ticks that suit the zoomed view; those of the program are kept
    aside and come back with the full view."""
    for key, axis in (('x', ax.xaxis), ('y', ax.yaxis)):
        mine = _mine(axis)
        if z[key] is not None:
            if mine:
                continue
            isLog = axis.get_scale() == 'log'
            loc = axis.get_major_locator()
            need = isLog
            if not need and isinstance(loc, FixedLocator):
                lo, hi = sorted(z[key])
                locs = np.asarray(loc.locs, float)
                need = int(np.sum((locs >= lo) & (locs <= hi))) < 2
            if not need:
                continue
            z['ticks'][key] = (loc, axis.get_major_formatter(),
                               axis.get_minor_locator(),
                               axis.get_minor_formatter())
            if isLog:
                axis.set_major_locator(ZoomLogLocator())
                if key == 'y':
                    axis.set_major_formatter(ZoomLogFormatter())
                    ot = axis.get_offset_text()
                    ot.set_size(max(axis.label.get_size() - 1.0, 6.0))
                else:
                    axis.set_major_formatter(FuncFormatter(_plain))
                axis.set_minor_locator(
                    LogLocator(base=10.0, subs=np.arange(2, 10)))
                axis.set_minor_formatter(NullFormatter())
            else:
                axis.set_major_locator(ZoomAutoLocator())
                axis.set_major_formatter(ScalarFormatter(useOffset=False))
        elif mine:
            _restoreTicks(axis, z['ticks'].pop(key, None))


def _restoreTicks(axis, saved):
    if saved is None:
        return
    axis.set_major_locator(saved[0])
    axis.set_major_formatter(saved[1])
    axis.set_minor_locator(saved[2])
    axis.set_minor_formatter(saved[3])


# ----------------------------------------------------------------------
# the kept view
# ----------------------------------------------------------------------
def _state(ax):
    z = ax.__dict__.get('_view')
    if z is None:
        lim = [tuple(ax.get_xlim()), tuple(ax.get_ylim())]
        z = {'x': None, 'y': None, 'home': list(lim), 'applied': list(lim),
             'scales': (ax.get_xscale(), ax.get_yscale()), 'ticks': {}}
        ax._view = z
    return z


def _drop(ax, restoreLimits):
    """Forget the user's view of ax. restoreLimits: put the full view
    back (False when the program has just set its own limits)."""
    z = ax.__dict__.pop('_view', None)
    if z is None:
        return
    for key, axis in (('x', ax.xaxis), ('y', ax.yaxis)):
        if _mine(axis):
            _restoreTicks(axis, z['ticks'].get(key))
    if restoreLimits:
        ax.set_xlim(z['home'][0], auto=None)
        ax.set_ylim(z['home'][1], auto=None)


def apply(cv):
    """Before a drawing: note what the program set since the last one
    (its limits are the full view), then put the user's view back."""
    live = []
    # 1. what did the program do (all axes first: setting the limits
    #    of one axes can move another that shares the range)
    for ax in cv.figure.axes:
        z = ax.__dict__.get('_view')
        if z is None:
            continue
        if not zoomable(ax) or \
                (ax.get_xscale(), ax.get_yscale()) != z['scales']:
            _drop(ax, False)
            continue
        cur = [tuple(ax.get_xlim()), tuple(ax.get_ylim())]
        gone = False
        for i in (0, 1):
            if _same(cur[i], z['applied'][i]):
                continue
            # the program set this range anew
            if not _same(cur[i], z['home'][i]):
                if i == 0:
                    gone = True         # another x-range: other content
                elif isMap(ax):
                    z['y'] = None       # another delay range of a map
            z['home'][i] = cur[i]
        if gone or (z['x'] is None and z['y'] is None):
            _drop(ax, False)
            continue
        live.append((ax, z))
    # 2. the user's view
    for ax, z in live:
        ax.set_xlim(z['x'] or z['home'][0], auto=None)
        ax.set_ylim(z['y'] or z['home'][1], auto=None)
        _ticks(ax, z)
    for ax, z in live:
        z['applied'] = [tuple(ax.get_xlim()), tuple(ax.get_ylim())]


def reset(cv, ax=None):
    """Back to the full view: of ax and the axes that share a range
    with it, or of every axes of the canvas."""
    if ax is None:
        todo = list(cv.figure.axes)
    else:
        todo = [ax] + partners(cv, ax, 'x') + partners(cv, ax, 'y')
    done = False
    for a in todo:
        if a.__dict__.get('_view') is not None:
            _drop(a, True)
            done = True
    return done


# ----------------------------------------------------------------------
# changing the view
# ----------------------------------------------------------------------
def _fwd(axis, v):
    return axis.get_transform().transform(np.asarray(v, float))


def _inv(axis, s):
    return axis.get_transform().inverted().transform(np.asarray(s, float))


def _setAxis(cv, ax, i, sLo, sHi, link=True):
    """Give axis i (0 x, 1 y) of ax the view [sLo, sHi] (in the
    coordinates of its scale), kept inside the full view; the full
    view itself means: not zoomed. Axes that share the range follow."""
    z = _state(ax)
    axis = ax.xaxis if i == 0 else ax.yaxis
    key = 'x' if i == 0 else 'y'
    home = z['home'][i]
    h = _fwd(axis, home)
    if not np.all(np.isfinite(h)) or h[0] == h[1]:
        return
    h0, h1 = sorted(h)
    lo, hi = sorted((float(sLo), float(sHi)))
    w = max(hi - lo, MIN_FRACTION * (h1 - h0))
    c = 0.5 * (lo + hi)
    lo, hi = c - 0.5 * w, c + 0.5 * w
    if w >= (h1 - h0) * (1 - 1e-9):
        z[key] = None
    else:
        if lo < h0:
            lo, hi = h0, h0 + w
        if hi > h1:
            lo, hi = h1 - w, h1
        d = _inv(axis, [lo, hi])
        if home[0] > home[1]:
            d = d[::-1]
        z[key] = (float(d[0]), float(d[1]))
    if link:
        for o in partners(cv, ax, key):
            zo = _state(o)
            ho = _fwd(o.xaxis if i == 0 else o.yaxis, zo['home'][i])
            if _same(tuple(sorted(ho)), (h0, h1)):
                zo[key] = z[key]


def _view(ax, i):
    """Present view of axis i in the coordinates of its scale."""
    axis = ax.xaxis if i == 0 else ax.yaxis
    z = ax.__dict__.get('_view')
    if z is not None:
        lim = z['x' if i == 0 else 'y'] or z['home'][i]
    else:
        lim = ax.get_xlim() if i == 0 else ax.get_ylim()
    return sorted(_fwd(axis, lim))


def _scaled(ax, x, y):
    """Pixel -> coordinates of the scales of ax (as they are drawn)."""
    fx = (x - ax.bbox.x0) / max(ax.bbox.width, 1e-9)
    fy = (y - ax.bbox.y0) / max(ax.bbox.height, 1e-9)
    out = []
    for i, f in ((0, fx), (1, fy)):
        axis = ax.xaxis if i == 0 else ax.yaxis
        lim = ax.get_xlim() if i == 0 else ax.get_ylim()
        a, b = _fwd(axis, lim)
        out.append(a + f * (b - a))
    return out


def zoomAt(cv, ax, x, y, factor, mode='xy'):
    """Zoom ax by factor (< 1: in) about the pixel (x, y)."""
    apply(cv)
    c = _scaled(ax, x, y)
    for i, key in ((0, 'x'), (1, 'y')):
        if key not in mode:
            continue
        lo, hi = _view(ax, i)
        _setAxis(cv, ax, i, c[i] + (lo - c[i]) * factor,
                 c[i] + (hi - c[i]) * factor)
    _settle(cv, ax)


def panStart(cv, ax):
    apply(cv)
    return {'ax': ax, 'view': [_view(ax, 0), _view(ax, 1)],
            'lim': [tuple(ax.get_xlim()), tuple(ax.get_ylim())],
            'box': (ax.bbox.width, ax.bbox.height)}


def panBy(cv, st, dx, dy):
    """Move the view of a pan started with panStart by (dx, dy) pixels.
    An axis that shows its full range stays where it is."""
    ax = st['ax']
    for i, d in ((0, dx), (1, dy)):
        lo, hi = st['view'][i]
        axis = ax.xaxis if i == 0 else ax.yaxis
        h = _fwd(axis, _home(ax)[i])
        if (hi - lo) >= abs(h[1] - h[0]) * (1 - 1e-9):
            continue
        a, b = _fwd(axis, st['lim'][i])
        s = (b - a) * d / max(st['box'][i], 1e-9)
        _setAxis(cv, ax, i, lo - s, hi - s)
    _settle(cv, ax)


def zoomTo(cv, ax, x0, y0, x1, y1, mode='xy'):
    """Zoom ax into the rectangle between two pixels."""
    apply(cv)
    b = ax.bbox
    x0, x1 = (min(max(v, b.x0), b.x1) for v in (x0, x1))
    y0, y1 = (min(max(v, b.y0), b.y1) for v in (y0, y1))
    p, q = _scaled(ax, x0, y0), _scaled(ax, x1, y1)
    for i, key in ((0, 'x'), (1, 'y')):
        if key in mode and p[i] != q[i]:
            _setAxis(cv, ax, i, p[i], q[i])
    _settle(cv, ax)


def _settle(cv, ax):
    """After a change: drop states that hold nothing, and show it."""
    for a in [ax] + partners(cv, ax, 'x') + partners(cv, ax, 'y'):
        z = a.__dict__.get('_view')
        if z is not None and z['x'] is None and z['y'] is None:
            _drop(a, True)
    apply(cv)
