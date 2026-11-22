"""Keeps the text of the plots readable whatever the size of the window
and whatever numbers are shown.

Before a canvas is drawn, fitCanvas() goes through its axes:
  - a title or an axis label that does not fit the room it has gets
    smaller letters and, if that is not enough, a shorter wording;
  - a title is moved aside from the 'x10^n' factor of the y-axis
    instead of on top of it;
  - the 'x10^n' factor of a color bar is written to the right of the
    bar's left edge, where no title can be;
  - axes stacked on top of (or next to) each other share one axis label;
  - a legend that is taller than its axes is set in several columns.
Legends made with legend() can be dragged with the mouse; the place
chosen by the user is kept when the plot is drawn again."""
import re

from PySide6 import QtCore

MIN_PT = 7.0            # smallest letters used to make a text fit
STEP = 0.5
LEG_PT = 8.0            # letters of a legend
LEG_MIN_PT = 7.0


# ----------------------------------------------------------------------
# shorter wordings
# ----------------------------------------------------------------------
def shorterLabel(text):
    """Shorter wordings of an axis label, best first."""
    out = []
    m = re.match(r'^(Amplitude|Residual) \((.*)\)$', text)
    if m:
        out.append(m.group(1))
    if text == 'Normalized amplitude':
        out += ['Norm. amplitude', 'Normalized']
    if text.startswith('Delay time ('):
        out.append('Delay (' + text[len('Delay time ('):])
    # (Auto Analysis window)
    if text == 'Noise (same unit as the data)':
        out.append('Noise')
    if text.startswith('RMS of the '):
        out += ['RMS of the fit', 'RMS']
    if text.endswith(' summed over the wavelengths (noise-weighted)'):
        base = text[:-len(' summed over the wavelengths (noise-weighted)')]
        out += [base + ' summed over the wavelengths', base + ', summed']
    if text.endswith(' (each scaled to 1)'):
        out.append(text[:-len(' (each scaled to 1)')] + ' (scaled)')
    return out


def setTitle(ax, text, *shorter, size=10.0, bold=False):
    """Title of ax; the shorter wordings are used, in this order, when
    the full one does not fit even with small letters."""
    ax.set_title(text, fontsize=size,
                 fontweight='bold' if bold else 'normal')
    ax._tfTitle = {'texts': [text] + [s for s in shorter if s],
                   'size': float(size), 'shown': None}


def setSpan(ax, xLeft, xRight):
    """The part of the canvas width (in points) the title and the
    x-label of ax may use. Without it: up to half way to the next
    axes."""
    ax._tfSpan = (float(xLeft), float(xRight))


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def isColorbar(ax):
    return getattr(ax, '_colorbar', None) is not None or \
        'outline' in ax.spines


def _info(ax, key, artist, alternatives):
    """Texts and letter size the artist is meant to have. If the
    program changed the artist since the last fit, its present state
    is taken as the intended one."""
    cur = artist.get_text()
    info = ax.__dict__.get(key)
    state = (cur, float(artist.get_fontsize()))
    size = state[1]
    if info is not None:
        shown = info.get('shown')
        if shown is None and info['texts'] and info['texts'][0] == cur:
            return info                           # set by setTitle
        if shown == state:
            return info                           # untouched since
        if shown is not None and shown[1] == state[1]:
            size = info['size']       # new text, letters not touched
    info = {'texts': [cur] + list(alternatives(cur)),
            'size': size, 'shown': None}
    setattr(ax, key, info)
    return info


def _choose(artist, info, room, measure):
    """Give the artist a text that fits into room (points): the full
    wording with letters up to one point smaller; failing that a
    shorter wording; failing that the smallest letters. Returns its
    length."""
    length = 0.0
    for floor in (max(MIN_PT, info['size'] - 1.0), MIN_PT):
        for txt in info['texts']:
            artist.set_text(txt)
            s = info['size']
            while True:
                artist.set_fontsize(s)
                length = measure()
                if length <= room or s <= floor:
                    break
                s = max(floor, s - STEP)
            if length <= room:
                break
        if length <= room:
            break
    info['shown'] = (artist.get_text(), float(artist.get_fontsize()))
    return length


def _offsetText(axis):
    """Bring the 'x10^n' text of an axis up to date (matplotlib does
    it only while drawing) and return it ('' if there is none)."""
    t = axis.get_offset_text()
    try:
        fmt = axis.get_major_formatter()
        fmt.set_locs(axis.get_majorticklocs())
        s = fmt.get_offset()
    except Exception:
        s = ''
    t.set_text(s)
    return s if t.get_visible() else ''


def _neighbors(ax, plots, k, W, H):
    """Default room of ax: (xLeft, xRight, yBottom, yTop) in points,
    reaching half way to the axes next to it."""
    b = ax.get_window_extent()
    xl, xr, yb, yt = 0.0, W, 0.0, H
    for o in plots:
        if o is ax:
            continue
        c = o.get_window_extent()
        if min(b.y1, c.y1) - max(b.y0, c.y0) > 0:        # same row
            if c.x1 <= b.x0 + 1:
                xl = max(xl, 0.5 * (c.x1 + b.x0) / k)
            elif c.x0 >= b.x1 - 1:
                xr = min(xr, 0.5 * (c.x0 + b.x1) / k)
        if min(b.x1, c.x1) - max(b.x0, c.x0) > 0:        # same column
            if c.y1 <= b.y0 + 1:
                yb = max(yb, 0.5 * (c.y1 + b.y0) / k)
            elif c.y0 >= b.y1 - 1:
                yt = min(yt, 0.5 * (c.y0 + b.y1) / k)
    return xl, xr, yb, yt


# ----------------------------------------------------------------------
# the fit
# ----------------------------------------------------------------------
def fitCanvas(cv):
    fig = cv.figure
    W, H = cv.sizePt()
    if W < 30 or H < 30:
        return
    r = cv.get_renderer()
    k = fig.dpi / 72.0                       # pixels per point
    axes = [a for a in fig.axes if a.get_visible()]
    plots = []
    for ax in axes:
        if isColorbar(ax):
            # 'x10^n' of a color bar: from the bar's left edge to the
            # right (matplotlib writes it to the left, over the map)
            t = ax.yaxis.get_offset_text()
            t.set_ha('left')
            t.set_x(0.0)
            ax.yaxis.OFFSETTEXTPAD = 7.0     # above the top tick label
        else:
            plots.append(ax)

    # axes that share one label
    yOwner, xOwner = {}, {}
    for grp, owner in ((getattr(cv, 'yGroups', []), yOwner),
                       (getattr(cv, 'xGroups', []), xOwner)):
        for g in grp:
            g = [a for a in g if a in plots]
            if len(g) > 1:
                for a in g:
                    owner[a] = g

    for ax in plots:
        b = ax.get_window_extent()
        aL, aR, aB, aT = b.x0 / k, b.x1 / k, b.y0 / k, b.y1 / k
        if aR - aL < 2 or aT - aB < 2:
            continue
        xl, xr, yb, yt = _neighbors(ax, plots, k, W, H)
        span = ax.__dict__.get('_tfSpan')
        if span is not None:
            xl, xr = span

        # ---- title ---------------------------------------------------
        t = ax.title
        if t.get_text() and t.get_visible():
            info = _info(ax, '_tfTitle', t, lambda s: [])
            left, right = xl + 2.0, xr - 2.0
            if _offsetText(ax.yaxis):
                wOff = ax.yaxis.get_offset_text() \
                    .get_window_extent(r).width / k
                left = max(left, aL + wOff + 5.0)
            w = _choose(t, info, right - left,
                        lambda: t.get_window_extent(r).width / k)
            c = 0.5 * (aL + aR)
            if w > right - left:
                c = left + 0.5 * w
            else:
                c = min(max(c, left + 0.5 * w), right - 0.5 * w)
            t.set_x((c - aL) / (aR - aL))

        # ---- y label -------------------------------------------------
        lab = ax.yaxis.label
        grp = yOwner.get(ax)
        if grp is not None and grp[0] is not ax and \
                lab.get_text() == grp[0].yaxis.label.get_text():
            lab.set_visible(False)
            ax._tfHidY = True
        else:
            if ax.__dict__.get('_tfHidY'):
                lab.set_visible(True)
                ax._tfHidY = False
            if lab.get_text() and lab.get_visible():
                lo, hi, cB, cT = yb, yt, aB, aT
                if grp is not None:
                    for o in grp[1:]:
                        if o.yaxis.label.get_text() != lab.get_text():
                            continue
                        c = o.get_window_extent()
                        ol = _neighbors(o, plots, k, W, H)
                        lo, hi = min(lo, ol[2]), max(hi, ol[3])
                        cB, cT = min(cB, c.y0 / k), max(cT, c.y1 / k)
                # not above the axes: that room belongs to the titles
                hi = min(hi, cT + 3.0)
                info = _info(ax, '_tfYLab', lab, shorterLabel)
                h = _choose(lab, info, hi - lo - 2.0,
                            lambda: lab.get_window_extent(r).height / k)
                c = 0.5 * (cB + cT)
                c = min(max(c, lo + 1.0 + 0.5 * h), hi - 1.0 - 0.5 * h) \
                    if h <= hi - lo - 2.0 else 0.5 * (lo + hi)
                lab.set_y((c - aB) / (aT - aB))

        # ---- x label -------------------------------------------------
        lab = ax.xaxis.label
        grp = xOwner.get(ax)
        if grp is not None and grp[0] is not ax and \
                lab.get_text() == grp[0].xaxis.label.get_text():
            lab.set_visible(False)
            ax._tfHidX = True
        else:
            if ax.__dict__.get('_tfHidX'):
                lab.set_visible(True)
                ax._tfHidX = False
            if lab.get_text() and lab.get_visible():
                lo, hi, cL, cR = xl, xr, aL, aR
                if grp is not None:
                    for o in grp[1:]:
                        if o.xaxis.label.get_text() != lab.get_text():
                            continue
                        c = o.get_window_extent()
                        ol = o.__dict__.get('_tfSpan') or \
                            _neighbors(o, plots, k, W, H)[:2]
                        lo, hi = min(lo, ol[0]), max(hi, ol[1])
                        cL, cR = min(cL, c.x0 / k), max(cR, c.x1 / k)
                info = _info(ax, '_tfXLab', lab, shorterLabel)
                w = _choose(lab, info, hi - lo - 2.0,
                            lambda: lab.get_window_extent(r).width / k)
                c = 0.5 * (cL + cR)
                c = min(max(c, lo + 1.0 + 0.5 * w), hi - 1.0 - 0.5 * w) \
                    if w <= hi - lo - 2.0 else 0.5 * (lo + hi)
                lab.set_x((c - aL) / (aR - aL))

        _fitLegend(ax, r, k)


# ----------------------------------------------------------------------
# legends
# ----------------------------------------------------------------------
def legend(ax, handles=None, labels=None, shorter=None, **kw):
    """Legend of ax that can be dragged with the mouse. A place chosen
    by the user is kept when the plot is drawn again; the number of
    columns follows the height of the axes. shorter: a second wording
    of the entries, used when the full one does not fit the axes."""
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    elif labels is None:
        labels = [h.get_label() for h in handles]
    if not handles:
        return None
    kw.setdefault('loc', 'best')
    kw.setdefault('fontsize', LEG_PT)
    kw.setdefault('borderpad', 0.35)
    kw.setdefault('labelspacing', 0.3)
    kw.setdefault('handlelength', 1.8)
    kw.setdefault('handletextpad', 0.5)
    kw.setdefault('columnspacing', 1.0)
    kw.setdefault('borderaxespad', 0.4)
    old = ax.__dict__.get('_tfLeg')
    # columns and letters as found for the last legend of these axes
    # (fitCanvas looks again when the entries or the size change)
    if shorter is not None and len(shorter) != len(labels):
        shorter = None
    ax._tfLeg = {'handles': list(handles), 'labels': list(labels),
                 'shorter': list(shorter) if shorter else None,
                 'short': bool(old.get('short')) if old and shorter else False,
                 'kw': kw, 'key': None, 'leg': None,
                 'ncols': old['ncols'] if old else 1,
                 'size': old['size'] if old else float(kw['fontsize'])}
    return _makeLegend(ax, ax._tfLeg['ncols'], ax._tfLeg['size'])


def _makeLegend(ax, ncols, size):
    sp = ax._tfLeg
    old = ax.get_legend()
    if old is not None:
        try:
            old.set_draggable(False)
        except Exception:
            pass
    kw = dict(sp['kw'])
    kw['ncols'] = max(1, min(int(ncols), len(sp['handles'])))
    kw['fontsize'] = size
    pos = ax.__dict__.get('_tfLegPos')
    if pos is not None:
        kw['loc'] = pos
        # (a legend that started beside the axes: the place chosen by
        # the user counts from the axes themselves)
        kw.pop('bbox_to_anchor', None)
    words = sp['shorter'] if sp.get('short') and sp.get('shorter') \
        else sp['labels']
    leg = ax.legend(sp['handles'], words, **kw)
    sp['leg'], sp['ncols'], sp['size'] = leg, kw['ncols'], size
    cv = ax.figure.canvas
    if getattr(cv, 'legendsHidden', False):
        leg.set_visible(False)
    try:
        d = leg.set_draggable(True, use_blit=True, update='loc')
        done = d.finalize_offset

        def finalize():
            # (called by matplotlib when the mouse button is released)
            try:        # moved, or only clicked?
                ox, oy = d.offsetbox._offset
                moved = abs(ox - d.offsetbox_x) > 1.5 or \
                    abs(oy - d.offsetbox_y) > 1.5
            except Exception:
                moved = True
            before = leg._loc
            done()
            if not moved:
                leg._loc = before         # a click alone changes nothing
            elif isinstance(leg._loc, tuple):
                ax._tfLegPos = _insidePicture(ax, leg)
                leg._loc = ax._tfLegPos
            _redrawLater(ax.figure.canvas)

        d.finalize_offset = finalize
    except Exception:
        pass
    return leg


def _insidePicture(ax, leg):
    """Place of the legend (lower left corner, as a fraction of the
    axes), moved back into the picture if it was dragged over its
    edge, where it could not be seen or taken hold of again."""
    fx, fy = float(leg._loc[0]), float(leg._loc[1])
    try:
        fig = ax.figure
        a = ax.get_window_extent()
        e = leg.get_window_extent(fig.canvas.get_renderer())
        x = min(max(a.x0 + fx * a.width, 0.0),
                max(fig.bbox.width - e.width, 0.0))
        y = min(max(a.y0 + fy * a.height, 0.0),
                max(fig.bbox.height - e.height, 0.0))
        fx, fy = (x - a.x0) / a.width, (y - a.y0) / a.height
    except Exception:
        pass
    return (fx, fy)


def _redrawLater(cv):
    """Draw the canvas again once the present mouse event has been
    dealt with. Not at once: when a drag ends, matplotlib still holds
    the legend as 'animated' (left out of a normal drawing) while it
    calls finalize; a redraw asked for at that moment is carried out by
    its very next blit, without the legend, and the legend is gone
    from the picture."""
    def go():
        try:
            cv.draw_idle()
        except Exception:               # the window was closed meanwhile
            pass
    QtCore.QTimer.singleShot(0, go)


def _fitLegend(ax, r, k):
    sp = ax.__dict__.get('_tfLeg')
    leg = ax.get_legend()
    if leg is not None and leg.get_animated() and not getattr(
            getattr(leg, '_draggable', None), 'got_artist', False):
        # 'animated' is for the time of a drag only; a legend left in
        # that state would not be drawn
        leg.set_animated(False)
    if sp is None or leg is None or leg is not sp['leg']:
        return
    b = ax.get_window_extent()
    aw, ah = b.width / k, b.height / k
    n = len(sp['handles'])
    key = (n, int(aw // 3), int(ah // 3))
    if sp['key'] == key:
        return
    sp['key'] = key
    if not leg.get_visible():
        return
    size0 = float(sp['kw']['fontsize'])
    sizes = [size0] + ([LEG_MIN_PT] if size0 > LEG_MIN_PT else [])
    best = None
    for short in ([False, True] if sp.get('shorter') else [False]):
        for size in sizes:
            for ncols in range(1, min(3, n) + 1):
                if (ncols, size, short) != (sp['ncols'], sp['size'],
                                            bool(sp.get('short'))):
                    sp['short'] = short
                    _makeLegend(ax, ncols, size)
                e = sp['leg'].get_window_extent(r)
                over = max(0.0, e.height / k - 0.88 * ah) / ah + \
                    max(0.0, e.width / k - 0.9 * aw) / aw
                if over == 0:
                    return
                if best is None or over < best[0]:
                    best = (over, ncols, size, short)
    # nothing fits: the form that sticks out least
    if best[1:] != (sp['ncols'], sp['size'], bool(sp.get('short'))):
        sp['short'] = best[3]
        _makeLegend(ax, best[1], best[2])


def setLegendText(ax, index, text):
    """Change one entry of the legend of ax (no new legend is made)."""
    sp = ax.__dict__.get('_tfLeg')
    leg = ax.get_legend()
    if leg is None:
        return
    if sp is not None and leg is sp['leg']:
        sp['labels'][index] = text
        if sp.get('short'):       # (the short wording is on show)
            return
    leg.get_texts()[index].set_text(text)


def legendAt(cv, ev):
    """The legend under the mouse event, or None."""
    for ax in cv.figure.axes:
        leg = ax.get_legend()
        if leg is not None and leg.get_visible() and ax.get_visible():
            try:
                if leg.contains(ev)[0]:
                    return leg
            except Exception:
                pass
    return None


def resetLegends(cv):
    """Forget the places chosen by the user."""
    for ax in cv.figure.axes:
        if ax.__dict__.get('_tfLegPos') is not None:
            ax._tfLegPos = None
            sp = ax.__dict__.get('_tfLeg')
            if sp is not None and ax.get_legend() is sp['leg']:
                _makeLegend(ax, sp['ncols'], sp['size'])
    cv.draw_idle()


def showLegends(cv, on):
    cv.legendsHidden = not on
    for ax in cv.figure.axes:
        leg = ax.get_legend()
        if leg is not None:
            leg.set_visible(bool(on))
            sp = ax.__dict__.get('_tfLeg')
            if sp is not None:
                sp['key'] = None
    cv.draw_idle()


def hasLegend(cv):
    return any(ax.get_legend() is not None for ax in cv.figure.axes)


def hasMovedLegend(cv):
    return any(ax.__dict__.get('_tfLegPos') is not None
               for ax in cv.figure.axes)
