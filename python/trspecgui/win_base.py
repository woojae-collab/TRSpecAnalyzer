"""Base of the separate windows (corrections, analyses)."""
import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from .textfit import setSpan, setTitle
from .widgets import (AX_BG, FONT_PT, BlitManager, FewTicks, PlotCanvas,
                      sciTicks, styleAxes)

GREEN = (0.30, 0.65, 0.35)
BLUE = (0.30, 0.55, 0.80)
PURPLE = (0.55, 0.30, 0.70)

# windows that rest on the working data: closed when the data change
ANALYSIS_KEYS = ('tf', 'cohcmp', 'coherence', 'ga', 'kinfit', 'svd', 'lda', 'mcr',
                 'target', 'auto', 'bandpos')


def fitSize(w, h):
    """Window size, reduced to what the screen offers (port of
    fitFigurePosition)."""
    try:
        g = QtGui.QGuiApplication.primaryScreen().availableGeometry()
        return min(w, g.width() - 40), min(h, g.height() - 100)
    except Exception:
        return w, h


class ChildWindow(QtWidgets.QWidget):
    """A separate window of the program. One window per KEY: opening
    it again brings the existing one to the front."""

    KEY = ''

    def __init__(self, app, title, w, h):
        super().__init__(None)
        self.app = app
        self.setWindowTitle(title)
        self.resize(*fitSize(w, h))
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        app.childWindows[self.KEY] = self

    def closeEvent(self, ev):
        if self.app.childWindows.get(self.KEY) is self:
            del self.app.childWindows[self.KEY]
        self.onClosed()
        ev.accept()

    def onClosed(self):
        pass

    def alert(self, message, title='TRSpecAnalyzer'):
        self.app.alert(message, title, parent=self)

    def safe(self, fcn, what):
        """Run fcn; an error is shown instead of being lost (port of
        the safeCb wrappers of the analysis windows)."""
        try:
            fcn()
        except Exception as err:
            if not isinstance(err, ValueError):     # (not a refused entry)
                import traceback
                traceback.print_exc()
            self.alert('"%s" failed:\n%s' % (what, err), 'Callback error')
            lb = getattr(self, 'lblStatus', None)
            if lb is not None:
                lb.setText('Error in %s.' % what)

    def inform(self, message, title='TRSpecAnalyzer'):
        self.app.inform(message, title, parent=self)


def openChild(app, cls, *args):
    """Show the window of class cls (bring it to the front if it is
    open already). Returns the window."""
    w = app.childWindows.get(cls.KEY)
    if w is not None:
        try:
            w.raise_()
            w.activateWindow()
            return w
        except RuntimeError:
            app.childWindows.pop(cls.KEY, None)
    w = cls(app, *args)
    w.show()
    return w


def absMax(A):
    """Largest finite |value| of A; 1 if there is none or it is 0."""
    A = np.asarray(A, float)
    fin = np.abs(A[np.isfinite(A)])
    cl = float(np.max(fin)) if fin.size else 1.0
    return cl if cl > 0 else 1.0


def mapOn(app, ax, cax, wl, t, data, cl=None, alpha=None):
    """Symmetric color map of data (n_wl x n_t) on ax, color bar on
    cax. Returns the mesh (None if there are too few points)."""
    ax.clear()
    ax.set_facecolor(AX_BG)
    styleAxes(ax)
    if cl is None:
        cl = absMax(data)
    h = app.plotHeatmapOn(ax, t, data, wl=np.asarray(wl, float))
    if cax is not None:
        cax.clear()
    if h is not None:
        h.set_cmap(app.getColormap())
        h.set_clim(-cl, cl)
        if alpha is not None:
            h.set_alpha(np.asarray(alpha, float))
        if cax is not None:
            cb = ax.figure.colorbar(h, cax=cax)
            cb.set_label(app.sigName(), fontsize=FONT_PT)
            cb.ax.tick_params(labelsize=FONT_PT, length=3)
            cb.solids.set_alpha(1.0)
            sciTicks(cb.ax, 'y')
    ax.set_xlabel(app.xAxLabel())
    ax.set_ylabel('Delay time (%s)' % app.tUnitAx())
    return h


def mapLayout(cv, ax, cax, L=52.0, R=74.0, T=24.0, B=36.0):
    """Place a map axes and its color bar on a canvas."""
    W, H = cv.sizePt()
    cv.place(ax, L, B, W - L - R, H - T - B)
    if cax is not None:
        cv.place(cax, W - R + 8.0, B, 10.0, H - T - B)
        setSpan(ax, 0.0, W - R + 6.0)     # the title ends before the bar
    else:
        setSpan(ax, 0.0, W)


class MapTriple:
    """Three maps side by side (data, fit, residual), each as one axes
    or as a log panel above a linear one, with a vertical marker that
    is moved by click and drag on any of them. Port of buildMapAxes /
    drawMapStack / the map mouse handling of the analysis windows."""

    def __init__(self, app, onPick):
        self.app = app
        self.onPick = onPick
        self.cv = PlotCanvas()
        self.blit = BlitManager(self.cv)
        self.groups = []
        self.nRows = 0
        self.lines = []
        self.dragging = False
        self.cv.layoutFcn = self._layout
        self.cv.mpl_connect('button_press_event', self._press)
        self.cv.mpl_connect('motion_notify_event', self._move)
        self.cv.mpl_connect('button_release_event', self._release)

    def allAxes(self):
        return [ax for g in self.groups for ax in g['axes']]

    def build(self, nRows):
        fig = self.cv.figure
        fig.clear()
        self.groups = []
        self.lines = []
        self.blit.setArtists([])
        for _ in range(3):
            axes = [fig.add_axes([0.1, 0.1, 0.2, 0.3]) for _r in range(nRows)]
            cax = fig.add_axes([0.9, 0.1, 0.01, 0.8])
            cax.set_visible(False)
            for ax in axes:
                ax.set_facecolor(AX_BG)
            self.groups.append({'axes': axes, 'cax': cax})
        self.nRows = nRows
        self.cv.yGroups = [g['axes'] for g in self.groups]
        # (data, fit and residual show the same region when one of
        # them is zoomed)
        self.cv.zoomLinks = [[g['axes'][r] for g in self.groups]
                             for r in range(nRows)]
        self._layout()

    def _layout(self):
        if not self.groups:
            return
        cv = self.cv
        W, H = cv.sizePt()
        # (R: the color bar, its numbers and its label)
        # (L as the plots below the maps in the Global Analysis
        # window, whose rows end where the last map ends)
        L, R, T, B, gapV, gapH = 54.0, 60.0, 20.0, 34.0, 6.0, 4.0
        wG = (W - 2 * gapH) / 3.0
        for i, g in enumerate(self.groups):
            x0 = i * (wG + gapH)
            w = wG - L - R
            n = len(g['axes'])
            hTot = H - T - B - gapV * (n - 1)
            h = hTot / n
            # a title ends before the color bar, whose 'x10^n' factor
            # is written above the bar
            xr = x0 + L + w + 3.0 if g['cax'].get_visible() else x0 + wG
            for r, ax in enumerate(g['axes']):
                cv.place(ax, x0 + L, B + (n - 1 - r) * (h + gapV), w, h)
                setSpan(ax, x0, xr)
            cv.place(g['cax'], x0 + L + w + 5.0, B, 8.0, H - T - B)

    def draw(self, t, wl, mats, clims, showCbar, titles, mode, thr, tmin,
             tmax, fitRange, marker):
        """mats: three matrices n_wl x n_t. mode 'linear' | 'log' |
        'split'; thr: split threshold; tmin, tmax: shown delay range;
        fitRange: None or (lo, hi), cuts the shown range. A title may
        be a tuple: the full wording, then shorter ones for a narrow
        window."""
        app = self.app
        nRows = 2 if mode == 'split' else 1
        if nRows != self.nRows or len(self.groups) != 3:
            self.build(nRows)
        t = np.asarray(t, float)
        thr = max(float(thr), np.finfo(float).eps)
        if not np.isfinite(tmin):
            tmin = float(np.min(t))
        if not np.isfinite(tmax):
            tmax = float(np.max(t))
        if tmax <= tmin:
            tmax = float(np.max(t))
        if fitRange is not None:
            tmin = max(tmin, fitRange[0])
            tmax = min(tmax, fitRange[1])
            if tmax <= tmin:
                tmin, tmax = fitRange
        if mode == 'linear':
            masks = [(t >= tmin) & (t <= tmax)]
            scales = ['linear']
            ylims = [(tmin, tmax)]
        elif mode == 'log':
            pos = (t > 0) & (t >= tmin) & (t <= tmax)
            masks = [pos]
            scales = ['log']
            ylims = [(float(np.min(t[pos])), float(np.max(t[pos])))
                     if pos.any() else (1.0, 10.0)]
        else:
            masks = [(t >= thr) & (t <= tmax), (t <= thr) & (t >= tmin)]
            scales = ['log', 'linear']
            ylims = [(thr, tmax), (tmin, thr)]
        cmap = app.getColormap()
        self.lines = []
        for g, C, clim, cbOn, title in zip(self.groups, mats, clims,
                                           showCbar, titles):
            n = len(g['axes'])
            mesh = None
            for k, ax in enumerate(g['axes']):
                ax.clear()
                ax.set_facecolor(AX_BG)
                styleAxes(ax)
                m = masks[k]
                ax.set_yscale(scales[k])
                if np.count_nonzero(m) >= 2:
                    h = app.plotHeatmapOn(ax, t[m], C[:, m], wl=wl)
                    if h is not None:
                        h.set_cmap(cmap)
                        h.set_clim(clim[0], clim[1])
                        mesh = h
                ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
                lo, hi = ylims[k]
                if np.isfinite(lo) and np.isfinite(hi) and lo < hi:
                    if scales[k] == 'log' and lo <= 0:
                        lo = np.finfo(float).eps
                    ax.set_ylim(lo, hi)
                if scales[k] == 'linear':
                    # (lower panel of the split view: the tick at the
                    # threshold is labeled on the panel above)
                    ax.yaxis.set_major_locator(FewTicks(
                        4, prune='upper' if n > 1 else None))
                if k == n - 1:
                    ax.set_xlabel(app.xAxLabel())
                else:
                    ax.tick_params(labelbottom=False)
                ax.set_ylabel('Delay (%s)' % app.tUnitAx())
                if k == 0:
                    ti = title if isinstance(title, (tuple, list)) \
                        else (title,)
                    setTitle(ax, *ti, size=FONT_PT, bold=True)
                self.lines.append(ax.axvline(marker, color='k',
                                             linewidth=1.0, alpha=0.7))
            cax = g['cax']
            cax.clear()
            if cbOn and mesh is not None:
                cax.set_visible(True)
                cb = self.cv.figure.colorbar(mesh, cax=cax)
                cb.set_label(app.sigName(), fontsize=FONT_PT - 1,
                             labelpad=2.0)
                cb.ax.tick_params(labelsize=FONT_PT - 1, length=2)
                sciTicks(cb.ax, 'y')
            else:
                cax.set_visible(False)
        self.blit.setArtists(self.lines)
        self._layout()
        self.cv.draw_idle()

    def setMarker(self, wl):
        for ln in self.lines:
            ln.set_xdata([wl, wl])
        self.blit.update()

    def _press(self, ev):
        if ev.button != 1 or ev.inaxes not in self.allAxes():
            return
        self.dragging = True
        self.dragAx = ev.inaxes
        self.onPick(float(ev.xdata))

    def _move(self, ev):
        if not self.dragging:
            return
        try:
            x, _ = self.dragAx.transData.inverted().transform((ev.x, ev.y))
        except Exception:
            return
        if np.isfinite(x):
            self.onPick(float(x))

    def _release(self, _ev=None):
        self.dragging = False

    def click(self, wl):
        """A click at wavelength wl (used by the tests)."""
        self.onPick(float(wl))
