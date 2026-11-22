"""Small building blocks of the windows: numeric edit field, colored
button, matplotlib canvas with a layout given in points, and the
helper that moves a crosshair without redrawing the map under it."""
import html
import math
import os
import sys

import matplotlib

matplotlib.use('QtAgg')
# The plots look the same on every computer: the font that comes with
# matplotlib (not one that happens to be installed), and titles that
# stay where the layout puts them (matplotlib would lift a title above
# a 'x10^n' factor, out of the picture).
matplotlib.rcParams['font.family'] = ['DejaVu Sans']
matplotlib.rcParams['mathtext.fontset'] = 'dejavusans'
matplotlib.rcParams['axes.titley'] = 1.0

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure                              # noqa: E402
from matplotlib.ticker import MaxNLocator                         # noqa: E402
from PySide6 import QtCore, QtGui, QtWidgets                      # noqa: E402

from . import textfit, zoom                                       # noqa: E402

# line colors of MATLAB's default color order ("lines")
LINES = [(0.0000, 0.4470, 0.7410), (0.8500, 0.3250, 0.0980),
         (0.9290, 0.6940, 0.1250), (0.4940, 0.1840, 0.5560),
         (0.4660, 0.6740, 0.1880), (0.3010, 0.7450, 0.9330),
         (0.6350, 0.0780, 0.1840)]

AX_BG = (0.93, 0.93, 0.93)       # axes color under the 2D maps
FONT_PT = 9                      # font size inside the plots


def lineColor(k):
    return LINES[k % len(LINES)]


def tip(text):
    """Tooltip text that wraps (Qt wraps rich text only)."""
    return '<qt>' + html.escape(text) + '</qt>'


def prepareEnvironment():
    """Settings that must be made before the QApplication exists."""
    if sys.platform == 'win32':
        # no dark-mode handling by Qt: this program is always light
        os.environ.setdefault('QT_QPA_PLATFORM', 'windows:darkmode=0')


def applyLightTheme(qapp):
    """Fusion style with a fixed light palette, so that the program
    looks the same whatever theme the system uses. The color scheme
    is set to 'light' as well: on a system in dark mode Fusion draws
    the frame of a checkbox in a color meant for dark windows, which
    cannot be seen on our light ones."""
    try:
        qapp.styleHints().setColorScheme(QtCore.Qt.ColorScheme.Light)
    except Exception:                 # Qt older than 6.8
        pass
    qapp.setStyle('Fusion')
    pal = QtGui.QPalette()
    C = QtGui.QColor
    R = QtGui.QPalette.ColorRole
    G = QtGui.QPalette.ColorGroup
    colors = {
        R.Window: (240, 240, 240), R.WindowText: (20, 20, 20),
        R.Base: (255, 255, 255), R.AlternateBase: (245, 245, 245),
        R.Text: (20, 20, 20), R.Button: (240, 240, 240),
        R.ButtonText: (20, 20, 20), R.BrightText: (255, 255, 255),
        R.ToolTipBase: (255, 255, 225), R.ToolTipText: (20, 20, 20),
        R.Highlight: (0, 120, 215), R.HighlightedText: (255, 255, 255),
        R.PlaceholderText: (130, 130, 130),
        R.Light: (255, 255, 255), R.Midlight: (202, 202, 202),
        R.Mid: (184, 184, 184), R.Dark: (159, 159, 159),
        R.Shadow: (118, 118, 118),
        R.Link: (0, 0, 255), R.LinkVisited: (128, 0, 128),
    }
    if hasattr(R, 'Accent'):
        colors[R.Accent] = (0, 120, 215)
    # every color is given, for every state: nothing is left to the
    # system's own (possibly dark) colors
    for role, rgb in colors.items():
        for grp in (G.Active, G.Inactive, G.Disabled):
            pal.setColor(grp, role, C(*rgb))
    for role in (R.WindowText, R.Text, R.ButtonText):
        pal.setColor(G.Disabled, role, C(150, 150, 150))
    pal.setColor(G.Disabled, R.Highlight, C(170, 170, 170))
    qapp.setPalette(pal)


class NumEdit(QtWidgets.QLineEdit):
    """Numeric edit field: holds a float, shows it with a format, and
    reports a change only for a number inside the limits (anything
    else is put back). Stands in for MATLAB's numeric uieditfield."""

    valueChanged = QtCore.Signal(float)

    def __init__(self, value=0.0, fmt='%.4g', limits=(-math.inf, math.inf),
                 width=None, tooltip=None, parent=None, integer=False):
        super().__init__(parent)
        self._int = bool(integer)
        self._fmt = '%d' if integer else fmt
        self._lo, self._hi = limits
        self._value = float(value)
        self._shown = ''
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight |
                          QtCore.Qt.AlignmentFlag.AlignVCenter)
        if width:
            self.setFixedWidth(width)
        if tooltip:
            self.setToolTip(tip(tooltip))
        self.editingFinished.connect(self._commit)
        self._show()

    def value(self):
        return self._value

    def setValue(self, v):
        """Set the value from the program (no signal)."""
        self._value = float(v)
        self._show()

    def setLimits(self, lo, hi):
        self._lo, self._hi = float(lo), float(hi)

    def intValue(self):
        return int(round(self._value))

    def limits(self):
        return self._lo, self._hi

    def _show(self):
        v = self._value
        if math.isnan(v):
            s = 'NaN'
        elif math.isinf(v):
            s = 'Inf' if v > 0 else '-Inf'
        else:
            s = (self._fmt % v).strip()
        self._shown = s
        self.setText(s)
        self.setCursorPosition(0)

    def _commit(self):
        txt = self.text().strip()
        if txt == self._shown:
            return
        try:
            v = float(txt)
        except ValueError:
            v = math.nan
        if math.isnan(v) or v < self._lo or v > self._hi:
            if math.isfinite(self._lo) or math.isfinite(self._hi):
                msg = 'Value must be a number between %g and %g.' % (
                    self._lo, self._hi)
            else:
                msg = 'Value must be a number.'
            QtWidgets.QToolTip.showText(
                self.mapToGlobal(QtCore.QPoint(0, self.height())), msg, self)
            self._show()
            return
        if self._int:
            v = float(round(v))
        self._value = v
        self._show()
        self.valueChanged.emit(v)

    def enter(self, text):
        """Type a value and confirm it (used by the tests)."""
        self.setText(str(text))
        self._commit()


def colorButton(text, rgb, tooltip=None, width=None):
    """Button with a background color and white text."""
    b = QtWidgets.QPushButton(text)
    c = [int(round(255 * x)) for x in rgb]
    hov = [min(255, int(round(x * 1.12 + 8))) for x in c]
    prs = [int(round(x * 0.85)) for x in c]
    brd = [int(round(x * 0.75)) for x in c]
    b.setStyleSheet(
        'QPushButton { background-color: rgb(%d,%d,%d); color: white;'
        ' border: 1px solid rgb(%d,%d,%d); border-radius: 3px;'
        ' padding: 4px 8px; }'
        'QPushButton:hover { background-color: rgb(%d,%d,%d); }'
        'QPushButton:pressed { background-color: rgb(%d,%d,%d); }'
        'QPushButton:disabled { background-color: rgb(190,190,190);'
        ' color: rgb(235,235,235); border-color: rgb(170,170,170); }'
        % (c[0], c[1], c[2], brd[0], brd[1], brd[2],
           hov[0], hov[1], hov[2], prs[0], prs[1], prs[2]))
    if tooltip:
        b.setToolTip(tip(tooltip))
    if width:
        b.setMinimumWidth(width)     # grows if the text needs more
    return b


def button(text, slot=None, tooltip=None, width=None):
    b = QtWidgets.QPushButton(text)
    if slot is not None:
        b.clicked.connect(lambda *_: slot())
    if tooltip:
        b.setToolTip(tip(tooltip))
    if width:
        b.setMinimumWidth(width)     # grows if the text needs more
    return b


def label(text, right=False, color=None, italic=False):
    lb = QtWidgets.QLabel(text)
    if right:
        lb.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight |
                        QtCore.Qt.AlignmentFlag.AlignVCenter)
    css = ''
    if color is not None:
        css += 'color: rgb(%d,%d,%d);' % tuple(int(round(255 * x)) for x in color)
    if italic:
        css += 'font-style: italic;'
    if css:
        lb.setStyleSheet(css)
    return lb


class ElideLabel(QtWidgets.QLabel):
    """One line of text that ends in "..." when the room is too short
    for it; the whole text is then the tooltip. text() always gives
    the whole text."""

    def __init__(self, text='', right=False, color=None):
        super().__init__(text)
        self._color = None
        if color is not None:
            self._color = QtGui.QColor(*(int(round(255 * x)) for x in color))
        if right:
            self.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight |
                              QtCore.Qt.AlignmentFlag.AlignVCenter)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                           QtWidgets.QSizePolicy.Policy.Preferred)

    def minimumSizeHint(self):
        return QtCore.QSize(20, super().minimumSizeHint().height())

    def setText(self, text):
        super().setText(text)
        self._sync()
        self.update()

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self._sync()

    def shownText(self):
        """The text as it stands in the window (cut or whole)."""
        return self.fontMetrics().elidedText(
            self.text(), QtCore.Qt.TextElideMode.ElideRight,
            self.contentsRect().width())

    def _sync(self):
        tipText = self.text() if self.shownText() != self.text() else ''
        if tipText != self.toolTip():
            self.setToolTip(tipText)

    def paintEvent(self, _ev):
        p = QtGui.QPainter(self)
        if self._color is not None:
            p.setPen(self._color)
        else:
            p.setPen(self.palette().color(QtGui.QPalette.ColorRole.WindowText))
        flags = self.alignment().value | \
            QtCore.Qt.TextFlag.TextSingleLine.value
        p.drawText(self.contentsRect(), flags, self.shownText())
        p.end()


class FitLabel(QtWidgets.QLabel):
    """One line of text that never widens its window: of the wordings
    given to setTexts (the full one first), the first that fits the
    room is shown. The full wording is the tooltip while a shorter one
    is on show."""

    def __init__(self, color=None):
        super().__init__('')
        self._texts = ['']
        if color is not None:
            self.setStyleSheet('color: rgb(%d,%d,%d);' % tuple(
                int(round(255 * x)) for x in color))
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                           QtWidgets.QSizePolicy.Policy.Preferred)
        self.setMinimumWidth(20)

    def setTexts(self, *texts):
        self._texts = [t for t in texts if t is not None] or ['']
        self._choose()

    def fullText(self):
        return self._texts[0]

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self._choose()

    def _choose(self):
        room = self.contentsRect().width()
        fm = self.fontMetrics()
        txt = self._texts[-1]
        for t in self._texts:
            if fm.horizontalAdvance(t) <= room:
                txt = t
                break
        if txt != self.text():
            self.setText(txt)
        self.setToolTip('' if txt == self._texts[0] else self._texts[0])


def hbox(spacing=4, margins=(0, 0, 0, 0)):
    """(widget, layout) of one row of controls."""
    w = QtWidgets.QWidget()
    lay = QtWidgets.QHBoxLayout(w)
    lay.setContentsMargins(*margins)
    lay.setSpacing(spacing)
    return w, lay


class FlowLayout(QtWidgets.QLayout):
    """A row of controls that continues on the next line when the
    window is too narrow for it. addStretch() marks the place where
    the spare room of a line goes (what follows is pushed to the
    right). addFill() adds a last control that takes the room left on
    the last line, whatever it would like to have (a status text): it
    never makes a line of its own unless less than minWidth is left."""

    def __init__(self, parent=None, hspacing=4, vspacing=4):
        super().__init__(parent)
        self._items = []
        self._hs, self._vs = hspacing, vspacing
        self._stretchAt = None
        self._fillAt = None
        self._fillMin = 0
        self.setContentsMargins(0, 0, 0, 0)

    def addItem(self, item):
        self._items.append(item)

    def addStretch(self, *_):
        self._stretchAt = len(self._items)

    def addFill(self, widget, minWidth=160):
        self.addWidget(widget)
        self._fillAt = len(self._items) - 1
        self._fillMin = int(minWidth)

    def addPair(self, first, second, spacing=4):
        """Two controls that stay together (a label and its field)."""
        w = QtWidgets.QWidget()
        h = QtWidgets.QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(spacing)
        h.addWidget(first)
        h.addWidget(second)
        self.addWidget(w)
        return w

    def count(self):
        return len(self._items)

    def itemAt(self, i):
        return self._items[i] if 0 <= i < len(self._items) else None

    def takeAt(self, i):
        return self._items.pop(i) if 0 <= i < len(self._items) else None

    def expandingDirections(self):
        return QtCore.Qt.Orientation(0)

    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, width):
        return self._arrange(QtCore.QRect(0, 0, width, 0), False)

    def setGeometry(self, rect):
        super().setGeometry(rect)
        self._arrange(rect, True)

    def _hint(self, i, it):
        """Size an item asks for; the filler asks for its least width."""
        sz = it.sizeHint()
        if i == self._fillAt:
            return QtCore.QSize(self._fillMin, sz.height())
        return sz

    def sizeHint(self):
        m = self.contentsMargins()
        sz = [self._hint(i, it) for i, it in enumerate(self._items)]
        w = sum(q.width() for q in sz) + self._hs * max(len(sz) - 1, 0)
        h = max([q.height() for q in sz] or [0])
        return QtCore.QSize(w + m.left() + m.right(),
                            h + m.top() + m.bottom())

    def minimumSize(self):
        m = self.contentsMargins()
        sz = [self._hint(i, it) for i, it in enumerate(self._items)]
        return QtCore.QSize(
            max([q.width() for q in sz] or [0]) + m.left() + m.right(),
            max([q.height() for q in sz] or [0]) + m.top() + m.bottom())

    def _arrange(self, rect, apply):
        m = self.contentsMargins()
        x0 = rect.x() + m.left()
        width = rect.width() - m.left() - m.right()
        lines, cur, used = [], [], 0
        fill = None
        for i, it in enumerate(self._items):
            if it.widget() is not None and it.widget().isHidden():
                continue
            if i == self._fillAt:
                fill = (i, it)
                continue
            sz = it.sizeHint()
            need = sz.width() + (self._hs if cur else 0)
            if cur and used + need > width:
                lines.append((cur, used))
                cur, used, need = [], 0, sz.width()
            cur.append((i, it, sz))
            used += need
        if cur:
            lines.append((cur, used))
        if fill is not None:
            # the rest of the last line, or a line of its own
            fh = fill[1].sizeHint().height()
            left = width - (lines[-1][1] + self._hs) if lines else -1
            if left >= self._fillMin:
                line, used = lines[-1]
                lines[-1] = (line + [(fill[0], fill[1],
                                      QtCore.QSize(left, fh))], width)
            else:
                lines.append(([(fill[0], fill[1],
                                QtCore.QSize(max(width, 0), fh))], width))
        y = rect.y() + m.top()
        for line, used in lines:
            lh = max(sz.height() for _i, _it, sz in line)
            x = x0
            for i, it, sz in line:
                if i == self._stretchAt:
                    x += max(0, width - used)
                if apply:
                    it.setGeometry(QtCore.QRect(
                        x, y + (lh - sz.height()) // 2,
                        sz.width(), sz.height()))
                x += sz.width() + self._hs
            y += lh + self._vs
        if lines:
            y -= self._vs
        return y - rect.y() + m.bottom()


def flowbox(spacing=4):
    """(widget, layout) of a row of controls that wraps when the
    window is narrow."""
    w = QtWidgets.QWidget()
    lay = FlowLayout(w, spacing, 4)
    sp = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Policy.Preferred,
                               QtWidgets.QSizePolicy.Policy.Minimum)
    sp.setHeightForWidth(True)
    w.setSizePolicy(sp)
    return w, lay


class VScroll(QtWidgets.QScrollArea):
    """Lets its content scroll up and down when the window is too low
    for it. Never sideways: the area is as wide as the content needs
    (at least minWidth)."""

    def __init__(self, inner, minWidth=0, fixedWidth=False):
        super().__init__()
        self._minW = int(minWidth)
        self._fixedW = bool(fixedWidth)
        self.setWidget(inner)
        self.setWidgetResizable(True)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.viewport().setAutoFillBackground(False)
        inner.setAutoFillBackground(False)

    def _width(self):
        return max(self._minW, self.widget().minimumSizeHint().width()) \
            + self.verticalScrollBar().sizeHint().width()

    def minimumSizeHint(self):
        return QtCore.QSize(self._width(), 100)

    def sizeHint(self):
        s = self.widget().sizeHint()
        w = self._width() if self._fixedW else \
            max(s.width() + self.verticalScrollBar().sizeHint().width(),
                self._width())
        return QtCore.QSize(w, s.height())


def sidePanel(title, minWidth):
    """Titled panel for a column of controls: (box, layout). The box
    is as wide as its controls need with the font in use (at least
    minWidth) and lets them scroll when the window is too low."""
    box = panel(title)
    outer = QtWidgets.QVBoxLayout(box)
    outer.setContentsMargins(2, 12, 2, 4)
    inner = QtWidgets.QWidget()
    lay = QtWidgets.QVBoxLayout(inner)
    lay.setContentsMargins(6, 0, 4, 4)
    lay.setSpacing(6)
    outer.addWidget(VScroll(inner, minWidth - 4, fixedWidth=True))
    box.setSizePolicy(QtWidgets.QSizePolicy.Policy.Fixed,
                      QtWidgets.QSizePolicy.Policy.Expanding)
    return box, lay


def panel(title):
    """Titled panel (uipanel)."""
    g = QtWidgets.QGroupBox(title)
    g.setStyleSheet(
        'QGroupBox { font-weight: bold; border: 1px solid rgb(190,190,190);'
        ' border-radius: 3px; margin-top: 9px; }'
        'QGroupBox::title { subcontrol-origin: margin; left: 8px;'
        ' padding: 0 3px; }')
    return g


ZOOM_HELP = (
    'Mouse wheel: zoom in and out at the pointer.\n'
    '    Shift + wheel, or the wheel over the numbers of the x-axis: '
    'only the x-axis.\n'
    '    Ctrl + wheel, or the wheel over the numbers of the y-axis: '
    'only the y-axis.\n'
    'Drag with the right button: zoom into the rectangle.\n'
    'Drag with the left button: move the zoomed plot. Where a click '
    'already means something (a crosshair, a marker), press the wheel '
    'and drag instead.\n'
    'Double click, or "Reset view" in this menu: back to the full view. '
    'Turning the wheel back as far as it goes does the same.\n\n'
    'A zoomed plot stays zoomed when another delay or wavelength is '
    'chosen. It returns to the full view when other data, another '
    'unit, scale or display range are set.')


class PlotCanvas(FigureCanvasQTAgg):
    """Matplotlib canvas. The axes are placed by layoutFcn, which is
    called whenever the size changes; before every drawing the texts
    are made to fit (textfit.fitCanvas) and the view the user zoomed to
    is put back (zoom.apply). Legends can be dragged; the wheel zooms,
    a drag with the right button zooms into a rectangle, a drag with
    the left or the middle button moves a zoomed plot; a right click
    offers to copy or save the picture and to reset the view."""

    def __init__(self, parent=None):
        fig = Figure(facecolor='white')
        super().__init__(fig)
        if parent is not None:
            self.setParent(parent)
        self.layoutFcn = None
        self.xGroups = []          # axes side by side with one x-label
        self.yGroups = []          # axes above each other, one y-label
        self.zoomLinks = []        # axes that are zoomed together
        self.legendsHidden = False
        self.ownsLeft = False      # a window uses the left button itself
        self._handCursor = False
        self._drag = None          # a pan or a zoom rectangle under way
        self._band = None
        self.setMinimumSize(160, 110)
        self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Expanding,
                           QtWidgets.QSizePolicy.Policy.Expanding)
        # (the menu opens when the right button is let go without a
        # drag; Qt's own request would also come after a zoom drag)
        self.setContextMenuPolicy(
            QtCore.Qt.ContextMenuPolicy.PreventContextMenu)
        self.mpl_connect('motion_notify_event', self._hover)
        self.mpl_connect('button_press_event', self._press)
        self.mpl_connect('button_release_event', self._release)
        self._ready = True

    def mpl_connect(self, s, func):
        if s == 'button_press_event' and getattr(self, '_ready', False):
            self.ownsLeft = True
        return super().mpl_connect(s, func)

    def sizeHint(self):
        # (matplotlib reports the present size, which would make a
        # canvas that was large once claim that room for ever)
        return self.minimumSize()

    def minimumSizeHint(self):
        return self.minimumSize()

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self.relayout()

    def relayout(self):
        if self.layoutFcn is not None:
            try:
                self.layoutFcn()
            except Exception:
                if os.environ.get('TRSPEC_DEBUG'):
                    raise

    def draw(self):
        try:
            zoom.apply(self)
            textfit.fitCanvas(self)
        except Exception:
            if os.environ.get('TRSPEC_DEBUG'):
                raise
        super().draw()

    def sizePt(self):
        """Width and height of the figure in points."""
        w, h = self.figure.get_size_inches()
        return w * 72.0, h * 72.0

    def place(self, ax, x, y, w, h):
        """Put ax at (x, y, w, h), given in points from the lower left
        corner of the figure."""
        W, H = self.sizePt()
        if W <= 0 or H <= 0:
            return
        w = max(w, 4.0)
        h = max(h, 4.0)
        ax.set_position([x / W, y / H, w / W, h / H])

    def legendAt(self, ev):
        """The legend under a mouse event (None if there is none):
        a click there moves the legend and means nothing else."""
        return textfit.legendAt(self, ev)

    # ---- zoom and pan ----------------------------------------------
    def wheelEvent(self, event):
        d = event.angleDelta()
        steps = (d.y() or d.x()) / 120.0
        x, y = self.mouseEventCoords(event)
        if steps and self.wheelAt(x, y, steps, event.modifiers()):
            event.accept()
        else:
            event.ignore()          # (a panel around may scroll)

    def wheelAt(self, x, y, steps, mods=None):
        """Wheel turned by steps notches at the pixel (x, y) of the
        figure. True if a plot was zoomed."""
        ax, mode = zoom.axesAt(self, x, y)
        if ax is None:
            return False
        if mode == 'xy' and mods is not None:
            if mods & QtCore.Qt.KeyboardModifier.ShiftModifier:
                mode = 'x'
            elif mods & QtCore.Qt.KeyboardModifier.ControlModifier:
                mode = 'y'
        steps = max(-5.0, min(5.0, float(steps)))
        zoom.zoomAt(self, ax, x, y, zoom.STEP ** (-steps), mode)
        self.draw_idle()
        return True

    def resetView(self, ax=None):
        """Full view of ax (and the axes that share a range with it),
        or of all axes of this canvas."""
        if zoom.reset(self, ax):
            self.draw_idle()
            return True
        return False

    def _qpos(self, ev):
        g = getattr(ev, 'guiEvent', None)
        if g is not None and hasattr(g, 'position'):
            return g.position().toPoint()
        r = self.devicePixelRatioF() or 1.0
        return QtCore.QPoint(int(round(ev.x / r)),
                             int(round(self.height() - ev.y / r)))

    def _press(self, ev):
        if self._drag is not None:
            return
        ax, mode = zoom.axesAt(self, ev.x, ev.y)
        if ev.button == 3:
            self._drag = {'kind': 'box', 'x0': ev.x, 'y0': ev.y,
                          'ax': ax if mode == 'xy' else None,
                          'pos': self._qpos(ev), 'moved': False}
            return
        if ev.button == 1 and (self.ownsLeft or
                               textfit.legendAt(self, ev) is not None):
            return
        if ev.button not in (1, 2) or ax is None or mode != 'xy':
            return
        if ev.dblclick:
            if ev.button == 1:
                self.resetView(ax)
            return
        if zoom.isZoomed(ax):
            self._drag = {'kind': 'pan', 'x0': ev.x, 'y0': ev.y,
                          'button': ev.button,
                          'st': zoom.panStart(self, ax)}
            self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)

    def _dragMove(self, ev):
        d = self._drag
        if d['kind'] == 'pan':
            zoom.panBy(self, d['st'], ev.x - d['x0'], ev.y - d['y0'])
            self.draw_idle()
            return
        if d['ax'] is None:
            return
        if abs(ev.x - d['x0']) > 4 or abs(ev.y - d['y0']) > 4:
            d['moved'] = True
        if d['moved']:
            if self._band is None:
                self._band = QtWidgets.QRubberBand(
                    QtWidgets.QRubberBand.Shape.Rectangle, self)
            r = self.devicePixelRatioF() or 1.0
            b = d['ax'].bbox
            box = QtCore.QRect(
                QtCore.QPoint(int(round(b.x0 / r)),
                              int(round(self.height() - b.y1 / r))),
                QtCore.QPoint(int(round(b.x1 / r)),
                              int(round(self.height() - b.y0 / r))))
            rect = QtCore.QRect(d['pos'], self._qpos(ev)).normalized()
            self._band.setGeometry(rect.intersected(box))
            self._band.show()

    def _release(self, ev):
        d = self._drag
        if d is None:
            return
        if d['kind'] == 'pan':
            if ev.button != d['button']:
                return
            self._drag = None
            self.unsetCursor()
            self._handCursor = False
            return
        if ev.button != 3:
            return
        self._drag = None
        if self._band is not None:
            self._band.hide()
        if d['moved'] and d['ax'] is not None:
            mode = ('x' if abs(ev.x - d['x0']) >= 5 else '') + \
                ('y' if abs(ev.y - d['y0']) >= 5 else '')
            if mode:
                zoom.zoomTo(self, d['ax'], d['x0'], d['y0'], ev.x, ev.y,
                            mode)
                self.draw_idle()
        elif not d['moved']:
            self._menu(d['pos'])

    def _hover(self, ev):
        if self._drag is not None:
            self._dragMove(ev)
            return
        hand = ev.button is None and textfit.legendAt(self, ev) is not None
        if hand and not self._handCursor:
            self.setCursor(QtCore.Qt.CursorShape.OpenHandCursor)
            self._handCursor = True
        elif not hand and self._handCursor and ev.button is None:
            self.unsetCursor()
            self._handCursor = False

    # ---- right-click menu ------------------------------------------
    def menuActions(self):
        """(key, text, enabled) of the entries of the right-click menu;
        None stands for a separating line."""
        out = [('copy', 'Copy picture', True),
               ('save', 'Save picture as...', True)]
        if textfit.hasLegend(self):
            out.append(None)
            out.append(('legend', 'Show legend' if self.legendsHidden
                        else 'Hide legend', True))
            if textfit.hasMovedLegend(self):
                out.append(('legendAuto',
                            'Legend back to its automatic place', True))
        out.append(None)
        out.append(('reset', 'Reset view', zoom.anyZoomed(self)))
        out.append(('zoomHelp', 'Zoom with the mouse...', True))
        return out

    def menuDo(self, key):
        if key == 'copy':
            QtWidgets.QApplication.clipboard().setPixmap(self.grab())
        elif key == 'save':
            fn, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, 'Save picture', 'figure.png',
                'PNG image (*.png);;PDF file (*.pdf);;SVG file (*.svg)')
            if fn:
                try:
                    self.figure.savefig(fn, dpi=200)
                except Exception as err:
                    QtWidgets.QMessageBox.warning(
                        self, 'Save picture', 'Could not save:\n%s' % err)
        elif key == 'legend':
            textfit.showLegends(self, self.legendsHidden)
        elif key == 'legendAuto':
            textfit.resetLegends(self)
        elif key == 'reset':
            self.resetView()
        elif key == 'zoomHelp':
            QtWidgets.QMessageBox.information(self, 'Zoom and pan',
                                              ZOOM_HELP)

    def _menu(self, pos):
        m = QtWidgets.QMenu(self)
        acts = {}
        for item in self.menuActions():
            if item is None:
                m.addSeparator()
                continue
            a = m.addAction(item[1])
            a.setEnabled(item[2])
            acts[a] = item[0]
        act = m.exec(self.mapToGlobal(pos))
        if act is not None:
            self.menuDo(acts[act])


class BlitManager:
    """Moves a few artists (a crosshair) over a picture that is kept
    as it is: the picture is drawn once, stored, and put back before
    the artists are drawn at their new place."""

    def __init__(self, canvas):
        self.canvas = canvas
        self.artists = []
        self._bg = None
        canvas.mpl_connect('draw_event', self._onDraw)

    def setArtists(self, artists):
        self.artists = list(artists)
        for a in self.artists:
            a.set_animated(True)
        self._bg = None

    def _onDraw(self, _ev):
        cv = self.canvas
        self._bg = cv.copy_from_bbox(cv.figure.bbox)
        self._drawArtists()

    def _drawArtists(self):
        fig = self.canvas.figure
        for a in self.artists:
            try:
                fig.draw_artist(a)
            except Exception:
                pass

    def update(self):
        cv = self.canvas
        if self._bg is None:
            cv.draw_idle()
            return
        cv.restore_region(self._bg)
        self._drawArtists()
        cv.blit(cv.figure.bbox)


def styleAxes(ax, grid=False):
    """Font size and tick style used in all plots."""
    ax.tick_params(labelsize=FONT_PT, length=3)
    ax.xaxis.label.set_size(FONT_PT)
    ax.yaxis.label.set_size(FONT_PT)
    ax.xaxis.get_offset_text().set_size(FONT_PT - 1)
    ax.yaxis.get_offset_text().set_size(FONT_PT - 1)
    if grid:
        ax.grid(True, color='0.15', alpha=0.15, linewidth=0.6)
        ax.set_axisbelow(True)


class FewTicks(MaxNLocator):
    """At most nmax intervals between the ticks of a linear axis, and
    fewer when the axis is too short for that many labels."""

    def __init__(self, nmax=4, prune=None):
        super().__init__(nbins='auto', prune=prune, min_n_ticks=1)
        self._nmax = int(nmax)

    def _raw_ticks(self, vmin, vmax):
        n = self._nmax
        if self.axis is not None:
            try:
                n = int(min(max(self.axis.get_tick_space(), 1), n))
            except Exception:
                pass
        self._nbins = n
        try:
            return super()._raw_ticks(vmin, vmax)
        finally:
            self._nbins = 'auto'


def sciTicks(ax, axis='y'):
    """Tick labels of a linear axis with a common factor (x 10^n)
    when the numbers are small or large, so that they stay short."""
    try:
        ax.ticklabel_format(axis=axis, style='sci', scilimits=(-2, 4),
                            useMathText=True)
    except (AttributeError, ValueError):
        pass
    for a in (ax.xaxis, ax.yaxis):
        a.get_offset_text().set_size(FONT_PT - 1)


def alert(parent, message, title='TRSpecAnalyzer'):
    QtWidgets.QMessageBox.warning(parent, title, message)


def inform(parent, message, title='TRSpecAnalyzer'):
    QtWidgets.QMessageBox.information(parent, title, message)


def confirm(parent, message, title='Confirm', yes='Yes', no='Cancel'):
    box = QtWidgets.QMessageBox(parent)
    box.setIcon(QtWidgets.QMessageBox.Icon.Question)
    box.setWindowTitle(title)
    box.setText(message)
    bYes = box.addButton(yes, QtWidgets.QMessageBox.ButtonRole.AcceptRole)
    box.addButton(no, QtWidgets.QMessageBox.ButtonRole.RejectRole)
    box.exec()
    return box.clickedButton() is bYes


class DataTable(QtWidgets.QTableWidget):
    """Small editable table (uitable): every column is 'num', 'bool'
    or 'text'. cellEdited(row, column) is sent after a change made by
    the user; a number that cannot be read is put back."""

    cellEdited = QtCore.Signal(int, int)

    def __init__(self, headers, kinds, editable=None, widths=None,
                 numFmt='%.6g', parent=None):
        super().__init__(0, len(headers), parent)
        self.kinds = list(kinds)
        self.editable = list(editable) if editable is not None \
            else [True] * len(headers)
        self.numFmt = numFmt
        self._filling = False
        self._vals = []
        self.setHorizontalHeaderLabels(list(headers))
        if widths:
            for c, w in enumerate(widths):
                if w:
                    self.setColumnWidth(c, w)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().setDefaultSectionSize(24)
        self.itemChanged.connect(self._onItem)

    def setData(self, rows):
        """rows: list of rows, each a list with one value per column."""
        F = QtCore.Qt.ItemFlag
        C = QtCore.Qt.CheckState
        self._filling = True
        try:
            self._vals = [list(r) for r in rows]
            self.setRowCount(len(rows))
            for i, r in enumerate(rows):
                for j, v in enumerate(r):
                    kind = self.kinds[j]
                    it = QtWidgets.QTableWidgetItem()
                    fl = F.ItemIsEnabled | F.ItemIsSelectable
                    if kind == 'bool':
                        if self.editable[j]:
                            fl |= F.ItemIsUserCheckable
                        it.setCheckState(C.Checked if v else C.Unchecked)
                    else:
                        if self.editable[j]:
                            fl |= F.ItemIsEditable
                        it.setText(self._txt(j, v))
                        if kind == 'num':
                            it.setTextAlignment(
                                QtCore.Qt.AlignmentFlag.AlignRight |
                                QtCore.Qt.AlignmentFlag.AlignVCenter)
                    it.setFlags(fl)
                    self.setItem(i, j, it)
        finally:
            self._filling = False

    def _txt(self, j, v):
        if self.kinds[j] == 'num':
            if v is None:
                return ''
            v = float(v)
            if math.isnan(v):
                return 'NaN'
            if math.isinf(v):
                return 'Inf' if v > 0 else '-Inf'
            return self.numFmt % v
        return '' if v is None else str(v)

    def data(self):
        return [list(r) for r in self._vals]

    def value(self, i, j):
        return self._vals[i][j]

    def setValue(self, i, j, v):
        """Set one cell from the program (no signal)."""
        self._vals[i][j] = v
        it = self.item(i, j)
        self._filling = True
        try:
            if self.kinds[j] == 'bool':
                it.setCheckState(QtCore.Qt.CheckState.Checked if v
                                 else QtCore.Qt.CheckState.Unchecked)
            else:
                it.setText(self._txt(j, v))
        finally:
            self._filling = False

    def _onItem(self, it):
        if self._filling:
            return
        i, j = it.row(), it.column()
        if i >= len(self._vals):
            return
        kind = self.kinds[j]
        if kind == 'bool':
            new = it.checkState() == QtCore.Qt.CheckState.Checked
            if new == bool(self._vals[i][j]):
                return
            self._vals[i][j] = new
        elif kind == 'num':
            txt = it.text().strip()
            try:
                new = float(txt)
            except ValueError:
                new = math.nan
            if math.isnan(new) and txt.lower() != 'nan':
                self.setValue(i, j, self._vals[i][j])       # put back
                return
            self._vals[i][j] = new
            self.setValue(i, j, new)
        else:
            self._vals[i][j] = it.text()
        self.cellEdited.emit(i, j)

    def edit_(self, i, j, v):
        """Change a cell as the user would (used by the tests)."""
        it = self.item(i, j)
        if self.kinds[j] == 'bool':
            it.setCheckState(QtCore.Qt.CheckState.Checked if v
                             else QtCore.Qt.CheckState.Unchecked)
        else:
            it.setText(str(v))


def monoFont(size=None):
    f = QtGui.QFontDatabase.systemFont(
        QtGui.QFontDatabase.SystemFont.FixedFont)
    if size:
        f.setPointSize(size)
    return f


def textArea(lines=None):
    """Read-only text box with a fixed-pitch font (uitextarea)."""
    t = QtWidgets.QPlainTextEdit()
    t.setReadOnly(True)
    t.setFont(monoFont())
    t.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
    if lines:
        t.setPlainText('\n'.join(lines))
    return t


def runBusy(app, parent, message, func):
    """Run func() while a 'please wait' box is shown; returns its
    result or raises its error. The work runs in a second thread so
    that the windows stay alive. Headless: plain call."""
    if getattr(app, 'headless', False):
        return func()
    box = {}

    class Worker(QtCore.QThread):
        def run(self):
            try:
                box['out'] = func()
            except BaseException as err:        # noqa: B036
                box['err'] = err

    dlg = QtWidgets.QProgressDialog(message, None, 0, 0, parent)
    dlg.setWindowTitle('Please wait')
    dlg.setWindowModality(QtCore.Qt.WindowModality.ApplicationModal)
    dlg.setCancelButton(None)
    dlg.setMinimumDuration(300)
    dlg.setWindowFlag(QtCore.Qt.WindowType.WindowCloseButtonHint, False)
    th = Worker()
    loop = QtCore.QEventLoop()
    th.finished.connect(loop.quit)
    th.start()
    dlg.show()
    loop.exec()
    th.wait()
    dlg.close()
    dlg.deleteLater()
    if 'err' in box:
        raise box['err']
    return box.get('out')


class _Relay(QtCore.QObject):
    """Takes the progress reports of the worker thread over to the
    thread of the windows."""

    def __init__(self, fcn):
        super().__init__()
        self.fcn = fcn

    @QtCore.Slot(float, str)
    def step(self, fr, msg):
        self.fcn(fr, msg)


def runProgress(app, parent, title, message, func, onStep=None):
    """Run func(prog) while a box with a progress bar is shown;
    prog(fraction, message) moves the bar and changes its text (and
    calls onStep, in the thread of the windows). Returns the result of
    func or raises its error. Headless: plain call."""
    if getattr(app, 'headless', False):
        def prog(fr, msg):
            if onStep is not None:
                onStep(float(fr), str(msg))
        return func(prog)
    box = {}

    class Worker(QtCore.QThread):
        step = QtCore.Signal(float, str)

        def run(self):
            try:
                box['out'] = func(
                    lambda fr, msg: self.step.emit(float(fr), str(msg)))
            except BaseException as err:        # noqa: B036
                box['err'] = err

    dlg = QtWidgets.QProgressDialog(message, None, 0, 1000, parent)
    dlg.setWindowTitle(title)
    dlg.setWindowModality(QtCore.Qt.WindowModality.ApplicationModal)
    dlg.setCancelButton(None)
    dlg.setAutoClose(False)
    dlg.setAutoReset(False)
    dlg.setMinimumDuration(0)
    dlg.setMinimumWidth(380)
    dlg.setWindowFlag(QtCore.Qt.WindowType.WindowCloseButtonHint, False)

    def show(fr, msg):
        dlg.setValue(int(round(1000 * max(0.0, min(1.0, fr)))))
        dlg.setLabelText(msg)
        if onStep is not None:
            onStep(fr, msg)

    relay = _Relay(show)
    th = Worker()
    th.step.connect(relay.step, QtCore.Qt.ConnectionType.QueuedConnection)
    loop = QtCore.QEventLoop()
    th.finished.connect(loop.quit)
    dlg.setValue(20)
    dlg.show()
    th.start()
    loop.exec()
    th.wait()
    QtCore.QCoreApplication.processEvents()     # reports still waiting
    dlg.close()
    dlg.deleteLater()
    if 'err' in box:
        raise box['err']
    return box.get('out')
