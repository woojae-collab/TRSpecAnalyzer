"""'Check axes orientation' dialog shown when one file is loaded.
Port of openSingleFileCheck."""
import numpy as np
from PySide6 import QtCore, QtWidgets

from .state import detect_axis_swap
from .textfit import setSpan
from .widgets import (AX_BG, FONT_PT, PlotCanvas, button, colorButton, label,
                      sciTicks, styleAxes, tip)


def orient(wl, t, A, transpose):
    """(wl, t, A) as parsed, or with the two axes interchanged."""
    wl = np.asarray(wl, float).ravel()
    t = np.asarray(t, float).ravel()
    A = np.asarray(A, float)
    if not transpose:
        return wl.copy(), t.copy(), A.copy()
    wlC, tC, AC = t.copy(), wl.copy(), A.T.copy()
    iw = np.argsort(wlC, kind='stable')
    it = np.argsort(tC, kind='stable')
    return wlC[iw], tC[it], AC[iw, :][:, it]


class SingleFileCheck(QtWidgets.QDialog):
    def __init__(self, app, wl, t, A, fname):
        super().__init__(app.fig)
        self.app = app
        self.wl, self.t, self.A = wl, t, A
        self.fname = fname
        self.canceled = True
        self.swapGuess = detect_axis_swap(wl, t)
        self.setWindowTitle('Check axes orientation')
        from .win_base import fitSize
        self.resize(*fitSize(820, 640))
        self.setModal(True)
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(12, 12, 12, 12)
        v.setSpacing(8)
        self.info = label('', color=(0.2, 0.2, 0.2))
        self.info.setWordWrap(True)
        self.info.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        v.addWidget(self.info)
        self.cv = PlotCanvas()
        v.addWidget(self.cv, 1)
        self.ax = self.cv.figure.add_axes([0.1, 0.1, 0.75, 0.8])
        self.cax = self.cv.figure.add_axes([0.9, 0.1, 0.02, 0.8])
        self.cv.layoutFcn = self._layout
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(10)
        row.addWidget(label('If wavelength ends up on the y-axis or delay on '
                            'the x-axis, use:', italic=True,
                            color=(0.35, 0.35, 0.35)))
        self.cbT = QtWidgets.QCheckBox('Transpose axes')
        self.cbT.setChecked(bool(self.swapGuess))
        self.cbT.setToolTip(tip(
            'Swap the parsed wavelength and delay axes before loading. '
            'Useful when the source file was saved with rows = delay and '
            'columns = wavelength instead of the expected layout.'))
        self.cbT.toggled.connect(lambda *_: self.drawPreview())
        row.addWidget(self.cbT)
        row.addWidget(label(fname, italic=True, color=(0.4, 0.4, 0.4)), 1)
        v.addLayout(row)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)
        row.addStretch(1)
        self.btnLoad = colorButton('Load', (0.30, 0.55, 0.80), width=140)
        self.btnLoad.clicked.connect(lambda *_: self.doLoad())
        self.btnLoad.setDefault(True)
        row.addWidget(self.btnLoad)
        row.addWidget(button('Cancel', self.reject, width=140))
        v.addLayout(row)
        self.drawPreview()

    def _layout(self):
        W, H = self.cv.sizePt()
        L, R, T, B = 52.0, 74.0, 22.0, 36.0
        self.cv.place(self.ax, L, B, W - L - R, H - T - B)
        self.cv.place(self.cax, W - R + 8.0, B, 10.0, H - T - B)
        setSpan(self.ax, 0.0, W - R + 6.0)

    def currentOrientation(self):
        return orient(self.wl, self.t, self.A, self.cbT.isChecked())

    def drawPreview(self):
        app = self.app
        wlC, tC, AC = self.currentOrientation()
        msg = ('File: %s\n\n'
               'With current orientation:\n'
               '  x-axis = %s: %d points, %.3g ... %.3g\n'
               '  y-axis = Delay time: %d points, %.3g ... %.3g' % (
                   self.fname,
                   'Wavenumber' if app.specUnit == 'cm-1' else 'Wavelength',
                   wlC.size, np.min(wlC), np.max(wlC),
                   tC.size, np.min(tC), np.max(tC)))
        if self.swapGuess and not self.cbT.isChecked():
            msg += ('\n\nWARNING: Heuristic suggests these axes may be '
                    'swapped. Tick "Transpose axes" if the map below looks '
                    'wrong.')
        elif self.swapGuess and self.cbT.isChecked():
            msg += ('\n\nTranspose applied (heuristic had suggested a swap). '
                    'Verify the map below looks right.')
        self.info.setText(msg)
        ax = self.ax
        ax.clear()
        ax.set_facecolor(AX_BG)
        styleAxes(ax)
        fin = np.abs(AC[np.isfinite(AC)])
        cl = float(np.max(fin)) if fin.size else 1.0
        if cl == 0:
            cl = 1.0
        h = app.plotHeatmapOn(ax, tC, AC, wl=wlC)
        self.cax.clear()
        if h is not None:
            h.set_cmap(app.getColormap())
            h.set_clim(-cl, cl)
            cb = self.cv.figure.colorbar(h, cax=self.cax)
            cb.set_label('ΔA or Intensity', fontsize=FONT_PT)
            cb.ax.tick_params(labelsize=FONT_PT, length=3)
            sciTicks(cb.ax, 'y')
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('Delay time (%s)' % app.tUnitAx())
        ax.grid(True, color='0.15', alpha=0.15, linewidth=0.6)
        ax.set_title('Preview of loaded data', fontsize=FONT_PT + 1)
        self.cv.relayout()
        self.cv.draw_idle()

    def doLoad(self):
        self.canceled = False
        self.accept()


def openSingleFileCheck(app, wl, t, A, fname):
    """Show the dialog; returns (wl, t, A, canceled)."""
    dlg = SingleFileCheck(app, wl, t, A, fname)
    dlg.exec()
    if dlg.canceled:
        return wl, t, A, True
    wl2, t2, A2 = dlg.currentOrientation()
    return wl2, t2, A2, False
