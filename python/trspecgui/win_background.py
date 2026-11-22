"""Background Correction window. Port of openBackgroundWindow."""
import warnings

import numpy as np
from matplotlib import colormaps
from PySide6 import QtWidgets

from .textfit import legend
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label, sciTicks, styleAxes)
from .win_base import GREEN, Z_LABEL, ChildWindow


def meanSpectrum(A, N):
    """Mean of the first N columns, NaN left out."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        return np.nanmean(A[:, :N], axis=1)


class BackgroundWindow(ChildWindow):
    KEY = 'bg'

    def __init__(self, app):
        super().__init__(app, 'Background Correction', 720, 520)
        self.applyBusy = False
        gl = QtWidgets.QGridLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setColumnStretch(2, 1)
        gl.setRowStretch(1, 1)
        gl.addWidget(label('# of initial delay points to average:',
                           right=True), 0, 0)
        n = app.delay.size
        defaultN = app.bgN if app.bgApplied else min(5, n)
        self.edN = NumEdit(defaultN, limits=(1, n), width=120, integer=True)
        gl.addWidget(self.edN, 0, 1)
        self.lblInfo = label('', color=(0.25, 0.25, 0.25))
        gl.addWidget(self.lblInfo, 0, 2)
        self.cv = PlotCanvas()
        gl.addWidget(self.cv, 1, 0, 1, 3)
        self.ax = self.cv.figure.add_axes([0.1, 0.1, 0.8, 0.8])
        self.cv.layoutFcn = self._layout
        row = QtWidgets.QHBoxLayout()
        row.addStretch(1)
        self.btnApply = colorButton('Apply && Close', GREEN)
        self.btnApply.clicked.connect(lambda *_: self.applyBG())
        row.addWidget(self.btnApply)
        row.addWidget(button('Cancel', self.close))
        gl.addLayout(row, 2, 0, 1, 3)
        self.edN.valueChanged.connect(lambda *_: self.drawPreview())
        self.drawPreview()

    def _layout(self):
        W, H = self.cv.sizePt()
        self.cv.place(self.ax, 54.0, 36.0, W - 54.0 - 16.0, H - 36.0 - 24.0)

    def drawPreview(self):
        app = self.app
        N = self.edN.intValue()
        if N < 1 or N > app.delay.size:
            return
        ax = self.ax
        ax.clear()
        styleAxes(ax, grid=True)
        cmap = colormaps['viridis'].resampled(max(N, 2))
        for k in range(N):
            c = cmap(k)
            ax.plot(app.wavelength, app.deltaA_raw[:, k], '-',
                    color=(c[0], c[1], c[2], 0.6), linewidth=0.7)
        bg = meanSpectrum(app.deltaA_raw, N)
        hMean, = ax.plot(app.wavelength, bg, 'k-', linewidth=2,
                         label='Mean of first %d point(s)' % N)
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        legend(ax, handles=[hMean], fontsize=FONT_PT)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel(Z_LABEL)
        ax.set_title('Early-time spectra and their mean (background)',
                     fontsize=FONT_PT + 1)
        sciTicks(ax, 'y')
        self.lblInfo.setText('Using t = %.4g ... %.4g %s' % (
            app.delay[0], app.delay[N - 1], app.tUnitTxt()))
        self.cv.draw_idle()

    def applyBG(self):
        if self.applyBusy:
            return
        app = self.app
        N = self.edN.intValue()
        if N < 1 or N > app.delay.size:
            return
        self.applyBusy = True
        try:
            app.bgSpectrum = meanSpectrum(app.deltaA_raw, N)
            app.bgN = N
            app.bgApplied = True
            app.updateAll()
        finally:
            self.applyBusy = False
        self.close()
