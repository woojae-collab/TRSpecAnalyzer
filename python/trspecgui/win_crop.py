"""Crop Data window. Port of openCropWindow."""
import numpy as np
from matplotlib.patheffects import withStroke
from PySide6 import QtWidgets

from .state import ui_str
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label)
from .win_base import GREEN, ChildWindow, absMax, mapLayout, mapOn


class CropWindow(ChildWindow):
    KEY = 'crop'

    def __init__(self, app):
        super().__init__(app, 'Crop Data', 820, 670)
        self.applyBusy = False
        self.wlFull = (float(np.min(app.wavelength)),
                       float(np.max(app.wavelength)))
        self.tFull = (float(np.min(app.delay)), float(np.max(app.delay)))
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(10, 10, 10, 10)
        v.setSpacing(8)
        g = QtWidgets.QGridLayout()
        g.setHorizontalSpacing(5)
        g.setVerticalSpacing(3)
        g.setColumnStretch(6, 1)
        g.addWidget(label('%s range (%s):' % (app.xSymTxt(), app.xUnitTxt()),
                          right=True), 0, 0)
        self.edWLmin = NumEdit(self.wlFull[0], width=90)
        self.edWLmax = NumEdit(self.wlFull[1], width=90)
        g.addWidget(self.edWLmin, 0, 1)
        g.addWidget(label('to'), 0, 2)
        g.addWidget(self.edWLmax, 0, 3)
        g.addWidget(button(ui_str('Full \\lambda range'), self.setFullWL,
                           width=130), 0, 5)
        g.addWidget(label('t range (%s):' % app.tUnitTxt(), right=True), 1, 0)
        self.edTmin = NumEdit(self.tFull[0], width=90)
        self.edTmax = NumEdit(self.tFull[1], width=90)
        g.addWidget(self.edTmin, 1, 1)
        g.addWidget(label('to'), 1, 2)
        g.addWidget(self.edTmax, 1, 3)
        g.addWidget(button('Full t range', self.setFullT, width=130), 1, 5)
        v.addLayout(g)
        self.cv = PlotCanvas()
        v.addWidget(self.cv, 1)
        self.ax = self.cv.figure.add_axes([0.1, 0.1, 0.75, 0.8])
        self.cax = self.cv.figure.add_axes([0.9, 0.1, 0.02, 0.8])
        self.cv.layoutFcn = lambda: mapLayout(self.cv, self.ax, self.cax)
        row = QtWidgets.QHBoxLayout()
        row.addStretch(1)
        row.addWidget(button('Cancel', self.close))
        self.btnApply = colorButton('Apply Crop', GREEN)
        self.btnApply.clicked.connect(lambda *_: self.doApply())
        row.addWidget(self.btnApply)
        v.addLayout(row)
        for ed in (self.edWLmin, self.edWLmax, self.edTmin, self.edTmax):
            ed.valueChanged.connect(lambda *_: self.drawPreview())
        self.drawPreview()

    def ranges(self):
        wl1, wl2 = self.edWLmin.value(), self.edWLmax.value()
        t1, t2 = self.edTmin.value(), self.edTmax.value()
        if wl1 > wl2:
            wl1, wl2 = wl2, wl1
        if t1 > t2:
            t1, t2 = t2, t1
        return wl1, wl2, t1, t2

    def drawPreview(self):
        app = self.app
        wl1, wl2, t1, t2 = self.ranges()
        wl, t = app.wavelength, app.delay
        C = app.deltaA
        W, T = np.meshgrid(wl, t)
        inside = (W >= wl1) & (W <= wl2) & (T >= t1) & (T <= t2)
        alpha = np.where(inside, 1.0, 0.2)[:-1, :-1]
        mapOn(app, self.ax, self.cax, wl, t, C, absMax(C), alpha=alpha)
        ax = self.ax
        pe = [withStroke(linewidth=2.6, foreground=(0, 0, 0, 0.55))]
        for x in (wl1, wl2):
            ax.axvline(x, color='w', linewidth=1.2, path_effects=pe)
        for y in (t1, t2):
            ax.axhline(y, color='w', linewidth=1.2, path_effects=pe)
        nWL = int(np.count_nonzero((wl >= wl1) & (wl <= wl2)))
        nT = int(np.count_nonzero((t >= t1) & (t <= t2)))
        ax.set_title('Kept: %s = [%.1f, %.1f] %s  (%d pts)  |  '
                     't = [%.3g, %.3g] %s  (%d pts)' % (
                         app.xSym(), wl1, wl2, app.xUnitAx(), nWL,
                         t1, t2, app.tUnitAx(), nT), fontsize=FONT_PT + 1)
        self.cv.relayout()
        self.cv.draw_idle()

    def setFullWL(self):
        self.edWLmin.setValue(self.wlFull[0])
        self.edWLmax.setValue(self.wlFull[1])
        self.drawPreview()

    def setFullT(self):
        self.edTmin.setValue(self.tFull[0])
        self.edTmax.setValue(self.tFull[1])
        self.drawPreview()

    def doApply(self):
        if self.applyBusy:
            return
        self.applyBusy = True
        try:
            self.app.applyCropByRange(
                self.edWLmin.value(), self.edWLmax.value(),
                self.edTmin.value(), self.edTmax.value())
        finally:
            self.applyBusy = False
        self.close()
