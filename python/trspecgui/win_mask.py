"""Mask Wavelength Regions window. Port of openMaskWindow."""
import numpy as np
from matplotlib.patches import Rectangle
from PySide6 import QtCore, QtWidgets

from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label)
from .win_base import BLUE, ChildWindow, absMax, mapLayout, mapOn


class MaskWindow(ChildWindow):
    KEY = 'mask'

    def __init__(self, app):
        super().__init__(app, 'Mask Wavelength Regions', 960, 680)
        self.pickStage = 1
        u = app.xUnitTxt()
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(12, 12, 12, 12)
        v.setSpacing(8)
        info = label(
            'Mask wavelength regions to exclude from analysis. Click once '
            'on the 2D map to pick the start wavelength (From), click again '
            'for the end (To). Then choose a mode and press Add region. '
            '"Set to NaN" removes the region from fits and shows it as a '
            'transparent gap; "Set to 0" forces values to 0 while keeping '
            'them in the wavelength axis (not recommended for fitting).',
            italic=True, color=(0.3, 0.3, 0.3))
        info.setWordWrap(True)
        v.addWidget(info)
        self.cv = PlotCanvas()
        v.addWidget(self.cv, 1)
        self.ax = self.cv.figure.add_axes([0.1, 0.1, 0.75, 0.8])
        self.cax = self.cv.figure.add_axes([0.9, 0.1, 0.02, 0.8])
        self.cv.layoutFcn = lambda: mapLayout(self.cv, self.ax, self.cax)
        self.cv.mpl_connect('button_press_event', self.onMapClick)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('From (%s):' % u, right=True))
        self.edWL1 = NumEdit(0, width=90)
        row.addWidget(self.edWL1)
        row.addWidget(label('To (%s):' % u, right=True))
        self.edWL2 = NumEdit(0, width=90)
        row.addWidget(self.edWL2)
        row.addWidget(label('Mode:', right=True))
        self.ddMode = QtWidgets.QComboBox()
        self.ddMode.addItem('Set to NaN (mask out)', 'nan')
        self.ddMode.addItem('Set to 0', 'zero')
        self.ddMode.setMinimumWidth(170)
        row.addWidget(self.ddMode)
        b = colorButton('Add region', BLUE, width=120)
        b.clicked.connect(lambda *_: self.doAdd())
        row.addWidget(b)
        row.addStretch(1)
        self.lblPick = label('Click map to pick From / To', italic=True,
                             color=(0.45, 0.45, 0.45))
        row.addWidget(self.lblPick)
        v.addLayout(row)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)
        self.tbl = QtWidgets.QTableWidget(0, 3)
        self.tbl.setHorizontalHeaderLabels(
            ['From (%s)' % u, 'To (%s)' % u, 'Mode'])
        self.tbl.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tbl.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.tbl.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.tbl.setColumnWidth(0, 100)
        self.tbl.setColumnWidth(1, 100)
        self.tbl.setColumnWidth(2, 200)
        self.tbl.setFixedHeight(180)
        row.addWidget(self.tbl, 1)
        col = QtWidgets.QVBoxLayout()
        col.setSpacing(6)
        col.addWidget(button('Remove selected', self.doRemove, width=150))
        col.addWidget(button('Clear all', self.doClear, width=150))
        col.addStretch(1)
        row.addLayout(col)
        v.addLayout(row)
        self.edWL1.valueChanged.connect(lambda *_: self.updatePickStage())
        self.edWL2.valueChanged.connect(lambda *_: self.updatePickStage())
        self.drawMap()
        self.refreshTable()

    def drawMap(self):
        app = self.app
        wl, t = app.wavelength, app.delay
        mapOn(app, self.ax, self.cax, wl, t, app.deltaA, absMax(app.deltaA))
        ax = self.ax
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        ax.set_ylim(float(np.min(t)), float(np.max(t)))
        ax.grid(True, color='0.15', alpha=0.15, linewidth=0.6)
        ax.set_title('Masked regions preview  (click to pick wavelengths)',
                     fontsize=FONT_PT + 1)
        y0, y1 = float(np.min(t)), float(np.max(t))
        for w1, w2, mode in app.maskedRegions:
            c = (0.2, 0.5, 0.8) if mode == 'zero' else (0.8, 0.2, 0.2)
            ax.add_patch(Rectangle((w1, y0), w2 - w1, y1 - y0,
                                   facecolor=c + (0.18,), edgecolor=c,
                                   linewidth=1.2))
        self.cv.relayout()
        self.cv.draw_idle()

    def refreshTable(self):
        regs = self.app.maskedRegions
        self.tbl.setRowCount(len(regs))
        for k, (w1, w2, mode) in enumerate(regs):
            vals = ['%.6g' % w1, '%.6g' % w2,
                    'Set to NaN' if mode == 'nan' else 'Set to 0']
            for j, s in enumerate(vals):
                it = QtWidgets.QTableWidgetItem(s)
                if j < 2:
                    it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight |
                                        QtCore.Qt.AlignmentFlag.AlignVCenter)
                self.tbl.setItem(k, j, it)
        self.tbl.clearSelection()

    def onMapClick(self, ev):
        if ev.inaxes is not self.ax or ev.button != 1 or ev.xdata is None:
            return
        self.pickWavelength(float(ev.xdata))

    def pickWavelength(self, wlClick):
        if self.pickStage == 1:
            self.edWL1.setValue(wlClick)
            self.pickStage = 2
            self.lblPick.setText('Next click sets the To (end) wavelength')
        else:
            self.edWL2.setValue(wlClick)
            self.pickStage = 1
            self.lblPick.setText(
                'Next click restarts picking from the From edge')

    def updatePickStage(self):
        self.pickStage = 1
        self.lblPick.setText('Click map to pick From / To')

    def doAdd(self):
        app = self.app
        wl1, wl2 = self.edWL1.value(), self.edWL2.value()
        if wl1 > wl2:
            wl1, wl2 = wl2, wl1
        if not (np.isfinite(wl1) and np.isfinite(wl2)) or wl1 >= wl2:
            self.alert('Invalid wavelength range. From must be < To.',
                       'Input error')
            return
        wlLo = float(np.min(app.wavelength))
        wlHi = float(np.max(app.wavelength))
        if wl2 < wlLo or wl1 > wlHi:
            self.alert('Region [%.1f, %.1f] %s is outside data range '
                       '[%.1f, %.1f] %s.' % (wl1, wl2, app.xUnitTxt(), wlLo,
                                             wlHi, app.xUnitTxt()),
                       'Input error')
            return
        app.maskedRegions.append((float(wl1), float(wl2),
                                  self.ddMode.currentData()))
        self.refreshTable()
        app.updateAll()
        self.drawMap()
        self.pickStage = 1
        self.lblPick.setText('Added region %d.  Click map to pick the next '
                             'From / To.' % len(app.maskedRegions))

    def selectedRow(self):
        rows = self.tbl.selectionModel().selectedRows()
        return rows[0].row() if rows else None

    def doRemove(self):
        app = self.app
        r = self.selectedRow()
        if r is None or r >= len(app.maskedRegions):
            self.alert('Select a row in the table first.', 'Nothing selected')
            return
        del app.maskedRegions[r]
        self.refreshTable()
        app.updateAll()
        self.drawMap()

    def doClear(self):
        app = self.app
        if not app.maskedRegions:
            return
        app.maskedRegions = []
        self.refreshTable()
        app.updateAll()
        self.drawMap()
