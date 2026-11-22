"""SVD Analysis window. Port of openSvdWindow."""
import os

import numpy as np
from PySide6 import QtWidgets

from .state import ui_str, write_data_file, write_matrix
from .textfit import legend
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, label, lineColor,
                      panel, sciTicks, sidePanel, styleAxes)
from .win_base import ChildWindow, mapLayout, mapOn


class SvdWindow(ChildWindow):
    KEY = 'svd'

    def __init__(self, app):
        super().__init__(app, 'SVD Analysis', 1280, 820)
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)
        left, lv = sidePanel('SVD setup', 330)
        lv.setSpacing(8)
        gl.addWidget(left)

        def numRow(text, ed):
            row = QtWidgets.QHBoxLayout()
            row.setSpacing(6)
            row.addWidget(label(text, right=True))
            row.addWidget(ed)
            row.addStretch(1)
            lv.addLayout(row)

        self.edN = NumEdit(3, limits=(1, 50), width=90, integer=True)
        numRow('Components to keep (N):', self.edN)
        self.edShow = NumEdit(3, limits=(1, 12), width=90, integer=True)
        numRow('Basis vectors to show:', self.edShow)
        self.ddTScale = QtWidgets.QComboBox()
        self.ddTScale.addItems(['Linear', 'Log'])
        self.ddTScale.setCurrentText('Log')
        self.ddTScale.setMinimumWidth(90)
        numRow('Temporal axis scale:', self.ddTScale)
        self.infoLbl = label('', color=(0.3, 0.3, 0.3))
        self.infoLbl.setWordWrap(True)
        lv.addWidget(self.infoLbl)
        sigPanel = panel('Singular values (log scale)')
        lv.addWidget(sigPanel, 1)
        sl = QtWidgets.QVBoxLayout(sigPanel)
        sl.setContentsMargins(4, 10, 4, 4)
        self.cvSig = PlotCanvas()
        sl.addWidget(self.cvSig)
        self.axSig = self.cvSig.figure.add_axes([0.2, 0.15, 0.75, 0.8])
        self.cvSig.layoutFcn = lambda: self._lay(self.cvSig, self.axSig, 12.0)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(button('Export U / S / V', self.exportAll))
        row.addWidget(button('Export reconstruction',
                             self.exportReconstruction))
        row.addStretch(1)
        lv.addLayout(row)

        right = panel('Bases and reconstruction')
        gl.addWidget(right, 1)
        rg = QtWidgets.QGridLayout(right)
        rg.setContentsMargins(10, 14, 10, 10)
        rg.setHorizontalSpacing(10)
        rg.setVerticalSpacing(8)
        rg.setRowStretch(0, 10)
        rg.setRowStretch(1, 10)
        rg.setRowStretch(2, 12)

        def sub(title, r, c, span=1, isMap=False):
            p = panel(title)
            rg.addWidget(p, r, c, 1, span)
            pl = QtWidgets.QVBoxLayout(p)
            pl.setContentsMargins(4, 10, 4, 4)
            cv = PlotCanvas()
            pl.addWidget(cv)
            ax = cv.figure.add_axes([0.15, 0.15, 0.7, 0.75])
            cax = None
            if isMap:
                cax = cv.figure.add_axes([0.9, 0.15, 0.02, 0.75])
                cv.layoutFcn = lambda: mapLayout(cv, ax, cax, T=20.0, R=58.0)
            else:
                cv.layoutFcn = lambda: self._lay(cv, ax, 10.0)
            return cv, ax, cax

        self.cvU, self.axU, _ = sub(
            'Spectral basis  U(:, k)  vs wavelength', 0, 0)
        self.cvV, self.axV, _ = sub('Temporal basis  V(:, k)  vs delay', 0, 1)
        self.cvExp, self.axExp, self.caxExp = sub('Original data', 1, 0,
                                                  isMap=True)
        self.cvRec, self.axRec, self.caxRec = sub(
            'Reconstruction (N components)', 1, 1, isMap=True)
        self.cvRes, self.axRes, self.caxRes = sub(
            'Residual = original - reconstruction', 2, 0, span=2, isMap=True)

        D = app.deltaA
        self.wl = app.wavelength.copy()
        self.t = app.delay.copy()
        nanMask = ~np.isfinite(D)
        self.nNaN = int(np.count_nonzero(nanMask))
        self.Dz = np.where(nanMask, 0.0, D)
        self.shape = D.shape
        self.ok = True
        try:
            U, s, Vh = np.linalg.svd(self.Dz, full_matrices=False)
        except Exception as err:
            self.ok = False
            self.alert('SVD failed:\n%s\n\nTry reducing data size (crop) or '
                       'check for non-finite values.' % err, 'SVD error')
            return
        self.U, self.sigma, self.V = U, s, Vh.T
        self.nMax = s.size
        sv = app.svdSaved
        if isinstance(sv, dict):
            try:
                self.edN.setValue(min(int(sv['N']), self.nMax))
                self.edShow.setValue(min(int(sv['show']), self.nMax))
                self.ddTScale.setCurrentText(sv['tscale'])
            except Exception:
                pass
        self.edN.valueChanged.connect(lambda *_: self.refreshAll())
        self.edShow.valueChanged.connect(lambda *_: self.refreshAll())
        self.ddTScale.currentTextChanged.connect(lambda *_: self.refreshAll())
        self.refreshAll()

    @staticmethod
    def _lay(cv, ax, top):
        W, H = cv.sizePt()
        cv.place(ax, 52.0, 32.0, W - 52.0 - 16.0, H - 32.0 - top)

    def nKeep(self):
        return min(self.edN.intValue(), self.nMax)

    def reconstruction(self, nKeep):
        return (self.U[:, :nKeep] * self.sigma[:nKeep]) @ self.V[:, :nKeep].T

    def refreshAll(self):
        app = self.app
        sigma, wl, t = self.sigma, self.wl, self.t
        nKeep = self.nKeep()
        nShow = min(self.edShow.intValue(), self.nMax)
        app.svdSaved = {'N': nKeep, 'show': nShow,
                        'tscale': self.ddTScale.currentText()}
        ax = self.axSig
        ax.clear()
        styleAxes(ax, grid=True)
        kAll = np.arange(1, min(self.nMax, max(20, 2 * nShow)) + 1)
        ax.semilogy(kAll, sigma[kAll - 1], 'o-', markersize=5,
                    linewidth=1.2, markerfacecolor='none')
        ax.semilogy(np.arange(1, nShow + 1), sigma[:nShow], 'o',
                    color=(0.20, 0.45, 0.80), markersize=6)
        ax.axvline(nKeep + 0.5, color='r', linestyle='--', linewidth=1.0)
        ax.text(nKeep + 0.5, 0.97, ' N=%d' % nKeep, color='r',
                transform=ax.get_xaxis_transform(), ha='left', va='top',
                fontsize=FONT_PT)
        ax.set_xlabel('Component index')
        ax.set_ylabel('$\\sigma_i$')
        self.cvSig.draw_idle()

        ax = self.axU
        ax.clear()
        styleAxes(ax, grid=True)
        for k in range(nShow):
            ax.plot(wl, self.U[:, k], '-', linewidth=1.2, color=lineColor(k),
                    label='k=%d  ($\\sigma$=%.3g)' % (k + 1, sigma[k]))
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('$U_k$')
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        legend(ax, fontsize=8)
        self.cvU.draw_idle()

        scale = self.ddTScale.currentText().lower()
        ax = self.axV
        ax.clear()
        styleAxes(ax, grid=True)
        for k in range(nShow):
            ax.plot(t, self.V[:, k], '-', linewidth=1.2, color=lineColor(k),
                    label='k=%d' % (k + 1))
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel('$V_k$')
        self._applyScale(ax, 'x', scale)
        legend(ax, fontsize=8)
        self.cvV.draw_idle()

        Dn = self.reconstruction(nKeep)
        Rn = self.Dz - Dn
        cl = float(np.max(np.abs(self.Dz))) if self.Dz.size else 1.0
        if not np.isfinite(cl) or cl == 0:
            cl = 1.0
        clr = float(np.max(np.abs(Rn))) if Rn.size else cl
        if not np.isfinite(clr) or clr == 0:
            clr = cl
        for cv, ax, cax, C, c in ((self.cvExp, self.axExp, self.caxExp,
                                   self.Dz, cl),
                                  (self.cvRec, self.axRec, self.caxRec,
                                   Dn, cl),
                                  (self.cvRes, self.axRes, self.caxRes,
                                   Rn, clr)):
            mapOn(app, ax, cax, wl, t, C, c)
            cax.yaxis.label.set_visible(False)
            ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
            self._applyScale(ax, 'y', scale)
            cv.relayout()
            cv.draw_idle()
        varExpl = float(np.sum(sigma[:nKeep] ** 2) / np.sum(sigma ** 2) * 100)
        self.infoLbl.setText(ui_str(
            'Data: %d wavelengths x %d delays\n'
            'NaN entries replaced with 0: %d\n'
            'N=%d components -> %.4f %% variance explained\n'
            '\\sigma_1 / \\sigma_N = %.3g / %.3g\n'
            'Tip: use the knee of the \\sigma plot as the starting '
            'Components value in Global Analysis.' % (
                self.shape[0], self.shape[1], self.nNaN, nKeep, varExpl,
                sigma[0], sigma[nKeep - 1])))

    def _applyScale(self, ax, which, scale):
        t = self.t
        pos = t[t > 0]
        lo, hi = float(np.min(t)), float(np.max(t))
        isLog = scale == 'log' and pos.size > 0
        if isLog:
            lo = float(np.min(pos))
        if which == 'x':
            ax.set_xscale('log' if isLog else 'linear')
            if hi > lo:
                ax.set_xlim(lo, hi)
        else:
            ax.set_yscale('log' if isLog else 'linear')
            if hi > lo:
                ax.set_ylim(lo, hi)

    # ------------------------------------------------------------------
    def exportAll(self):
        app = self.app
        fn = app.putFile('CSV files (*.csv)',
                         'Choose base filename for U / S / V exports',
                         'svd_result.csv', parent=self)
        if not fn:
            return
        outDir = os.path.dirname(fn)
        base = os.path.splitext(os.path.basename(fn))[0]
        nKeep = self.nKeep()
        try:
            head = np.concatenate([[np.nan], np.arange(1, nKeep + 1)])
            write_matrix(os.path.join(outDir, base + '_U.csv'), np.vstack(
                [head, np.column_stack([self.wl, self.U[:, :nKeep]])]))
            write_matrix(os.path.join(outDir, base + '_V.csv'), np.vstack(
                [head, np.column_stack([self.t, self.V[:, :nKeep]])]))
            write_matrix(os.path.join(outDir, base + '_sigma.csv'),
                         np.column_stack([np.arange(1, self.nMax + 1),
                                          self.sigma]))
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')
            return
        self.inform('Saved:\n  %s_U.csv\n  %s_V.csv\n  %s_sigma.csv\n\n'
                    'Folder:\n%s' % (base, base, base, outDir),
                    'Export complete')

    def exportReconstruction(self):
        app = self.app
        fn = app.putFile('CSV files (*.csv)', 'Save reconstruction',
                         'svd_reconstruction.csv', parent=self)
        if not fn:
            return
        try:
            write_data_file(fn, self.wl, self.t,
                            self.reconstruction(self.nKeep()), ',')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')
            return
        self.inform('Saved: ' + fn, 'Export complete')
