"""Single-Wavelength Kinetic Fit window ("Fit trace..."). Port of
openKineticFitWindow."""
import os

import numpy as np
from PySide6 import QtCore, QtWidgets

from trspec.kinetics import irf_mode_norm
from trspec.single import fit_single_trace

from .state import ui_str
from .textfit import legend, setTitle
from .widgets import (FONT_PT, DataTable, NumEdit, PlotCanvas, button,
                      colorButton, label, panel, sciTicks, sidePanel,
                      styleAxes, textArea, tip)
from .win_base import GREEN, ChildWindow



def adjustVec(v0, N, pad):
    v = list(v0)
    return v + [pad] * (N - len(v)) if len(v) < N else v[:N]


def numText(v):
    """A number as writetable prints it."""
    if np.isnan(v):
        return 'NaN'
    if np.isinf(v):
        return 'Inf' if v > 0 else '-Inf'
    return '%.15g' % v


def writeTable(path, header, columns, delim=','):
    """Columns of numbers with a header line (port of the writetable
    calls of the export functions)."""
    cols = [np.asarray(c, float).ravel() for c in columns]
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(delim.join(header) + '\n')
        for i in range(cols[0].size):
            fh.write(delim.join(numText(c[i]) for c in cols) + '\n')


class KineticFitWindow(ChildWindow):
    KEY = 'kinfit'

    def __init__(self, app):
        super().__init__(app, 'Single-Wavelength Kinetic Fit', 1180, 720)
        u = app.tUnitTxt()
        self.lastFit = self._emptyFit()
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(8)

        left, lv = sidePanel('Fit setup', 350)
        gl.addWidget(left)

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('%s (%s):' % (app.xSymTxt(), app.xUnitTxt())))
        self.edWL = NumEdit(app.selWL, width=90,
                            limits=(float(np.min(app.wavelength)),
                                    float(np.max(app.wavelength))))
        row.addWidget(self.edWL)
        row.addWidget(button('Use crosshair', self.useCrosshair))
        row.addStretch(1)
        lv.addLayout(row)

        # a band instead of one wavelength
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        self.cbBand = QtWidgets.QCheckBox('Band, to:')
        self.cbBand.setChecked(bool(app.kfBandOn))
        self.cbBand.setToolTip(tip(
            'Fit the kinetics of a band instead of one wavelength: from '
            'the value above to the value given here.'))
        row.addWidget(self.cbBand)
        self.edB2 = NumEdit(self.bandEndStart(), width=80,
                            limits=(float(np.min(app.wavelength)),
                                    float(np.max(app.wavelength))))
        row.addWidget(self.edB2)
        self.ddBand = QtWidgets.QComboBox()
        self.ddBand.addItem('Integral', 'integral')
        self.ddBand.addItem('Mean', 'mean')
        self.ddBand.setCurrentIndex(
            max(self.ddBand.findData(app.kfBandMode), 0))
        self.ddBand.setToolTip(tip(
            'Integral: integral of the signal over the spectral axis of '
            'the band (trapezoid). Mean: mean over the wavelengths of the '
            'band, in the unit of the signal. On an evenly spaced axis '
            'both give the same time constants. Wavelengths without data '
            '(masked) are left out; a delay at which one of the others '
            'has no value stays empty.'))
        row.addWidget(self.ddBand)
        row.addStretch(1)
        lv.addLayout(row)
        self.edB2.setEnabled(self.cbBand.isChecked())
        self.ddBand.setEnabled(self.cbBand.isChecked())

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('Components:'))
        self.ddN = QtWidgets.QComboBox()
        self.ddN.addItems(['1', '2', '3', '4', '5'])
        self.ddN.setCurrentText(str(int(app.kfNComp)))
        self.ddN.setMinimumWidth(60)
        row.addWidget(self.ddN)
        self.cbInf = QtWidgets.QCheckBox(
            ui_str('Include \\tau = \\infty'))
        self.cbInf.setChecked(bool(app.kfHasInf))
        row.addWidget(self.cbInf)
        row.addStretch(1)
        lv.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('Trace scale:'))
        self.ddScale = QtWidgets.QComboBox()
        self.ddScale.addItems(['Linear', 'Log'])
        self.ddScale.setCurrentText('Log')
        self.ddScale.setMinimumWidth(90)
        row.addWidget(self.ddScale)
        row.addStretch(1)
        lv.addLayout(row)

        lv.addWidget(label('Initial time constants (%s, positive):' % u,
                           italic=True, color=(0.3, 0.3, 0.3)))
        self.tauTable = DataTable(
            ui_str(['\\tau_{init} (%s)' % u, '\\tau fixed', '\\beta_{init}',
                    'Stretched']),
            ['num', 'bool', 'num', 'bool'], widths=[85, 70, 70, 70])
        self.tauTable.setFixedHeight(150)
        lv.addWidget(self.tauTable)

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('IRF t_0 (%s):' % u))
        self.edT0 = NumEdit(app.kfT0, width=80)
        row.addWidget(self.edT0)
        self.cbT0fix = QtWidgets.QCheckBox('Fixed')
        self.cbT0fix.setChecked(bool(app.kfT0Fixed))
        row.addWidget(self.cbT0fix)
        row.addStretch(1)
        lv.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('IRF FWHM (%s):' % u))
        self.edFW = NumEdit(app.kfFWHM, width=80, limits=(1e-6, np.inf))
        row.addWidget(self.edFW)
        self.cbFWfix = QtWidgets.QCheckBox('Fixed')
        self.cbFWfix.setChecked(bool(app.kfFWHMFixed))
        row.addWidget(self.cbFWfix)
        row.addStretch(1)
        lv.addLayout(row)

        app.kfIRFTreat = irf_mode_norm(app.kfIRFTreat)
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('IRF model:'))
        self.ddIRF = QtWidgets.QComboBox()
        self.ddIRF.addItem('None (pure exponentials)', 'none')
        self.ddIRF.addItem('Numerical conv.', 'numerical')
        self.ddIRF.setCurrentIndex(self.ddIRF.findData(app.kfIRFTreat))
        self.ddIRF.setMinimumWidth(165)
        self.ddIRF.setToolTip(tip(
            'Numerical conv.: the model is convolved with the Gaussian IRF '
            '(t0, FWHM). Standard exp components use the closed-form erfcx '
            'expression, stretched components a numerical quadrature at '
            'each delay.\nNone: no IRF convolution -- fit a plain sum of '
            'exponentials exp(-(t-t0)/tau)*H(t-t0). Combine with "Fit t '
            'range" to fit only the post-IRF part of the trace; t0 and '
            'FWHM are excluded from the fit.'))
        row.addWidget(self.ddIRF)
        row.addStretch(1)
        lv.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(4)
        row.addWidget(label('Fit t range (%s):' % u))
        self.edFitTmin = NumEdit(float(np.min(app.delay)), width=62,
                                 tooltip='Lower bound of the fitted delay '
                                         'range.')
        row.addWidget(self.edFitTmin)
        row.addWidget(label('~'))
        self.edFitTmax = NumEdit(float(np.max(app.delay)), width=62,
                                 tooltip='Upper bound of the fitted delay '
                                         'range.')
        row.addWidget(self.edFitTmax)
        row.addWidget(button('Full', self.setFitTFull,
                             tooltip='Reset to the full delay range.'))
        row.addStretch(1)
        lv.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(5)
        self.btnRun = colorButton('Run Fit', GREEN)
        self.btnRun.clicked.connect(lambda *_: self.doRun())
        row.addWidget(self.btnRun, 1)
        b = button('Reset', self.doReset)
        row.addWidget(b, 1)
        b = button('Export...', self.doExport)
        row.addWidget(b, 1)
        lv.addLayout(row)

        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        lv.addWidget(self.lblStatus)
        self.txtResult = textArea(['Run a fit to see results.'])
        self.txtResult.setMinimumHeight(110)
        lv.addWidget(self.txtResult, 1)

        right = panel('Trace and fit')
        gl.addWidget(right, 1)
        rv = QtWidgets.QVBoxLayout(right)
        rv.setContentsMargins(6, 10, 6, 6)
        self.cv = PlotCanvas()
        rv.addWidget(self.cv)
        fig = self.cv.figure
        self.axTrace = fig.add_axes([0.1, 0.45, 0.85, 0.5])
        self.axRes = fig.add_axes([0.1, 0.1, 0.85, 0.25], sharex=self.axTrace)
        self.cv.layoutFcn = self._layout

        self.fillTable()
        self.onIRFModeChange()
        lo, hi = self.clampFitRange(app.kfFitTMin, app.kfFitTMax)
        self.edFitTmin.setValue(lo)
        self.edFitTmax.setValue(hi)
        self.drawTraceOnly()
        sv = app.kfSaved
        if isinstance(sv, dict) and sv.get('lastFit') is not None and \
                len(sv['lastFit'].get('fitVec', [])) and \
                len(sv['lastFit']['t']) == app.delay.size:
            try:
                self.lastFit = dict(sv['lastFit'])
                self.edWL.setValue(self.lastFit['wlActual'])
                fi = self.fitInfo()
                self.cbBand.setChecked(bool(fi['band']))
                if fi['band']:
                    self.edB2.setValue(fi['w2'])
                    self.ddBand.setCurrentIndex(
                        max(self.ddBand.findData(fi['mode']), 0))
                self.edB2.setEnabled(self.cbBand.isChecked())
                self.ddBand.setEnabled(self.cbBand.isChecked())
                self.drawWithFit()
                self.txtResult.setPlainText('\n'.join(sv['lines']))
                self.lblStatus.setText('%s  [restored]' % sv['status'])
            except Exception:
                self.lblStatus.setText('Stored result could not be restored '
                                       '-- run again.')
        self.ddN.currentTextChanged.connect(lambda *_: self.onNCompChange())
        self.ddScale.currentTextChanged.connect(
            lambda *_: self.refreshAxesScale(True))
        self.edWL.valueChanged.connect(lambda *_: self.drawTraceOnly())
        self.cbBand.toggled.connect(lambda *_: self.onBandChange())
        self.edB2.valueChanged.connect(lambda *_: self.onBandChange())
        self.ddBand.currentIndexChanged.connect(
            lambda *_: self.onBandChange())
        self.ddIRF.currentIndexChanged.connect(
            lambda *_: self.onIRFModeChange())

    # ------------------------------------------------------------------
    @staticmethod
    def _emptyFit():
        return {'wlActual': np.nan, 'iw': -1, 'data': np.zeros(0),
                't': np.zeros(0), 'mask': np.zeros(0, bool),
                'tau': np.zeros(0), 'beta': np.zeros(0), 'stretchOn': [],
                't0': np.nan, 'fwhm': np.nan, 'hasInf': False,
                'A': np.zeros(0), 'fitVec': np.zeros(0), 'res': np.zeros(0),
                'rms': np.nan}

    def _layout(self):
        cv = self.cv
        W, H = cv.sizePt()
        L, R, T, B, gap = 58.0, 16.0, 22.0, 34.0, 54.0
        avail = max(H - T - B - gap, 20.0)
        hRes = avail / 3.0
        cv.place(self.axRes, L, B, W - L - R, hRes)
        cv.place(self.axTrace, L, B + hRes + gap, W - L - R, avail - hRes)

    def clampFitRange(self, lo, hi):
        app = self.app
        dLo, dHi = float(np.min(app.delay)), float(np.max(app.delay))
        if lo is None or hi is None or not np.isfinite(lo) or \
                not np.isfinite(hi):
            return dLo, dHi
        if lo > hi:
            lo, hi = hi, lo
        lo, hi = max(lo, dLo), min(hi, dHi)
        if not (hi > lo):
            return dLo, dHi
        return float(lo), float(hi)

    def fitRangeToStore(self, lo, hi):
        app = self.app
        dLo, dHi = float(np.min(app.delay)), float(np.max(app.delay))
        span = dHi - dLo
        a, b = min(lo, hi), max(lo, hi)
        if not np.isfinite(a) or not np.isfinite(b) or span <= 0 or \
                ((a - dLo) <= 1e-9 * span and (dHi - b) <= 1e-9 * span):
            return None, None
        return float(a), float(b)

    # ------------------------------------------------------------------
    def fillTable(self):
        app = self.app
        N = int(app.kfNComp)
        app.kfTauInit = adjustVec(app.kfTauInit, N, 1.0)
        app.kfTauFixed = adjustVec(app.kfTauFixed, N, False)
        app.kfBetaInit = adjustVec(app.kfBetaInit, N, 1.0)
        app.kfBetaFixed = adjustVec(app.kfBetaFixed, N, True)
        app.kfStretchOn = adjustVec(app.kfStretchOn, N, False)
        self.tauTable.setData(
            [[float(app.kfTauInit[i]), bool(app.kfTauFixed[i]),
              float(app.kfBetaInit[i]), bool(app.kfStretchOn[i])]
             for i in range(N)])

    def readTable(self):
        app = self.app
        d = self.tauTable.data()
        app.kfTauInit = [float(r[0]) for r in d]
        app.kfTauFixed = [bool(r[1]) for r in d]
        app.kfBetaInit = [float(r[2]) for r in d]
        app.kfStretchOn = [bool(r[3]) for r in d]
        app.kfBetaFixed = [not bool(r[3]) for r in d]

    def onNCompChange(self):
        self.readTable()
        self.app.kfNComp = int(self.ddN.currentText())
        self.fillTable()

    def onIRFModeChange(self):
        isNone = self.ddIRF.currentData() == 'none'
        self.edFW.setEnabled(not isNone)
        self.cbFWfix.setEnabled(not isNone)
        if isNone:
            self.cbT0fix.setChecked(True)
        self.cbT0fix.setEnabled(not isNone)

    def setFitTFull(self):
        self.edFitTmin.setValue(float(np.min(self.app.delay)))
        self.edFitTmax.setValue(float(np.max(self.app.delay)))

    def getFitTRangeIfRestricted(self):
        app = self.app
        lo = min(self.edFitTmin.value(), self.edFitTmax.value())
        hi = max(self.edFitTmin.value(), self.edFitTmax.value())
        if not np.isfinite(lo) or not np.isfinite(hi):
            return None
        dLo, dHi = float(np.min(app.delay)), float(np.max(app.delay))
        span = dHi - dLo
        if span <= 0:
            return None
        if (lo - dLo) > 1e-9 * span or (dHi - hi) > 1e-9 * span:
            return [lo, hi]
        return None

    def useCrosshair(self):
        self.edWL.setValue(self.app.selWL)
        self.drawTraceOnly()

    def bandEndStart(self):
        """Other end of the band when the window opens: the value used
        last, or 5 % of the spectral axis beside the wavelength."""
        app = self.app
        lo, hi = float(np.min(app.wavelength)), float(np.max(app.wavelength))
        v = app.kfBandTo
        if v is not None and np.isfinite(v) and lo <= v <= hi:
            return float(v)
        step = 0.05 * (hi - lo)
        v = app.selWL + step
        if v > hi:
            v = app.selWL - step
        return float(min(max(v, lo), hi))

    def onBandChange(self):
        app = self.app
        on = self.cbBand.isChecked()
        self.edB2.setEnabled(on)
        self.ddBand.setEnabled(on)
        app.kfBandOn = bool(on)
        app.kfBandTo = float(self.edB2.value())
        app.kfBandMode = self.ddBand.currentData()
        self.drawTraceOnly()

    def traceNow(self):
        """(d, info): the trace the controls ask for."""
        return self.app.kinTrace(self.edWL.value(), self.cbBand.isChecked(),
                                 self.edB2.value(), self.ddBand.currentData())

    def fitInfo(self):
        """Description of the trace of the last fit (a result kept by
        an earlier version has none: one wavelength)."""
        fi = self.lastFit.get('info')
        if not isinstance(fi, dict) or 'band' not in fi:
            fi = self.app.kinTrace(self.lastFit['wlActual'], False)[1]
        return fi

    def refreshAxesScale(self, redraw=False):
        app = self.app
        isLog = self.ddScale.currentText() == 'Log' and \
            bool(np.any(app.delay > 0))
        if isLog:
            pos = app.delay[app.delay > 0]
            xl = [float(np.min(pos)), float(np.max(app.delay))]
        else:
            xl = [float(np.min(app.delay)), float(np.max(app.delay))]
        self.axTrace.set_xscale('log' if isLog else 'linear')
        fr = self.getFitTRangeIfRestricted()
        if fr is not None:
            if isLog:
                fr[0] = max(fr[0], xl[0])
            if fr[0] < fr[1]:
                xl = fr
        if xl[1] > xl[0]:
            self.axTrace.set_xlim(xl[0], xl[1])
        if redraw:
            self.cv.draw_idle()

    def _styleBoth(self, info):
        app = self.app
        for ax in (self.axTrace, self.axRes):
            styleAxes(ax, grid=True)
            ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
            sciTicks(ax, 'y')
        yl = app.kinTraceLabel(info, 'ylabel')
        self.axTrace.set_ylabel(yl)
        self.axRes.set_ylabel('Residual (%s)' % yl)

    def drawTraceOnly(self):
        app = self.app
        wlReq = self.edWL.value()
        d, info = self.traceNow()
        ax = self.axTrace
        ax.clear()
        self.axRes.clear()
        if not info['err']:
            ax.plot(app.delay, d, 'k.', markersize=5.5, label='Data')
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        self._styleBoth(info)
        if info['err']:
            setTitle(ax, info['err'], 'No trace', 'No trace')
        elif info['band']:
            at = app.kinTraceLabel(info, 'at')
            setTitle(ax, 'Trace of the ' + at, 'Trace of the band', 'Band')
            legend(ax, fontsize=FONT_PT)
        else:
            at = app.kinTraceLabel(info, 'at')
            setTitle(ax, 'Trace at %s  (requested %.2f %s)' % (
                at, wlReq, app.xUnitAx()), 'Trace at ' + at, at)
            legend(ax, fontsize=FONT_PT)
        self.axRes.set_title('Residual (no fit yet)', fontsize=FONT_PT + 1)
        self.refreshAxesScale()
        self.cv.draw_idle()

    def drawWithFit(self):
        app = self.app
        lf = self.lastFit
        ax = self.axTrace
        ax.clear()
        self.axRes.clear()
        ax.plot(lf['t'], lf['data'], 'k.', markersize=5, label='Data')
        ax.plot(lf['t'], lf['fitVec'], 'r-', linewidth=1.6, label='Fit')
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        self.axRes.plot(lf['t'], lf['res'], '-', color=(0.45, 0.45, 0.45),
                        linewidth=0.9)
        self.axRes.axhline(0, color='k', linestyle=':', linewidth=0.8)
        fi = self.fitInfo()
        self._styleBoth(fi)
        at = app.kinTraceLabel(fi, 'at')
        if fi['band']:
            setTitle(ax, 'Fit of the %s  (RMS = %.3g)' % (at, lf['rms']),
                     'Fit of the band  (RMS = %.3g)' % lf['rms'], 'Band')
        else:
            setTitle(ax, 'Fit at %s  (RMS = %.3g)' % (at, lf['rms']),
                     'Fit at ' + at, at)
        legend(ax, fontsize=FONT_PT)
        self.axRes.set_title('Residual (RMS = %.3g)' % lf['rms'],
                             fontsize=FONT_PT + 1)
        self.refreshAxesScale()
        self.cv.draw_idle()

    # ------------------------------------------------------------------
    def doReset(self):
        app = self.app
        app.kfNComp = 3
        app.kfTauInit = [1.0, 10.0, 100.0]
        app.kfTauFixed = [False, False, False]
        app.kfBetaInit = [1.0, 1.0, 1.0]
        app.kfBetaFixed = [True, True, True]
        app.kfStretchOn = [False, False, False]
        app.kfHasInf = False
        app.kfT0, app.kfT0Fixed = 0.0, True
        app.kfFWHM, app.kfFWHMFixed = 0.15, True
        app.kfIRFTreat = 'numerical'
        b = QtCore.QSignalBlocker(self.ddN)
        self.ddN.setCurrentText('3')
        del b
        self.cbInf.setChecked(False)
        self.edT0.setValue(0.0)
        self.cbT0fix.setChecked(True)
        self.edFW.setValue(0.15)
        self.cbFWfix.setChecked(True)
        self.ddIRF.setCurrentIndex(self.ddIRF.findData('numerical'))
        wasBand = self.cbBand.isChecked()
        b1 = QtCore.QSignalBlocker(self.cbBand)
        b2 = QtCore.QSignalBlocker(self.ddBand)
        self.cbBand.setChecked(False)
        self.ddBand.setCurrentIndex(self.ddBand.findData('integral'))
        del b1, b2
        self.edB2.setEnabled(False)
        self.ddBand.setEnabled(False)
        app.kfBandOn = False
        app.kfBandMode = 'integral'
        app.kfFitTMin = None
        app.kfFitTMax = None
        self.setFitTFull()
        self.fillTable()
        self.onIRFModeChange()
        if wasBand:
            self.drawTraceOnly()
        self.lblStatus.setText('Reset to defaults.')
        app.kfSaved = None

    def doRun(self):
        app = self.app
        self.readTable()
        app.kfHasInf = self.cbInf.isChecked()
        app.kfT0 = self.edT0.value()
        app.kfT0Fixed = self.cbT0fix.isChecked()
        app.kfFWHM = self.edFW.value()
        app.kfFWHMFixed = self.cbFWfix.isChecked()
        app.kfIRFTreat = self.ddIRF.currentData()
        if app.kfIRFTreat == 'none':
            app.kfT0Fixed = True
            app.kfFWHMFixed = True
        tA, tB = self.edFitTmin.value(), self.edFitTmax.value()
        app.kfFitTMin, app.kfFitTMax = self.fitRangeToStore(tA, tB)
        tSelChk = (app.delay >= min(tA, tB)) & (app.delay <= max(tA, tB))
        if np.count_nonzero(tSelChk) < app.kfNComp + 2:
            self.alert('Fit t range contains too few delay points. Widen '
                       'the range.', 'Invalid range')
            return
        if any(v <= 0 for v in app.kfTauInit):
            self.alert('All tau values must be positive.', 'Invalid input')
            return
        if any(v <= 0 or v > 2 for v in app.kfBetaInit):
            self.alert('beta_init must be in (0, 2]. Use 1 for ordinary '
                       'exponential.', 'Invalid input')
            return
        if app.kfFWHM <= 0 and app.kfIRFTreat != 'none':
            self.alert('IRF FWHM must be positive.', 'Invalid input')
            return
        d, info = self.traceNow()
        if info['err']:
            self.alert(info['err'], 'No data')
            return
        wlActual, iw = info['wl'], info['iw']
        if not np.any(np.isfinite(d)):
            self.alert('There is no delay at which every wavelength of the '
                       'band holds data.' if info['band'] else
                       'Trace at this wavelength is all NaN.', 'No data')
            return
        self.lblStatus.setText('Fitting...')
        QtWidgets.QApplication.setOverrideCursor(
            QtCore.Qt.CursorShape.WaitCursor)
        QtWidgets.QApplication.processEvents()
        tRange = None
        if app.kfFitTMin is not None:
            tRange = [app.kfFitTMin, app.kfFitTMax]
        try:
            r = fit_single_trace(
                app.delay, d, app.kfTauInit, app.kfTauFixed, app.kfBetaInit,
                app.kfBetaFixed, app.kfStretchOn, app.kfT0, app.kfT0Fixed,
                app.kfFWHM, app.kfFWHMFixed, app.kfHasInf, app.kfIRFTreat,
                tRange)
        except Exception as err:
            QtWidgets.QApplication.restoreOverrideCursor()
            self.alert('Fit failed:\n%s' % err, 'Fit Error')
            self.lblStatus.setText('Fit failed.')
            self.lastFit = self._emptyFit()
            return
        QtWidgets.QApplication.restoreOverrideCursor()
        lf = self.lastFit = {
            'wlActual': wlActual, 'iw': iw, 'info': dict(info),
            'data': np.array(d, float),
            't': app.delay.copy(), 'mask': np.asarray(r.mask, bool),
            'tau': np.asarray(r.tau, float),
            'beta': np.asarray(r.beta, float),
            'stretchOn': list(app.kfStretchOn), 't0': float(r.t0),
            'fwhm': float(r.fwhm), 'hasInf': bool(app.kfHasInf),
            'A': np.asarray(r.amp, float).ravel(),
            'fitVec': np.asarray(r.fit, float),
            'res': np.asarray(r.residual, float), 'rms': float(r.rms)}
        self.drawWithFit()
        lines = self.resultLines(r, lf)
        self.txtResult.setPlainText('\n'.join(lines))
        self.lblStatus.setText('Done. RMS = %.3g' % r.rms)
        app.kfSaved = {'lastFit': lf, 'status': self.lblStatus.text(),
                       'lines': lines}

    def resultLines(self, r, lf):
        app = self.app
        u = app.tUnitTxt()
        tauOut, betaOut = lf['tau'], lf['beta']
        L = list(app.kinTraceLabel(self.fitInfo(), 'lines'))
        L.append('Iterations : %d' % r.iters)
        L.append('Final   RMS: %.4g' % r.rms)
        if r.exit_flag <= 0:
            L += ['', '** NOT CONVERGED **']
        L += ['', 'Time constants & beta:']
        for i in range(tauOut.size):
            tag = ''
            if app.kfStretchOn[i]:
                tag = ' [stretched, beta=%.3f]' % betaOut[i]
                if np.isfinite(r.beta_std[i]):
                    tag = ' [stretched, beta=%.3f+/-%.3f]' % (
                        betaOut[i], r.beta_std[i])
            if app.kfTauFixed[i]:
                tag += ' (tau fixed)'
            tag += undetNote(r.tau_std, i)
            if np.isfinite(r.tau_std[i]):
                L.append('  tau_%d = %10.4g +/- %.2g %s%s' % (
                    i + 1, tauOut[i], r.tau_std[i], u, tag))
            else:
                L.append('  tau_%d = %10.4g %s%s' % (i + 1, tauOut[i], u, tag))
        if app.kfHasInf:
            L.append('  tau_inf  (offset)')
        L += ['', 'Amplitudes:']
        nMain = tauOut.size
        for i in range(nMain):
            L.append('  A_%d   = %12.4g' % (i + 1, lf['A'][i]))
        if app.kfHasInf:
            L.append('  A_inf = %12.4g' % lf['A'][nMain])
        L.append('')
        if np.isfinite(r.t0_std):
            L.append('IRF t0   = %.4g +/- %.2g %s' % (r.t0, r.t0_std, u))
        else:
            L.append('IRF t0   = %.4g %s%s' % (
                r.t0, u, (' (fixed)' if app.kfT0Fixed else '')
                + undetNote(r.t0_std, 0)))
        if r.fwhm_at_limit:
            L.append('IRF FWHM = %.4g %s  [not resolved]' % (r.fwhm, u))
        elif getattr(r, 'fwhm_too_wide', False):
            L.append('IRF FWHM = %.4g %s  [not determined by the data]'
                     % (r.fwhm, u))
        elif np.isfinite(r.fwhm_std):
            L.append('IRF FWHM = %.4g +/- %.2g %s' % (r.fwhm, r.fwhm_std, u))
        else:
            L.append('IRF FWHM = %.4g %s%s' % (
                r.fwhm, u, (' (fixed)' if app.kfFWHMFixed else '')
                + undetNote(r.fwhm_std, 0)))
        L += irfLimitNote(app, r)
        L.append('IRF model: ' + ('none (pure exponentials)'
                                  if app.kfIRFTreat == 'none'
                                  else 'numerical conv.'))
        return L

    # ------------------------------------------------------------------
    def doExport(self):
        app = self.app
        lf = self.lastFit
        if not len(lf['fitVec']):
            self.alert('Run a fit first.', 'No results')
            return
        fi = self.fitInfo()
        full = app.putFile(
            'CSV (*.csv);;DAT (tab) (*.dat)',
            app.kinTraceLabel(fi, 'export'),
            app.kinTraceLabel(fi, 'file') + '.csv',
            parent=self)
        if not full:
            return
        dlm = ',' if full.lower().endswith('.csv') else '\t'
        try:
            writeTable(full, ['delay_%s' % app.tUnitHdr(), 'data_dA',
                              'fit_dA', 'residual_dA'],
                       [lf['t'], lf['data'], lf['fitVec'], lf['res']], dlm)
            paramPath = os.path.splitext(full)[0] + '_params.txt'
            uh = app.tUnitHdr()
            with open(paramPath, 'w', encoding='utf-8', newline='\n') as fh:
                fh.write('Single-trace kinetic fit\n')
                fh.write(app.kinTraceLabel(fi, 'param') + '\n')
                fh.write('RMS       : %.6g\n\n' % lf['rms'])
                fh.write('Time constants:\n')
                for i in range(len(lf['tau'])):
                    if i < len(lf['stretchOn']) and lf['stretchOn'][i]:
                        fh.write('  tau_%d = %.6g %s   beta_%d = %.4f  '
                                 '[stretched]\n' % (i + 1, lf['tau'][i], uh,
                                                    i + 1, lf['beta'][i]))
                    else:
                        fh.write('  tau_%d = %.6g %s\n' % (
                            i + 1, lf['tau'][i], uh))
                fh.write('\nAmplitudes:\n')
                nMain = len(lf['tau'])
                for i in range(nMain):
                    fh.write('  A_%d   = %.6g\n' % (i + 1, lf['A'][i]))
                if lf['hasInf']:
                    fh.write('  A_inf = %.6g\n' % lf['A'][nMain])
                fh.write('\nIRF t0   = %.6g %s\n' % (lf['t0'], uh))
                fh.write('IRF FWHM = %.6g %s\n' % (lf['fwhm'], uh))
            self.inform('Saved trace to:\n  %s\n\nParameters in:\n  %s' % (
                full, paramPath), 'Exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')


def undetNote(std, k):
    """Remark for a parameter whose uncertainty is infinite."""
    try:
        v = np.atleast_1d(np.asarray(std, float))
        if v.size > k and np.isinf(v[k]):
            return '  [not determined by the data]'
    except Exception:
        pass
    return ''


def irfLimitNote(app, r):
    """Lines of the result box about a fitted IRF width that the data
    do not resolve or that is wider than half the range of the delays,
    or that came out of the fit in two steps. Port of irfLimitNote."""
    L = []
    if getattr(r, 'fwhm_at_limit', False):
        lim = getattr(r, 'fwhm_lo', np.nan)
        L.append('  FWHM is close to, or below, half the')
        L.append('  delay step at time zero (%.3g %s):' % (lim, app.tUnitTxt()))
        L.append('  the delays do not resolve such a')
        L.append('  width. Fix it at a known value.')
    if getattr(r, 'fwhm_too_wide', False):
        lim = getattr(r, 'fwhm_hi', np.nan)
        # (three digits; a whole number from 100 on)
        limTxt = '%d' % int(np.floor(lim + 0.5)) if lim >= 100 else '%.3g' % lim
        L.append('  FWHM is above half the range of the')
        L.append('  delays (%s %s): the data do not' % (limTxt, app.tUnitTxt()))
        L.append('  determine the IRF width. Fix it at a')
        L.append('  known value.')
    if not L and getattr(r, 'irf_retry', 0) == 1:
        L.append('  The direct fit lost the IRF width;')
        L.append('  this is the result of a fit in two')
        L.append('  steps (IRF fixed, then free).')
    return L
