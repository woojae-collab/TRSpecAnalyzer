"""Lifetime Density Analysis window. Port of openLdaWindow."""
import numpy as np
from PySide6 import QtWidgets

from trspec.lda import compute_lcurve, compute_lda

from .state import ui_str
from .textfit import legend, setTitle
from .widgets import (AX_BG, FONT_PT, BlitManager, NumEdit, PlotCanvas,
                      VScroll, button, colorButton, label, panel, runBusy,
                      sciTicks, sidePanel, styleAxes, textArea, tip)
from .win_base import ChildWindow, absMax, mapLayout
from .win_kinfit import writeTable

EPS = float(np.finfo(float).eps)
RUN_COLOR = (0.20, 0.55, 0.40)
BLUE = (0.0, 0.0, 1.0)
IDLE_TEXT = ['Click "Run LDA" to compute.']
T_MAP = 'Amplitude density A(λ, τ)'
T_1D = 'Click on the map to see A(τ)'
T_SPEC = 'Spectrum at clicked τ'
T_KIN = 'Reconstruction at clicked λ'
T_L = 'L-curve (run "L-curve..." to populate)'


def emptyResult():
    return {'Amap': None, 'Drec': None, 'tauGrid': None, 'alpha': np.nan,
            'regType': '', 'rms': np.nan}


class LdaWindow(ChildWindow):
    KEY = 'lda'

    def __init__(self, app):
        super().__init__(app, 'Lifetime Density Analysis', 1380, 900)
        u = app.tUnitTxt()
        self.res = emptyResult()
        self.curLambda = float(app.selWL)
        self.curTau = np.nan
        self.xhV = None
        self.xhH = None
        self.dragging = False
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: setup ==================================
        left, lv = sidePanel('LDA setup', 360)
        gl.addWidget(left)

        def bold(text):
            lb = label(text)
            lb.setStyleSheet('font-weight: bold;')
            lv.addWidget(lb)

        def pair(text, field, stretch=True):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(6)
            r.addWidget(label(text))
            r.addWidget(field, 1 if stretch else 0)
            if not stretch:
                r.addStretch(1)
            lv.addLayout(r)

        bold('Lifetime grid (log-spaced):')
        self.edTauMin = NumEdit(max(0.5 * app.gaFWHM, EPS),
                                limits=(EPS, np.inf))
        pair(ui_str('\\tau_{min} (%s):' % u), self.edTauMin)
        tmaxD = float(np.max(app.delay))
        self.tspan = tmaxD - float(np.min(app.delay))
        if self.tspan <= 0:
            self.tspan = 1000.0
        self.edTauMax = NumEdit(self._tauMaxDefault(), limits=(EPS, np.inf))
        pair(ui_str('\\tau_{max} (%s):' % u), self.edTauMax)
        self.edK = NumEdit(100, limits=(10, 500), integer=True)
        pair('Grid points (K):', self.edK)
        bold('IRF parameters:')
        t0Init = app.gaResultT0 if app.gaResultT0 is not None else app.gaT0
        self.edT0 = NumEdit(t0Init)
        pair('t_0 (%s):' % u, self.edT0)
        fwInit = app.gaResultFWHM if app.gaResultFWHM is not None \
            else app.gaFWHM
        self.edFW = NumEdit(max(fwInit, EPS), limits=(EPS, np.inf))
        pair('FWHM (%s):' % u, self.edFW)
        self.cbNoIRF = QtWidgets.QCheckBox('No IRF (pure exponentials)')
        self.cbNoIRF.setToolTip(tip(
            'Basis becomes exp(-(t-t0)/tau)*H(t-t0) without IRF '
            'convolution; FWHM is ignored. Intended for data (or a Fit t '
            'range) starting after the IRF.'))
        lv.addWidget(self.cbNoIRF)
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(4)
        r.addWidget(label('Fit t range (%s):' % u))
        self.edFitTmin = NumEdit(
            float(np.min(app.delay)), width=58,
            tooltip='Delay points outside this range are excluded from the '
                    'LDA solve (reconstruction is still shown on the full '
                    'axis).\nNOTE: when restricting to t >= T, tau values '
                    'much shorter than T are no longer constrained by the '
                    'data -- raise tau_min accordingly and do not interpret '
                    'the fast-tau part of the map.')
        r.addWidget(self.edFitTmin)
        r.addWidget(label('~'))
        self.edFitTmax = NumEdit(
            tmaxD, width=58,
            tooltip='Upper bound of the fitted delay range.')
        r.addWidget(self.edFitTmax)
        r.addWidget(button('Full', self.setFitTFull,
                           tooltip='Reset to the full delay range.'))
        r.addStretch(1)
        lv.addLayout(r)
        bold('Regularization:')
        self.ddReg = QtWidgets.QComboBox()
        self.ddReg.addItem(ui_str('L2 (||A||^2)'), 'l2')
        self.ddReg.addItem(ui_str('1st-deriv (smooth in log-\\tau)'),
                           'l2deriv')
        self.ddReg.setCurrentIndex(1)
        pair('Type:', self.ddReg)
        self.edAlpha = NumEdit(1e-2, limits=(0.0, np.inf), width=100)
        pair(ui_str('\\alpha:'), self.edAlpha, stretch=False)
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        self.btnRun = colorButton('Run LDA', RUN_COLOR)
        self.btnRun.setStyleSheet(self.btnRun.styleSheet()
                                  + 'QPushButton { font-weight: bold; }')
        self.btnLcurve = button(
            'L-curve...',
            tooltip='Sweep alpha from 1e-6 to 1e2 and find the corner of '
                    'the L-curve. The corner alpha is then loaded into the '
                    'alpha field.')
        self.btnReset = button('Reset')
        for b in (self.btnRun, self.btnLcurve, self.btnReset):
            r.addWidget(b, 1)
        lv.addLayout(r)
        bold('Export:')
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        self.btnExpAmap = button(ui_str('A(\\lambda,\\tau)'))
        self.btnExpRec = button('Recon.')
        self.btnExp1D = button(ui_str('A(\\lambda*)'))
        for b in (self.btnExpAmap, self.btnExpRec, self.btnExp1D):
            b.setEnabled(False)
            r.addWidget(b, 1)
        lv.addLayout(r)
        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        lv.addWidget(self.lblStatus)
        self.txtResult = textArea(IDLE_TEXT)
        self.txtResult.setLineWrapMode(      # (sentences, not columns)
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.txtResult.setMinimumHeight(110)
        lv.addWidget(self.txtResult, 1)

        # ================= right: results ===============================
        right = panel('Results')
        gl.addWidget(right, 1)
        ro = QtWidgets.QVBoxLayout(right)
        ro.setContentsMargins(2, 12, 2, 4)
        rin = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(rin)
        grid.setContentsMargins(4, 0, 4, 2)
        grid.setSpacing(8)
        ro.addWidget(VScroll(rin))
        self.cvMap = PlotCanvas()
        self.cv1D = PlotCanvas()
        self.cvSpec = PlotCanvas()
        self.cvKin = PlotCanvas()
        self.cvL = PlotCanvas()
        self.cvMap.setMinimumHeight(240)
        self.cv1D.setMinimumHeight(240)
        for cv in (self.cvSpec, self.cvKin, self.cvL):
            cv.setMinimumHeight(165)
        grid.addWidget(self.cvMap, 0, 0)
        grid.addWidget(self.cv1D, 0, 1)
        grid.addWidget(self.cvSpec, 1, 0)
        grid.addWidget(self.cvKin, 1, 1)
        grid.addWidget(self.cvL, 2, 0, 1, 2)
        grid.setRowStretch(0, 3)
        grid.setRowStretch(1, 2)
        grid.setRowStretch(2, 2)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        self.axMap = self.cvMap.figure.add_subplot(111)
        self.caxMap = self.cvMap.figure.add_axes([0.9, 0.1, 0.02, 0.8])
        self.ax1D = self.cv1D.figure.add_subplot(111)
        self.axSpec = self.cvSpec.figure.add_subplot(111)
        self.axKin = self.cvKin.figure.add_subplot(111)
        self.axL = self.cvL.figure.add_subplot(111)

        def lay1(cv, ax, R=22.0):
            # (R: about as wide as the blank left of the label of the
            # y-axis. A(lambda) under the map has the wavelength axis
            # of the map and therefore its width: R = the room of the
            # color bar)
            W, H = cv.sizePt()
            cv.place(ax, 58.0, 34.0, W - 58.0 - R, H - 34.0 - 22.0)

        self.cvMap.layoutFcn = self._layMap
        self.cv1D.layoutFcn = lambda: lay1(self.cv1D, self.ax1D)
        self.cvSpec.layoutFcn = lambda: lay1(self.cvSpec, self.axSpec, 78.0)
        self.cvKin.layoutFcn = lambda: lay1(self.cvKin, self.axKin)
        self.cvL.layoutFcn = lambda: lay1(self.cvL, self.axL)
        self.blMap = BlitManager(self.cvMap)
        self.cvMap.mpl_connect('button_press_event', self.onMapPress)
        self.cvMap.mpl_connect('motion_notify_event', self.onMapMotion)
        self.cvMap.mpl_connect('button_release_event',
                               lambda ev: self.onMapUp())
        self._emptyPlots()

        # ---- what the controls do --------------------------------------
        self.cbNoIRF.toggled.connect(lambda *_: self.onNoIRFChange())
        self.btnRun.clicked.connect(
            lambda *_: self.safe(self.doRun, 'Run LDA'))
        self.btnLcurve.clicked.connect(
            lambda *_: self.safe(self.doLcurve, 'L-curve'))
        self.btnReset.clicked.connect(
            lambda *_: self.safe(self.doReset, 'Reset'))
        self.btnExpAmap.clicked.connect(
            lambda *_: self.safe(self.exportAmap, 'Export A-map'))
        self.btnExpRec.clicked.connect(
            lambda *_: self.safe(self.exportRec, 'Export recon'))
        self.btnExp1D.clicked.connect(
            lambda *_: self.safe(self.export1D, 'Export 1D'))

        # ---- result of an earlier visit --------------------------------
        sv = app.ldaSaved
        if isinstance(sv, dict) and isinstance(sv.get('result'), dict) and \
                sv['result'].get('Amap') is not None and \
                np.shape(sv['result']['Amap'])[0] == app.wavelength.size:
            try:
                st = sv['settings']
                self.edTauMin.setValue(st['tauMin'])
                self.edTauMax.setValue(st['tauMax'])
                self.edK.setValue(st['K'])
                self.edT0.setValue(st['t0'])
                self.edFW.setValue(st['fwhm'])
                self.ddReg.setCurrentIndex(0 if st['reg'] == 'l2' else 1)
                self.edAlpha.setValue(st['alpha'])
                if 'noIRF' in st:
                    self.cbNoIRF.setChecked(bool(st['noIRF']))
                    self.onNoIRFChange()
                lo, hi = st.get('fitTmin'), st.get('fitTmax')
                if lo is not None and np.isfinite(lo):
                    self.edFitTmin.setValue(max(float(np.min(app.delay)),
                                                lo))
                if hi is not None and np.isfinite(hi):
                    self.edFitTmax.setValue(min(tmaxD, hi))
                r0 = sv['result']
                self.res = {'Amap': np.asarray(r0['Amap'], float),
                            'Drec': np.asarray(r0['Drec'], float),
                            'tauGrid': np.asarray(r0['tauGrid'],
                                                  float).ravel(),
                            'alpha': float(r0['alpha']),
                            'regType': r0['regType'],
                            'rms': float(r0['rms'])}
                cl, ct = sv['curLambda'], sv['curTau']
                self.curLambda = np.nan if cl is None else float(cl)
                self.curTau = np.nan if ct is None else float(ct)
                self._setExport(True)
                self.plotResults()
                self.lblStatus.setText('%s  [restored]' % sv['status'])
                self.txtResult.setPlainText('\n'.join(sv['lines']))
            except Exception:
                self.res = emptyResult()
                self._setExport(False)
                self._emptyPlots()
                self.lblStatus.setText(
                    'Stored result could not be restored -- run again.')

    # ------------------------------------------------------------------
    def _tauMaxDefault(self):
        v = 3 * float(np.max(self.app.delay))
        if not np.isfinite(v) or v <= 0:
            v = 5 * self.tspan
        return v

    def _layMap(self):
        mapLayout(self.cvMap, self.axMap, self.caxMap, L=58.0, R=78.0,
                  T=22.0, B=34.0)

    def _setExport(self, on):
        for b in (self.btnExpAmap, self.btnExpRec, self.btnExp1D):
            b.setEnabled(bool(on))

    def _tauLabel(self):
        return 'Lifetime τ (%s)' % self.app.tUnitAx()

    def _emptyPlots(self):
        """The five plots before anything was computed."""
        app = self.app
        for ax in (self.axMap, self.ax1D, self.axSpec, self.axKin, self.axL):
            ax.clear()
        self.caxMap.clear()
        self.caxMap.set_visible(False)
        styleAxes(self.axMap)
        self.axMap.set_yscale('log')
        self.axMap.set_xlabel(app.xAxLabel())
        self.axMap.set_ylabel(self._tauLabel())
        setTitle(self.axMap, T_MAP)
        styleAxes(self.ax1D, grid=True)
        self.ax1D.set_xscale('log')
        self.ax1D.set_xlabel(self._tauLabel())
        self.ax1D.set_ylabel('Amplitude')
        setTitle(self.ax1D, T_1D)
        styleAxes(self.axSpec, grid=True)
        self.axSpec.set_xlabel(app.xAxLabel())
        self.axSpec.set_ylabel('Amplitude')
        setTitle(self.axSpec, T_SPEC)
        styleAxes(self.axKin, grid=True)
        self.axKin.set_xlabel('Delay time (%s)' % app.tUnitAx())
        self.axKin.set_ylabel(self.app.sigName())
        setTitle(self.axKin, T_KIN)
        self._emptyL()
        self.xhV = self.xhH = None
        self.blMap.setArtists([])
        self._layMap()
        for cv in (self.cvMap, self.cv1D, self.cvSpec, self.cvKin, self.cvL):
            cv.draw_idle()

    def _emptyL(self):
        ax = self.axL
        ax.clear()
        styleAxes(ax, grid=True)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('||residual||')
        ax.set_ylabel('||solution||')
        setTitle(ax, T_L, 'L-curve')
        self.cvL.draw_idle()

    def buildGrid(self):
        return np.logspace(np.log10(self.edTauMin.value()),
                           np.log10(self.edTauMax.value()),
                           int(self.edK.value()))

    def irfMode(self):
        return 'none' if self.cbNoIRF.isChecked() else 'irf'

    # ------------------------------------------------------------------
    def doRun(self):
        app = self.app
        if self.edTauMin.value() <= 0 or \
                self.edTauMax.value() <= self.edTauMin.value():
            self.alert('Need 0 < tau_min < tau_max.', 'Invalid')
            return
        if self.edFW.value() <= 0:
            self.alert('IRF FWHM must be positive.', 'Invalid')
            return
        tauGrid = self.buildGrid()
        self.lblStatus.setText('Running LDA...')
        D, t = app.deltaA, app.delay
        t0, fw, al = self.edT0.value(), self.edFW.value(), \
            self.edAlpha.value()
        reg, irfM = self.ddReg.currentData(), self.irfMode()
        tR = [self.edFitTmin.value(), self.edFitTmax.value()]
        try:
            o = runBusy(app, self, 'Running LDA...', lambda: compute_lda(
                D, t, tauGrid, t0, fw, al, reg, irfM, tR))
        except Exception as err:
            self.alert('LDA failed:\n%s' % err, 'LDA error')
            self.lblStatus.setText('LDA failed.')
            return
        self.res = {'Amap': np.asarray(o.amap, float),
                    'Drec': np.asarray(o.drec, float), 'tauGrid': tauGrid,
                    'alpha': float(al), 'regType': reg, 'rms': float(o.rms)}
        self._setExport(True)
        self.lblStatus.setText('Done. RMS = %.3g (rel. %.2g%%)'
                               % (o.rms, 100 * o.ratio))
        self.txtResult.setPlainText('\n'.join([
            'K = %d log-spaced taus' % o.K,
            'tau range: %.3g .. %.3g %s' % (tauGrid[0], tauGrid[-1],
                                            app.tUnitTxt()),
            'alpha = %.3g  (%s)' % (al, reg),
            'RMS residual: %.4g' % o.rms,
            'Relative residual: %.3f%%' % (100 * o.ratio),
            '',
            'Click on the heatmap to inspect a particular wavelength.']))
        self.plotResults()
        app.ldaSaved = {
            'result': dict(self.res),
            'curLambda': float(self.curLambda), 'curTau': float(self.curTau),
            'status': self.lblStatus.text(),
            'lines': self.txtResult.toPlainText().split('\n'),
            'settings': {
                'tauMin': self.edTauMin.value(),
                'tauMax': self.edTauMax.value(),
                'K': int(self.edK.value()),
                't0': self.edT0.value(), 'fwhm': self.edFW.value(),
                'reg': reg, 'alpha': float(al),
                'noIRF': self.cbNoIRF.isChecked(),
                'fitTmin': self.edFitTmin.value(),
                'fitTmax': self.edFitTmax.value()}}

    def plotResults(self):
        app = self.app
        r = self.res
        if r['Amap'] is None:
            return
        A = r['Amap']
        cl = absMax(A)
        wl = np.asarray(app.wavelength, float).ravel()
        tauG = np.asarray(r['tauGrid'], float).ravel()
        ax, cax = self.axMap, self.caxMap
        ax.clear()
        cax.clear()
        styleAxes(ax)
        ax.set_facecolor(AX_BG)
        h = None
        if wl.size > 1 and tauG.size > 1:
            # (each cell takes the color of its lower left corner, as
            # MATLAB's surface with flat faces does)
            h = ax.pcolormesh(wl, tauG, np.ma.masked_invalid(A.T[:-1, :-1]),
                              shading='flat', cmap=app.getColormap(),
                              vmin=-cl, vmax=cl, rasterized=True)
        ax.set_yscale('log')
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        ax.set_ylim(float(np.min(tauG)), float(np.max(tauG)))
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel(self._tauLabel())
        setTitle(ax, 'Amplitude density A(λ, τ)   α=%.3g' % r['alpha'],
                 'A(λ, τ)   α=%.3g' % r['alpha'])
        cax.set_visible(h is not None)
        if h is not None:
            cb = self.cvMap.figure.colorbar(h, cax=cax)
            cb.set_label(self.app.sigName() + '  amplitude', fontsize=FONT_PT)
            cb.ax.tick_params(labelsize=FONT_PT, length=3)
            sciTicks(cb.ax, 'y')
        if not np.isfinite(self.curLambda):
            self.curLambda = 0.5 * (float(np.min(wl)) + float(np.max(wl)))
        if not np.isfinite(self.curTau):
            self.curTau = float(np.sqrt(np.min(tauG) * np.max(tauG)))
        kw = dict(color='k', linewidth=1.0, alpha=0.7)
        self.xhV = ax.axvline(self.curLambda, **kw)
        self.xhH = ax.axhline(self.curTau, **kw)
        self.blMap.setArtists([self.xhV, self.xhH])
        self._layMap()
        self.cvMap.draw_idle()
        self.update1DandKin()

    # ---- the point picked on the map ---------------------------------
    def onMapPress(self, ev):
        if ev.button != 1 or ev.inaxes is not self.axMap or \
                self.res['Amap'] is None:
            return
        self.dragging = True
        self.applyPoint(ev.xdata, ev.ydata)

    def onMapMotion(self, ev):
        if not self.dragging or ev.x is None or ev.y is None:
            return
        try:
            x, y = self.axMap.transData.inverted().transform((ev.x, ev.y))
        except Exception:
            return
        self.applyPoint(x, y)

    def onMapUp(self):
        self.dragging = False

    def applyPoint(self, x, y):
        if x is None or y is None or not np.isfinite(x) or \
                not np.isfinite(y) or self.res['Amap'] is None:
            return
        wl = self.app.wavelength
        tg = self.res['tauGrid']
        self.curLambda = float(min(max(x, np.min(wl)), np.max(wl)))
        self.curTau = float(min(max(y, np.min(tg)), np.max(tg)))
        if self.xhV is not None:
            self.xhV.set_xdata([self.curLambda, self.curLambda])
        if self.xhH is not None:
            self.xhH.set_ydata([self.curTau, self.curTau])
        self.blMap.update()
        self.update1DandKin()

    def update1DandKin(self):
        app = self.app
        r = self.res
        if r['Amap'] is None:
            return
        wl = np.asarray(app.wavelength, float).ravel()
        tg = r['tauGrid']
        iw = int(np.argmin(np.abs(wl - self.curLambda)))
        wlAct = wl[iw]
        ax = self.ax1D
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(tg, r['Amap'][iw, :], '-', color=BLUE, linewidth=1.4)
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        if np.isfinite(self.curTau):
            ax.axvline(self.curTau, color='k', linestyle='--',
                       linewidth=0.8, alpha=0.6)
        ax.set_xlabel(self._tauLabel())
        ax.set_ylabel('Amplitude')
        ax.set_xscale('log')
        ax.set_xlim(float(np.min(tg)), float(np.max(tg)))
        setTitle(ax, 'A(τ) at %s = %.2f %s' % (app.xSym(), wlAct,
                                               app.xUnitAx()))
        sciTicks(ax, 'y')
        self.cv1D.draw_idle()

        ax = self.axSpec
        ax.clear()
        styleAxes(ax, grid=True)
        if np.isfinite(self.curTau):
            kTau = int(np.argmin(np.abs(tg - self.curTau)))
            ax.plot(wl, r['Amap'][:, kTau], '-', color=BLUE, linewidth=1.4)
            ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
            if np.isfinite(self.curLambda):
                ax.axvline(self.curLambda, color='k', linestyle='--',
                           linewidth=0.8, alpha=0.6)
            ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
            setTitle(ax, 'A(λ) at τ = %.3g %s' % (tg[kTau], app.tUnitTxt()))
            sciTicks(ax, 'y')
        else:
            setTitle(ax, T_SPEC)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('Amplitude')
        self.cvSpec.draw_idle()

        ax = self.axKin
        ax.clear()
        styleAxes(ax, grid=True)
        t = np.asarray(app.delay, float).ravel()
        ax.plot(t, app.deltaA[iw, :], 'k.', markersize=4, label='Data')
        ax.plot(t, r['Drec'][iw, :], 'r-', linewidth=1.4,
                label='LDA reconstr.')
        ax.plot(t, app.deltaA[iw, :] - r['Drec'][iw, :], '-',
                color=(0.5, 0.5, 0.5), linewidth=0.8, label='Residual')
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel(self.app.sigName())
        if np.any(t > 0):
            ax.set_xscale('log')
            xl = [float(np.min(t[t > 0])), float(np.max(t))]
            fr = self.getFitTRangeIfRestricted()
            if fr is not None:
                fr[0] = max(fr[0], xl[0])
                if fr[0] < fr[1]:
                    xl = fr
            if xl[1] > xl[0]:
                ax.set_xlim(xl[0], xl[1])
        setTitle(ax, 'Reconstruction at %s = %.2f %s' % (
            app.xSym(), wlAct, app.xUnitAx()))
        sciTicks(ax, 'y')
        legend(ax, fontsize=8)
        self.cvKin.draw_idle()

    # ------------------------------------------------------------------
    def doLcurve(self):
        app = self.app
        tauGrid = self.buildGrid()
        self.lblStatus.setText('Computing L-curve (sweeping alpha)...')
        D, t = app.deltaA, app.delay
        t0, fw = self.edT0.value(), self.edFW.value()
        reg, irfM = self.ddReg.currentData(), self.irfMode()
        tR = [self.edFitTmin.value(), self.edFitTmax.value()]
        try:
            aList, rN, sN, aCorner = runBusy(
                app, self, 'Computing L-curve...', lambda: compute_lcurve(
                    D, t, tauGrid, t0, fw, reg, None, irfM, tR))
        except Exception as err:
            self.alert('L-curve failed:\n%s' % err, 'Error')
            self.lblStatus.setText('L-curve failed.')
            return
        ax = self.axL
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(rN, sN, '-o', color=BLUE, markersize=4, linewidth=1.0)
        iC = int(np.argmin(np.abs(aList - aCorner)))
        ax.plot([rN[iC]], [sN[iC]], 's', markerfacecolor='none',
                markeredgecolor='r', markersize=11, markeredgewidth=2.0)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('||residual||')
        ax.set_ylabel('||solution||')
        setTitle(ax, 'L-curve  (corner α = %.3g)' % aCorner)
        self.cvL.draw_idle()
        self.lcurve = {'alpha': aList, 'res': rN, 'sol': sN,
                       'corner': float(aCorner)}
        self.edAlpha.setValue(float(aCorner))
        self.lblStatus.setText(
            'L-curve done. Suggested alpha = %.3g (loaded).' % aCorner)

    def doReset(self):
        app = self.app
        self.edTauMin.setValue(max(0.5 * app.gaFWHM, EPS))
        self.edTauMax.setValue(self._tauMaxDefault())
        self.edK.setValue(100)
        self.edT0.setValue(0.0)
        self.edFW.setValue(max(app.gaFWHM, EPS))
        self.ddReg.setCurrentIndex(1)
        self.edAlpha.setValue(1e-2)
        self.res = emptyResult()
        self.curTau = np.nan
        self._setExport(False)
        self._emptyPlots()
        self.lblStatus.setText('Reset.')
        self.txtResult.setPlainText('\n'.join(IDLE_TEXT))
        self.cbNoIRF.setChecked(False)
        self.onNoIRFChange()
        self.setFitTFull()
        app.ldaSaved = None

    def onNoIRFChange(self):
        self.edFW.setEnabled(not self.cbNoIRF.isChecked())

    def setFitTFull(self):
        self.edFitTmin.setValue(float(np.min(self.app.delay)))
        self.edFitTmax.setValue(float(np.max(self.app.delay)))

    def getFitTRangeIfRestricted(self):
        """[lo, hi] of the fit range, or None if it is the full axis."""
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

    # ---- exports -----------------------------------------------------
    def exportAmap(self):
        app, r = self.app, self.res
        if r['Amap'] is None:
            return
        hdr = [app.xHdr()] + ['tau=%.4g_%s' % (x, app.tUnitHdr())
                              for x in r['tauGrid']]
        cols = [app.wavelength] + [r['Amap'][:, k]
                                   for k in range(r['Amap'].shape[1])]
        self.doExportTable(cols, hdr, 'LDA_Amap', 'Export A(lambda, tau)')

    def exportRec(self):
        app, r = self.app, self.res
        if r['Drec'] is None:
            return
        hdr = [app.xHdr()] + ['t=%.4g_%s' % (x, app.tUnitHdr())
                              for x in np.ravel(app.delay)]
        cols = [app.wavelength] + [r['Drec'][:, j]
                                   for j in range(r['Drec'].shape[1])]
        self.doExportTable(cols, hdr, 'LDA_recon', 'Export reconstruction')

    def export1D(self):
        app, r = self.app, self.res
        if r['Amap'] is None:
            return
        iw = int(np.argmin(np.abs(app.wavelength - self.curLambda)))
        self.doExportTable(
            [r['tauGrid'], r['Amap'][iw, :]],
            ['tau_%s' % app.tUnitHdr(), 'amplitude'],
            'LDA_Atau_%.0f%s' % (app.wavelength[iw], app.xUnitTxt()),
            'Export A(tau)')

    def doExportTable(self, cols, hdr, defaultName, dlgTitle):
        full = self.app.putFile('CSV (*.csv);;DAT (tab) (*.dat)', dlgTitle,
                                defaultName + '.csv', parent=self)
        if not full:
            return
        dlm = ',' if full.lower().endswith('.csv') else '\t'
        try:
            writeTable(full, hdr, cols, dlm)
            self.inform('Saved: %s' % full, 'Exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Error')
