"""Global Analysis window. Port of openGlobalAnalysisWindow."""
import os

import numpy as np
from PySide6 import QtCore, QtWidgets

from trspec.global_analysis import fit_global_analysis
from trspec.kinetics import build_ga_basis_ext, ca_options, irf_mode_norm
from trspec.sequential import (compute_eads_from_dads,
                               compute_sas_sequential, suggest_n_components)

from .state import ui_str, write_data_file, write_matrix
from .textfit import legend, setTitle
from .widgets import (FONT_PT, DataTable, NumEdit, PlotCanvas, VScroll,
                      button, colorButton, flowbox, label, lineColor, panel,
                      runBusy, sciTicks, sidePanel, styleAxes, textArea, tip)
from .win_base import GREEN, Z_LABEL, ChildWindow, MapTriple, absMax
from .win_kinfit import irfLimitNote, undetNote, writeTable

AMP_LABEL = 'Amplitude (ΔA or Intensity)'
RES_LABEL = 'Residual (ΔA or Intensity)'
CA_NAMES = ['Gaussian', '1st derivative', '2nd derivative']


def irfTip():
    return ('Numerical conv.: the model is convolved with the Gaussian IRF '
            '(t0, FWHM). Standard exp components use the closed-form erfcx '
            'expression, stretched components a numerical quadrature at '
            'each delay.\nNone: no IRF convolution -- fit a plain sum of '
            'exponentials exp(-(t-t0)/tau)*H(t-t0). Use after cropping the '
            'data to start beyond the IRF; t0 and FWHM are excluded from '
            'the fit (t0 acts as a fixed onset).')


CA_TIP = ('Signal that exists only while pump and probe overlap (cross-phase '
          'modulation, two-photon absorption, stimulated Raman of the '
          'solvent). It is fitted with up to three extra terms at time zero, '
          'each with a free spectrum: a Gaussian and its first and second '
          'time derivative.\nWithout these terms an artifact distorts the '
          'fitted IRF and can appear as a very fast component. With them, '
          'time constants below the IRF width are less well determined: '
          'switch them on only when the data show an artifact.\nThe artifact '
          'spectra are drawn dashed in the DADS plot and are not part of the '
          'DADS / EADS.')
CAW_TIP = ('Width (FWHM) of the artifact terms.\nSame as IRF: the IRF FWHM, '
           'so the artifact helps to determine the IRF.\nOwn: a separate '
           'width, fitted (kept within a factor of 5 of the IRF FWHM) or '
           'fixed at the value on the right. The IRF FWHM itself is then '
           'hardly determined by the fit: keep it fixed.')
SEQ_TIP = ('How a stretched step is read in the sequential model S1 -> S2 -> '
           '... (SAS and populations; the fit itself is not affected, no new '
           'fit is needed).\nClock from entry: the stretched decay of a '
           'state starts when the population enters that state (a '
           'distribution of lifetimes of the state).\nClock from excitation: '
           'all steps share one clock started by the pump at t0 -- a '
           'time-dependent rate k(t) = (beta/tau)*(t/tau)^(beta-1), as in '
           'dispersive or diffusion-limited kinetics.\nThe two differ only '
           'when a component other than the fastest one is stretched.')


def caLines(app, r, caOrder, caMode, caFwhm, caFwhmStd, caSpec, tauOut,
            fwhmOut, t0Out, fwFixed, t0Fixed, fwStart, t0Start):
    """Lines of a result box about the coherent-artifact terms (shared
    by the global and the target window)."""
    u = app.tUnitTxt()
    L = ['', 'Coherent artifact: %d term(s)' % caOrder]
    if caMode == 'irf':
        L.append('  width = IRF FWHM')
    elif caMode == 'fixed':
        L.append('  width = %.4g %s (fixed)' % (caFwhm, u))
    elif np.isfinite(caFwhmStd):
        L.append('  width = %.4g +/- %.2g %s' % (caFwhm, caFwhmStd, u))
    else:
        L.append('  width = %.4g %s%s' % (caFwhm, u, undetNote(caFwhmStd, 0)))
    for i in range(caOrder):
        col = np.where(np.isfinite(caSpec[:, i]), np.abs(caSpec[:, i]), -1.0)
        iMx = int(np.argmax(col))
        L.append('  %-14s largest %.3g at %.4g %s' % (
            CA_NAMES[i], caSpec[iMx, i], app.wavelength[iMx], app.xUnitTxt()))
    if np.any(np.asarray(tauOut) < fwhmOut):
        L += ['  A tau below the IRF FWHM is',
              '  hard to tell from the artifact.']
    if caMode != 'irf' and not fwFixed and caOrder >= 2:
        L += ['  With a width of their own the',
              '  derivative terms absorb a change of',
              '  the IRF FWHM: it is then hardly',
              '  determined. Fix it, or use', '  "Same as IRF".']
    if (not fwFixed and (fwhmOut > 3 * fwStart or fwhmOut < fwStart / 3)) or \
            (not t0Fixed and abs(t0Out - t0Start) > 3 * max(fwStart, fwhmOut)):
        L += ['  t0 or FWHM moved far from the start',
              '  value. With too few components the',
              '  artifact terms can leave time zero',
              '  and stand in for kinetics: check the',
              '  spectra, or fix the IRF.']
    return L


class GlobalAnalysisWindow(ChildWindow):
    KEY = 'ga'

    def __init__(self, app):
        super().__init__(app, 'Global Analysis', 1420, 1040)
        u = app.tUnitTxt()
        self.curKinWL = float(app.selWL)
        self.cbComps = []
        self.cbCAshow = None
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: fit setup ============================
        left, lv = sidePanel('Fit setup', 356)
        gl.addWidget(left)

        def hrow(spacing=6):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(spacing)
            lv.addLayout(r)
            return r

        row = hrow()
        row.addWidget(label('Components:'))
        self.ddN = QtWidgets.QComboBox()
        self.ddN.addItems(['1', '2', '3', '4', '5'])
        self.ddN.setCurrentText(str(int(app.gaNComp)))
        self.ddN.setMinimumWidth(60)
        row.addWidget(self.ddN)
        self.cbInf = QtWidgets.QCheckBox(ui_str('Include \\tau = \\infty'))
        self.cbInf.setChecked(bool(app.gaHasInf))
        row.addWidget(self.cbInf)
        row.addWidget(button(
            'SVD?', self.onSuggestN,
            tooltip='Estimate the number of components from the singular '
                    'values of the data (count above the noise floor) and '
                    'set Components accordingly. The ratios are printed in '
                    'the result box; the choice remains yours.'))
        row.addStretch(1)

        lv.addWidget(label('Initial time constants (%s, positive):' % u,
                           italic=True, color=(0.3, 0.3, 0.3)))
        self.tauTable = DataTable(
            ui_str(['\\tau_{init} (%s)' % u, '\\tau fixed', '\\beta_{init}',
                    'Stretched']),
            ['num', 'bool', 'num', 'bool'], widths=[85, 70, 70, 70])
        self.tauTable.setFixedHeight(150)
        self.tauTable.cellEdited.connect(
            lambda *_: self.refreshSeqRowState())
        lv.addWidget(self.tauTable)

        row = hrow()
        row.addWidget(label('IRF t_0 (%s):' % u))
        self.edT0 = NumEdit(app.gaT0, width=80)
        row.addWidget(self.edT0)
        self.cbT0fix = QtWidgets.QCheckBox('Fixed')
        self.cbT0fix.setChecked(bool(app.gaT0Fixed))
        row.addWidget(self.cbT0fix)
        row.addStretch(1)

        row = hrow()
        row.addWidget(label('IRF FWHM (%s):' % u))
        self.edFW = NumEdit(app.gaFWHM, width=80, limits=(1e-6, np.inf))
        row.addWidget(self.edFW)
        self.cbFWfix = QtWidgets.QCheckBox('Fixed')
        self.cbFWfix.setChecked(bool(app.gaFWHMFixed))
        row.addWidget(self.cbFWfix)
        row.addStretch(1)

        app.gaIRFTreat = irf_mode_norm(app.gaIRFTreat)
        row = hrow()
        row.addWidget(label('IRF model:'))
        self.ddIRF = QtWidgets.QComboBox()
        self.ddIRF.addItem('None (pure exponentials)', 'none')
        self.ddIRF.addItem('Numerical conv.', 'numerical')
        self.ddIRF.setCurrentIndex(self.ddIRF.findData(app.gaIRFTreat))
        self.ddIRF.setMinimumWidth(165)
        self.ddIRF.setToolTip(tip(irfTip()))
        row.addWidget(self.ddIRF)
        row.addStretch(1)

        app.gaCA = ca_options(app.gaCA)
        row = hrow()
        row.addWidget(label('Coh. artifact:'))
        self.ddCA = QtWidgets.QComboBox()
        for k, txt in enumerate(['None', '1 term: Gaussian',
                                 '2 terms: + 1st derivative',
                                 '3 terms: + 2nd derivative']):
            self.ddCA.addItem(txt, k)
        self.ddCA.setCurrentIndex(int(app.gaCA['order']))
        self.ddCA.setToolTip(tip(CA_TIP))
        row.addWidget(self.ddCA, 1)

        caW0 = app.gaCA['fwhm']
        if not (np.isfinite(caW0) and caW0 > 0):
            caW0 = app.gaFWHM
        row = hrow()
        self.lblCAW = label('Artifact width:')
        row.addWidget(self.lblCAW)
        self.ddCAW = QtWidgets.QComboBox()
        for txt, dat in (('Same as IRF', 'irf'), ('Own, fitted', 'free'),
                         ('Own, fixed', 'fixed')):
            self.ddCAW.addItem(txt, dat)
        self.ddCAW.setCurrentIndex(self.ddCAW.findData(app.gaCA['mode']))
        self.ddCAW.setToolTip(tip(CAW_TIP))
        row.addWidget(self.ddCAW, 1)
        self.edCAW = NumEdit(caW0, width=70, limits=(1e-6, np.inf),
                             tooltip='Own width of the artifact terms (FWHM, '
                                     '%s): starting value or fixed value.' % u)
        row.addWidget(self.edCAW)

        row = hrow(4)
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

        row = hrow(5)
        self.btnRun = colorButton('Run Fit', GREEN)
        self.btnRun.clicked.connect(lambda *_: self.doRun())
        row.addWidget(self.btnRun, 1)
        row.addWidget(button('Reset', self.doReset), 1)

        row = hrow(3)
        row.addWidget(label('Export:', right=True))
        for txt, fn in (('DADS', self.exportDADS), ('EADS', self.exportEADS),
                        ('Fit', lambda: self.exportMatrix2D(
                            'fit', 'Global fit', 'fit')),
                        ('Resid.', lambda: self.exportMatrix2D(
                            'residual', 'Residual', 'residual')),
                        (ui_str('Kin \\lambda'), self.exportKinetics),
                        ('Popul.', self.exportPopulation)):
            b = button(txt, fn)
            b.setMinimumWidth(30)
            b.setStyleSheet('padding: 3px 2px;')
            row.addWidget(b, 1)
            if txt == 'EADS':
                self.btnExpEADS = b
            if txt == 'Popul.':
                b.setToolTip(tip('Export the sequential-model population '
                                 'kinetics c_i(t) shown in the bottom-right '
                                 'panel.'))

        if app.gaSeqModel not in ('age', 't0'):
            app.gaSeqModel = 'age'
        row = hrow()
        self.lblSeq = label('Stretched step:')
        self.lblSeq.setToolTip(tip(SEQ_TIP))
        row.addWidget(self.lblSeq)
        self.ddSeq = QtWidgets.QComboBox()
        self.ddSeq.addItem('Clock from entry into state', 'age')
        self.ddSeq.addItem('Clock from excitation (t0)', 't0')
        self.ddSeq.setCurrentIndex(self.ddSeq.findData(app.gaSeqModel))
        self.ddSeq.setToolTip(tip(SEQ_TIP))
        row.addWidget(self.ddSeq, 1)

        row = hrow()
        row.addWidget(label('Map scale:', right=True))
        self.ddMapScale = QtWidgets.QComboBox()
        for txt, dat in (('Linear', 'linear'), ('Log', 'log'),
                         ('Split', 'split')):
            self.ddMapScale.addItem(txt, dat)
        self.ddMapScale.setCurrentIndex(1)
        self.ddMapScale.setMinimumWidth(90)
        row.addWidget(self.ddMapScale)
        row.addWidget(label('Thr. (%s):' % u, right=True))
        self.edMapThr = NumEdit(app.splitThreshold, width=55,
                                limits=(np.finfo(float).eps, np.inf))
        row.addWidget(self.edMapThr)
        row.addStretch(1)

        row = hrow()
        lb = label('t min (%s):' % u, right=True)
        lb.setToolTip(tip('Lower bound of the delay axis on the 2D maps. '
                          'Applies to linear mode (whole axis) and split '
                          'mode (linear lower panel only).'))
        row.addWidget(lb)
        self.edMapTmin = NumEdit(-1, width=70)
        row.addWidget(self.edMapTmin)
        lb = label('t max (%s):' % u, right=True)
        lb.setToolTip(tip('Upper bound of the delay axis on the 2D maps. '
                          'Default: max delay in the data.'))
        row.addWidget(lb)
        self.edMapTmax = NumEdit(float(np.max(app.delay)), width=70)
        row.addWidget(self.edMapTmax)
        row.addStretch(1)

        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        lv.addWidget(self.lblStatus)
        self.txtResult = textArea(['Results will appear here after fitting.'])
        self.txtResult.setMinimumHeight(110)
        lv.addWidget(self.txtResult, 1)

        # ================= right: results =============================
        # (the plots keep a least height; a window too low for all of
        # them lets them scroll)
        right = panel('Results')
        gl.addWidget(right, 1)
        ro = QtWidgets.QVBoxLayout(right)
        ro.setContentsMargins(2, 12, 2, 4)
        rin = QtWidgets.QWidget()
        rv = QtWidgets.QVBoxLayout(rin)
        rv.setContentsMargins(4, 0, 4, 2)
        rv.setSpacing(6)
        ro.addWidget(VScroll(rin))
        self.maps = MapTriple(app, self.applyGaPoint)
        self.maps.cv.setMinimumHeight(180)
        rv.addWidget(self.maps.cv, 100)
        self.cvDADS = PlotCanvas()
        self.cvDADS.setMinimumHeight(160)
        rv.addWidget(self.cvDADS, 85)
        fig = self.cvDADS.figure
        self.axDADS = fig.add_axes([0.07, 0.15, 0.4, 0.75])
        self.axDADSnorm = fig.add_axes([0.57, 0.15, 0.4, 0.75])
        self.cvDADS.layoutFcn = lambda: self._lay2(
            self.cvDADS, self.axDADS, self.axDADSnorm)
        self.cvEADS = PlotCanvas()
        self.cvEADS.setMinimumHeight(160)
        rv.addWidget(self.cvEADS, 85)
        fig = self.cvEADS.figure
        self.axEADS = fig.add_axes([0.07, 0.15, 0.4, 0.75])
        self.axEADSnorm = fig.add_axes([0.57, 0.15, 0.4, 0.75])
        self.cvEADS.layoutFcn = lambda: self._lay2(
            self.cvEADS, self.axEADS, self.axEADSnorm)
        # the two rows of controls wrap when the window is narrow
        self.cbHost, self.cbLay = flowbox(10)
        self.cbLay.setContentsMargins(0, 2, 0, 2)
        self.cbHost.setMinimumHeight(24)
        rv.addWidget(self.cbHost)

        kcHost, kc = flowbox(8)
        self.edKinWL = NumEdit(app.selWL, fmt='%.2f', width=90)
        kc.addPair(label('%s (%s):' % (app.xSymTxt(), app.xUnitTxt()),
                         right=True), self.edKinWL)
        kc.addWidget(button('Use main selection', self.syncKinWLFromMain))
        self.ddKinScale = QtWidgets.QComboBox()
        self.ddKinScale.addItems(['Log', 'Linear'])
        self.ddKinScale.setMinimumWidth(80)
        kc.addPair(label('Scale:', right=True), self.ddKinScale)
        kc.addWidget(label('(or click on the experimental map above)',
                           color=(0.35, 0.35, 0.35), italic=True))
        rv.addWidget(kcHost)
        self.cvKin = PlotCanvas()
        self.cvKin.setMinimumHeight(170)
        rv.addWidget(self.cvKin, 95)
        fig = self.cvKin.figure
        self.axKinFit = fig.add_axes([0.06, 0.15, 0.35, 0.75])
        self.axKinRes = fig.add_axes([0.48, 0.15, 0.22, 0.75])
        self.axPop = fig.add_axes([0.76, 0.15, 0.22, 0.75])
        self.cvKin.layoutFcn = self._layKin

        self._initAxes()
        self.fillTable()
        self.refreshSeqRowState()
        self.updateEadsPanelTitle()
        self.onIRFModeChange()
        lo, hi = self.clampFitRange(app.gaFitTMin, app.gaFitTMax)
        self.edFitTmin.setValue(lo)
        self.edFitTmax.setValue(hi)

        self.ddN.currentTextChanged.connect(lambda *_: self.onNCompChange())
        self.ddIRF.currentIndexChanged.connect(
            lambda *_: self.onIRFModeChange())
        self.ddCA.currentIndexChanged.connect(lambda *_: self.syncCAEnable())
        self.ddCAW.currentIndexChanged.connect(lambda *_: self.syncCAEnable())
        self.ddSeq.currentIndexChanged.connect(
            lambda *_: self.onSeqModelChange())
        self.ddMapScale.currentIndexChanged.connect(
            lambda *_: self.onMapScaleChange())
        for ed in (self.edMapThr, self.edMapTmin, self.edMapTmax):
            ed.valueChanged.connect(lambda *_: self.onMapScaleChange())
        self.edKinWL.valueChanged.connect(self.onKinWLChange)
        self.ddKinScale.currentTextChanged.connect(
            lambda *_: self.onKinScaleChange())

        gaOK = app.gaResultDADS is not None and \
            np.shape(app.gaResultDADS)[0] == app.wavelength.size and \
            (app.gaResultFit is None or
             np.shape(app.gaResultFit) == (app.wavelength.size,
                                           app.delay.size))
        if app.gaResultDADS is not None and not gaOK:
            app.clearGaResults()
        if gaOK:
            self.lblStatus.setText('Restored previous fit. RMS = %.3g'
                                   % app.gaResultRMS)
            self.txtResult.setPlainText(
                'Previous fit restored (RMS = %.3g).  Run again to refit.'
                % app.gaResultRMS)
            self.plotResults()
            self.updateEadsPanelTitle()
        self.refreshSeqRowState()

    # ------------------------------------------------------------------
    # layout helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _lay2(cv, a1, a2):
        W, H = cv.sizePt()
        L, R, T, B, gap = 54.0, 16.0, 20.0, 32.0, 56.0
        w = (W - L - R - gap) / 2.0
        cv.place(a1, L, B, w, H - T - B)
        cv.place(a2, L + w + gap, B, w, H - T - B)

    def _layKin(self):
        cv = self.cvKin
        W, H = cv.sizePt()
        L, R, T, B, gap = 54.0, 16.0, 20.0, 32.0, 56.0
        w = (W - L - R - 2 * gap) / 7.0
        h = H - T - B
        cv.place(self.axKinFit, L, B, 3 * w, h)
        cv.place(self.axKinRes, L + 3 * w + gap, B, 2 * w, h)
        cv.place(self.axPop, L + 5 * w + 2 * gap, B, 2 * w, h)

    def _initAxes(self):
        app = self.app
        tl = 'Delay time (%s)' % app.tUnitAx()
        for ax, xl, yl, ti in (
                (self.axDADS, app.xAxLabel(), AMP_LABEL,
                 'Decay-associated difference spectra'),
                (self.axDADSnorm, app.xAxLabel(), 'Normalized amplitude',
                 'DADS (peak-normalized)'),
                (self.axEADS, app.xAxLabel(), AMP_LABEL, ''),
                (self.axEADSnorm, app.xAxLabel(), 'Normalized amplitude', ''),
                (self.axKinFit, tl, Z_LABEL,
                 'Kinetics fit at selected %s' % app.xSym()),
                (self.axKinRes, tl, RES_LABEL, 'Fit residual'),
                (self.axPop, tl, 'Population',
                 'EADS populations $c_i(t)$')):
            styleAxes(ax, grid=True)
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            setTitle(ax, ti, {'Decay-associated difference spectra': 'DADS',
                              'EADS populations $c_i(t)$': 'Populations'
                              }.get(ti, ''))
        self.maps.build(1)
        for g, ti in zip(self.maps.groups, (
                (ui_str('Experimental  (click to pick \\lambda)'),
                 'Experimental'), ('Global fit',), ('Residual',))):
            ax = g['axes'][0]
            styleAxes(ax)
            setTitle(ax, *ti, size=FONT_PT, bold=True)
            ax.set_xlabel(app.xAxLabel())
            ax.set_ylabel('Delay (%s)' % app.tUnitAx())

    # ------------------------------------------------------------------
    # small helpers (shared with the single-trace window in MATLAB)
    # ------------------------------------------------------------------
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

    def getFitTRangeIfRestricted(self):
        app = self.app
        a, b = app.gaFitTMin, app.gaFitTMax
        if a is None or b is None or not np.isfinite(a) or not np.isfinite(b):
            return None
        lo, hi = min(a, b), max(a, b)
        dLo, dHi = float(np.min(app.delay)), float(np.max(app.delay))
        span = dHi - dLo
        if span <= 0:
            return None
        if (lo - dLo) > 1e-9 * span or (dHi - hi) > 1e-9 * span:
            return [lo, hi]
        return None

    def resStretch(self):
        app = self.app
        if app.gaResultTau is None:
            return np.asarray(app.gaStretchOn, bool)
        n = len(app.gaResultTau)
        so = app.gaResultStretchOn
        if so is None or len(so) != n:
            so = app.gaStretchOn
        if so is None or len(so) != n:
            so = np.zeros(n, bool)
        return np.asarray(so, bool)

    def resIRF(self):
        m = self.app.gaResultIRF or self.app.gaIRFTreat
        return irf_mode_norm(m)

    def seqModelTxt(self):
        return 'clock from t0' if self.app.gaSeqModel == 't0' \
            else 'clock from entry'

    # ------------------------------------------------------------------
    # table, controls
    # ------------------------------------------------------------------
    def fillTable(self):
        app = self.app
        N = int(app.gaNComp)
        tau = list(app.gaTauInit)
        fix = list(app.gaTauFixed)
        if len(tau) < N:
            last = tau[-1] if tau else 1.0
            extra = N - len(tau)
            tau += [last * 10.0 ** k for k in range(1, extra + 1)]
            fix += [False] * extra
        tau, fix = tau[:N], (fix + [False] * N)[:N]
        bet = (list(app.gaBetaInit) + [1.0] * N)[:N]
        son = (list(app.gaStretchOn) + [False] * N)[:N]
        app.gaTauInit, app.gaTauFixed = tau, fix
        app.gaBetaInit, app.gaStretchOn = bet, son
        self.tauTable.setData([[float(tau[i]), bool(fix[i]), float(bet[i]),
                                bool(son[i])] for i in range(N)])

    def readTable(self):
        app = self.app
        d = self.tauTable.data()
        if not d:
            return
        app.gaTauInit = [float(r[0]) for r in d]
        app.gaTauFixed = [bool(r[1]) for r in d]
        app.gaBetaInit = [float(r[2]) for r in d]
        app.gaStretchOn = [bool(r[3]) for r in d]

    def onSuggestN(self):
        app = self.app
        nSug, _, txt = suggest_n_components(app.deltaA)
        if not txt:
            return
        nSug = max(1, min(nSug, 5))
        self.ddN.setCurrentText(str(nSug))
        self.onNCompChange()
        self.lblStatus.setText('SVD suggests N = %d' % nSug)
        self.txtResult.setPlainText('\n'.join([
            'SVD rank estimate: N = %d' % nSug, txt, '',
            'Components set to this value; adjust and press Run Fit.']))

    def onNCompChange(self):
        self.readTable()
        self.app.gaNComp = int(self.ddN.currentText())
        self.fillTable()
        self.refreshSeqRowState()

    def refreshSeqRowState(self):
        app = self.app
        anyStretch = any(bool(r[3]) for r in self.tauTable.data())
        if not anyStretch and app.gaResultTau is not None:
            anyStretch = bool(np.any(self.resStretch()))
        self.ddSeq.setEnabled(anyStretch)
        self.lblSeq.setEnabled(anyStretch)

    def onMapScaleChange(self):
        b = QtCore.QSignalBlocker(self.ddKinScale)
        self.ddKinScale.setCurrentText(
            'Linear' if self.ddMapScale.currentData() == 'linear' else 'Log')
        del b
        if self.app.gaResultFit is None:
            return
        self.plotMaps()
        self.onKinScaleChange()

    def onKinScaleChange(self):
        self.plotKinetics()
        self.plotPopulation()

    def setFitTFull(self):
        self.edFitTmin.setValue(float(np.min(self.app.delay)))
        self.edFitTmax.setValue(float(np.max(self.app.delay)))

    def onIRFModeChange(self):
        isNone = self.ddIRF.currentData() == 'none'
        self.edFW.setEnabled(not isNone)
        self.cbFWfix.setEnabled(not isNone)
        if isNone:
            self.cbT0fix.setChecked(True)
        self.cbT0fix.setEnabled(not isNone)
        self.syncCAEnable()

    def syncCAEnable(self):
        hasIRF = self.ddIRF.currentData() != 'none'
        useCA = hasIRF and self.ddCA.currentData() > 0
        self.ddCA.setEnabled(hasIRF)
        self.ddCAW.setEnabled(useCA)
        self.lblCAW.setEnabled(useCA)
        self.edCAW.setEnabled(useCA and self.ddCAW.currentData() != 'irf')

    def doReset(self):
        app = self.app
        app.gaNComp = 3
        app.gaTauInit = [1.0, 10.0, 100.0]
        app.gaTauFixed = [False, False, False]
        app.gaBetaInit = [1.0, 1.0, 1.0]
        app.gaStretchOn = [False, False, False]
        app.gaIRFTreat = 'numerical'
        app.gaHasInf = False
        app.gaT0, app.gaT0Fixed = 0.0, True
        app.gaFWHM, app.gaFWHMFixed = 0.15, True
        b = QtCore.QSignalBlocker(self.ddN)
        self.ddN.setCurrentText('3')
        del b
        self.cbInf.setChecked(False)
        self.edT0.setValue(0.0)
        self.cbT0fix.setChecked(True)
        self.edFW.setValue(0.15)
        self.cbFWfix.setChecked(True)
        self.ddIRF.setCurrentIndex(self.ddIRF.findData('numerical'))
        app.gaCA = ca_options(None)
        self.ddCA.setCurrentIndex(0)
        self.ddCAW.setCurrentIndex(0)
        self.edCAW.setValue(app.gaFWHM)
        if self.ddSeq.currentData() != 'age':
            self.ddSeq.setCurrentIndex(0)       # calls onSeqModelChange
        app.gaFitTMin = None
        app.gaFitTMax = None
        self.setFitTFull()
        self.fillTable()
        self.refreshSeqRowState()
        self.onIRFModeChange()
        self.updateEadsPanelTitle()
        self.axPop.clear()
        styleAxes(self.axPop, grid=True)
        self.axPop.set_title('EADS populations $c_i(t)$', fontsize=FONT_PT + 1)
        self.cvKin.draw_idle()
        self.lblStatus.setText('Reset to defaults.')

    # ------------------------------------------------------------------
    # fit
    # ------------------------------------------------------------------
    def doRun(self):
        app = self.app
        self.readTable()
        app.gaHasInf = self.cbInf.isChecked()
        app.gaT0 = self.edT0.value()
        app.gaT0Fixed = self.cbT0fix.isChecked()
        app.gaFWHM = self.edFW.value()
        app.gaFWHMFixed = self.cbFWfix.isChecked()
        app.gaIRFTreat = self.ddIRF.currentData()
        app.gaSeqModel = self.ddSeq.currentData()
        app.gaCA = {'order': int(self.ddCA.currentData()),
                    'mode': self.ddCAW.currentData(),
                    'fwhm': float(self.edCAW.value())}
        if app.gaIRFTreat == 'none':
            app.gaT0Fixed = True
            app.gaFWHMFixed = True
        tA, tB = self.edFitTmin.value(), self.edFitTmax.value()
        app.gaFitTMin, app.gaFitTMax = self.fitRangeToStore(tA, tB)
        tSelChk = (app.delay >= min(tA, tB)) & (app.delay <= max(tA, tB))
        if np.count_nonzero(tSelChk) < app.gaNComp + 2:
            self.alert('Fit t range contains too few delay points. Widen '
                       'the range.', 'Invalid range')
            return
        if any(v <= 0 for v in app.gaTauInit):
            self.alert('All tau values must be positive.', 'Invalid input')
            return
        if any(v <= 0 or v > 2 for v in app.gaBetaInit):
            self.alert('beta_init must be in (0, 2]. Use 1 for ordinary '
                       'exponential.', 'Invalid input')
            return
        if app.gaFWHM <= 0:
            self.alert('IRF FWHM must be positive.', 'Invalid input')
            return
        self.lblStatus.setText('Fitting...')
        tRange = None
        if app.gaFitTMin is not None:
            tRange = [app.gaFitTMin, app.gaFitTMax]
        D, t = app.deltaA, app.delay
        args = (D, t, list(app.gaTauInit), app.gaT0, app.gaFWHM,
                list(app.gaTauFixed), app.gaT0Fixed, app.gaFWHMFixed,
                app.gaHasInf, list(app.gaBetaInit), list(app.gaStretchOn),
                app.gaIRFTreat, tRange, dict(app.gaCA))
        try:
            r = runBusy(app, self, 'Fitting...',
                        lambda: fit_global_analysis(*args))
        except Exception as err:
            self.alert('Fit failed:\n%s' % err, 'Fit error')
            self.lblStatus.setText('Fit failed.')
            app.clearGaResults()
            self.refreshSeqRowState()
            return
        tauOut, betaOut = np.asarray(r.tau, float), np.asarray(r.beta, float)
        app.gaResultTau = tauOut
        app.gaResultBeta = betaOut
        app.gaResultT0 = float(r.t0)
        app.gaResultFWHM = float(r.fwhm)
        app.gaResultHasInf = bool(app.gaHasInf)
        app.gaResultDADS = np.asarray(r.dads, float)
        app.gaResultFit = np.asarray(r.fit, float)
        app.gaResultRMS = float(r.rms)
        app.gaResultTauStd = np.asarray(r.tau_std, float)
        app.gaResultStretchOn = list(app.gaStretchOn)
        app.gaResultIRF = app.gaIRFTreat
        app.gaResultTRange = tRange
        if r.ca_order > 0:
            app.gaResultCA = {'order': int(r.ca_order),
                              'spec': np.asarray(r.ca_spec, float),
                              'fwhm': float(r.ca_fwhm),
                              'fwhmStd': float(r.ca_fwhm_std),
                              'mode': r.ca_mode}
        else:
            app.gaResultCA = None
        useStretched = any(app.gaStretchOn)
        sasErr = ''
        try:
            if useStretched:
                app.gaResultEADS = compute_sas_sequential(
                    app.gaDataNoCA(), t, tauOut, betaOut, app.gaStretchOn,
                    app.gaHasInf, r.t0, r.fwhm, app.gaIRFTreat, tRange,
                    app.gaSeqModel)[0]
            else:
                app.gaResultEADS = compute_eads_from_dads(
                    app.gaResultDADS, tauOut, app.gaHasInf)[0]
        except Exception as err:
            app.gaResultEADS = None
            sasErr = str(err)
        if app.gaResultEADS is not None and \
                not np.any(np.isfinite(app.gaResultEADS)):
            sasErr = 'all-NaN result'
            app.gaResultEADS = None
        lines = self.resultLines(r, useStretched, sasErr)
        self.txtResult.setPlainText('\n'.join(lines))
        self.lblStatus.setText('Done. RMS = %.3g' % r.rms)
        self.plotResults()
        self.updateEadsPanelTitle()
        self.refreshSeqRowState()

    def resultLines(self, r, useStretched, sasErr):
        app = self.app
        u = app.tUnitTxt()
        tauOut, betaOut = r.tau, r.beta
        L = ['Fit converged (%d iters)' % r.iters,
             'Initial RMS : %.4g' % r.initial_rms,
             'Final   RMS : %.4g' % r.rms]
        if r.rms >= 0.99 * r.initial_rms and r.iters > 2:
            L += ['', '** NOT CONVERGED - try different',
                  '   initial tau values, or fix IRF/t0.']
        L += ['', 'Time constants:']
        for i in range(len(tauOut)):
            tag = ''
            if app.gaStretchOn[i]:
                tag = ' [stretched, beta=%.3f]' % betaOut[i]
                if np.isfinite(r.beta_std[i]):
                    tag = ' [stretched, beta=%.3f+/-%.3f]' % (
                        betaOut[i], r.beta_std[i])
            if app.gaTauFixed[i]:
                tag += ' (tau fixed)'
            tag += undetNote(r.tau_std, i)
            if np.isfinite(r.tau_std[i]):
                L.append('  tau_%d = %10.4g +/- %.2g %s%s' % (
                    i + 1, tauOut[i], r.tau_std[i], u, tag))
            else:
                L.append('  tau_%d = %10.4g %s%s' % (i + 1, tauOut[i], u, tag))
        if app.gaHasInf:
            L.append('  tau_inf  (offset component)')
        L.append('')
        if np.isfinite(r.t0_std):
            L.append('IRF t0   = %.4g +/- %.2g %s' % (r.t0, r.t0_std, u))
        else:
            L.append('IRF t0   = %.4g %s%s' % (
                r.t0, u, (' (fixed)' if app.gaT0Fixed else '')
                + undetNote(r.t0_std, 0)))
        if r.fwhm_at_limit:
            L.append('IRF FWHM = %.4g %s  [not resolved]' % (r.fwhm, u))
        elif np.isfinite(r.fwhm_std):
            L.append('IRF FWHM = %.4g +/- %.2g %s' % (r.fwhm, r.fwhm_std, u))
        else:
            L.append('IRF FWHM = %.4g %s%s' % (
                r.fwhm, u, (' (fixed)' if app.gaFWHMFixed else '')
                + undetNote(r.fwhm_std, 0)))
        L += irfLimitNote(app, r)
        L.append('IRF model: ' + ('none (pure exponentials)'
                                  if app.gaIRFTreat == 'none'
                                  else 'numerical conv.'))
        if r.ca_order > 0:
            L += caLines(app, r, r.ca_order, r.ca_mode, r.ca_fwhm,
                         r.ca_fwhm_std, np.asarray(r.ca_spec, float), tauOut,
                         r.fwhm, r.t0, app.gaFWHMFixed, app.gaT0Fixed,
                         app.gaFWHM, app.gaT0)
        if useStretched:
            L += ['', 'Right-panel "EADS" -> "SAS"',
                  '(sequential w/ stretched; clock',
                  'set by "Stretched step").']
            if sasErr:
                L += ['', '*** SAS computation FAILED ***', sasErr]
        return L

    def onSeqModelChange(self):
        app = self.app
        new = self.ddSeq.currentData()
        if app.gaSeqModel == new:
            return
        app.gaSeqModel = new
        if app.gaResultTau is None or not np.any(self.resStretch()):
            return
        self.lblStatus.setText('Recomputing SAS...')
        QtWidgets.QApplication.processEvents()
        try:
            E = compute_sas_sequential(
                app.gaDataNoCA(), app.delay, app.gaResultTau,
                app.gaResultBeta, self.resStretch(), app.gaResultHasInf,
                app.gaResultT0, app.gaResultFWHM, self.resIRF(),
                app.gaResultTRange, app.gaSeqModel)[0]
            if not np.any(np.isfinite(E)):
                raise ValueError('all-NaN result')
            app.gaResultEADS = E
            self.lblStatus.setText('SAS recomputed (%s).' % self.seqModelTxt())
        except Exception as err:
            app.gaResultEADS = None
            self.lblStatus.setText('SAS computation failed: %s' % err)
        self.plotEADS()
        self.plotPopulation()
        self.updateEadsPanelTitle()

    # ------------------------------------------------------------------
    # populations
    # ------------------------------------------------------------------
    def computePopulation(self):
        app = self.app
        if app.gaResultTau is None:
            return None, None, 0
        popT = app.delay.copy()
        popC = None
        try:
            if np.any(self.resStretch()):
                popC = compute_sas_sequential(
                    app.gaDataNoCA(), app.delay, app.gaResultTau,
                    app.gaResultBeta, self.resStretch(), app.gaResultHasInf,
                    app.gaResultT0, app.gaResultFWHM, self.resIRF(),
                    app.gaResultTRange, app.gaSeqModel)[1]
            else:
                _, tauSorted, Bm = compute_eads_from_dads(
                    app.gaResultDADS, app.gaResultTau, app.gaResultHasInf)
                tauCols = np.concatenate([tauSorted, [np.inf]]) \
                    if app.gaResultHasInf else np.asarray(tauSorted, float)
                n = tauCols.size
                basis = build_ga_basis_ext(
                    popT, tauCols, np.ones(n), np.zeros(n, bool),
                    app.gaResultT0, app.gaResultFWHM, False, self.resIRF())
                popC = (Bm @ basis.T).T
        except Exception:
            popC = None
        if popC is None or np.size(popC) == 0:
            return None, popT, 0
        popC = np.asarray(popC, float)
        return popC, popT, popC.shape[1]

    def plotPopulation(self):
        app = self.app
        ax = self.axPop
        ax.clear()
        styleAxes(ax, grid=True)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel('Population')
        popC, popT, nSp = self.computePopulation()
        if popC is None:
            if app.gaResultTau is None:
                setTitle(ax, 'EADS populations $c_i(t)$', 'Populations')
            else:
                setTitle(ax, 'Populations $c_i(t)$: not available',
                         'Populations: not available')
            self.cvKin.draw_idle()
            return
        for j in range(nSp):
            ax.plot(popT, popC[:, j], '-', linewidth=1.4, color=lineColor(j),
                    label='S%d' % (j + 1))
        xl = [float(np.min(popT)), float(np.max(popT))]
        isLog = False
        if self.ddKinScale.currentText() == 'Log':
            tPos = popT[popT > 0]
            if tPos.size:
                isLog = True
                xl = [float(np.min(tPos)), float(np.max(popT))]
        ax.set_xscale('log' if isLog else 'linear')
        fr = self.getFitTRangeIfRestricted()
        if fr is not None:
            if isLog:
                fr[0] = max(fr[0], xl[0])
            if fr[0] < fr[1]:
                xl = fr
        ax.set_xlim(xl[0], xl[1])
        legend(ax, frameon=False, fontsize=8)
        setTitle(ax, ('SAS' if np.any(self.resStretch()) else 'EADS')
                 + ' populations $c_i(t)$', 'Populations $c_i(t)$',
                 'Populations')
        self.cvKin.draw_idle()

    # ------------------------------------------------------------------
    # plots
    # ------------------------------------------------------------------
    def updateEadsPanelTitle(self):
        if np.any(self.resStretch()):
            setTitle(self.axEADS, 'SAS (sequential; stretched steps: %s)'
                     % self.seqModelTxt(), 'SAS (sequential)', 'SAS')
            setTitle(self.axEADSnorm, 'SAS (peak-normalized)',
                     'SAS (normalized)')
            self.btnExpEADS.setText('SAS')
        else:
            setTitle(self.axEADS, 'Evolution-associated difference spectra '
                     '(sequential model)', 'EADS (sequential model)', 'EADS')
            setTitle(self.axEADSnorm, 'EADS (peak-normalized)',
                     'EADS (normalized)')
            self.btnExpEADS.setText('EADS')
        self.cvEADS.draw_idle()

    def plotResults(self):
        if self.app.gaResultDADS is None:
            return
        self.rebuildCompCheckboxes()
        self.plotDADS()
        self.plotEADS()
        self.plotMaps()
        self.plotKinetics()
        self.plotPopulation()

    def onCompVis(self):
        self.plotDADS()
        self.plotEADS()

    def syncKinWLFromMain(self):
        self.curKinWL = float(self.app.selWL)
        self.edKinWL.setValue(self.curKinWL)
        self.maps.setMarker(self.curKinWL)
        self.plotKinetics()

    def onKinWLChange(self, val):
        self.curKinWL = float(val)
        self.maps.setMarker(self.curKinWL)
        self.plotKinetics()

    def rebuildCompCheckboxes(self):
        app = self.app
        while self.cbLay.count():
            it = self.cbLay.takeAt(0)
            wd = it.widget()
            if wd is not None:
                wd.hide()
                wd.setParent(None)
                wd.deleteLater()
        self.cbComps = []
        self.cbCAshow = None
        if app.gaResultDADS is None:
            return
        nMain = len(app.gaResultTau)
        self.cbLay.addWidget(label('Show components:', italic=True))
        for i in range(nMain):
            cb = QtWidgets.QCheckBox(ui_str('\\tau_%d = %.3g %s' % (
                i + 1, app.gaResultTau[i], app.tUnitAx())))
            cb.setChecked(True)
            cb.toggled.connect(lambda *_: self.onCompVis())
            self.cbLay.addWidget(cb)
            self.cbComps.append(cb)
        if app.gaResultHasInf:
            cb = QtWidgets.QCheckBox(ui_str('\\tau = \\infty'))
            cb.setChecked(True)
            cb.toggled.connect(lambda *_: self.onCompVis())
            self.cbLay.addWidget(cb)
            self.cbComps.append(cb)
        ca = app.gaResultCA
        if isinstance(ca, dict) and np.size(ca.get('spec', [])):
            cb = QtWidgets.QCheckBox('artifact')
            cb.setChecked(True)
            cb.setToolTip(tip('Show the spectra of the coherent-artifact '
                              'terms (dashed) in the DADS plot.'))
            cb.toggled.connect(lambda *_: self.onCompVis())
            self.cbLay.addWidget(cb)
            self.cbCAshow = cb

    def getCompVisible(self):
        nc = np.shape(self.app.gaResultDADS)[1]
        if len(self.cbComps) != nc:
            return [True] * nc
        return [cb.isChecked() for cb in self.cbComps]

    def _tauLegend(self, i, nMain, withValue=True):
        app = self.app
        if i < nMain:
            if withValue:
                return '$\\tau_%d$ = %.3g %s' % (i + 1, app.gaResultTau[i],
                                               app.tUnitAx())
            return '$\\tau_%d$' % (i + 1)
        return '$\\tau$ = $\\infty$'

    def plotDADS(self):
        app = self.app
        if app.gaResultDADS is None:
            return
        A = app.gaResultDADS
        wl = app.wavelength
        nc = A.shape[1]
        nMain = len(app.gaResultTau)
        vis = self.getCompVisible()
        ax = self.axDADS
        ax.clear()
        styleAxes(ax, grid=True)
        n = 0
        for i in range(nc):
            if not vis[i]:
                continue
            ax.plot(wl, A[:, i], '-', color=lineColor(i), linewidth=1.4,
                    label=self._tauLegend(i, nMain))
            n += 1
        ca = app.gaResultCA
        showCA = self.cbCAshow is None or self.cbCAshow.isChecked()
        if showCA and isinstance(ca, dict) and \
                np.shape(ca['spec'])[0] == wl.size:
            caNm = ['artifact: Gaussian', 'artifact: 1st deriv.',
                    'artifact: 2nd deriv.']
            caCol = [(0.30,) * 3, (0.50,) * 3, (0.68,) * 3]
            for i in range(min(3, np.shape(ca['spec'])[1])):
                ax.plot(wl, ca['spec'][:, i], '--', color=caCol[i],
                        linewidth=1.1, label=caNm[i])
                n += 1
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel(AMP_LABEL)
        setTitle(ax, 'Decay-associated difference spectra   (RMS = %.3g)'
                 % app.gaResultRMS, 'DADS   (RMS = %.3g)' % app.gaResultRMS,
                 'DADS')
        if n:
            legend(ax, fontsize=8)
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        sciTicks(ax, 'y')

        ax = self.axDADSnorm
        ax.clear()
        styleAxes(ax, grid=True)
        n = 0
        for i in range(nc):
            if not vis[i]:
                continue
            ax.plot(wl, A[:, i] / absMax(A[:, i]), '-', color=lineColor(i),
                    linewidth=1.4, label=self._tauLegend(i, nMain, False))
            n += 1
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('Normalized amplitude')
        setTitle(ax, 'DADS (peak-normalized)', 'DADS (normalized)')
        if n:
            legend(ax, fontsize=8)
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        ax.set_ylim(-1.1, 1.1)
        self.cvDADS.draw_idle()

    def plotEADS(self):
        app = self.app
        wl = app.wavelength
        for ax, yl in ((self.axEADS, AMP_LABEL),
                       (self.axEADSnorm, 'Normalized amplitude')):
            ax.clear()
            styleAxes(ax, grid=True)
            ax.set_xlabel(app.xAxLabel())
            ax.set_ylabel(yl)
        if app.gaResultEADS is None:
            if np.any(self.resStretch()):
                setTitle(self.axEADS,
                         'Species-associated spectra (sequential, with '
                         'stretched) - run fit to compute',
                         'SAS - run fit to compute', 'SAS')
                setTitle(self.axEADSnorm, 'SAS (peak-normalized)',
                         'SAS (normalized)')
            else:
                setTitle(self.axEADS,
                         'Evolution-associated difference spectra  '
                         '(sequential) - run fit to compute',
                         'EADS - run fit to compute', 'EADS')
                setTitle(self.axEADSnorm, 'EADS (peak-normalized)',
                         'EADS (normalized)')
            self.cvEADS.draw_idle()
            return
        E = np.asarray(app.gaResultEADS, float)
        nc = E.shape[1]
        tau = np.asarray(app.gaResultTau, float)
        nMain = tau.size
        sortIdx = np.argsort(tau, kind='stable')
        tauSorted = tau[sortIdx]
        permFit = list(sortIdx) + ([nMain] if app.gaResultHasInf else [])
        vis = self.getCompVisible()

        def shown(iS):
            fi = permFit[iS] if iS < len(permFit) else iS
            return fi >= len(vis) or vis[fi], fi

        ax = self.axEADS
        n = 0
        for iS in range(nc):
            ok, fi = shown(iS)
            if not ok:
                continue
            if iS < nMain:
                lg = '$S_%d$  $\\tau$ = %.3g %s' % (iS + 1, tauSorted[iS],
                                                  app.tUnitAx())
            else:
                lg = '$S_%d$  $\\tau$ = $\\infty$' % (iS + 1)
            ax.plot(wl, E[:, iS], '-', color=lineColor(fi), linewidth=1.4,
                    label=lg)
            n += 1
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        if n:
            legend(ax, fontsize=8)
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        visIdx = [shown(iS)[0] for iS in range(nc)]
        if any(visIdx):
            Es = E[:, visIdx]
            fin = Es[np.isfinite(Es)]
            if fin.size and np.min(fin) < np.max(fin):
                pad = 0.08 * (np.max(fin) - np.min(fin))
                ax.set_ylim(np.min(fin) - pad, np.max(fin) + pad)
        sciTicks(ax, 'y')
        ax = self.axEADSnorm
        n = 0
        for iS in range(nc):
            ok, fi = shown(iS)
            if not ok:
                continue
            ax.plot(wl, E[:, iS] / absMax(E[:, iS]), '-', color=lineColor(fi),
                    linewidth=1.4, label='$S_%d$' % (iS + 1))
            n += 1
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        if n:
            legend(ax, fontsize=8)
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        ax.set_ylim(-1.1, 1.1)
        self.updateEadsPanelTitle()

    def plotMaps(self):
        app = self.app
        if app.gaResultFit is None:
            return
        cl = absMax(app.deltaA)
        R = app.deltaA - app.gaResultFit
        fin = np.abs(R[np.isfinite(R)])
        clr = float(np.max(fin)) if fin.size else 0.0
        if not clr > 0:
            clr = cl / 10
        resRMS = float(np.sqrt(np.mean(R[np.isfinite(R)] ** 2))) \
            if fin.size else np.nan
        self.maps.draw(
            app.delay, app.wavelength, [app.deltaA, app.gaResultFit, R],
            [(-cl, cl), (-cl, cl), (-clr, clr)], [True, False, True],
            [(ui_str('Experimental  (click to pick \\lambda)'),
              'Experimental'), 'Global fit',
             ('Residual  (RMS = %.3g)' % resRMS, 'RMS = %.3g' % resRMS)],
            self.ddMapScale.currentData(), self.edMapThr.value(),
            self.edMapTmin.value(), self.edMapTmax.value(),
            self.getFitTRangeIfRestricted(), self.curKinWL)

    def applyGaPoint(self, wlClick):
        app = self.app
        self.curKinWL = float(min(max(wlClick, np.min(app.wavelength)),
                                  np.max(app.wavelength)))
        self.edKinWL.setValue(self.curKinWL)
        self.maps.setMarker(self.curKinWL)
        self.plotKinetics()

    def plotKinetics(self):
        app = self.app
        if app.gaResultFit is None:
            return
        iw = int(np.argmin(np.abs(app.wavelength - self.curKinWL)))
        wlActual = app.wavelength[iw]
        dataRow = app.deltaA[iw, :]
        fitRow = app.gaResultFit[iw, :]
        tl = 'Delay time (%s)' % app.tUnitAx()
        ax = self.axKinFit
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(app.delay, dataRow, 'k.', markersize=5, label='Data')
        ax.plot(app.delay, fitRow, 'r-', linewidth=1.6, label='Fit')
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel(tl)
        ax.set_ylabel(Z_LABEL)
        legend(ax, fontsize=8)
        at = '%s = %.2f %s' % (app.xSym(), wlActual, app.xUnitAx())
        setTitle(ax, 'Kinetics at  %s   (requested %.2f %s)' % (
            at, self.curKinWL, app.xUnitAx()), 'Kinetics at  ' + at, at)
        sciTicks(ax, 'y')
        ax2 = self.axKinRes
        ax2.clear()
        styleAxes(ax2, grid=True)
        ax2.plot(app.delay, dataRow - fitRow, '-', color=(0.45, 0.45, 0.45),
                 linewidth=0.9)
        ax2.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax2.set_xlabel(tl)
        ax2.set_ylabel(RES_LABEL)
        setTitle(ax2, 'Fit residual at ' + at, 'Fit residual')
        sciTicks(ax2, 'y')
        useLog = self.ddKinScale.currentText() == 'Log' and \
            bool(np.any(app.delay > 0))
        if useLog:
            xl = [float(np.min(app.delay[app.delay > 0])),
                  float(np.max(app.delay))]
        else:
            xl = [float(np.min(app.delay)), float(np.max(app.delay))]
        fr = self.getFitTRangeIfRestricted()
        if fr is not None:
            if useLog:
                fr[0] = max(fr[0], xl[0])
            if fr[0] < fr[1]:
                xl = fr
        for a in (ax, ax2):
            a.set_xscale('log' if useLog else 'linear')
            a.set_xlim(xl[0], xl[1])
        self.cvKin.draw_idle()

    # ------------------------------------------------------------------
    # exports
    # ------------------------------------------------------------------
    def _putTable(self, title, defName):
        full = self.app.putFile('CSV (*.csv);;DAT (tab) (*.dat)', title,
                                defName, parent=self)
        if not full:
            return None, None
        return full, (',' if full.lower().endswith('.csv') else '\t')

    def exportDADS(self):
        app = self.app
        if app.gaResultDADS is None:
            self.alert('Run the fit first.', 'No results')
            return
        nMain = len(app.gaResultTau)
        hdr = [app.xHdr()] + ['DADS_tau=%.4g_%s' % (app.gaResultTau[i],
                                                   app.tUnitHdr())
                              for i in range(nMain)]
        if app.gaResultHasInf:
            hdr.append('DADS_tau=inf')
        cols = [app.wavelength] + [app.gaResultDADS[:, i]
                                   for i in range(app.gaResultDADS.shape[1])]
        ca = app.gaResultCA
        if isinstance(ca, dict) and \
                np.shape(ca['spec'])[0] == app.wavelength.size:
            caHdr = ['artifact_Gaussian', 'artifact_1st_derivative',
                     'artifact_2nd_derivative']
            for i in range(min(3, np.shape(ca['spec'])[1])):
                cols.append(ca['spec'][:, i])
                hdr.append(caHdr[i])
        full, dlm = self._putTable('Export DADS', 'DADS.csv')
        if not full:
            return
        try:
            writeTable(full, hdr, cols, dlm)
            self.inform('Saved: %s' % full, 'DADS exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')

    def exportEADS(self):
        app = self.app
        if app.gaResultEADS is None:
            self.alert('Spectra not available. Either the fit has not been '
                       'run yet, or the computation failed (near-degenerate '
                       'lifetimes).', 'No results')
            return
        if np.any(self.resStretch()):
            lab, sci = 'SAS', 'Species-associated spectra'
        else:
            lab, sci = 'EADS', 'Evolution-associated DS'
        E = np.asarray(app.gaResultEADS, float)
        nDecay = len(app.gaResultTau)
        nTot = E.shape[1]
        tauSorted = np.sort(np.asarray(app.gaResultTau, float))
        hdr = [app.xHdr()] + ['%s_S%d_tau=%.4g_%s' % (
            lab, i + 1, tauSorted[i], app.tUnitHdr()) for i in range(nDecay)]
        if app.gaResultHasInf:
            hdr.append('%s_S%d_tau=inf' % (lab, nTot))
        full, dlm = self._putTable('Export %s' % lab, lab + '.csv')
        if not full:
            return
        try:
            writeTable(full, hdr, [app.wavelength] + [E[:, i]
                                                      for i in range(nTot)],
                       dlm)
            self.inform('Saved %s (%s, %d species).\n\nFile: %s' % (
                lab, sci, nTot, full), '%s exported' % lab)
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')

    def exportMatrix2D(self, which, lab, tag):
        app = self.app
        if app.gaResultFit is None:
            self.alert('Run the fit first.', 'No results')
            return
        M = app.gaResultFit if which == 'fit' \
            else app.deltaA - app.gaResultFit
        full = app.putFile(
            'CSV comma-delimited (*.csv);;DAT tab-delimited (*.dat);;'
            'TXT tab-delimited (*.txt)', 'Export %s (2D matrix)' % lab,
            'GA_%s.csv' % tag, parent=self)
        if not full:
            return
        dlm = ',' if full.lower().endswith('.csv') else '\t'
        try:
            write_data_file(full, app.wavelength, app.delay, M, dlm)
            self.inform('Saved %s matrix  (%d λ x %d t).\n\nFile: %s' % (
                lab, app.wavelength.size, app.delay.size, full),
                lab + ' exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')

    def exportKinetics(self):
        app = self.app
        if app.gaResultFit is None:
            self.alert('Run the fit first.', 'No results')
            return
        iw = int(np.argmin(np.abs(app.wavelength - self.curKinWL)))
        wlActual = app.wavelength[iw]
        dataCol = app.deltaA[iw, :]
        fitCol = app.gaResultFit[iw, :]
        full, dlm = self._putTable(
            'Export kinetics @ %.2f %s' % (wlActual, app.xUnitTxt()),
            'kinetics_%.1f%s.csv' % (wlActual, app.xUnitTxt()))
        if not full:
            return
        try:
            writeTable(full, ['delay_%s' % app.tUnitHdr(), 'data_dA',
                              'fit_dA', 'residual_dA'],
                       [app.delay, dataCol, fitCol, dataCol - fitCol], dlm)
            self.inform('Saved kinetics at %s = %.2f %s.\n\nColumns: delay, '
                        'data, fit, residual.\nFile: %s' % (
                            'nu' if app.specUnit == 'cm-1' else 'lambda',
                            wlActual, app.xUnitTxt(), full),
                        'Kinetics exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')

    def exportPopulation(self):
        app = self.app
        popC, popT, nSp = self.computePopulation()
        if popC is None:
            self.alert('Run a fit first.', 'No result')
            return
        full = app.putFile('CSV (*.csv)', 'Export population kinetics',
                           'ga_populations.csv', parent=self)
        if not full:
            return
        try:
            write_matrix(full, np.column_stack([popT, popC]),
                         header='t_' + app.tUnitTxt() + ''.join(
                             ',S%d' % (j + 1) for j in range(nSp)))
            self.lblStatus.setText('Populations saved: %s'
                                   % os.path.basename(full))
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')
