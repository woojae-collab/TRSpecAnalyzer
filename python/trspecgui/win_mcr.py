"""MCR-ALS window. Port of openMcrWindow."""
import time

import numpy as np
from PySide6 import QtWidgets

from trspec.mcr import compute_mcr

from .textfit import legend, setTitle
from .widgets import (BlitManager, NumEdit, PlotCanvas, VScroll, button,
                      colorButton, label, lineColor, panel, runBusy,
                      sciTicks, sidePanel, styleAxes, textArea, tip)
from .win_base import Z_LABEL, ChildWindow
from .win_kinfit import writeTable

EPS = float(np.finfo(float).eps)
RUN_COLOR = (0.85, 0.50, 0.20)
BLUE = (0.0, 0.0, 1.0)
IDLE_TEXT = ['Click "Run MCR-ALS" to compute.']
T_C = 'C profiles (concentration vs time)'
T_S = 'S spectra'
T_REC = 'Reconstruction at clicked λ'
T_CONV = 'Convergence'
INIT_ITEMS = [('SVD-based (default)', 'svd'), ('Random', 'random'),
              ('GA result: EADS / SAS', 'gaeads'), ('GA result: DADS', 'ga')]


def emptyResult():
    return {'C': None, 'S': None, 'Drec': None, 'lof': np.nan, 'iter': 0,
            'lofHistory': np.zeros(0), 'converged': False}


def sub(n):
    """n as subscript digits (C_{1} of the MATLAB legends)."""
    return ''.join(chr(8320 + int(ch)) for ch in str(int(n)))


def tf2str(b):
    return 'ON' if b else 'off'


class McrWindow(ChildWindow):
    KEY = 'mcr'

    def __init__(self, app):
        super().__init__(app, 'MCR-ALS', 1380, 880)
        self.res = emptyResult()
        self.curLambda = float(app.selWL)
        self.xhS = None
        self.draggingS = False
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: setup ==================================
        left, lv = sidePanel('MCR-ALS setup', 340)
        gl.addWidget(left)

        def bold(text):
            lb = label(text)
            lb.setStyleSheet('font-weight: bold;')
            lv.addWidget(lb)

        def hrow(spacing=6):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(spacing)
            lv.addLayout(r)
            return r

        r = hrow()
        r.addWidget(label('Components (N):'))
        self.edN = NumEdit(3, limits=(1, 30), integer=True)
        r.addWidget(self.edN, 1)
        bold('Constraints:')
        self.cbNNc = QtWidgets.QCheckBox(
            'Non-negativity on C (concentrations)')
        self.cbNNc.setChecked(True)
        self.cbUni = QtWidgets.QCheckBox(
            'Unimodality on C (single peak per profile)')
        self.cbNNs = QtWidgets.QCheckBox(
            'Non-negativity on S (spectra) - usually OFF for TA')
        for cb in (self.cbNNc, self.cbUni, self.cbNNs):
            lv.addWidget(cb)
        bold('Initialization:')
        r = hrow()
        r.addWidget(label('Method:'))
        self.ddInit = QtWidgets.QComboBox()
        for nm, key in INIT_ITEMS:
            self.ddInit.addItem(nm, key)
        self.ddInit.setToolTip(tip(
            'Where the iteration starts. MCR-ALS can have several '
            'solutions of the same quality, so the start can change the '
            'result: starts that agree indicate a well-determined '
            'solution.\nSVD-based: the leading singular vectors of the '
            'data. No assumptions, the same result every run.\nRandom: '
            'random positive profiles. Every run takes a new seed (shown '
            'in the result box) unless one is entered below.\nGA result: '
            'the spectra of the stored global fit as starting spectra '
            '(components the fit does not have are taken from the data). '
            'EADS / SAS start from the sequential picture -- populations '
            'that rise and decay -- and usually give species-like results. '
            'DADS start from the parallel picture (every component decays '
            'from time zero) and stay close to it.'))
        r.addWidget(self.ddInit, 1)
        r = hrow()
        self.lblSeed = label('Random seed:')
        r.addWidget(self.lblSeed)
        self.edSeed = NumEdit(
            0, limits=(0, 4294967295), integer=True, width=95,
            tooltip='Seed of the random start. 0: a new seed for every '
                    'run; the one used is shown in the result box. Enter '
                    'that number here to repeat a run exactly.')
        r.addWidget(self.edSeed)
        self.lblSeedHint = label('0 = new every run')
        r.addWidget(self.lblSeedHint, 1)
        bold('Convergence:')
        r = hrow()
        r.addWidget(label('Tol:'))
        self.edTol = NumEdit(1e-4, limits=(EPS, np.inf), width=80)
        r.addWidget(self.edTol)
        r.addWidget(label('Max iter:'))
        self.edMaxIt = NumEdit(200, limits=(1, 10000), integer=True,
                               width=60)
        r.addWidget(self.edMaxIt)
        r.addStretch(1)
        r = hrow()
        r.addWidget(label('Time scale:'))
        self.ddTimeScale = QtWidgets.QComboBox()
        self.ddTimeScale.addItem('Log', 'log')
        self.ddTimeScale.addItem('Linear', 'linear')
        self.ddTimeScale.setToolTip(tip(
            'Delay axis scale for the C(t) profiles and the '
            'reconstruction-at-lambda panel.'))
        r.addWidget(self.ddTimeScale, 1)
        r = hrow()
        self.btnRun = colorButton('Run MCR-ALS', RUN_COLOR)
        self.btnRun.setStyleSheet(self.btnRun.styleSheet()
                                  + 'QPushButton { font-weight: bold; }')
        self.btnReset = button('Reset')
        r.addWidget(self.btnRun, 1)
        r.addWidget(self.btnReset, 1)
        bold('Export:')
        r = hrow()
        self.btnExpC = button('C')
        self.btnExpS = button('S')
        self.btnExpRec = button('Recon.')
        self.btnExpConv = button('Conv.')
        for b in (self.btnExpC, self.btnExpS, self.btnExpRec,
                  self.btnExpConv):
            b.setEnabled(False)
            r.addWidget(b, 1)
        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        lv.addWidget(self.lblStatus)
        self.txtResult = textArea(IDLE_TEXT)
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
        self.cvC = PlotCanvas()
        self.cvS = PlotCanvas()
        self.cvRec = PlotCanvas()
        self.cvConv = PlotCanvas()
        for cv, (i, j) in ((self.cvC, (0, 0)), (self.cvS, (0, 1)),
                           (self.cvRec, (1, 0)), (self.cvConv, (1, 1))):
            cv.setMinimumHeight(210)
            grid.addWidget(cv, i, j)
        for k in (0, 1):
            grid.setRowStretch(k, 1)
            grid.setColumnStretch(k, 1)
        self.axC = self.cvC.figure.add_subplot(111)
        self.axS = self.cvS.figure.add_subplot(111)
        self.axRec = self.cvRec.figure.add_subplot(111)
        self.axConv = self.cvConv.figure.add_subplot(111)

        def lay1(cv, ax):
            W, H = cv.sizePt()
            cv.place(ax, 58.0, 34.0, W - 58.0 - 16.0, H - 34.0 - 22.0)

        self.cvC.layoutFcn = lambda: lay1(self.cvC, self.axC)
        self.cvS.layoutFcn = lambda: lay1(self.cvS, self.axS)
        self.cvRec.layoutFcn = lambda: lay1(self.cvRec, self.axRec)
        self.cvConv.layoutFcn = lambda: lay1(self.cvConv, self.axConv)
        self.blS = BlitManager(self.cvS)
        self.cvS.mpl_connect('button_press_event', self.onSPress)
        self.cvS.mpl_connect('motion_notify_event', self.onSMotion)
        self.cvS.mpl_connect('button_release_event',
                             lambda ev: self.onSUp())
        self._emptyPlots()

        # ---- what the controls do --------------------------------------
        self.btnRun.clicked.connect(
            lambda *_: self.safe(self.doRun, 'Run MCR-ALS'))
        self.btnReset.clicked.connect(
            lambda *_: self.safe(self.doReset, 'Reset'))
        self.ddTimeScale.currentIndexChanged.connect(
            lambda *_: self.onTimeScaleChange())
        self.ddInit.currentIndexChanged.connect(
            lambda *_: self.onInitChange())
        self.btnExpC.clicked.connect(
            lambda *_: self.safe(self.exportC, 'Export C'))
        self.btnExpS.clicked.connect(
            lambda *_: self.safe(self.exportS, 'Export S'))
        self.btnExpRec.clicked.connect(
            lambda *_: self.safe(self.exportRec, 'Export recon'))
        self.btnExpConv.clicked.connect(
            lambda *_: self.safe(self.exportConv, 'Export conv'))

        # ---- result of an earlier visit --------------------------------
        sv = app.mcrSaved
        if isinstance(sv, dict) and isinstance(sv.get('result'), dict) and \
                sv['result'].get('S') is not None and \
                np.shape(sv['result']['S'])[0] == app.wavelength.size and \
                np.shape(sv['result']['C'])[0] == app.delay.size:
            try:
                st = sv['settings']
                self.edN.setValue(st['N'])
                self.cbNNc.setChecked(bool(st['nnC']))
                self.cbUni.setChecked(bool(st['uniC']))
                self.cbNNs.setChecked(bool(st['nnS']))
                self._setCombo(self.ddInit, st['init'])
                if 'seed' in st:
                    self.edSeed.setValue(st['seed'])
                self.edTol.setValue(st['tol'])
                self.edMaxIt.setValue(st['maxIt'])
                self.ddTimeScale.blockSignals(True)
                self._setCombo(self.ddTimeScale, st['tscale'])
                self.ddTimeScale.blockSignals(False)
                r0 = sv['result']
                n = np.shape(r0['S'])[1] if np.ndim(r0['S']) == 2 else 1
                self.res = {
                    'C': np.asarray(r0['C'], float).reshape(-1, n),
                    'S': np.asarray(r0['S'], float).reshape(-1, n),
                    'Drec': np.asarray(r0['Drec'], float),
                    'lof': float(r0['lof']), 'iter': int(r0['iter']),
                    'lofHistory': np.atleast_1d(
                        np.asarray(r0['lofHistory'], float)),
                    'converged': bool(r0['converged'])}
                cl = sv['curLambda']
                self.curLambda = np.nan if cl is None else float(cl)
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
        self.onInitChange()

    # ------------------------------------------------------------------
    @staticmethod
    def _setCombo(dd, key):
        i = dd.findData(key)
        if i >= 0:
            dd.setCurrentIndex(i)

    def _setExport(self, on):
        for b in (self.btnExpC, self.btnExpS, self.btnExpRec,
                  self.btnExpConv):
            b.setEnabled(bool(on))

    def _emptyPlots(self):
        app = self.app
        for ax in (self.axC, self.axS, self.axRec, self.axConv):
            ax.clear()
            styleAxes(ax, grid=True)
        self.axC.set_xlabel('Delay time (%s)' % app.tUnitAx())
        self.axC.set_ylabel('Concentration')
        setTitle(self.axC, T_C, 'C profiles')
        self.axS.set_xlabel(app.xAxLabel())
        self.axS.set_ylabel('Pure spectrum (%s)' % Z_LABEL)
        setTitle(self.axS, T_S)
        self.axRec.set_xlabel('Delay time (%s)' % app.tUnitAx())
        self.axRec.set_ylabel(Z_LABEL)
        setTitle(self.axRec, T_REC)
        self.axConv.set_xlabel('ALS iteration')
        self.axConv.set_ylabel('LOF (%)')
        self.axConv.set_yscale('log')
        setTitle(self.axConv, T_CONV)
        self.xhS = None
        self.blS.setArtists([])
        for cv in (self.cvC, self.cvS, self.cvRec, self.cvConv):
            cv.draw_idle()

    def isLog(self, t):
        return self.ddTimeScale.currentData() == 'log' and bool(np.any(t > 0))

    def _timeAxis(self, ax, t):
        if self.isLog(t):
            ax.set_xscale('log')
            lo, hi = float(np.min(t[t > 0])), float(np.max(t))
        else:
            ax.set_xscale('linear')
            lo, hi = float(np.min(t)), float(np.max(t))
        if hi > lo:
            ax.set_xlim(lo, hi)

    # ------------------------------------------------------------------
    def doRun(self):
        app = self.app
        N = int(round(self.edN.value()))
        nnC, uniC, nnS = self.cbNNc.isChecked(), self.cbUni.isChecked(), \
            self.cbNNs.isChecked()
        tol = self.edTol.value()
        maxIt = int(round(self.edMaxIt.value()))
        mode = self.ddInit.currentData()
        seed = self.edSeed.value()
        Sinit = None
        initNote = ''       # why the SVD start was used instead
        gaKind = ''
        gaNote = ''
        if mode in ('ga', 'gaeads'):
            nWL = app.wavelength.size
            Sga = None
            E, Dd = app.gaResultEADS, app.gaResultDADS
            if mode == 'gaeads' and E is not None and np.size(E) and \
                    np.shape(E)[0] == nWL:
                Sga = np.asarray(E, float)
                so = app.gaResultStretchOn
                gaKind = 'SAS' if so is not None and np.any(so) else 'EADS'
            elif Dd is not None and np.size(Dd) and np.shape(Dd)[0] == nWL:
                Sga = np.asarray(Dd, float)
                gaKind = 'DADS'
                if mode == 'gaeads':
                    gaNote = '      (EADS not available)'
            if Sga is None:
                mode = 'svd'
                initNote = ' (no GA result)'
            else:
                Sga = Sga.reshape(nWL, -1)
                Sinit = Sga[:, :min(N, Sga.shape[1])]
                mode = 'custom'
        D = app.deltaA
        if D.shape[1] < N or D.shape[0] < N:
            self.alert('Need at least %d delays and %d wavelengths.'
                       % (N, N), 'Too few points')
            return
        self.lblStatus.setText('Running MCR-ALS...')
        t0 = time.time()
        o = runBusy(app, self, 'Running MCR-ALS...', lambda: compute_mcr(
            D, N, nnC, uniC, nnS, tol, maxIt, mode, None, Sinit, seed))
        elapsed = time.time() - t0
        self.res = {'C': np.asarray(o.C, float), 'S': np.asarray(o.S, float),
                    'Drec': np.asarray(o.S, float) @ np.asarray(o.C, float).T,
                    'lof': float(o.lof), 'iter': int(o.iters),
                    'lofHistory': np.asarray(o.lof_history, float),
                    'converged': bool(o.converged)}
        self._setExport(True)
        tag = 'converged' if o.converged else 'NOT converged'
        self.lblStatus.setText('Done. LOF = %.3f%% (%s, %d iters, %.1fs)'
                               % (o.lof, tag, o.iters, elapsed))
        L = ['Components: %d' % N,
             'LOF: %.4f%%' % o.lof,
             'Iterations: %d  (%s)' % (o.iters, tag),
             '',
             'Constraints:',
             '  NN(C):  %s' % tf2str(nnC),
             '  Uni(C): %s' % tf2str(uniC),
             '  NN(S):  %s' % tf2str(nnS),
             '']
        if mode == 'random':
            L.append('Init: random (seed %d)' % o.seed)
        elif mode == 'custom':
            if o.n_s_init is None or np.isnan(o.n_s_init):
                L.append('Init: SVD (GA spectra unusable)')
            else:
                L.append('Init: GA result (%d %s)' % (o.n_s_init, gaKind))
                if gaNote:
                    L.append(gaNote)
                if o.n_s_init < N:
                    L.append('      + %d from the data (SVD)'
                             % (N - o.n_s_init))
        else:
            L.append('Init: SVD' + initNote)
        L += ['', 'Click on either plot to see',
              'reconstruction at that wavelength.']
        self.txtResult.setPlainText('\n'.join(L))
        self.plotResults()
        app.mcrSaved = {
            'result': dict(self.res),
            'curLambda': float(self.curLambda),
            'status': self.lblStatus.text(),
            'lines': list(L),
            'settings': {
                'N': int(self.edN.value()), 'nnC': nnC, 'uniC': uniC,
                'nnS': nnS, 'init': self.ddInit.currentData(),
                'seed': int(self.edSeed.value()), 'tol': float(tol),
                'maxIt': int(self.edMaxIt.value()),
                'tscale': self.ddTimeScale.currentData()}}

    def plotResults(self):
        app = self.app
        r = self.res
        if r['C'] is None:
            return
        t = np.asarray(app.delay, float).ravel()
        wl = np.asarray(app.wavelength, float).ravel()
        N = r['C'].shape[1]
        ax = self.axC
        ax.clear()
        styleAxes(ax, grid=True)
        for j in range(N):
            ax.plot(t, r['C'][:, j], '-', color=lineColor(j), linewidth=1.4,
                    label='C' + sub(j + 1))
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel('Concentration')
        setTitle(ax, 'C profiles')
        self._timeAxis(ax, t)
        sciTicks(ax, 'y')
        legend(ax, fontsize=8)
        self.cvC.draw_idle()

        ax = self.axS
        ax.clear()
        styleAxes(ax, grid=True)
        for j in range(N):
            ax.plot(wl, r['S'][:, j], '-', color=lineColor(j), linewidth=1.4,
                    label='S' + sub(j + 1))
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        if np.isfinite(self.curLambda):
            self.xhS = ax.axvline(self.curLambda, color='k', linewidth=1.0,
                                  alpha=0.7)
            self.blS.setArtists([self.xhS])
        else:
            self.xhS = None
            self.blS.setArtists([])
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('Pure spectrum (%s)' % Z_LABEL)
        setTitle(ax, T_S)
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        sciTicks(ax, 'y')
        legend(ax, fontsize=8)
        self.cvS.draw_idle()

        ax = self.axConv
        ax.clear()
        styleAxes(ax, grid=True)
        h = r['lofHistory']
        if h.size > 1:
            ax.plot(np.arange(1, h.size + 1), h, '-o', color=BLUE,
                    markersize=3, linewidth=1.0)
        ax.set_xlabel('ALS iteration')
        ax.set_ylabel('LOF (%)')
        setTitle(ax, 'Convergence (final LOF = %.3f%%)' % r['lof'],
                 'Final LOF = %.3f%%' % r['lof'])
        ax.set_yscale('log')
        self.cvConv.draw_idle()
        self.updateRec()

    # ---- wavelength picked on the spectra panel ----------------------
    def onSPress(self, ev):
        if ev.button != 1 or ev.inaxes is not self.axS or \
                self.res['C'] is None:
            return
        if self.cvS.legendAt(ev) is not None:        # (moving the legend)
            return
        self.draggingS = True
        self.applyPoint(ev.xdata)

    def onSMotion(self, ev):
        if not self.draggingS or ev.x is None or ev.y is None:
            return
        try:
            x, _ = self.axS.transData.inverted().transform((ev.x, ev.y))
        except Exception:
            return
        self.applyPoint(x)

    def onSUp(self):
        self.draggingS = False

    def applyPoint(self, x):
        if x is None or not np.isfinite(x):
            return
        wl = self.app.wavelength
        self.curLambda = float(min(max(x, np.min(wl)), np.max(wl)))
        if self.xhS is not None:
            self.xhS.set_xdata([self.curLambda, self.curLambda])
            self.blS.update()
        self.updateRec()

    def updateRec(self):
        app = self.app
        r = self.res
        if r['Drec'] is None:
            return
        wl = np.asarray(app.wavelength, float).ravel()
        t = np.asarray(app.delay, float).ravel()
        iw = int(np.argmin(np.abs(wl - self.curLambda)))
        ax = self.axRec
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(t, app.deltaA[iw, :], 'k.', markersize=4, label='Data')
        ax.plot(t, r['Drec'][iw, :], 'r-', linewidth=1.4,
                label='MCR reconstr.')
        ax.plot(t, app.deltaA[iw, :] - r['Drec'][iw, :], '-',
                color=(0.5, 0.5, 0.5), linewidth=0.8, label='Residual')
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel(Z_LABEL)
        setTitle(ax, 'Reconstruction at %s = %.2f %s' % (
            app.xSym(), wl[iw], app.xUnitAx()))
        self._timeAxis(ax, t)
        sciTicks(ax, 'y')
        legend(ax, fontsize=8)
        self.cvRec.draw_idle()

    def onTimeScaleChange(self):
        if self.res['C'] is None:
            return
        self.plotResults()

    def onInitChange(self):
        rnd = self.ddInit.currentData() == 'random'
        self.edSeed.setEnabled(rnd)
        c1 = (0.10, 0.10, 0.10) if rnd else (0.60, 0.60, 0.60)
        c2 = (0.45, 0.45, 0.45) if rnd else (0.70, 0.70, 0.70)
        for lb, c in ((self.lblSeed, c1), (self.lblSeedHint, c2)):
            lb.setStyleSheet('color: rgb(%d,%d,%d);' % tuple(
                int(round(255 * x)) for x in c))

    def doReset(self):
        self.edN.setValue(3)
        self.cbNNc.setChecked(True)
        self.cbUni.setChecked(False)
        self.cbNNs.setChecked(False)
        self._setCombo(self.ddInit, 'svd')
        self.edSeed.setValue(0)
        self.onInitChange()
        self.edTol.setValue(1e-4)
        self.edMaxIt.setValue(200)
        self.res = emptyResult()
        self.ddTimeScale.blockSignals(True)
        self._setCombo(self.ddTimeScale, 'log')
        self.ddTimeScale.blockSignals(False)
        self._setExport(False)
        self._emptyPlots()
        self.lblStatus.setText('Reset.')
        self.txtResult.setPlainText('\n'.join(IDLE_TEXT))
        self.app.mcrSaved = None

    # ---- exports -----------------------------------------------------
    def exportC(self):
        app, r = self.app, self.res
        if r['C'] is None:
            return
        n = r['C'].shape[1]
        self.exportMat([app.delay] + [r['C'][:, j] for j in range(n)],
                       ['delay_' + app.tUnitHdr()]
                       + ['C%d' % (j + 1) for j in range(n)],
                       'MCR_C', 'Export C profiles')

    def exportS(self):
        app, r = self.app, self.res
        if r['S'] is None:
            return
        n = r['S'].shape[1]
        self.exportMat([app.wavelength] + [r['S'][:, j] for j in range(n)],
                       [app.xHdr()] + ['S%d' % (j + 1) for j in range(n)],
                       'MCR_S', 'Export S spectra')

    def exportRec(self):
        app, r = self.app, self.res
        if r['Drec'] is None:
            return
        hdr = [app.xHdr()] + ['t=%.4g_%s' % (x, app.tUnitHdr())
                              for x in np.ravel(app.delay)]
        self.exportMat([app.wavelength] + [r['Drec'][:, j] for j in
                                           range(r['Drec'].shape[1])],
                       hdr, 'MCR_recon', 'Export reconstruction')

    def exportConv(self):
        h = self.res['lofHistory']
        if self.res['C'] is None:
            return
        self.exportMat([np.arange(1, h.size + 1), h], ['iter', 'LOF_pct'],
                       'MCR_conv', 'Export convergence')

    def exportMat(self, cols, hdr, defName, dlgTitle):
        full = self.app.putFile('CSV (*.csv);;DAT (tab) (*.dat)', dlgTitle,
                                defName + '.csv', parent=self)
        if not full:
            return
        dlm = ',' if full.lower().endswith('.csv') else '\t'
        try:
            writeTable(full, hdr, cols, dlm)
            self.inform('Saved: %s' % full, 'Exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Error')
