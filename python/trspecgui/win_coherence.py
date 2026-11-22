"""Vibrational Coherence Analysis window and its Time-Frequency window.
Ports of openCoherenceWindow and openTimeFrequencyWindow."""
import re
import time

import numpy as np
from PySide6 import QtWidgets

from trspec.coherence import compute_coherence, compute_time_frequency
from trspec.single import fit_single_trace

from .state import ui_str
from .textfit import setTitle
from .widgets import (AX_BG, FONT_PT, BlitManager, NumEdit, PlotCanvas,
                      VScroll, button, colorButton, label, panel, runBusy,
                      sciTicks, sidePanel, styleAxes, textArea, tip)
from .win_base import ChildWindow, mapLayout, openChild
from .win_kinfit import numText, writeTable

RUN_COLOR = (0.30, 0.65, 0.75)
RUN_COLOR_TF = (0.25, 0.55, 0.75)
BLUE = (0.0, 0.0, 1.0)
GREEN = (0.20, 0.55, 0.40)
PURPLE = (0.55, 0.30, 0.70)
IDLE_TEXT = ['Click "Run Coherence" to compute the FFT map.']
T_MAP = '|FFT|$^2$ map'
T_PS = 'Power spectrum at clicked λ'
T_AVG = 'Average over λ range'
T_PROF = 'Spectral profile at clicked ν'
T_RES = 'Residual time trace at clicked λ'
P_LABEL = '|FFT|$^2$'
TF_IDLE = ['Click "Run" to compute the ν-t map.']
TF_MAP = 'ν-t map'
TF_DECAY = 'Power vs time at selected ν'
TF_SPEC = 'Spectrum at selected t'
TF_TRACE = 'Input trace (analysis window)'


def freqAxLabel(unit):
    u = str(unit).lower()
    if u == 'cm-1':
        return 'Wavenumber (cm$^{-1}$)'
    if u == 'thz':
        return 'Frequency (THz)'
    return 'Frequency (Hz)'


def freqUnitTxt(unit):
    """Unit of the frequency axis as plain text (buttons, text box)."""
    u = str(unit).lower()
    return 'cm-1' if u == 'cm-1' else ('THz' if u == 'thz' else 'Hz')


def freqUnitAx(unit):
    """... and for a plot title."""
    u = str(unit).lower()
    return 'cm$^{-1}$' if u == 'cm-1' else ('THz' if u == 'thz' else 'Hz')


def validName(s):
    """Column name as matlab.lang.makeValidName(..., 'ReplacementStyle',
    'underscore') gives it."""
    s = re.sub(r'[^A-Za-z0-9_]', '_', str(s))
    if not re.match(r'[A-Za-z]', s):
        s = 'x' + s
    return s


def nanmaxOr(a, default=np.nan):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(np.max(a)) if a.size else default


def setCombo(dd, key):
    i = dd.findData(key)
    if i >= 0:
        dd.setCurrentIndex(i)


def lay1(cv, ax, L=58.0, B=34.0, R=16.0, T=22.0):
    W, H = cv.sizePt()
    cv.place(ax, L, B, W - L - R, H - B - T)


def layStack(cv, top, bot, gap, L=58.0, B=34.0, R=16.0, T=22.0):
    """Two axes above each other on one canvas."""
    W, H = cv.sizePt()
    h = max((H - T - B - gap) / 2.0, 4.0)
    cv.place(bot, L, B, W - L - R, h)
    cv.place(top, L, B + h + gap, W - L - R, h)


def exportTable(win, cols, hdr, defName, dlgTitle):
    """Table with sanitized column names (exportTbl of both windows)."""
    full = win.app.putFile('CSV (*.csv);;DAT (tab) (*.dat)', dlgTitle,
                           defName + '.csv', parent=win)
    if not full:
        return
    dlm = ',' if full.lower().endswith('.csv') else '\t'
    try:
        writeTable(full, [validName(h) for h in hdr], cols, dlm)
        win.inform('Saved: %s' % full, 'Exported')
    except Exception as err:
        win.alert('Export failed:\n%s' % err, 'Error')


def writeMatrix(path, first, firstFmt, head, rows, dlm):
    """Header line (first, then one number per column) and the rows."""
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(first + ''.join(dlm + (firstFmt % v) for v in head) + '\n')
        for r in rows:
            fh.write(dlm.join(numText(v) for v in r) + '\n')


def emptyCoh():
    return {'Pmap': None, 'freq': None, 'Rresid': None, 'tWin': None,
            'RinWin': None, 'info': {}}


class CoherenceWindow(ChildWindow):
    KEY = 'coherence'

    def __init__(self, app):
        super().__init__(app, 'Vibrational Coherence Analysis', 1380, 900)
        u = app.tUnitTxt()
        self.res = emptyCoh()
        self.curLambda = float(app.selWL)
        self.curFreq = np.nan
        self.xhV = self.xhH = self.yhAvgLo = self.yhAvgHi = None
        self.dragging = False
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: setup ==================================
        left, lv = sidePanel('Coherence setup', 360)
        lv.setSpacing(5)
        gl.addWidget(left)

        def bold(text):
            lb = label(text)
            lb.setStyleSheet('font-weight: bold;')
            lv.addWidget(lb)

        def pair(text, field):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(6)
            lb = label(text)
            r.addWidget(lb)
            r.addWidget(field, 1)
            lv.addLayout(r)
            return lb

        def pair2(t1, f1, t2, f2):
            """Two fields in one row (the MATLAB window has them below
            each other; a row saves the height a small screen lacks)."""
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(5)
            r.addWidget(label(t1))
            r.addWidget(f1, 1)
            r.addWidget(label(t2))
            r.addWidget(f2, 1)
            lv.addLayout(r)

        bold('Residual source:')
        self.ddResid = QtWidgets.QComboBox()
        self.ddResid.addItem('Use GA residual', 'ga')
        self.ddResid.addItem('Per-pixel exp fit', 'pixfit')
        pair('Method:', self.ddResid)
        self.edNexp = NumEdit(
            2, limits=(1, 5), integer=True,
            tooltip='Number of exponentials fitted at every wavelength. '
                    'Used by "Per-pixel exp fit" only.')
        self.lblNexp = pair('Per-pixel N exp:', self.edNexp)
        bold('Time window:')
        self.edTMin = NumEdit(self.defaultTMin())
        self.edTMax = NumEdit(self.defaultTMax())
        pair2('t_min (%s):' % u, self.edTMin, 't_max (%s):' % u, self.edTMax)
        bold('FFT options:')
        self.ddApod = QtWidgets.QComboBox()
        for nm, key in (('Rectangular', 'rect'), ('Hann', 'hann'),
                        ('Hamming', 'hamming'), ('Blackman', 'blackman')):
            self.ddApod.addItem(nm, key)
        self.ddApod.setCurrentIndex(1)
        pair('Apodization:', self.ddApod)
        self.edZpad = NumEdit(2, limits=(1, 16), integer=True)
        pair('Zero-pad factor:', self.edZpad)
        self.ddDet = QtWidgets.QComboBox()
        for nm, key in (('None', 'none'), ('Mean', 'mean'),
                        ('Linear', 'linear')):
            self.ddDet.addItem(nm, key)
        self.ddDet.setCurrentIndex(2)
        pair('Detrend:', self.ddDet)
        bold('Output:')
        self.ddFreqU = QtWidgets.QComboBox()
        self.ddFreqU.addItem(ui_str('Wavenumber (cm^{-1})'), 'cm-1')
        self.ddFreqU.addItem('THz', 'thz')
        pair('Frequency unit:', self.ddFreqU)
        self.ddNorm = QtWidgets.QComboBox()
        for nm, key in (('Peak (global)', 'peak'),
                        ('Per-wavelength', 'perwl'), ('None', 'none')):
            self.ddNorm.addItem(nm, key)
        pair('Normalize:', self.ddNorm)
        self.ddYScale = QtWidgets.QComboBox()
        self.ddYScale.addItem('Linear', 'linear')
        self.ddYScale.addItem('Log', 'log')
        self.ddYScale.setToolTip(tip(
            'Apply log scale to both top plots: the |FFT|² colormap on the '
            'heatmap (z-axis) and the y-axis of the power spectrum at λ. '
            'Log highlights weak modes.'))
        pair('Heatmap z / 1D y scale:', self.ddYScale)
        bold(ui_str('Average power over \\lambda range:'))
        wl = np.asarray(app.wavelength, float)
        self.wlLo0, self.wlHi0 = float(np.min(wl)), float(np.max(wl))
        lim = (self.wlLo0, max(self.wlHi0, self.wlLo0))
        self.edAvgWLmin = NumEdit(self.wlLo0, limits=lim)
        self.edAvgWLmax = NumEdit(self.wlHi0, limits=lim)
        xs, xu = app.xSymTxt(), app.xUnitTxt()
        pair2('%s_min (%s):' % (xs, xu), self.edAvgWLmin,
              '%s_max (%s):' % (xs, xu), self.edAvgWLmax)
        bold(ui_str('Display \\nu range:'))
        self.edDispVmin = NumEdit(0.0)
        self.edDispVmax = NumEdit(
            0.0, tooltip='Visualization window on the wavenumber axis. '
                         'ν_max = 0 means "auto" (full Nyquist range of the '
                         'current FFT).')
        pair2(ui_str('\\nu_min:'), self.edDispVmin,
              ui_str('\\nu_max:'), self.edDispVmax)
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        self.btnRun = colorButton('Run Coherence', RUN_COLOR)
        self.btnRun.setStyleSheet(self.btnRun.styleSheet()
                                  + 'QPushButton { font-weight: bold; }')
        self.btnReset = button('Reset')
        r.addWidget(self.btnRun, 1)
        r.addWidget(self.btnReset, 1)
        lv.addLayout(r)
        self.btnTF = button(
            'Time-frequency (STFT / wavelet)...',
            tooltip='Open the nu-t analysis of this residual: when does '
                    'each mode live and how fast does it dephase. '
                    'Available after a run.')
        self.btnTF.setEnabled(False)
        lv.addWidget(self.btnTF)
        bold('Export:')
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(5)
        self.btnExpMap = button(ui_str('|FFT|^2 map'))
        self.btnExpSpec = button(ui_str('P(\\nu)@\\lambda'))
        self.btnExpAvg = button(
            ui_str('Avg P(\\nu)'),
            tooltip='Export the power spectrum averaged over the current '
                    '"Average over λ range" window.')
        self.btnExpProf = button(ui_str('A(\\lambda)@\\nu'))
        for b in (self.btnExpMap, self.btnExpSpec, self.btnExpAvg,
                  self.btnExpProf):
            b.setEnabled(False)
            r.addWidget(b, 1)
        lv.addLayout(r)
        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        lv.addWidget(self.lblStatus)
        self.txtResult = textArea(IDLE_TEXT)
        self.txtResult.setMinimumHeight(100)
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
        self.cvP = PlotCanvas()
        self.cvProf = PlotCanvas()
        self.cvRes = PlotCanvas()
        self.cvMap.setMinimumHeight(250)
        self.cvP.setMinimumHeight(250)
        self.cvProf.setMinimumHeight(215)
        self.cvRes.setMinimumHeight(215)
        grid.addWidget(self.cvMap, 0, 0)
        grid.addWidget(self.cvP, 0, 1)
        grid.addWidget(self.cvProf, 1, 0)
        grid.addWidget(self.cvRes, 1, 1)
        grid.setRowStretch(0, 13)
        grid.setRowStretch(1, 10)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        self.axMap = self.cvMap.figure.add_subplot(111)
        self.caxMap = self.cvMap.figure.add_axes([0.9, 0.1, 0.02, 0.8])
        self.axPSpec = self.cvP.figure.add_subplot(211)
        self.axPSpecAvg = self.cvP.figure.add_subplot(212)
        self.axProf = self.cvProf.figure.add_subplot(211)
        self.axProfPh = self.cvProf.figure.add_subplot(212)
        self.axRes = self.cvRes.figure.add_subplot(111)
        self.cvMap.layoutFcn = self._layMap
        self.cvP.layoutFcn = lambda: layStack(
            self.cvP, self.axPSpec, self.axPSpecAvg, 26.0)
        self.cvProf.layoutFcn = lambda: layStack(
            self.cvProf, self.axProf, self.axProfPh, 10.0)
        self.cvRes.layoutFcn = lambda: lay1(self.cvRes, self.axRes)
        self.blMap = BlitManager(self.cvMap)
        self.cvMap.mpl_connect('button_press_event', self.onMapPress)
        self.cvMap.mpl_connect('motion_notify_event', self.onMapMotion)
        self.cvMap.mpl_connect('button_release_event',
                               lambda ev: self.onMapUp())
        self._emptyPlots()

        # ---- what the controls do --------------------------------------
        self.btnRun.clicked.connect(
            lambda *_: self.safe(self.doRun, 'Run Coherence'))
        self.btnReset.clicked.connect(
            lambda *_: self.safe(self.doReset, 'Reset'))
        self.btnTF.clicked.connect(
            lambda *_: self.safe(app.openTimeFrequencyWindow,
                                 'Time-frequency'))
        self.ddResid.currentIndexChanged.connect(
            lambda *_: self.safe(self.syncNexpEnable, 'Method'))
        self.ddYScale.currentIndexChanged.connect(
            lambda *_: self._user(self.onYScaleChange))
        self.edAvgWLmin.valueChanged.connect(
            lambda *_: self._user(self.onAvgRangeChange))
        self.edAvgWLmax.valueChanged.connect(
            lambda *_: self._user(self.onAvgRangeChange))
        self.edDispVmin.valueChanged.connect(
            lambda *_: self._user(self.onDispRangeChange))
        self.edDispVmax.valueChanged.connect(
            lambda *_: self._user(self.onDispRangeChange))
        self.btnExpMap.clicked.connect(
            lambda *_: self.safe(self.exportMap, 'Export map'))
        self.btnExpSpec.clicked.connect(
            lambda *_: self.safe(self.exportSpec, 'Export spec'))
        self.btnExpAvg.clicked.connect(
            lambda *_: self.safe(self.exportAvg, 'Export avg spec'))
        self.btnExpProf.clicked.connect(
            lambda *_: self.safe(self.exportProf, 'Export profile'))

        # ---- result of an earlier visit --------------------------------
        self._quiet = True
        sv = app.cohSaved
        if isinstance(sv, dict) and isinstance(sv.get('result'), dict) and \
                sv['result'].get('Pmap') is not None and \
                np.shape(sv['result']['Pmap'])[0] == app.wavelength.size:
            try:
                st = sv['settings']
                setCombo(self.ddResid, st['resid'])
                self.edNexp.setValue(st['nexp'])
                self.edTMin.setValue(st['tmin'])
                self.edTMax.setValue(st['tmax'])
                setCombo(self.ddApod, st['apod'])
                self.edZpad.setValue(st['zpad'])
                setCombo(self.ddDet, st['detrend'])
                setCombo(self.ddFreqU, st['frequ'])
                setCombo(self.ddNorm, st['norm'])
                setCombo(self.ddYScale, st['yscale'])
                self.edAvgWLmin.setValue(st['avgLo'])
                self.edAvgWLmax.setValue(st['avgHi'])
                self.edDispVmin.setValue(st['dispLo'])
                self.edDispVmax.setValue(st['dispHi'])
                r0 = sv['result']
                info = dict(r0['info'])
                if info.get('Phimap') is not None:
                    info['Phimap'] = np.asarray(info['Phimap'], float)
                self.res = {
                    'Pmap': np.asarray(r0['Pmap'], float),
                    'freq': np.asarray(r0['freq'], float).ravel(),
                    'Rresid': np.asarray(r0['Rresid'], float),
                    'tWin': np.asarray(r0['tWin'], float).ravel(),
                    'RinWin': np.asarray(r0['RinWin'], float),
                    'info': info}
                cl, cf = sv['curLambda'], sv['curFreq']
                self.curLambda = np.nan if cl is None else float(cl)
                self.curFreq = np.nan if cf is None else float(cf)
                self._setExport(True)
                self.plotResults()
                self.lblStatus.setText('%s  [restored]' % sv['status'])
                self.txtResult.setPlainText('\n'.join(sv['lines']))
            except Exception:
                self.res = emptyCoh()
                self._setExport(False)
                self._emptyPlots()
                self.lblStatus.setText(
                    'Stored result could not be restored -- run again.')
        self._quiet = False
        self.syncNexpEnable()

    # ------------------------------------------------------------------
    def _user(self, fcn):
        if not self._quiet:
            self.safe(fcn, 'Display')

    def _layMap(self):
        mapLayout(self.cvMap, self.axMap, self.caxMap, L=58.0, R=74.0,
                  T=22.0, B=34.0)

    def _setExport(self, on):
        for b in (self.btnExpMap, self.btnExpSpec, self.btnExpAvg,
                  self.btnExpProf, self.btnTF):
            b.setEnabled(bool(on))

    def _emptyPlots(self):
        app = self.app
        for ax in (self.axMap, self.axPSpec, self.axPSpecAvg, self.axProf,
                   self.axProfPh, self.axRes):
            ax.clear()
        self.caxMap.clear()
        self.caxMap.set_visible(False)
        styleAxes(self.axMap)
        self.axMap.set_xlabel(freqAxLabel('cm-1'))
        self.axMap.set_ylabel(app.xAxLabel())
        setTitle(self.axMap, T_MAP)
        for ax in (self.axPSpec, self.axPSpecAvg, self.axProf, self.axProfPh,
                   self.axRes):
            styleAxes(ax, grid=True)
        self.axPSpec.set_ylabel(P_LABEL)
        self.axPSpec.tick_params(labelbottom=False)
        self.axProf.tick_params(labelbottom=False)
        setTitle(self.axPSpec, T_PS)
        self.axPSpecAvg.set_xlabel(freqAxLabel('cm-1'))
        self.axPSpecAvg.set_ylabel(P_LABEL)
        setTitle(self.axPSpecAvg, T_AVG)
        self.axProf.set_ylabel(P_LABEL)
        setTitle(self.axProf, T_PROF)
        self.axProfPh.set_xlabel(app.xAxLabel())
        self.axProfPh.set_ylabel('Phase (rad)')
        self.axRes.set_xlabel('Delay time (%s)' % app.tUnitAx())
        self.axRes.set_ylabel('Residual')
        setTitle(self.axRes, T_RES)
        self.xhV = self.xhH = self.yhAvgLo = self.yhAvgHi = None
        self.blMap.setArtists([])
        self._layMap()
        for cv in (self.cvMap, self.cvP, self.cvProf, self.cvRes):
            cv.draw_idle()

    def syncNexpEnable(self):
        on = self.ddResid.currentData() == 'pixfit'
        self.edNexp.setEnabled(on)
        self.lblNexp.setEnabled(on)

    def defaultTMin(self):
        app = self.app
        fw = app.gaResultFWHM
        if fw is None or not np.isfinite(fw):
            fw = app.gaFWHM
        sigma = fw / (2 * np.sqrt(2 * np.log(2))) \
            if fw is not None and np.isfinite(fw) else 0.0
        base = {'fs': 100.0, 'ps': 0.3, 'ns': 3e-4, 'us': 3e-7}.get(
            app.timeUnit, 0.3)
        return max(base, 3 * sigma)

    def defaultTMax(self):
        app = self.app
        base = {'fs': 5000.0, 'ps': 5.0, 'ns': 0.005, 'us': 5e-6}.get(
            app.timeUnit, 5.0)
        return min(base, float(np.max(app.delay)))

    # ---- display range -----------------------------------------------
    def getDispLims(self):
        fr = self.res['freq']
        if fr is None or fr.size == 0:
            return None
        fLo, fHi = float(np.min(fr)), float(np.max(fr))
        lo, hi = self.edDispVmin.value(), self.edDispVmax.value()
        if not np.isfinite(hi) or hi <= 0:
            hi = fHi
        if not np.isfinite(lo):
            lo = fLo
        if lo > hi:
            lo, hi = hi, lo
        lo, hi = max(lo, fLo), min(hi, fHi)
        if lo >= hi:
            return None
        return [lo, hi]

    def getXlimForFreq(self, freq):
        lims = self.getDispLims()
        if lims is None:
            lims = [float(np.min(freq)), float(np.max(freq))]
        return lims

    def onYScaleChange(self):
        if self.res['Pmap'] is not None:
            self.plotResults()

    def onAvgRangeChange(self):
        if self.res['Pmap'] is None:
            return
        if self.yhAvgLo is not None:
            v = self.edAvgWLmin.value()
            self.yhAvgLo.set_ydata([v, v])
        if self.yhAvgHi is not None:
            v = self.edAvgWLmax.value()
            self.yhAvgHi.set_ydata([v, v])
        self.blMap.update()
        self.update1Ds()

    def onDispRangeChange(self):
        if self.res['Pmap'] is None:
            return
        self.plotResults()

    # ---- residual ----------------------------------------------------
    def computeResidual(self, method, Nexp):
        """(R, note): the oscillating part of the data."""
        app = self.app
        if method == 'ga':
            if app.gaResultFit is None or np.size(app.gaResultFit) == 0:
                raise ValueError(
                    'No GA fit available. Run Global Analysis first, or '
                    'switch to "Per-pixel exp fit".')
            return app.deltaA - app.gaResultFit, ''
        if method == 'pixfit':
            R, nFail, nTried = self.perPixelResidual(Nexp)
            note = ''
            if nFail > 0:
                note = '%d of %d rows: exp fit failed (mean-subtracted)' % (
                    nFail, nTried)
            return R, note
        raise ValueError('Unknown residual method: %s' % method)

    def perPixelResidual(self, Nexp):
        app = self.app
        D = app.deltaA
        nWl, nT = D.shape
        R = np.zeros((nWl, nT))
        t = app.delay
        t0in = app.gaResultT0 if app.gaResultT0 is not None else app.gaT0
        fwin = app.gaResultFWHM if app.gaResultFWHM is not None \
            else app.gaFWHM
        tauInit = np.logspace(np.log10(max(0.1, fwin)),
                              np.log10(float(np.max(t)) / 2), Nexp)
        nFail = nTried = 0
        lastMsg = ''
        for i in range(nWl):
            y = D[i, :]
            if not np.any(np.isfinite(y)):
                continue
            nTried += 1
            try:
                o = fit_single_trace(
                    t, y, tauInit, np.zeros(Nexp, bool), np.ones(Nexp),
                    np.ones(Nexp, bool), np.zeros(Nexp, bool), t0in, True,
                    fwin, True, False, 'numerical')
                R[i, :] = o.residual
            except Exception as err:
                R[i, :] = y - np.nanmean(y)
                nFail += 1
                lastMsg = str(err)
        if nTried > 0 and nFail == nTried:
            raise ValueError(
                'Per-pixel exp fit failed on every row (%s). Use the GA '
                'residual instead, or change the number of exponentials.'
                % lastMsg)
        return R, nFail, nTried

    # ------------------------------------------------------------------
    def doRun(self):
        app = self.app
        if self.edTMin.value() >= self.edTMax.value():
            self.alert('Need t_min < t_max.', 'Invalid')
            return
        method = self.ddResid.currentData()
        Nexp = int(round(self.edNexp.value()))
        tMin, tMax = self.edTMin.value(), self.edTMax.value()
        apod, zpad = self.ddApod.currentData(), \
            int(round(self.edZpad.value()))
        det, fu, nrm = self.ddDet.currentData(), self.ddFreqU.currentData(), \
            self.ddNorm.currentData()
        tU = app.timeUnit
        if method == 'pixfit':
            self.lblStatus.setText('Per-pixel exp fit on %d rows...'
                                   % app.deltaA.shape[0])
        else:
            self.lblStatus.setText('Computing residual...')
        box = {}

        def work():
            R, note = self.computeResidual(method, Nexp)
            t0 = time.time()
            o = compute_coherence(R, app.delay, tMin, tMax, apod, zpad, det,
                                  fu, nrm, tU)
            box['el'] = time.time() - t0
            return R, note, o

        Rresid, residNote, o = runBusy(app, self, 'Computing FFT map...',
                                       work)
        elapsed = box.get('el', 0.0)
        tA = np.asarray(app.delay, float).ravel()
        mask = (tA >= tMin) & (tA <= tMax)
        freq = np.asarray(o.freq, float).ravel()
        Pmap = np.asarray(o.pmap, float)
        info = {'dtUni': float(o.dt), 'Nfft': int(o.n_fft),
                'tUni': np.asarray(o.t_uni, float), 'Nwin': int(o.n_win),
                'apod': apod, 'zeroPad': int(o.zero_pad), 'detrend': det,
                'freqUnit': fu, 'normMode': nrm,
                'Phimap': np.asarray(o.phase, float)}
        self.res = {'Pmap': Pmap, 'freq': freq,
                    'Rresid': np.asarray(Rresid, float), 'tWin': tA[mask],
                    'RinWin': np.asarray(Rresid, float)[:, mask],
                    'info': info}
        self._setExport(True)
        app.tfSaved = None
        w = app.childWindows.get('tf')
        if w is not None:
            try:
                w.close()
            except Exception:
                pass
        self.lblStatus.setText('Done (%d freqs, %d wl, %.2fs).'
                               % (freq.size, Pmap.shape[0], elapsed))
        L = ['Residual: %s' % method]
        if residNote:
            L.append(residNote)
        L += ['Window: %.4g .. %.4g %s' % (tMin, tMax, app.tUnitTxt()),
              'dt (uniform): %.4g %s' % (info['dtUni'], app.tUnitTxt()),
              'Nfft: %d  (input %d)' % (info['Nfft'], info['Nwin']),
              'Apod: %s   Zero-pad: %dx' % (info['apod'], info['zeroPad']),
              'Detrend: %s' % info['detrend'],
              'Freq range: %.3g .. %.3g %s' % (freq[0], freq[-1],
                                               info['freqUnit']),
              '', 'Click on the heatmap to inspect',
              'a particular (lambda, nu).']
        self.txtResult.setPlainText('\n'.join(L))
        self.plotResults()
        app.cohSaved = {
            'result': dict(self.res),
            'curLambda': float(self.curLambda),
            'curFreq': float(self.curFreq),
            'status': self.lblStatus.text(), 'lines': list(L),
            'settings': {
                'resid': method, 'nexp': int(self.edNexp.value()),
                'tmin': float(tMin), 'tmax': float(tMax), 'apod': apod,
                'zpad': int(self.edZpad.value()), 'detrend': det,
                'frequ': fu, 'norm': nrm,
                'yscale': self.ddYScale.currentData(),
                'avgLo': self.edAvgWLmin.value(),
                'avgHi': self.edAvgWLmax.value(),
                'dispLo': self.edDispVmin.value(),
                'dispHi': self.edDispVmax.value()}}

    def plotResults(self):
        app = self.app
        r = self.res
        if r['Pmap'] is None:
            return
        wl = np.asarray(app.wavelength, float).ravel()
        freq, Pmap = r['freq'], r['Pmap']
        lims = self.getDispLims()
        Psub = Pmap
        if lims is not None:
            fm = (freq >= lims[0]) & (freq <= lims[1])
            if np.any(fm):
                Psub = Pmap[:, fm]
        if self.ddYScale.currentData() == 'log':
            pk = nanmaxOr(Psub)
            if not np.isfinite(pk) or pk <= 0:
                pk = 1.0
            floorVal = pk * 1e-4
            Cdat = np.log10(np.maximum(Pmap, floorVal))
            cLo, cHi = np.log10(floorVal), np.log10(pk)
            cbLabel = 'log$_{10}$|FFT|$^2$'
        else:
            Cdat = Pmap
            cbLabel = P_LABEL
            cl = nanmaxOr(Psub)
            if not np.isfinite(cl) or cl <= 0:
                cl = 1.0
            cLo, cHi = 0.0, cl
        ax, cax = self.axMap, self.caxMap
        ax.clear()
        cax.clear()
        styleAxes(ax)
        ax.set_facecolor(AX_BG)
        h = None
        if wl.size > 1 and freq.size > 1:
            h = ax.pcolormesh(freq, wl, np.ma.masked_invalid(Cdat[:-1, :-1]),
                              shading='flat', cmap='turbo', vmin=cLo,
                              vmax=cHi, rasterized=True)
        xl = self.getXlimForFreq(freq)
        ax.set_xlim(xl[0], xl[1])
        ax.set_ylim(float(np.min(wl)), float(np.max(wl)))
        ax.set_xlabel(freqAxLabel(r['info']['freqUnit']))
        ax.set_ylabel(app.xAxLabel())
        setTitle(ax, T_MAP)
        cax.set_visible(h is not None)
        if h is not None:
            cb = self.cvMap.figure.colorbar(h, cax=cax)
            cb.set_label(cbLabel, fontsize=FONT_PT)
            cb.ax.tick_params(labelsize=FONT_PT, length=3)
            if self.ddYScale.currentData() != 'log':
                sciTicks(cb.ax, 'y')
        if not np.isfinite(self.curLambda):
            self.curLambda = 0.5 * (float(np.min(wl)) + float(np.max(wl)))
        if not np.isfinite(self.curFreq):
            col = np.nanmax(np.where(np.isfinite(Pmap), Pmap, -np.inf),
                            axis=0)
            self.curFreq = float(freq[int(np.argmax(col))])
        kw = dict(color='k', linewidth=1.0, alpha=0.7)
        self.xhV = ax.axvline(self.curFreq, **kw)
        self.xhH = ax.axhline(self.curLambda, **kw)
        kd = dict(color=(0.30, 0.30, 0.30), linewidth=1.2, alpha=0.8,
                  linestyle=':')
        self.yhAvgLo = ax.axhline(self.edAvgWLmin.value(), **kd)
        self.yhAvgHi = ax.axhline(self.edAvgWLmax.value(), **kd)
        self.blMap.setArtists([self.xhV, self.xhH, self.yhAvgLo,
                               self.yhAvgHi])
        self._layMap()
        self.cvMap.draw_idle()
        self.update1Ds()

    # ---- the point picked on the map ---------------------------------
    def onMapPress(self, ev):
        if ev.button != 1 or ev.inaxes is not self.axMap or \
                self.res['Pmap'] is None:
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
        fr = self.res['freq']
        if fr is None or fr.size == 0 or x is None or y is None or \
                not np.isfinite(x) or not np.isfinite(y):
            return
        wl = self.app.wavelength
        self.curFreq = float(min(max(x, np.min(fr)), np.max(fr)))
        self.curLambda = float(min(max(y, np.min(wl)), np.max(wl)))
        if self.xhV is not None:
            self.xhV.set_xdata([self.curFreq, self.curFreq])
        if self.xhH is not None:
            self.xhH.set_ydata([self.curLambda, self.curLambda])
        self.blMap.update()
        self.update1Ds()

    def avgRange(self):
        lo = min(self.edAvgWLmin.value(), self.edAvgWLmax.value())
        hi = max(self.edAvgWLmin.value(), self.edAvgWLmax.value())
        wl = np.asarray(self.app.wavelength, float).ravel()
        return lo, hi, (wl >= lo) & (wl <= hi)

    def avgSpectrum(self):
        lo, hi, inR = self.avgRange()
        P = self.res['Pmap']
        if np.any(inR):
            with np.errstate(invalid='ignore'):
                Pavg = np.nanmean(P[inR, :], axis=0)
        else:
            Pavg = np.full(P.shape[1], np.nan)
        return Pavg, lo, hi, inR

    def update1Ds(self):
        app = self.app
        r = self.res
        if r['Pmap'] is None:
            return
        wl = np.asarray(app.wavelength, float).ravel()
        freq, Pmap = r['freq'], r['Pmap']
        iw = int(np.argmin(np.abs(wl - self.curLambda)))
        ifq = int(np.argmin(np.abs(freq - self.curFreq)))
        wlAct, fqAct = wl[iw], freq[ifq]
        ysc = self.ddYScale.currentData()
        xl = self.getXlimForFreq(freq)
        dash = dict(color='k', linestyle='--', linewidth=0.8, alpha=0.6)

        ax = self.axPSpec
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(freq, Pmap[iw, :], '-', color=BLUE, linewidth=1.4)
        ax.axvline(fqAct, **dash)
        ax.set_xlim(xl[0], xl[1])
        ax.set_ylabel(P_LABEL)
        setTitle(ax, 'Power spectrum at %s = %.2f %s' % (
            app.xSym(), wlAct, app.xUnitAx()))
        ax.set_yscale(ysc)
        ax.tick_params(labelbottom=False)
        if ysc != 'log':
            sciTicks(ax, 'y')

        ax = self.axPSpecAvg
        ax.clear()
        styleAxes(ax, grid=True)
        Pavg, lo, hi, inR = self.avgSpectrum()
        ax.plot(freq, Pavg, '-', color=GREEN, linewidth=1.4)
        ax.axvline(fqAct, **dash)
        ax.set_xlim(xl[0], xl[1])
        ax.set_xlabel(freqAxLabel(r['info']['freqUnit']))
        ax.set_ylabel(P_LABEL)
        n = int(np.count_nonzero(inR))
        setTitle(ax, 'Average over %s = %.0f .. %.0f %s (N=%d)' % (
            app.xSym(), lo, hi, app.xUnitAx(), n),
            'Average, %.0f .. %.0f %s (N=%d)' % (lo, hi, app.xUnitAx(), n))
        ax.set_yscale(ysc)
        if ysc != 'log':
            sciTicks(ax, 'y')
        self.cvP.draw_idle()

        ax = self.axProf
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(wl, Pmap[:, ifq], '-', color=BLUE, linewidth=1.4)
        ax.axvline(wlAct, **dash)
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        ax.set_ylabel(P_LABEL)
        setTitle(ax, 'Spectral profile at ν = %.1f %s' % (
            fqAct, r['info']['freqUnit']))
        ax.tick_params(labelbottom=False)
        sciTicks(ax, 'y')

        ax = self.axProfPh
        ax.clear()
        styleAxes(ax, grid=True)
        Phi = r['info'].get('Phimap')
        if Phi is not None and np.size(Phi):
            ph = np.array(Phi[:, ifq], float)
            pk = nanmaxOr(Pmap[:, ifq])
            if np.isfinite(pk) and pk > 0:
                ph[Pmap[:, ifq] < 0.05 * pk] = np.nan
            ax.plot(wl, ph, '-', color=PURPLE, linewidth=1.4)
            ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
            ax.axvline(wlAct, **dash)
            ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
            ax.set_ylim(-np.pi, np.pi)
            ax.set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
            ax.set_yticklabels(['-π', '-π/2', '0', 'π/2', 'π'])
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('Phase (rad)')
        self.cvProf.draw_idle()

        ax = self.axRes
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(r['tWin'], r['RinWin'][iw, :], 'k-', linewidth=1.0)
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel('Residual')
        setTitle(ax, 'Residual at %s = %.2f %s' % (
            app.xSym(), wlAct, app.xUnitAx()))
        sciTicks(ax, 'y')
        self.cvRes.draw_idle()

    def doReset(self):
        app = self.app
        self._quiet = True
        try:
            self.edTMin.setValue(self.defaultTMin())
            self.edTMax.setValue(self.defaultTMax())
            setCombo(self.ddResid, 'ga')
            self.edNexp.setValue(2)
            self.syncNexpEnable()
            setCombo(self.ddApod, 'hann')
            self.edZpad.setValue(2)
            setCombo(self.ddDet, 'linear')
            setCombo(self.ddFreqU, 'cm-1')
            setCombo(self.ddNorm, 'peak')
            setCombo(self.ddYScale, 'linear')
            self.edAvgWLmin.setValue(self.wlLo0)
            self.edAvgWLmax.setValue(self.wlHi0)
            self.edDispVmin.setValue(0.0)
            self.edDispVmax.setValue(0.0)
        finally:
            self._quiet = False
        self.res = emptyCoh()
        self.curFreq = np.nan
        self._setExport(False)
        app.tfSaved = None
        w = app.childWindows.get('tf')
        if w is not None:
            try:
                w.close()
            except Exception:
                pass
        self._emptyPlots()
        self.lblStatus.setText('Reset.')
        self.txtResult.setPlainText('\n'.join(IDLE_TEXT))
        app.cohSaved = None

    # ---- exports -----------------------------------------------------
    def exportMap(self):
        app, r = self.app, self.res
        if r['Pmap'] is None:
            return
        full = app.putFile('CSV (*.csv);;DAT (tab) (*.dat)',
                           'Export |FFT|² map', 'coherence_map.csv',
                           parent=self)
        if not full:
            return
        dlm = ',' if full.lower().endswith('.csv') else '\t'
        try:
            rows = np.column_stack([r['freq'], r['Pmap'].T])
            writeMatrix(full, 'wavenumber_%s' % r['info']['freqUnit'],
                        '%.6g', np.ravel(app.wavelength), rows, dlm)
            self.inform('Saved: %s' % full, 'Exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Error')

    def exportSpec(self):
        app, r = self.app, self.res
        if r['Pmap'] is None:
            return
        wl = np.asarray(app.wavelength, float).ravel()
        iw = int(np.argmin(np.abs(wl - self.curLambda)))
        exportTable(self, [r['freq'], r['Pmap'][iw, :]],
                    ['nu_%s' % r['info']['freqUnit'], 'power'],
                    'coherence_Pnu_%.0f%s' % (wl[iw], app.xUnitTxt()),
                    'Export power spectrum')

    def exportAvg(self):
        app, r = self.app, self.res
        if r['Pmap'] is None:
            return
        Pavg, lo, hi, inR = self.avgSpectrum()
        if not np.any(inR):
            self.alert('No wavelengths fall inside the current averaging '
                       'range.', 'Empty range')
            return
        exportTable(self, [r['freq'], Pavg],
                    ['nu_%s' % r['info']['freqUnit'], 'power_avg'],
                    'coherence_avgPnu_%.0f-%.0f%s' % (lo, hi,
                                                      app.xUnitTxt()),
                    'Export averaged power spectrum')

    def exportProf(self):
        app, r = self.app, self.res
        if r['Pmap'] is None:
            return
        freq = r['freq']
        ifq = int(np.argmin(np.abs(freq - self.curFreq)))
        cols = [app.wavelength, r['Pmap'][:, ifq]]
        hdr = [app.xHdr(), 'power']
        Phi = r['info'].get('Phimap')
        if Phi is not None and np.size(Phi):
            cols.append(Phi[:, ifq])
            hdr.append('phase_rad')
        exportTable(self, cols, hdr, 'coherence_profile_%.0f%s' % (
            freq[ifq], r['info']['freqUnit']), 'Export spectral profile')


# ======================================================================
def emptyTf():
    return {'Ptf': None, 'freq': None, 'tc': None, 'info': {},
            'trace': None, 'srcLabel': '', 'srcTxt': ''}


def openTimeFrequency(app):
    """Open the Time-Frequency window on the stored coherence residual
    (the checks of openTimeFrequencyWindow)."""
    sv = app.cohSaved
    if not isinstance(sv, dict) or not isinstance(sv.get('result'), dict) \
            or sv['result'].get('Rresid') is None or \
            np.size(sv['result']['Rresid']) == 0:
        app.alert('Run the Coherence analysis first -- its residual is the '
                  'input of the time-frequency analysis.', 'No residual')
        return None
    R = np.asarray(sv['result']['Rresid'], float)
    if R.ndim != 2 or R.shape[0] != app.wavelength.size or \
            R.shape[1] != app.delay.size:
        app.alert('The stored coherence residual no longer matches the '
                  'current data (crop / reload). Re-run the Coherence '
                  'analysis.', 'Stale residual')
        return None
    return openChild(app, TimeFrequencyWindow)


class TimeFrequencyWindow(ChildWindow):
    KEY = 'tf'

    def __init__(self, app):
        super().__init__(app, 'Time-Frequency Analysis (STFT / Wavelet)',
                         1380, 880)
        u = app.tUnitTxt()
        sv = app.cohSaved
        self.R = np.asarray(sv['result']['Rresid'], float)
        cs = self.cs = dict(sv['settings'])
        self.freqUnit = cs['frequ']
        self.res = emptyTf()
        self.curT = np.nan
        self.curNu = np.nan
        self.xhV = self.xhH = None
        self.dragging = False
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: setup ==================================
        left, lv = sidePanel('Time-frequency setup', 345)
        lv.setSpacing(5)
        gl.addWidget(left)

        def bold(text):
            lb = label(text)
            lb.setStyleSheet('font-weight: bold;')
            lv.addWidget(lb)

        def pair(text, field):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(6)
            lb = label(text)
            r.addWidget(lb)
            r.addWidget(field, 1)
            lv.addLayout(r)
            return lb

        def rng(t1, f1, t2, f2):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(4)
            r.addWidget(label(t1))
            r.addWidget(f1, 1)
            r.addWidget(label(t2))
            r.addWidget(f2, 1)
            lv.addLayout(r)

        wl = np.asarray(app.wavelength, float)
        wlLo0, wlHi0 = float(np.min(wl)), float(np.max(wl))
        lim = (wlLo0, max(wlHi0, wlLo0))
        bold('Input trace (coherence residual):')
        self.ddSource = QtWidgets.QComboBox()
        self.ddSource.addItem(ui_str('Residual at \\lambda'), 'single')
        self.ddSource.addItem(ui_str('\\lambda-range average'), 'avg')
        pair('Source:', self.ddSource)
        self.edLambda = NumEdit(
            self._lam0(), limits=lim,
            tooltip='Wavelength whose residual trace is analyzed.')
        pair('%s (%s):' % (app.xSymTxt(), app.xUnitTxt()), self.edLambda)
        self.edAvgLo = NumEdit(cs['avgLo'], limits=lim)
        self.edAvgHi = NumEdit(cs['avgHi'], limits=lim)
        rng(ui_str('Avg \\lambda:'), self.edAvgLo, '~', self.edAvgHi)
        bold('Analysis window:')
        self.edTMin = NumEdit(cs['tmin'])
        self.edTMax = NumEdit(cs['tmax'])
        rng('t_min (%s):' % u, self.edTMin, 't_max (%s):' % u, self.edTMax)
        self.ddDet = QtWidgets.QComboBox()
        for nm, key in (('None', 'none'), ('Mean', 'mean'),
                        ('Linear', 'linear')):
            self.ddDet.addItem(nm, key)
        setCombo(self.ddDet, cs['detrend'])
        pair('Detrend:', self.ddDet)
        bold('Method:')
        self.ddMethod = QtWidgets.QComboBox()
        self.ddMethod.addItem('STFT (sliding-window FFT)', 'stft')
        self.ddMethod.addItem('Wavelet (Morlet CWT)', 'cwt')
        self.ddMethod.setToolTip(tip(
            'STFT: fixed window -> uniform time / frequency resolution.\n'
            'Wavelet: window scales with 1/nu -> frequency resolution at '
            'low nu, time resolution at high nu.'))
        pair('Transform:', self.ddMethod)
        self.winDef = (cs['tmax'] - cs['tmin']) / 4.0
        self.edWin = NumEdit(
            self.winDef,
            tooltip='Sliding-window length. Frequency resolution ~ 1 / '
                    'window; time resolution ~ window. Longer = sharper '
                    'peaks, blurrier in time.')
        self.lblWin = pair('STFT window (%s):' % u, self.edWin)
        self.edHop = NumEdit(
            0.0, tooltip='Step between window centers. 0 = every sample.')
        self.lblHop = pair('STFT hop (%s):' % u, self.edHop)
        self.ddApod = QtWidgets.QComboBox()
        for nm, key in (('Hann', 'hann'), ('Hamming', 'hamming'),
                        ('Blackman', 'blackman'), ('Rectangular', 'rect')):
            self.ddApod.addItem(nm, key)
        setCombo(self.ddApod, cs['apod'])
        self.lblApod = pair('STFT apodization:', self.ddApod)
        self.edZpad = NumEdit(cs['zpad'], limits=(1, 32), integer=True)
        self.lblZpad = pair('STFT zero-pad factor:', self.edZpad)
        self.edOm0 = NumEdit(
            6, limits=(2, 40),
            tooltip='Morlet time-frequency trade-off. Higher = sharper '
                    'frequency, blurrier time (6 is the usual default; 4-5 '
                    'for fast dephasing, 10-12 for close modes).')
        self.lblOm0 = pair(ui_str('Wavelet \\omega_0:'), self.edOm0)
        self.edNf = NumEdit(200, limits=(16, 2000), integer=True)
        self.lblNf = pair(ui_str('Wavelet \\nu grid points:'), self.edNf)
        bold('Display (%s):' % freqUnitTxt(self.freqUnit))
        self.edDispVmin = NumEdit(cs['dispLo'])
        self.edDispVmax = NumEdit(
            cs['dispHi'],
            tooltip='Visualization window on the frequency axis; also '
                    'bounds the wavelet frequency grid at Run. 0 = auto '
                    '(up to Nyquist).')
        rng(ui_str('\\nu range:'), self.edDispVmin, '~', self.edDispVmax)
        self.ddZScale = QtWidgets.QComboBox()
        self.ddZScale.addItem('Linear', 'linear')
        self.ddZScale.addItem('Log', 'log')
        pair('Map z / 1D y scale:', self.ddZScale)
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        self.btnRun = colorButton('Run', RUN_COLOR_TF)
        self.btnRun.setStyleSheet(self.btnRun.styleSheet()
                                  + 'QPushButton { font-weight: bold; }')
        self.btnReset = button('Reset')
        r.addWidget(self.btnRun, 1)
        r.addWidget(self.btnReset, 1)
        lv.addLayout(r)
        bold('Export:')
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        self.btnExpMap = button(ui_str('\\nu-t map'))
        self.btnExpDec = button(
            ui_str('P(t)@\\nu'),
            tooltip='Time decay of the power at the selected frequency.')
        self.btnExpSpc = button(
            ui_str('P(\\nu)@t'),
            tooltip='Instantaneous spectrum at the selected time.')
        for b in (self.btnExpMap, self.btnExpDec, self.btnExpSpc):
            b.setEnabled(False)
            r.addWidget(b, 1)
        lv.addLayout(r)
        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        lv.addWidget(self.lblStatus)
        self.txtResult = textArea(TF_IDLE)
        self.txtResult.setMinimumHeight(100)
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
        self.cvDecay = PlotCanvas()
        self.cvSpecT = PlotCanvas()
        self.cvTrace = PlotCanvas()
        self.cvMap.setMinimumHeight(250)
        self.cvDecay.setMinimumHeight(250)
        self.cvSpecT.setMinimumHeight(180)
        self.cvTrace.setMinimumHeight(180)
        grid.addWidget(self.cvMap, 0, 0)
        grid.addWidget(self.cvDecay, 0, 1)
        grid.addWidget(self.cvSpecT, 1, 0)
        grid.addWidget(self.cvTrace, 1, 1)
        grid.setRowStretch(0, 3)
        grid.setRowStretch(1, 2)
        grid.setColumnStretch(0, 3)
        grid.setColumnStretch(1, 2)
        self.axMap = self.cvMap.figure.add_subplot(111)
        self.caxMap = self.cvMap.figure.add_axes([0.9, 0.1, 0.02, 0.8])
        self.axDecay = self.cvDecay.figure.add_subplot(111)
        self.axSpecT = self.cvSpecT.figure.add_subplot(111)
        self.axTrace = self.cvTrace.figure.add_subplot(111)
        self.cvMap.layoutFcn = self._layMap
        self.cvDecay.layoutFcn = lambda: lay1(self.cvDecay, self.axDecay)
        self.cvSpecT.layoutFcn = lambda: lay1(self.cvSpecT, self.axSpecT)
        self.cvTrace.layoutFcn = lambda: lay1(self.cvTrace, self.axTrace)
        self.blMap = BlitManager(self.cvMap)
        self.cvMap.mpl_connect('button_press_event', self.onMapPress)
        self.cvMap.mpl_connect('motion_notify_event', self.onMapMotion)
        self.cvMap.mpl_connect('button_release_event',
                               lambda ev: self.onMapUp())
        self._emptyPlots()

        # ---- what the controls do --------------------------------------
        self._quiet = False
        self.btnRun.clicked.connect(lambda *_: self.safe(self.doRun, 'Run'))
        self.btnReset.clicked.connect(
            lambda *_: self.safe(self.doReset, 'Reset'))
        self.btnExpMap.clicked.connect(
            lambda *_: self.safe(self.exportMap, 'Export map'))
        self.btnExpDec.clicked.connect(
            lambda *_: self.safe(self.exportDecay, 'Export decay'))
        self.btnExpSpc.clicked.connect(
            lambda *_: self.safe(self.exportSpecT, 'Export spectrum'))
        self.ddSource.currentIndexChanged.connect(
            lambda *_: self.onSourceChange())
        self.ddMethod.currentIndexChanged.connect(
            lambda *_: self.onMethodChange())
        self.edDispVmin.valueChanged.connect(
            lambda *_: self._user(self.replot))
        self.edDispVmax.valueChanged.connect(
            lambda *_: self._user(self.replot))
        self.ddZScale.currentIndexChanged.connect(
            lambda *_: self._user(self.replot))
        self.onSourceChange()
        self.onMethodChange()

        # ---- result of an earlier visit --------------------------------
        ts = app.tfSaved
        if isinstance(ts, dict) and isinstance(ts.get('result'), dict) and \
                ts['result'].get('Ptf') is not None and \
                np.size(ts['result'].get('trace')) == app.delay.size:
            self._quiet = True
            try:
                st = ts['settings']
                setCombo(self.ddSource, st['source'])
                self.edLambda.setValue(st['lambda'])
                self.edAvgLo.setValue(st['avgLo'])
                self.edAvgHi.setValue(st['avgHi'])
                self.edTMin.setValue(st['tmin'])
                self.edTMax.setValue(st['tmax'])
                setCombo(self.ddDet, st['detrend'])
                setCombo(self.ddMethod, st['method'])
                self.edWin.setValue(st['winLen'])
                self.edHop.setValue(st['hop'])
                setCombo(self.ddApod, st['apod'])
                self.edZpad.setValue(st['zpad'])
                self.edOm0.setValue(st['omega0'])
                self.edNf.setValue(st['nFreq'])
                self.edDispVmin.setValue(st['dispLo'])
                self.edDispVmax.setValue(st['dispHi'])
                setCombo(self.ddZScale, st['zscale'])
                r0 = ts['result']
                info = {k: (np.asarray(v, float) if isinstance(
                    v, (list, tuple, np.ndarray)) else v)
                    for k, v in dict(r0['info']).items()}
                self.res = {
                    'Ptf': np.asarray(r0['Ptf'], float),
                    'freq': np.asarray(r0['freq'], float).ravel(),
                    'tc': np.asarray(r0['tc'], float).ravel(),
                    'info': info,
                    'trace': np.asarray(r0['trace'], float).ravel(),
                    'srcLabel': r0['srcLabel'],
                    'srcTxt': r0.get('srcTxt', r0['srcLabel'])}
                self.curT = float(ts['curT'])
                self.curNu = float(ts['curNu'])
                self.onSourceChange()
                self.onMethodChange()
                self.setExports(True)
                self.plotResults()
                self.lblStatus.setText('%s  [restored]' % ts['status'])
                self.txtResult.setPlainText('\n'.join(ts['lines']))
            except Exception:
                self.res = emptyTf()
                self.setExports(False)
                self._emptyPlots()
                self.lblStatus.setText(
                    'Stored result could not be restored -- run again.')
            self._quiet = False

    # ------------------------------------------------------------------
    def _user(self, fcn):
        if not self._quiet:
            self.safe(fcn, 'Display')

    def _lam0(self):
        app = self.app
        lam = app.cohSaved.get('curLambda') if isinstance(
            app.cohSaved, dict) else None
        if lam is None or not np.isfinite(lam):
            n = app.wavelength.size
            lam = app.wavelength[max(int(np.floor(n / 2 + 0.5)) - 1, 0)]
        return float(lam)

    def _layMap(self):
        mapLayout(self.cvMap, self.axMap, self.caxMap, L=58.0, R=74.0,
                  T=22.0, B=34.0)

    def _emptyPlots(self):
        for ax in (self.axMap, self.axDecay, self.axSpecT, self.axTrace):
            ax.clear()
            styleAxes(ax)
        self.caxMap.clear()
        self.caxMap.set_visible(False)
        setTitle(self.axMap, TF_MAP)
        setTitle(self.axDecay, TF_DECAY, 'Power vs time')
        setTitle(self.axSpecT, TF_SPEC)
        setTitle(self.axTrace, TF_TRACE, 'Input trace')
        self.xhV = self.xhH = None
        self.blMap.setArtists([])
        self._layMap()
        for cv in (self.cvMap, self.cvDecay, self.cvSpecT, self.cvTrace):
            cv.draw_idle()

    def onSourceChange(self):
        single = self.ddSource.currentData() == 'single'
        self.edLambda.setEnabled(single)
        self.edAvgLo.setEnabled(not single)
        self.edAvgHi.setEnabled(not single)

    def onMethodChange(self):
        stft = self.ddMethod.currentData() == 'stft'

        def col(lb, on):
            c = (0.10, 0.10, 0.10) if on else (0.60, 0.60, 0.60)
            lb.setStyleSheet('color: rgb(%d,%d,%d);' % tuple(
                int(round(255 * x)) for x in c))

        for w in (self.edWin, self.edHop, self.ddApod, self.edZpad):
            w.setEnabled(stft)
        for lb in (self.lblWin, self.lblHop, self.lblApod, self.lblZpad):
            col(lb, stft)
        for w in (self.edOm0, self.edNf):
            w.setEnabled(not stft)
        for lb in (self.lblOm0, self.lblNf):
            col(lb, not stft)

    def setExports(self, on):
        for b in (self.btnExpMap, self.btnExpDec, self.btnExpSpc):
            b.setEnabled(bool(on))

    def getDispLims(self):
        """([lo, hi], ok) of the frequency axis shown."""
        fr = self.res['freq']
        if fr is None or fr.size == 0:
            return None, False
        fLo, fHi = float(np.min(fr)), float(np.max(fr))
        lo, hi = self.edDispVmin.value(), self.edDispVmax.value()
        if not np.isfinite(hi) or hi <= 0:
            hi = fHi
        if not np.isfinite(lo) or lo < 0:
            lo = fLo
        if lo > hi:
            lo, hi = hi, lo
        lo, hi = max(lo, fLo), min(hi, fHi)
        if lo >= hi:
            return None, False
        return [lo, hi], True

    def _lims(self):
        lims, ok = self.getDispLims()
        if not ok:
            fr = self.res['freq']
            lims = [float(np.min(fr)), float(np.max(fr))]
        return lims

    # ------------------------------------------------------------------
    def doRun(self):
        app = self.app
        wl = np.asarray(app.wavelength, float).ravel()
        if self.ddSource.currentData() == 'single':
            iw = int(np.argmin(np.abs(wl - self.edLambda.value())))
            y = self.R[iw, :]
            self.edLambda.setValue(float(wl[iw]))
            srcLabel = '%s = %.1f %s' % (app.xSym(), wl[iw], app.xUnitAx())
            srcTxt = '%s = %.1f %s' % (app.xSym(), wl[iw], app.xUnitTxt())
        else:
            lo = min(self.edAvgLo.value(), self.edAvgHi.value())
            hi = max(self.edAvgLo.value(), self.edAvgHi.value())
            m = (wl >= lo) & (wl <= hi)
            if not np.any(m):
                self.alert('No wavelengths inside the averaging range.',
                           'Empty range')
                return
            with np.errstate(invalid='ignore'):
                y = np.nanmean(self.R[m, :], axis=0)
            n = int(np.count_nonzero(m))
            srcLabel = 'avg %s = %.0f..%.0f %s (N=%d)' % (
                app.xSym(), lo, hi, app.xUnitAx(), n)
            srcTxt = 'avg %s = %.0f..%.0f %s (N=%d)' % (
                app.xSym(), lo, hi, app.xUnitTxt(), n)
        tMin, tMax = self.edTMin.value(), self.edTMax.value()
        if tMax <= tMin:
            self.alert('t_max must exceed t_min.', 'Invalid window')
            return
        det, meth = self.ddDet.currentData(), self.ddMethod.currentData()
        fu = self.freqUnit
        winLen, hop = self.edWin.value(), self.edHop.value()
        apod, zpad = self.ddApod.currentData(), self.edZpad.value()
        om0, nF = self.edOm0.value(), self.edNf.value()
        fLo, fHi = self.edDispVmin.value(), self.edDispVmax.value()
        tU = app.timeUnit
        self.lblStatus.setText('Computing...')
        t0 = time.time()
        try:
            o = runBusy(app, self, 'Computing...',
                        lambda: compute_time_frequency(
                            y, app.delay, tMin, tMax, det, meth, fu, 'peak',
                            winLen, hop, apod, zpad, om0, nF, fLo, fHi, tU))
        except Exception as err:
            self.alert('Time-frequency analysis failed:\n%s' % err, 'Error')
            self.lblStatus.setText('Failed.')
            return
        elapsed = time.time() - t0
        Ptf = np.asarray(o.ptf, float)
        freq = np.asarray(o.freq, float).ravel()
        tc = np.asarray(o.tc, float).ravel()
        info = {'method': o.method, 'dtUni': float(o.dt),
                'tUni': np.asarray(o.t_uni, float),
                'yu': np.asarray(o.yu, float),
                'coiFreq': np.asarray(o.coi_freq, float),
                'coiT': np.asarray(o.coi_t, float)}
        if o.method == 'stft':
            info.update(nWin=int(o.n_win), winLen=float(o.win_len),
                        hopN=int(o.hop_n), Nfft=int(o.n_fft), apod=apod,
                        zeroPad=int(o.zero_pad))
        else:
            info.update(omega0=float(o.omega0), nFreq=int(o.n_freq),
                        scales=np.asarray(o.scales, float))
        self.res = {'Ptf': Ptf, 'freq': freq, 'tc': tc, 'info': info,
                    'trace': np.asarray(y, float), 'srcLabel': srcLabel,
                    'srcTxt': srcTxt}
        # the strongest point outside the edge-affected part
        Pm = Ptf.copy()
        Pm[freq[:, None] < info['coiFreq'][None, :]] = -np.inf
        k = int(np.argmax(Pm))
        if not np.isfinite(Pm.flat[k]):
            k = int(np.argmax(Ptf))
        iF, iT = np.unravel_index(k, Ptf.shape)
        self.curNu, self.curT = float(freq[iF]), float(tc[iT])
        self.setExports(True)
        self.plotResults()
        L = ['Source  : %s' % srcTxt,
             'Window  : %.4g .. %.4g %s' % (tMin, tMax, app.tUnitTxt()),
             'dt      : %.4g %s (%d pts)' % (info['dtUni'], app.tUnitTxt(),
                                             info['tUni'].size),
             'Detrend : %s' % det]
        if o.method == 'stft':
            L += ['Method  : STFT',
                  '  window = %.4g %s (%d pts, %s)' % (
                      info['winLen'], app.tUnitTxt(), info['nWin'],
                      info['apod']),
                  '  hop    = %d pt, zero-pad x%d' % (info['hopN'],
                                                      info['zeroPad']),
                  '  d(nu)  = %.3g %s (FFT bin)' % (freq[1] - freq[0], fu),
                  '  edge-affected: t < %.4g, t > %.4g' % (
                      info['coiT'][0], info['coiT'][1])]
        else:
            L += ['Method  : Morlet CWT',
                  '  omega0 = %.3g, %d freq pts' % (info['omega0'],
                                                    info['nFreq']),
                  '  nu grid: %.4g .. %.4g %s' % (float(np.min(freq)),
                                                  float(np.max(freq)), fu),
                  '  Dashed line on map = cone of influence',
                  '  (below it: edge-affected, unreliable).']
        L += ['', 'Click + drag on the map to move the', 'cursor (nu, t).']
        self.txtResult.setPlainText('\n'.join(L))
        self.lblStatus.setText('Done (%d x %d, %.2fs).' % (
            Ptf.shape[0], Ptf.shape[1], elapsed))
        app.tfSaved = {
            'result': dict(self.res), 'curT': self.curT,
            'curNu': self.curNu, 'status': self.lblStatus.text(),
            'lines': list(L),
            'settings': {
                'source': self.ddSource.currentData(),
                'lambda': self.edLambda.value(),
                'avgLo': self.edAvgLo.value(),
                'avgHi': self.edAvgHi.value(), 'tmin': float(tMin),
                'tmax': float(tMax), 'detrend': det, 'method': meth,
                'winLen': float(winLen), 'hop': float(hop), 'apod': apod,
                'zpad': int(zpad), 'omega0': float(om0), 'nFreq': int(nF),
                'dispLo': float(fLo), 'dispHi': float(fHi),
                'zscale': self.ddZScale.currentData()}}

    def doReset(self):
        cs = self.cs
        self._quiet = True
        try:
            setCombo(self.ddSource, 'single')
            self.edLambda.setValue(self._lam0())
            self.edAvgLo.setValue(cs['avgLo'])
            self.edAvgHi.setValue(cs['avgHi'])
            self.edTMin.setValue(cs['tmin'])
            self.edTMax.setValue(cs['tmax'])
            setCombo(self.ddDet, cs['detrend'])
            setCombo(self.ddMethod, 'stft')
            self.edWin.setValue(self.winDef)
            self.edHop.setValue(0.0)
            setCombo(self.ddApod, cs['apod'])
            self.edZpad.setValue(cs['zpad'])
            self.edOm0.setValue(6)
            self.edNf.setValue(200)
            self.edDispVmin.setValue(cs['dispLo'])
            self.edDispVmax.setValue(cs['dispHi'])
            setCombo(self.ddZScale, 'linear')
        finally:
            self._quiet = False
        self.onSourceChange()
        self.onMethodChange()
        self.res = emptyTf()
        self.curT = self.curNu = np.nan
        self.setExports(False)
        self._emptyPlots()
        self.lblStatus.setText('Reset.')
        self.txtResult.setPlainText('\n'.join(TF_IDLE))
        self.app.tfSaved = None

    def replot(self):
        if self.res['Ptf'] is not None:
            self.plotResults()

    def plotResults(self):
        app = self.app
        r = self.res
        if r['Ptf'] is None:
            return
        P, freq, tc, info = r['Ptf'], r['freq'], r['tc'], r['info']
        lims = self._lims()
        fm = (freq >= lims[0]) & (freq <= lims[1])
        if not np.any(fm):
            fm = np.ones(freq.size, bool)
        Psub = P[fm, :]
        if self.ddZScale.currentData() == 'log':
            pk = nanmaxOr(Psub)
            if not np.isfinite(pk) or pk <= 0:
                pk = 1.0
            floorVal = pk * 1e-4
            Cdat = np.log10(np.maximum(P, floorVal))
            cLo, cHi = np.log10(floorVal), np.log10(pk)
            cbLabel = 'log$_{10}$ power'
        else:
            Cdat = P
            cHi = nanmaxOr(Psub)
            if not np.isfinite(cHi) or cHi <= 0:
                cHi = 1.0
            cLo, cbLabel = 0.0, 'Power (norm.)'
        ax, cax = self.axMap, self.caxMap
        ax.clear()
        cax.clear()
        styleAxes(ax)
        ax.set_facecolor(AX_BG)
        h = None
        if tc.size > 1 and freq.size > 1:
            h = ax.pcolormesh(tc, freq, np.ma.masked_invalid(Cdat[:-1, :-1]),
                              shading='flat', cmap='turbo', vmin=cLo,
                              vmax=cHi, rasterized=True)
        if info['method'] == 'cwt':
            cf = np.array(info['coiFreq'], float)
            cf[cf > lims[1]] = np.nan
            ax.plot(tc, cf, '--', color='w', linewidth=1.2)
        else:
            for e in np.ravel(info['coiT']):
                ax.axvline(e, color='w', linestyle='--', linewidth=1.2)
        self.xhV = ax.axvline(self.curT, color='k', linewidth=1.0)
        self.xhH = ax.axhline(self.curNu, color='k', linewidth=1.0)
        self.blMap.setArtists([self.xhV, self.xhH])
        ax.set_xlim(float(np.min(tc)), float(np.max(tc)))
        ax.set_ylim(lims[0], lims[1])
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel(freqAxLabel(self.freqUnit))
        if info['method'] == 'stft':
            head = 'STFT  (window %.3g %s)' % (info['winLen'], app.tUnitAx())
            short = 'STFT'
        else:
            head = 'Morlet CWT  (ω$_0$ = %.3g)' % info['omega0']
            short = 'Morlet CWT'
        setTitle(ax, '%s  --  %s' % (head, r['srcLabel']), head, short)
        cax.set_visible(h is not None)
        if h is not None:
            cb = self.cvMap.figure.colorbar(h, cax=cax)
            cb.set_label(cbLabel, fontsize=FONT_PT)
            cb.ax.tick_params(labelsize=FONT_PT, length=3)
        self._layMap()
        self.cvMap.draw_idle()
        self.update1Ds()

    def update1Ds(self):
        app = self.app
        r = self.res
        if r['Ptf'] is None:
            return
        P, freq, tc, info = r['Ptf'], r['freq'], r['tc'], r['info']
        useLog = self.ddZScale.currentData() == 'log'
        lims = self._lims()
        coi = np.asarray(info['coiFreq'], float)
        iF = int(np.argmin(np.abs(freq - self.curNu)))
        rowP = P[iF, :]
        reliable = coi <= freq[iF]
        rOK = np.where(reliable, rowP, np.nan)
        rNO = np.where(reliable, np.nan, rowP)
        dash = dict(color='k', linestyle='--', linewidth=0.8)
        ax = self.axDecay
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(tc, rNO, '-', color=(0.65, 0.65, 0.65), linewidth=1.2)
        ax.plot(tc, rOK, '-', color=(0.10, 0.30, 0.80), linewidth=1.6)
        ax.axvline(self.curT, **dash)
        ax.set_xlim(float(np.min(tc)), float(np.max(tc)))
        ax.set_yscale('log' if useLog else 'linear')
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel('Power')
        fa = 'P(t) at ν = %.1f %s' % (freq[iF], freqUnitAx(self.freqUnit))
        setTitle(ax, fa + '   (gray = edge-affected)', fa)
        self.cvDecay.draw_idle()

        iT = int(np.argmin(np.abs(tc - self.curT)))
        ax = self.axSpecT
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(freq, P[:, iT], '-', color=(0.10, 0.55, 0.35), linewidth=1.6)
        ax.axvline(self.curNu, **dash)
        cfT = coi[iT]
        if np.isfinite(cfT) and lims[0] < cfT < lims[1]:
            ax.axvline(cfT, color=(0.5, 0.5, 0.5), linestyle='--',
                       linewidth=0.9)
            ax.annotate('COI', xy=(cfT, 1.0), xycoords=('data', 'axes fraction'),
                        xytext=(3, -3), textcoords='offset points',
                        ha='left', va='top', fontsize=FONT_PT - 1,
                        color=(0.4, 0.4, 0.4))
        ax.set_xlim(lims[0], lims[1])
        ax.set_yscale('log' if useLog else 'linear')
        ax.set_xlabel(freqAxLabel(self.freqUnit))
        ax.set_ylabel('Power')
        ta = 'P(ν) at t = %.4g %s' % (tc[iT], app.tUnitAx())
        if np.isfinite(cfT):
            setTitle(ax, ta)
        else:
            setTitle(ax, ta + '  (edge-affected)', ta)
        self.cvSpecT.draw_idle()

        tU, yu = np.asarray(info['tUni'], float), np.asarray(info['yu'],
                                                             float)
        ax = self.axTrace
        ax.clear()
        styleAxes(ax, grid=True)
        yr = [float(np.min(yu)), float(np.max(yu))]
        if not (yr[1] > yr[0]):
            yr = [yr[0] - 1, yr[1] + 1]
        pad = 0.08 * (yr[1] - yr[0])
        yr = [yr[0] - pad, yr[1] + pad]
        if info['method'] == 'stft':
            halfSpan = info['winLen'] / 2.0
        else:
            halfSpan = np.sqrt(2) * float(np.interp(
                self.curNu, freq, np.asarray(info['scales'], float)))
        tA = max(tU[0], self.curT - halfSpan)
        tB = min(tU[-1], self.curT + halfSpan)
        if tB > tA:
            ax.fill([tA, tB, tB, tA], [yr[0], yr[0], yr[1], yr[1]],
                    facecolor=(0.95, 0.75, 0.30), alpha=0.25,
                    edgecolor='none')
        ax.plot(tU, yu, 'k-', linewidth=1.0)
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.axvline(self.curT, **dash)
        ax.set_xlim(float(tU[0]), float(tU[-1]))
        ax.set_ylim(yr[0], yr[1])
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel('Residual')
        if info['method'] == 'stft':
            setTitle(ax, 'Input trace  (shaded = STFT window at cursor)',
                     'Input trace  (shaded = STFT window)', 'Input trace')
        else:
            setTitle(ax, 'Input trace  (shaded = wavelet e-folding width at '
                         'cursor ν)',
                     'Input trace  (shaded = wavelet width)', 'Input trace')
        sciTicks(ax, 'y')
        self.cvTrace.draw_idle()

    # ---- the cursor on the map ---------------------------------------
    def onMapPress(self, ev):
        if ev.button != 1 or ev.inaxes is not self.axMap or \
                self.res['Ptf'] is None:
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
        r = self.res
        if r['Ptf'] is None or x is None or y is None or \
                not np.isfinite(x) or not np.isfinite(y):
            return
        tc, fr = r['tc'], r['freq']
        self.curT = float(min(max(x, np.min(tc)), np.max(tc)))
        self.curNu = float(min(max(y, np.min(fr)), np.max(fr)))
        if self.xhV is not None:
            self.xhV.set_xdata([self.curT, self.curT])
        if self.xhH is not None:
            self.xhH.set_ydata([self.curNu, self.curNu])
        self.blMap.update()
        self.update1Ds()
        if isinstance(self.app.tfSaved, dict):
            self.app.tfSaved['curT'] = self.curT
            self.app.tfSaved['curNu'] = self.curNu

    # ---- exports -----------------------------------------------------
    def exportMap(self):
        r = self.res
        if r['Ptf'] is None:
            return
        full = self.app.putFile('CSV (*.csv)', 'Export nu-t map',
                                'timefreq_map.csv', parent=self)
        if not full:
            return
        try:
            rows = np.column_stack([r['freq'], r['Ptf']])
            writeMatrix(full, 'nu_%s' % self.freqUnit, '%.10g', r['tc'],
                        rows, ',')
            self.inform('Saved: %s' % full, 'Exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Error')

    def exportDecay(self):
        r = self.res
        if r['Ptf'] is None:
            return
        iF = int(np.argmin(np.abs(r['freq'] - self.curNu)))
        rel = (np.asarray(r['info']['coiFreq'], float)
               <= r['freq'][iF]).astype(float)
        exportTable(self, [r['tc'], r['Ptf'][iF, :], rel],
                    ['t_%s' % self.app.tUnitTxt(), 'power', 'reliable'],
                    'timefreq_Pt_%.0f' % r['freq'][iF],
                    'Export decay at nu')

    def exportSpecT(self):
        r = self.res
        if r['Ptf'] is None:
            return
        iT = int(np.argmin(np.abs(r['tc'] - self.curT)))
        exportTable(self, [r['freq'], r['Ptf'][:, iT]],
                    ['nu_%s' % self.freqUnit, 'power'],
                    'timefreq_Pnu_t%.4g' % r['tc'][iT],
                    'Export spectrum at t')
