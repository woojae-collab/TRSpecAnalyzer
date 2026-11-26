"""Band position window: where a band lies at every delay (spectral
shift, dynamic Stokes shift), and the fit of that time course. Port of
openBandPositionWindow."""
import numpy as np
from PySide6 import QtWidgets

from trspec._util import Result
from trspec.bandpos import (band_model_curve, band_parabola,
                            fit_band_shift)

from .state import ui_str
from .textfit import legend
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label, panel, sidePanel, styleAxes, textArea, tip)
from .win_base import ChildWindow
from .win_kinfit import writeTable

RUN_COLOR = (0.45, 0.40, 0.70)
C_DATA = (0.10, 0.35, 0.75)
C_FIT = (0.85, 0.20, 0.15)
C_POS = (0.20, 0.50, 0.30)
IDLE_TEXT = ['Press "Run" to follow the band.']
T_SPEC = 'Band at the delay chosen'
T_POS = 'Position of the band'
T_C = 'C(t)'
T_W = 'Width of the band (FWHM)'
METHOD_LABELS = (('moment', 'First moment'),
                 ('parabola', 'Maximum (parabola)'),
                 ('lognormal', 'Log-normal fit'))


class BandPositionWindow(ChildWindow):
    KEY = 'bandpos'

    def __init__(self, app):
        super().__init__(app, 'Band Position', 1320, 820)
        u = app.tUnitTxt()
        xu = app.xUnitTxt()
        self.T = None               # positions (Result) of the last run
        self.F = None               # fit of the positions
        self.info = None
        self.S = None               # the setup the results belong to
        sv = app.bpSaved if isinstance(app.bpSaved, dict) else None
        S0 = app.bandSetupDefault()
        if sv and isinstance(sv.get('set'), dict):
            S0.update({k: sv['set'][k] for k in S0 if k in sv['set']})
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: setup ==================================
        left, lv = sidePanel('Band position setup', 380)
        gl.addWidget(left)

        def bold(text):
            lb = label(text)
            lb.setStyleSheet('font-weight: bold;')
            lv.addWidget(lb)

        def row(*items):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(6)
            for it in items:
                if isinstance(it, str):
                    r.addWidget(label(it))
                else:
                    r.addWidget(it)
            r.addStretch(1)
            lv.addLayout(r)

        bold('Band:')
        self.edX1 = NumEdit(S0['x1'], fmt='%.6g', width=80)
        self.edX2 = NumEdit(S0['x2'], fmt='%.6g', width=80)
        self.edX1.setToolTip(tip(
            'The part of the spectral axis that holds the band at all '
            'delays. Wavelengths without data (masked) are left out.'))
        row('from', self.edX1, 'to', self.edX2, xu)
        self.ddSign = QtWidgets.QComboBox()
        self.ddSign.addItem('Positive', 'pos')
        self.ddSign.addItem('Negative', 'neg')
        self.ddSign.setCurrentIndex(max(self.ddSign.findData(S0['sign']), 0))
        self.ddSign.setToolTip(tip(
            'A negative band (bleach, stimulated emission in dA) is turned '
            'over before its position is taken.'))
        row('The band is:', self.ddSign)
        bold('Position:')
        self.ddUnit = QtWidgets.QComboBox()
        self.ddUnit.addItem('Wavenumber (cm-1)', 'cm-1')
        self.ddUnit.addItem('Wavelength (nm)', 'nm')
        self.ddUnit.setCurrentIndex(max(self.ddUnit.findData(S0['posUnit']),
                                        0))
        row('given in:', self.ddUnit)
        self.cbDens = QtWidgets.QCheckBox('')
        self.cbDens.setChecked(bool(S0['density']))
        self.cbDens.setToolTip(tip(
            'For a signal that is a density along the spectral axis, as an '
            'emission intensity per nm: on the other axis the same signal '
            'is the one per unit of that axis, which is the signal times '
            'the square of the axis value / 1e7. Leave it off for dA and '
            'dT/T, which are no densities.'))
        lv.addWidget(self.cbDens)
        self.ddMethod = QtWidgets.QComboBox()
        for key, name in METHOD_LABELS:
            self.ddMethod.addItem(name, key)
        self.ddMethod.setCurrentIndex(
            max(self.ddMethod.findData(S0['method']), 0))
        self.ddMethod.setToolTip(tip(
            'First moment: the mean of the axis weighted with the signal, '
            'over the whole band; no assumption about the shape, but the '
            'band has to lie inside the range at all delays. Maximum '
            '(parabola): the top of a parabola through the highest points. '
            'Log-normal fit: the maximum of a log-normal band fitted to '
            'the whole range (height, position, width, asymmetry).'))
        row('taken as:', self.ddMethod)
        self.edFrac = NumEdit(S0['frac'], fmt='%.3g', limits=(0.05, 0.99),
                              width=60)
        self.lblFrac = label('Parabola: points from')
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(self.lblFrac)
        r.addWidget(self.edFrac)
        self.lblFrac2 = label('of the maximum')
        r.addWidget(self.lblFrac2)
        r.addStretch(1)
        lv.addLayout(r)
        self.edMin = NumEdit(S0['minPct'], fmt='%.3g', limits=(0, 100),
                             width=60)
        self.edMin.setToolTip(tip(
            'A delay gets a position only when the largest value of its '
            'band is at least this share of the largest value of the band '
            'at all delays: before the signal has risen there is no band.'))
        row('No position below', self.edMin, '% of the largest signal')

        bold('Fit of the position (sum of exponentials):')
        self.edT1 = NumEdit(S0['tFrom'], fmt='%.6g', width=80)
        self.edT2 = NumEdit(S0['tTo'], fmt='%.6g', width=80)
        self.edT1.setToolTip(tip(
            'Delays of the fit. Start after the instrument response: while '
            'the signal rises, the spectrum is a mixture of earlier and '
            'later times.'))
        row('t from', self.edT1, 'to', self.edT2, u)
        self.ddN = QtWidgets.QComboBox()
        self.ddN.addItems(['1', '2', '3', '4'])
        self.ddN.setCurrentText(str(int(S0['nExp'])))
        row('Exponentials:', self.ddN)
        self.edTau = QtWidgets.QLineEdit()
        self.edTau.setText(self.tauText(S0['tauStart']))
        self.edTau.setToolTip(tip(
            'Starting time constants, separated by spaces or commas, one '
            'per exponential. Empty: three automatic starts are tried and '
            'the best fit is kept.'))
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label(ui_str('Start \\tau (%s):' % u)))
        r.addWidget(self.edTau, 1)
        lv.addLayout(r)
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        self.btnRun = colorButton('Run', RUN_COLOR)
        self.btnRun.setStyleSheet(self.btnRun.styleSheet()
                                  + 'QPushButton { font-weight: bold; }')
        self.btnReset = button('Reset')
        self.btnExport = button('Export...')
        self.btnExport.setEnabled(False)
        for b in (self.btnRun, self.btnReset, self.btnExport):
            r.addWidget(b, 1)
        lv.addLayout(r)
        bold('Display:')
        self.edTS = NumEdit(S0['tSpec'], fmt='%.6g', width=80)
        self.cbLog = QtWidgets.QCheckBox('Log delay axis')
        self.cbLog.setChecked(bool(S0['logT']))
        row('Band at t:', self.edTS, u, self.cbLog)
        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        self.lblStatus.setWordWrap(True)
        lv.addWidget(self.lblStatus)
        self.txtResult = textArea(IDLE_TEXT)
        self.txtResult.setLineWrapMode(      # (sentences, not columns)
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.txtResult.setMinimumHeight(150)
        lv.addWidget(self.txtResult, 1)

        # ================= right: plots =================================
        right = panel('Results')
        gl.addWidget(right, 1)
        rg = QtWidgets.QGridLayout(right)
        rg.setContentsMargins(10, 14, 10, 10)
        rg.setHorizontalSpacing(10)
        rg.setVerticalSpacing(8)

        def sub(r0, c0):
            cv = PlotCanvas()
            rg.addWidget(cv, r0, c0)
            ax = cv.figure.add_axes([0.15, 0.15, 0.8, 0.75])
            cv.layoutFcn = lambda: self._lay(cv, ax)
            return cv, ax

        self.cvS, self.axS = sub(0, 0)
        self.cvP, self.axP = sub(0, 1)
        self.cvC, self.axC = sub(1, 0)
        self.cvW, self.axW = sub(1, 1)

        self.btnRun.clicked.connect(lambda *_: self.safe(self.doRun, 'Run'))
        self.btnReset.clicked.connect(lambda *_: self.doReset())
        self.btnExport.clicked.connect(lambda *_: self.doExport())
        self.ddUnit.currentIndexChanged.connect(lambda *_: self.syncControls())
        self.ddMethod.currentIndexChanged.connect(
            lambda *_: self.syncControls())
        self.edTS.valueChanged.connect(lambda *_: self.onDisplay())
        self.cbLog.toggled.connect(lambda *_: self.onDisplay())
        self.syncControls()

        if sv and isinstance(sv.get('res'), dict) and np.size(
                sv['res'].get('pos', [])) == np.size(app.delay):
            try:
                self.S = dict(S0)
                self.info = dict(sv['info'])
                self.T = self.asResult(sv['res'])
                self.F = self.asResult(sv['fit']) if isinstance(
                    sv.get('fit'), dict) else None
                self.txtResult.setPlainText('\n'.join(sv['lines']))
                self.lblStatus.setText('%s  [restored]' % sv['status'])
                self.btnExport.setEnabled(True)
            except Exception:
                self.T = self.F = self.info = self.S = None
        self.drawAll()

    # ------------------------------------------------------------------
    @staticmethod
    def _lay(cv, ax):
        W, H = cv.sizePt()
        cv.place(ax, 62.0, 34.0, W - 62.0 - 14.0, H - 34.0 - 24.0)

    @staticmethod
    def asResult(d):
        """A stored result (dictionary, as a session file gives it) with
        its arrays as arrays."""
        R = Result()
        for k, v in d.items():
            if isinstance(v, (list, tuple, np.ndarray)):
                a = np.asarray(v)
                if k in ('ok', 'mask'):
                    a = a.astype(bool)
                elif a.dtype != bool:
                    a = a.astype(float)
                R[k] = a
            else:
                R[k] = v
        return R

    @staticmethod
    def tauText(v):
        v = np.asarray(v, float).ravel() if v is not None else np.zeros(0)
        return ' '.join('%.4g' % x for x in v)

    def tauStart(self):
        """Starting time constants as entered: an array (empty: the
        automatic starts). Raises ValueError for text that is no list of
        numbers."""
        txt = self.edTau.text().replace(',', ' ').split()
        try:
            return np.array([float(x) for x in txt], float)
        except ValueError:
            raise ValueError('The starting time constants must be numbers, '
                             'separated by spaces or commas.') from None

    def setup(self):
        """The setup as entered."""
        return {'x1': self.edX1.value(), 'x2': self.edX2.value(),
                'sign': self.ddSign.currentData(),
                'posUnit': self.ddUnit.currentData(),
                'density': bool(self.cbDens.isChecked()),
                'method': self.ddMethod.currentData(),
                'frac': self.edFrac.value(), 'minPct': self.edMin.value(),
                'tFrom': self.edT1.value(), 'tTo': self.edT2.value(),
                'nExp': int(self.ddN.currentText()),
                'tauStart': [float(x) for x in self.tauStart()],
                'tSpec': self.edTS.value(),
                'logT': bool(self.cbLog.isChecked())}

    def syncControls(self):
        app = self.app
        pu = self.ddUnit.currentData()
        # (the box names the conversion to the other axis; it can be set
        # only while the position is given there)
        other = pu != app.specUnit
        self.cbDens.setEnabled(other)
        self.cbDens.setText(ui_str(
            'Signal per %s: convert to a signal per %s'
            % (app.xUnitTxt(),
               'nm' if app.specUnit == 'cm-1' else 'cm-1')))
        para = self.ddMethod.currentData() == 'parabola'
        for w in (self.lblFrac, self.edFrac, self.lblFrac2):
            w.setEnabled(para)

    # ------------------------------------------------------------------
    def doRun(self):
        app = self.app
        S = self.setup()
        if S['tauStart'] and len(S['tauStart']) != S['nExp']:
            self.alert('Give %d starting time constant%s, or leave the field '
                       'empty.' % (S['nExp'], '' if S['nExp'] == 1 else 's'),
                       'Invalid')
            return
        self.lblStatus.setText('Following the band...')
        T, info = app.bandTrack(S)
        if T is None:
            self.lblStatus.setText(info['err'])
            return
        nOk = int(np.count_nonzero(T.ok))
        if nOk == 0:
            self.T, self.F, self.info, self.S = T, None, info, S
            self.lblStatus.setText(
                'No position found at any delay: check the range and the '
                'sign of the band.')
            self.txtResult.setPlainText('\n'.join(
                app.bandShiftLines(S, T, info, None)))
            self.btnExport.setEnabled(False)
            app.bpSaved = None
            self.drawAll()
            return
        F = fit_band_shift(app.delay, T.pos, [S['tFrom'], S['tTo']],
                           S['nExp'], S['tauStart'] or None)
        self.T, self.F, self.info, self.S = T, F, info, S
        lines = app.bandShiftLines(S, T, info, F)
        self.txtResult.setPlainText('\n'.join(lines))
        if F.ok:
            self.lblStatus.setText(
                'Done. Position at %d of %d delays; fit RMS %.3g %s.'
                % (nOk, T.ok.size, F.rms, S['posUnit']))
        else:
            self.lblStatus.setText(
                'Position at %d of %d delays. No fit: %s'
                % (nOk, T.ok.size, F.msg))
        self.btnExport.setEnabled(True)
        app.bpSaved = {'set': dict(S), 'info': dict(info),
                       'res': {k: (np.array(v) if isinstance(v, np.ndarray)
                                   else v) for k, v in T.items()},
                       'fit': {k: (np.array(v) if isinstance(v, np.ndarray)
                                   else v) for k, v in F.items()},
                       'lines': list(lines),
                       'status': self.lblStatus.text()}
        self.drawAll()

    def doReset(self):
        app = self.app
        app.bpSaved = None
        S0 = app.bandSetupDefault()
        self.edX1.setValue(S0['x1'])
        self.edX2.setValue(S0['x2'])
        self.ddSign.setCurrentIndex(max(self.ddSign.findData(S0['sign']), 0))
        self.ddUnit.setCurrentIndex(max(self.ddUnit.findData(S0['posUnit']),
                                        0))
        self.cbDens.setChecked(bool(S0['density']))
        self.ddMethod.setCurrentIndex(
            max(self.ddMethod.findData(S0['method']), 0))
        self.edFrac.setValue(S0['frac'])
        self.edMin.setValue(S0['minPct'])
        self.edT1.setValue(S0['tFrom'])
        self.edT2.setValue(S0['tTo'])
        self.ddN.setCurrentText(str(int(S0['nExp'])))
        self.edTau.setText('')
        self.T = self.F = self.info = self.S = None
        self.btnExport.setEnabled(False)
        self.lblStatus.setText('Reset.')
        self.txtResult.setPlainText('\n'.join(IDLE_TEXT))
        self.syncControls()
        self.drawAll()

    def onDisplay(self):
        """The delay of the spectrum shown or the scale of the delay
        axis changed: redraw, and keep the choice with the result."""
        app = self.app
        if isinstance(app.bpSaved, dict) and isinstance(
                app.bpSaved.get('set'), dict):
            app.bpSaved['set']['tSpec'] = self.edTS.value()
            app.bpSaved['set']['logT'] = bool(self.cbLog.isChecked())
        self.drawAll()

    # ------------------------------------------------------------------
    def curves(self):
        """What the plots show: delays, positions, fitted curve on a
        fine grid, C(t). (None, ...) without a result."""
        app = self.app
        T, F = self.T, self.F
        if T is None:
            return None
        t = np.asarray(app.delay, float).ravel()
        ok = np.asarray(T.ok, bool)
        isLog = bool(self.cbLog.isChecked() and np.any(t[ok] > 0))
        # (on a logarithmic axis only delays above 0 are shown)
        shown = ok & (t > 0) if isLog else ok
        out = Result(t=t, ok=ok, shown=shown, isLog=isLog,
                     pos=np.asarray(T.pos, float), tf=None, yf=None, c=None,
                     cf=None)
        if F is not None and F.ok:
            tm = t[np.asarray(F.mask, bool)]
            lo, hi = float(tm[0]), float(tm[-1])
            if isLog and lo > 0:
                tf = np.exp(np.linspace(np.log(lo), np.log(hi), 300))
            else:
                tf = np.linspace(lo, hi, 300)
                if isLog:
                    tf = tf[tf > 0]
            tau = np.asarray(F.tau, float)
            amp = np.asarray(F.A, float)
            yf = F.v_inf + np.exp(-tf.reshape(-1, 1)
                                  / tau.reshape(1, -1)) @ amp
            out.tf, out.yf = tf, yf
            # (a time constant at the limit of the search: the position
            # at long delays and C(t) mean nothing)
            if F.shift != 0 and not np.any(np.asarray(F.tau_lim, bool)):
                out.c = (out.pos - F.v_inf) / F.shift
                out.cf = (yf - F.v_inf) / F.shift
        return out

    def _tAxis(self, ax, t):
        """Scale and limits of a delay axis for the delays t shown."""
        t = np.asarray(t, float)
        t = t[np.isfinite(t)]
        if self.cbLog.isChecked() and np.any(t > 0):
            ax.set_xscale('log')
            tp = t[t > 0]
            if tp.size and np.max(tp) > np.min(tp):
                ax.set_xlim(float(np.min(tp)) / 1.3, float(np.max(tp)) * 1.3)
        elif t.size and np.max(t) > np.min(t):
            pad = 0.03 * float(np.max(t) - np.min(t))
            ax.set_xlim(float(np.min(t)) - pad, float(np.max(t)) + pad)

    def drawAll(self):
        app = self.app
        u = app.tUnitAx()
        S = self.S if self.S is not None else None
        pu = (S or {}).get('posUnit', self.ddUnit.currentData())
        puAx = 'cm$^{-1}$' if pu == 'cm-1' else 'nm'
        posName = 'Wavenumber' if pu == 'cm-1' else 'Wavelength'
        for ax, title in ((self.axS, T_SPEC), (self.axP, T_POS),
                          (self.axC, T_C), (self.axW, T_W)):
            ax.clear()
            styleAxes(ax, grid=True)
            ax.set_title(title, fontsize=FONT_PT + 1)
        self.axS.set_xlabel('%s (%s)' % (posName, puAx))
        self.axS.set_ylabel('Signal')
        for ax in (self.axP, self.axC, self.axW):
            ax.set_xlabel('Delay time (%s)' % u)
        self.axP.set_ylabel('Position (%s)' % puAx)
        self.axC.set_ylabel('C(t)')
        self.axW.set_ylabel('FWHM (%s)' % puAx)
        Cv = self.curves()
        if Cv is None or S is None:
            for cv in (self.cvS, self.cvP, self.cvC, self.cvW):
                cv.draw_idle()
            return
        T, F = self.T, self.F
        t, ok = Cv.t, Cv.ok

        # ---- the band at one delay ----
        v, Y, info = app.bandSpectra(S)
        if not info['err']:
            j = int(np.argmin(np.abs(t - self.edTS.value())))
            ax = self.axS
            hs = [ax.plot(v, Y[:, j], '.-', color=C_DATA, markersize=4,
                          linewidth=0.9, label='data')[0]]
            ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
            if ok[j]:
                if S['method'] == 'parabola':
                    _, _, i1, i2, _ = band_parabola(v, Y[:, j], S['frac'])
                    vv = np.linspace(v[i1], v[i2], 60)
                    hs.append(ax.plot(
                        vv, band_model_curve(vv, 'parabola', T.par[:, j]),
                        '-', color=C_FIT, linewidth=1.6,
                        label='parabola')[0])
                elif S['method'] == 'lognormal':
                    vv = np.linspace(v[0], v[-1], 400)
                    hs.append(ax.plot(
                        vv, band_model_curve(vv, 'lognormal', T.par[:, j]),
                        '-', color=C_FIT, linewidth=1.4,
                        label='log-normal band')[0])
                hs.append(ax.axvline(
                    T.pos[j], color=C_POS, linestyle='--', linewidth=1.3,
                    label='position: %.6g' % T.pos[j]))
                ax.set_title('Band at t = %.4g %s' % (t[j], app.tUnitAx()),
                             fontsize=FONT_PT + 1)
            else:
                ax.set_title('Band at t = %.4g %s: no position'
                             % (t[j], app.tUnitAx()), fontsize=FONT_PT + 1)
            ax.set_xlim(float(v[0]), float(v[-1]))
            legend(ax, handles=hs, fontsize=FONT_PT)

        # ---- position along the delays ----
        ax = self.axP
        shown = Cv.shown
        # (as in the MATLAB window, where plot() of empty data gives no
        # handle: a curve only when it has points, a legend only when
        # there is a curve)
        hp = []
        if np.any(shown):
            hp.append(ax.plot(t[shown], Cv.pos[shown], 'o', color=C_DATA,
                              markersize=3.5, label='position')[0])
        hasFine = Cv.tf is not None and np.size(Cv.tf) > 0
        if hasFine:
            hp.append(ax.plot(Cv.tf, Cv.yf, '-', color=C_FIT, linewidth=1.4,
                              label='fit')[0])
        if Cv.c is not None:
            hp.append(ax.axhline(F.v_inf, color=(0.4, 0.4, 0.4),
                                 linestyle='--', linewidth=1.0,
                                 label='at long delays'))
        self._tAxis(ax, t[shown])
        if hp:
            legend(ax, handles=hp, fontsize=FONT_PT)

        # ---- C(t) ----
        ax = self.axC
        if Cv.c is not None:
            m = shown & np.asarray(F.mask, bool)
            hc = []
            if np.any(m):
                hc.append(ax.plot(t[m], Cv.c[m], 'o', color=C_DATA,
                                  markersize=3.5, label='data')[0])
            if hasFine:
                hc.append(ax.plot(Cv.tf, Cv.cf, '-', color=C_FIT,
                                  linewidth=1.4, label='fit')[0])
            ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
            self._tAxis(ax, t[m])
            if hc:
                legend(ax, handles=hc, fontsize=FONT_PT)
        elif Cv.tf is not None and np.any(np.asarray(F.tau_lim, bool)):
            ax.set_title('C(t): a time constant is not determined',
                         fontsize=FONT_PT + 1)
        else:
            ax.set_title('C(t): needs the fit of the position',
                         fontsize=FONT_PT + 1)

        # ---- width ----
        ax = self.axW
        w = np.asarray(T.width, float)
        mw = shown & np.isfinite(w)
        if np.any(mw):
            ax.plot(t[mw], w[mw], 'o', color=C_DATA, markersize=3.5)
            self._tAxis(ax, t[mw])
        for cv in (self.cvS, self.cvP, self.cvC, self.cvW):
            cv.draw_idle()

    # ------------------------------------------------------------------
    def doExport(self):
        if self.T is None or self.S is None:
            return
        app = self.app
        full = app.putFile('CSV (*.csv);;DAT (tab) (*.dat)',
                           'Save the position of the band',
                           '%s_band_position.csv' % (app.dataStem() or 'data'),
                           parent=self)
        if not full:
            return
        dlm = ',' if full.lower().endswith('.csv') else '\t'
        try:
            hdr, cols = app.bandExportTable(self.S, self.T, self.F)
            writeTable(full, hdr, cols, dlm)
            self.inform('Saved: %s' % full, 'Exported')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Error')
