"""Auto Analysis window. Port of openAutoWindow.

The last run is kept in app.autoSaved (keepState): a window opened again
shows it (restoreState) as long as the data of the main window are the
ones it describes. A change of the data drops it, a session stores it.
"""
import os
import time

import numpy as np
from matplotlib.ticker import FuncFormatter
from PySide6 import QtCore, QtWidgets

from trspec.auto import (auto_done_text, auto_preprocess,
                         auto_report_lines, auto_select_model,
                         auto_summarize)
from trspec.kinetics import ca_options
from trspec._util import Result

from .textfit import legend, setSpan, setTitle
from .widgets import (AX_BG, FitLabel, NumEdit, PlotCanvas, button,
                      colorButton, label, lineColor, runProgress, sciTicks,
                      sidePanel, styleAxes, textArea, tip)
from .win_base import ChildWindow, mapLayout, mapOn

RUN_COLOR = (0.20, 0.50, 0.45)
IDLE = 'Press "Run"'
HELP_TEXT = [
    'Run analyses the data as they are shown in the main window and '
    'changes nothing.',
    '',
    'Report: what was found, with the numbers behind every decision, and '
    'a list of things to check by eye.',
    '',
    'Corrections: the rise of every wavelength with the fitted time zero, '
    'the data after the corrections, the background and the noise of '
    'every channel.',
    '',
    'Components: singular values, RMS against the number of components, '
    'kinetics at a few wavelengths (data and fit), DADS and EADS of the '
    'fit in use (dashed: coherent-artifact terms).',
    '',
    'Residual: what that fit leaves, and its spectrum (oscillations).',
    '',
    '"Components used" switches the report and the plots to a fit with '
    'another number of components.',
    '',
    'The report describes the data. It does not decide which states and '
    'reaction steps are behind the time constants.']
CA_NAMES = ('artifact: Gaussian', 'artifact: 1st derivative',
            'artifact: 2nd derivative')
CA_COLORS = ((0.15,) * 3, (0.45,) * 3, (0.70,) * 3)
KIN_LEGEND_PT = 86.0        # room beside the kinetics for their legend


def lay1(cv, ax, L=58.0, B=36.0, R=16.0, T=24.0):
    # (B and T as the maps: the plots of a row begin and end together)
    W, H = cv.sizePt()
    cv.place(ax, L, B, W - L - R, H - B - T)
    setSpan(ax, 0.0, W)


def logY(ax):
    """Logarithmic y-axis with short tick labels (0.03 instead of
    3 x 10^-2, which would push the axis label off the picture)."""
    ax.set_yscale('log')

    def minor(v, _pos):
        lo, hi = ax.get_ylim()
        if not (lo > 0 and hi > 0 and v > 0):
            return ''
        dec = abs(np.log10(hi / lo))
        if dec > 1.5:
            return ''
        m = v / 10.0 ** np.floor(np.log10(v))
        keep = (2, 3, 5) if dec > 0.8 else (2, 3, 4, 5, 6, 8)
        if abs(m - round(m)) < 1e-6 and int(round(m)) in keep:
            return '%g' % v
        return ''

    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _pos: '%g' % v))
    ax.yaxis.set_minor_formatter(FuncFormatter(minor))


def asResult(obj):
    """obj with every dictionary in it as a Result (a kept run that came
    out of a session file holds plain dictionaries)."""
    if isinstance(obj, dict):
        return Result((k, asResult(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return [asResult(v) for v in obj]
    if isinstance(obj, tuple):
        return tuple(asResult(v) for v in obj)
    return obj


def scanNumbers(text):
    """Numbers at the start of text, separated by spaces, commas or
    semicolons; reading stops at the first entry that is no number
    (sscanf '%f')."""
    out = []
    for part in text.replace(',', ' ').replace(';', ' ').split():
        try:
            out.append(float(part))
        except ValueError:
            break
    return out


class AutoWindow(ChildWindow):
    KEY = 'auto'

    def __init__(self, app):
        super().__init__(app, 'Auto Analysis', 1380, 860)
        self.P = None
        self.S = None
        self.L = []
        self.sumErr = ''
        self.applied = False
        # data the result describes: those of the run; after "Apply" the
        # corrected data
        self.srcWl = self.srcT = self.srcD = None
        self.srcDesc = ''
        self.runD = None        # data given to the run (the plot "before
                                # the corrections")
        self.runSet = None      # options of the run
        self.kMode = 'lin'
        self.kinIdx = []
        self.maps = {}          # what the three maps show: (wl, t, Z)
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: what to do =============================
        left, lv = sidePanel('Automatic analysis', 340)
        gl.addWidget(left)

        def bold(text):
            lb = label(text)
            lb.setStyleSheet('font-weight: bold;')
            lv.addWidget(lb)

        def hrow(indent=0):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(6)
            r.setContentsMargins(indent, 0, 0, 0)
            lv.addLayout(r)
            return r

        def check(text, on, tooltip):
            cb = QtWidgets.QCheckBox(text)
            cb.setChecked(on)
            cb.setToolTip(tip(tooltip))
            lv.addWidget(cb)
            return cb

        bold('1. Corrections to propose')
        self.cbBg = check(
            'Background (delays before the signal)', True,
            'Finds the delays that lie before the signal and subtracts '
            'their mean spectrum. Skipped when that mean is zero within '
            'the noise, or when fewer than 3 delays lie safely before the '
            'signal.')
        self.cbMask = check(
            'Mask channels with excess noise', True,
            'Channels whose noise is more than 3 times that of their '
            'neighbors (pump scatter, dead pixels, the ends of the probe '
            'spectrum) can be set to NaN and left out of every fit. Which '
            'of them are dropped is set by "Mask rule". Dropped channels '
            'at an end of the spectrum are cut away instead when "Crop" '
            'is ticked.')
        r = hrow(18)
        self.lblMask = label('Mask rule:')
        r.addWidget(self.lblMask)
        self.ddMask = QtWidgets.QComboBox()
        self.ddMask.addItem('Only where the signal is lost', 'weak')
        self.ddMask.addItem('Every noisy channel', 'all')
        self.ddMask.setToolTip(tip(
            'Which of the noisy channels are dropped.\n'
            'Only where the signal is lost: a noisy channel is dropped '
            'when its signal is less than 10 times its noise, or when its '
            'noise is more than 10 times that of the neighbors. Channels '
            'that are noisier than their neighbors but show a clear '
            'signal stay in the data (the pump region with a strong '
            'bleach, the ends of the probe spectrum); the report lists '
            'them.\n'
            'Every noisy channel: every channel more than 3 times noisier '
            'than its neighbors is dropped, together with the neighbors '
            'that are still 1.5 times noisier.'))
        r.addWidget(self.ddMask, 1)
        self.cbChirp = check(
            'Chirp (time zero along the spectrum)', True,
            'Locates the steepest rise of every wavelength, fits the '
            'chirp curve of the Chirp Corr window through them and shifts '
            'every wavelength so that its time zero is 0. If time zero is '
            'the same everywhere but not at 0, the delay axis is shifted '
            'instead.')
        self.cbCrop = check(
            'Crop what the corrections leave empty', True,
            'Removes the first and last delays that the chirp correction '
            'leaves without data at some wavelengths, and noisy channels '
            'at the ends of the spectrum.')
        bold('2. Summary of the kinetics')
        r = hrow()
        r.addWidget(label('Largest number of components:'))
        self.edMax = NumEdit(
            5, limits=(1, 5), integer=True, width=60,
            tooltip='Global fits are run with 1, 2, ... components until '
                    'one more lowers the RMS by less than 5 % or brings a '
                    'time constant the data do not determine. This is the '
                    'largest number tried.')
        r.addWidget(self.edMax)
        r.addStretch(1)
        r = hrow()
        r.addWidget(label('Coherent artifact:'))
        self.ddCA = QtWidgets.QComboBox()
        for nm, key in (('Automatic', 'auto'), ('Off', 'off'),
                        ('On (3 terms)', 'on')):
            self.ddCA.addItem(nm, key)
        self.ddCA.setToolTip(tip(
            'Terms for a signal confined to time zero (cross-phase '
            'modulation, two-photon absorption, Raman of the solvent): a '
            'Gaussian of the IRF width and its first two derivatives, each '
            'with a free spectrum, as in the Global Analysis window. '
            'Automatic: used when the fits without them need a time '
            'constant below a third of the IRF width, or when the terms '
            'lower the RMS by 15 % or more. Off: never. On: always. The '
            'report says what was done and why.'))
        r.addWidget(self.ddCA, 1)
        self.cbAll = check(
            'Fit every number up to the largest', False,
            'Normally the fits stop at the first number of components '
            'that is rejected. Ticked: all numbers up to the largest are '
            'fitted, so that each of them can be picked under "Components '
            'used". The automatic choice stays the same. Takes longer.')
        self.btnRun = colorButton(
            'Run', RUN_COLOR,
            tooltip='Analyze the data as they are shown in the main '
                    'window. Nothing is changed.')
        self.btnRun.setStyleSheet(self.btnRun.styleSheet()
                                  + 'QPushButton { font-weight: bold; }')
        self.btnRun.setMinimumHeight(32)
        lv.addWidget(self.btnRun)
        bold('3. Use the result')
        r = hrow()
        r.addWidget(label('Components used:'))
        self.ddN = QtWidgets.QComboBox()
        self.ddN.addItem('(run first)')
        self.ddN.setEnabled(False)
        self.ddN.setToolTip(tip(
            'Number of components of the fit that section 3 of the '
            'report, the Components and Residual tabs and "Send model to '
            'Global Analysis" use. After a run it is the automatic choice. '
            'Every number in the table of fits can be picked instead, also '
            'one that the automatic choice rejected; the report then says '
            'so and why it was rejected. Nothing is fitted again.'))
        r.addWidget(self.ddN, 1)
        self.btnApply = button(
            'Apply corrections to the data',
            tooltip='Write the proposed corrections into the working data '
                    'of the main window. Like Crop, this makes the '
                    'corrected data the new baseline; "Reset all" in the '
                    'main window returns to the loaded data.')
        self.btnGA = button(
            'Send model to Global Analysis',
            tooltip='Open the Global Analysis window with the number of '
                    'components in use ("Components used"), the fitted '
                    'values as starting values and the same coherent-'
                    'artifact terms. A stored global-fit result is '
                    'discarded. Press "Run Fit" there.')
        self.btnSave = button('Save report...')
        for b in (self.btnApply, self.btnGA, self.btnSave):
            b.setEnabled(False)
            b.setMinimumHeight(28)
            lv.addWidget(b)
        # (three lines: the line after a run can name steps that were not
        # determined, and on a small screen with a scaled display the
        # longest one does not fit into two)
        self.lblStatus = label('Ready.', color=(0.2, 0.2, 0.2))
        self.lblStatus.setWordWrap(True)
        self.lblStatus.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft
                                    | QtCore.Qt.AlignmentFlag.AlignTop)
        self.lblStatus.setFixedHeight(
            3 * self.lblStatus.fontMetrics().lineSpacing() + 4)
        lv.addWidget(self.lblStatus)
        self.txtHelp = QtWidgets.QPlainTextEdit()
        self.txtHelp.setReadOnly(True)
        self.txtHelp.setLineWrapMode(
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.txtHelp.setPlainText('\n'.join(HELP_TEXT))
        self.txtHelp.setMinimumHeight(70)
        lv.addWidget(self.txtHelp, 1)

        # ================= right: report and pictures ===================
        self.tabs = QtWidgets.QTabWidget()
        gl.addWidget(self.tabs, 1)

        # ---- Report ----------------------------------------------------
        tabRep = QtWidgets.QWidget()
        g = QtWidgets.QVBoxLayout(tabRep)
        g.setContentsMargins(6, 6, 6, 6)
        self.txtRep = textArea(['Press "Run".'])
        self.txtRep.setLineWrapMode(      # (the long lines of section 5)
            QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
        g.addWidget(self.txtRep)
        self.tabs.addTab(tabRep, 'Report')

        def canvas(nAxes=1):
            cv = PlotCanvas()
            cv.setMinimumHeight(150)
            axs = [cv.figure.add_axes([0.1 + 0.4 * k, 0.15, 0.35, 0.7])
                   for k in range(nAxes)]
            return (cv, *axs)

        def mapCanvas():
            cv, ax = canvas()
            cax = cv.figure.add_axes([0.9, 0.15, 0.02, 0.7])
            cv.layoutFcn = lambda: mapLayout(
                cv, ax, cax if cax.get_visible() else None, L=58.0, R=74.0,
                T=24.0, B=36.0)
            return cv, ax, cax

        def plain(cv, ax, L=58.0, R=16.0):
            cv.layoutFcn = lambda: lay1(cv, ax, L=L, R=R)

        # ---- Corrections -----------------------------------------------
        tabCor = QtWidgets.QWidget()
        g = QtWidgets.QGridLayout(tabCor)
        g.setContentsMargins(8, 8, 8, 8)
        g.setSpacing(8)
        self.cvRaw, self.axRaw, self.caxRaw = mapCanvas()
        self.cvCor, self.axCor, self.caxCor = mapCanvas()
        self.cvBg, self.axBg = canvas()
        self.cvNz, self.axNz = canvas()
        # (under the maps: the room of their color bars stays free)
        plain(self.cvBg, self.axBg, R=74.0)
        plain(self.cvNz, self.axNz, R=74.0)
        g.addWidget(self.cvRaw, 0, 0)
        g.addWidget(self.cvCor, 0, 1)
        g.addWidget(self.cvBg, 1, 0)
        g.addWidget(self.cvNz, 1, 1)
        g.setRowStretch(0, 13)
        g.setRowStretch(1, 10)
        g.setColumnStretch(0, 1)
        g.setColumnStretch(1, 1)
        self.tabs.addTab(tabCor, 'Corrections')

        # ---- Components ------------------------------------------------
        tabCmp = QtWidgets.QWidget()
        g = QtWidgets.QGridLayout(tabCmp)
        g.setContentsMargins(8, 8, 8, 8)
        g.setSpacing(8)
        self.cvSV, self.axSV = canvas()
        self.cvRms, self.axRms = canvas()
        plain(self.cvSV, self.axSV)
        plain(self.cvRms, self.axRms)
        tl = QtWidgets.QHBoxLayout()
        tl.setSpacing(8)
        tl.addWidget(self.cvSV, 1)
        tl.addWidget(self.cvRms, 1)
        tlv = QtWidgets.QVBoxLayout()
        tlv.setSpacing(4)
        tlv.addLayout(tl, 1)
        g.addLayout(tlv, 0, 0)
        kin = QtWidgets.QVBoxLayout()
        kin.setSpacing(4)
        kc = QtWidgets.QHBoxLayout()
        kc.setSpacing(6)
        kc.addWidget(label('Kinetics at (%s):' % app.xUnitTxt()))
        self.edKin = QtWidgets.QLineEdit()
        self.edKin.setFixedWidth(170)
        self.edKin.setToolTip(tip(
            'Positions on the spectral axis whose kinetics are shown, '
            'data and fit: up to 6 numbers separated by spaces or commas, '
            'for example  520 600 680. The nearest channel of the data is '
            'taken. Empty: chosen automatically (the channel with the '
            'largest signal in each fifth of the spectrum).'))
        kc.addWidget(self.edKin)
        self.lblKin = FitLabel(color=(0.35, 0.35, 0.35))
        kc.addWidget(self.lblKin, 1)
        kin.addLayout(kc)
        # (the two plots left of the kinetics begin as far down as the
        # kinetics plot, which has this row of controls above it)
        # (the row and the spacing under it; no spacing follows a spacer)
        tlv.insertSpacing(0, self.edKin.sizeHint().height() + kin.spacing())
        self.cvKin, self.axKinL, self.axKinR = canvas(2)
        self.cvKin.layoutFcn = self._layKin
        kin.addWidget(self.cvKin, 1)
        g.addLayout(kin, 0, 1)
        self.cvDADS, self.axDADS = canvas()
        self.cvEADS, self.axEADS = canvas()
        # (the two are equally wide; the EADS do not leave the room of
        # the legend of the kinetics above them free)
        plain(self.cvDADS, self.axDADS)
        plain(self.cvEADS, self.axEADS)
        g.addWidget(self.cvDADS, 1, 0)
        g.addWidget(self.cvEADS, 1, 1)
        for k in (0, 1):
            g.setRowStretch(k, 1)
            g.setColumnStretch(k, 1)
        self.tabs.addTab(tabCmp, 'Components')

        # ---- Residual --------------------------------------------------
        tabRes = QtWidgets.QWidget()
        g = QtWidgets.QGridLayout(tabRes)
        g.setContentsMargins(8, 8, 8, 8)
        g.setSpacing(8)
        self.cvRes, self.axRes, self.caxRes = mapCanvas()
        self.cvPow, self.axPow = canvas()
        self.cvProf, self.axProf = canvas()
        # (beside the map: the room of its color bar, so that all plot
        # boxes are equally wide)
        # (R: about as wide as the blank left of the label of the y-axis)
        plain(self.cvPow, self.axPow, L=62.0, R=22.0)
        plain(self.cvProf, self.axProf, L=62.0, R=22.0)
        g.addWidget(self.cvRes, 0, 0, 2, 1)
        g.addWidget(self.cvPow, 0, 1)
        g.addWidget(self.cvProf, 1, 1)
        for k in (0, 1):
            g.setRowStretch(k, 1)
            g.setColumnStretch(k, 1)
        self.tabs.addTab(tabRes, 'Residual')

        # ---- what the controls do --------------------------------------
        self.cbMask.toggled.connect(self._maskEnable)
        self.btnRun.clicked.connect(lambda *_: self.safe(self.doRun, 'Run'))
        self.btnApply.clicked.connect(
            lambda *_: self.safe(self.doApply, 'Apply corrections'))
        self.btnGA.clicked.connect(
            lambda *_: self.safe(self.doSendGA, 'Send model'))
        self.btnSave.clicked.connect(
            lambda *_: self.safe(self.doSave, 'Save report'))
        self.ddN.activated.connect(
            lambda *_: self.safe(self.onPickN, 'Components used'))
        self.edKin.editingFinished.connect(
            lambda: self.safe(self.onKinEdit, 'Kinetics'))
        self.emptyAxes()
        self.restoreState()     # the last run, if there is one

    # ------------------------------------------------------------------
    def safe(self, fcn, what):
        """As in the other windows; in addition "Run" is switched on
        again after an error."""
        super().safe(fcn, what)
        self.btnRun.setEnabled(True)

    def _maskEnable(self, on):
        self.ddMask.setEnabled(bool(on))
        self.lblMask.setEnabled(bool(on))

    def _canvases(self):
        return (self.cvRaw, self.cvCor, self.cvBg, self.cvNz, self.cvSV,
                self.cvRms, self.cvKin, self.cvDADS, self.cvEADS,
                self.cvRes, self.cvPow, self.cvProf)

    def _maps(self):
        return ((self.cvRaw, self.axRaw, self.caxRaw),
                (self.cvCor, self.axCor, self.caxCor),
                (self.cvRes, self.axRes, self.caxRes))

    def _clear(self, ax, grid=True):
        ax.clear()
        ax.set_facecolor(AX_BG)
        ax.__dict__.pop('_tfLeg', None)
        styleAxes(ax, grid=grid)

    def _mapOff(self, cv, ax, cax):
        """Empty map axes without color bar."""
        self.maps.pop(ax, None)
        self._clear(ax, grid=False)
        cax.clear()
        cax.set_visible(False)
        cv.relayout()

    def emptyAxes(self):
        for cv, ax, cax in self._maps():
            self._mapOff(cv, ax, cax)
            setTitle(ax, IDLE)
        for ax in (self.axBg, self.axNz, self.axSV, self.axRms, self.axDADS,
                   self.axEADS, self.axPow, self.axProf):
            self._clear(ax, grid=False)
            setTitle(ax, IDLE)
        self.drawKinetics()
        for cv in self._canvases():
            cv.draw_idle()

    # ------------------------------------------------------------------
    # run
    # ------------------------------------------------------------------
    def options(self):
        opt = {'doBg': self.cbBg.isChecked(),
               'doMask': self.cbMask.isChecked(),
               'doChirp': self.cbChirp.isChecked(),
               'doCrop': self.cbCrop.isChecked(),
               'maxComp': int(round(self.edMax.value())),
               'ca': self.ddCA.currentData(),
               'tryAll': self.cbAll.isChecked()}
        if self.ddMask.currentData() == 'all':
            opt['maskSN'] = np.inf
        return opt

    def doRun(self):
        app = self.app
        opt = self.options()
        # (the options as they are set in the window: kept with the
        # result)
        self.runSet = {k: v for k, v in opt.items() if k != 'maskSN'}
        self.runSet['maskRule'] = self.ddMask.currentData()
        self.srcWl = np.array(app.wavelength, float).ravel()
        self.srcT = np.array(app.delay, float).ravel()
        self.srcD = np.array(app.deltaA, float)
        self.runD = self.srcD
        self.P = None
        self.S = None
        self.L = []
        self.sumErr = ''
        self.applied = False
        app.autoSaved = None        # a new run replaces the kept one
        for b in (self.btnApply, self.btnGA, self.btnSave, self.btnRun):
            b.setEnabled(False)
        self.ddN.clear()
        self.ddN.addItem('(run first)')
        self.ddN.setEnabled(False)
        tu, xu = app.tUnitHdr(), app.xUnitTxt()
        wl, t, D = self.srcWl, self.srcT, self.srcD

        def work(prog):
            prog(0.02, 'Corrections...')
            P = auto_preprocess(wl, t, D, opt, tu, xu)
            prog(0.12, 'Singular values...')
            S, err = None, ''
            try:
                S = auto_summarize(P, opt, prog)
            except Exception as errS:
                if os.environ.get('TRSPEC_DEBUG') and \
                        not isinstance(errS, ValueError):
                    import traceback
                    traceback.print_exc()
                err = str(errS)
            return P, S, err

        try:
            self.P, self.S, self.sumErr = runProgress(
                app, self, 'Auto Analysis', 'Corrections...', work,
                lambda fr, msg: self.lblStatus.setText(msg))
        except Exception as errP:
            self.txtRep.setPlainText('\n'.join(
                ['The analysis could not be run:', str(errP)]))
            self.emptyAxes()
            raise
        finally:
            self.btnRun.setEnabled(True)
        P = self.P
        self.srcDesc = app.dataSourceDesc or 'current data'
        self.refreshReport()
        self.fillPickN()
        self.drawCorrections()
        self.drawComponents()
        self.drawResidual()
        self.btnApply.setEnabled(bool(P.changed))
        self.btnGA.setEnabled(self.S is not None)
        self.btnSave.setEnabled(True)
        self.tabs.setCurrentIndex(0)
        self.lblStatus.setText(auto_done_text(P))
        self.keepState()

    def keepState(self):
        """Keep the run in app.autoSaved, so that the window shows it
        again after it was closed and a saved session holds it. Called
        after every step that changes what the window shows of a run."""
        if self.P is None:
            return
        self.app.autoSaved = {
            'P': self.P, 'S': self.S, 'sumErr': self.sumErr,
            'applied': bool(self.applied), 'srcDesc': self.srcDesc,
            'srcWl': self.srcWl, 'srcT': self.srcT, 'srcD': self.srcD,
            'runD': self.runD, 'set': dict(self.runSet),
            'kin': self.edKin.text()}

    def restoreState(self):
        """Show the kept run again if it still describes the data of the
        main window. Every change of those data drops it
        (invalidateChildResults); the comparison made here is a second
        check."""
        app = self.app
        st = getattr(app, 'autoSaved', None)
        if not isinstance(st, dict) or st.get('P') is None:
            return
        try:
            if not self.keptFor(st['srcWl'], st['srcT'], st['srcD']):
                app.autoSaved = None
                return
            # the options of that run
            o = dict(st['set'])
            self.cbBg.setChecked(bool(o['doBg']))
            self.cbMask.setChecked(bool(o['doMask']))
            self.ddMask.setCurrentIndex(
                max(0, self.ddMask.findData(str(o['maskRule']))))
            self._maskEnable(self.cbMask.isChecked())
            self.cbChirp.setChecked(bool(o['doChirp']))
            self.cbCrop.setChecked(bool(o['doCrop']))
            self.edMax.setValue(int(o['maxComp']))
            self.ddCA.setCurrentIndex(
                max(0, self.ddCA.findData(str(o['ca']))))
            self.cbAll.setChecked(bool(o['tryAll']))
            self.runSet = o
            # the result
            self.P = P = asResult(st['P'])
            self.S = asResult(st['S'])
            self.sumErr = str(st['sumErr'] or '')
            self.applied = bool(st['applied'])
            self.srcDesc = str(st['srcDesc'] or '')
            # (the data as they are now: see keptFor)
            self.srcWl = np.array(app.wavelength, float).ravel()
            self.srcT = np.array(app.delay, float).ravel()
            self.srcD = np.array(app.deltaA, float)
            self.runD = np.array(st['runD'], float)
            self.edKin.setText(str(st['kin'] or ''))
            self.refreshReport()
            self.fillPickN()
            self.drawCorrections()
            self.drawComponents()
            self.drawResidual()
            self.btnApply.setEnabled(bool(P.changed) and not self.applied)
            self.btnGA.setEnabled(self.S is not None)
            self.btnSave.setEnabled(True)
            if self.applied:
                msg = 'Corrections applied to the main data.'
            else:
                msg = auto_done_text(P)
            self.lblStatus.setText('%s  [restored]' % msg)
        except Exception:
            # A kept run that cannot be shown: back to the empty window.
            if os.environ.get('TRSPEC_DEBUG'):
                import traceback
                traceback.print_exc()
            app.autoSaved = None
            self.P = None
            self.S = None
            self.L = []
            self.sumErr = ''
            self.applied = False
            self.srcWl = self.srcT = self.srcD = None
            self.txtRep.setPlainText('Press "Run".')
            self.ddN.clear()
            self.ddN.addItem('(run first)')
            self.ddN.setEnabled(False)
            for b in (self.btnApply, self.btnGA, self.btnSave):
                b.setEnabled(False)
            self.emptyAxes()
            self.lblStatus.setText(
                'Stored result could not be restored -- run again.')

    def keptFor(self, wl0, t0, D0):
        """Are the data of the main window those a kept run describes?
        The axes have to agree exactly, the values within rounding:
        "Load Session" builds the working data again from the raw data
        and the corrections, and another version of the libraries or
        another computer may do that with a difference in the last
        digit."""
        app = self.app
        wl0 = np.asarray(wl0, float).ravel()
        t0 = np.asarray(t0, float).ravel()
        D0 = np.asarray(D0, float)
        wl = np.asarray(app.wavelength, float).ravel()
        t = np.asarray(app.delay, float).ravel()
        D = np.asarray(app.deltaA, float)
        if not (wl.shape == wl0.shape and np.array_equal(wl, wl0)
                and t.shape == t0.shape and np.array_equal(t, t0)
                and D.shape == D0.shape
                and np.array_equal(np.isnan(D), np.isnan(D0))):
            return False
        fin = ~np.isnan(D0)
        if not fin.any():
            return True                 # no numbers to compare
        with np.errstate(invalid='ignore'):
            sc = float(np.max(np.abs(D0[fin])))
            dv = np.abs(D[fin] - D0[fin])
            dv = dv[~np.isnan(dv)]
        if dv.size == 0 or np.isnan(sc):
            return True
        return bool(float(np.max(dv)) <= 1e-9 * sc)

    def onKinEdit(self):
        """Other wavelengths typed into "Kinetics at": they are kept with
        the run."""
        self.drawKinetics()
        self.keepState()

    def refreshReport(self):
        L = auto_report_lines(self.P, self.S, self.srcDesc, self.app.Version)
        if self.sumErr:
            L += ['', 'The summary of the kinetics failed: %s' % self.sumErr]
        if self.applied:
            L += ['', 'The corrections were applied to the data of the main '
                      'window.']
        self.L = L
        self.txtRep.setPlainText('\n'.join(L))

    def fillPickN(self):
        S = self.S
        self.ddN.clear()
        if S is None:
            self.ddN.addItem('(no summary)')
            self.ddN.setEnabled(False)
            return
        for q in range(1, len(S.tabN) + 1):
            if q == S.nStar:
                self.ddN.addItem('%d  (automatic choice)' % q)
            elif S.tabWhy[q - 1]:
                self.ddN.addItem('%d  (rejected)' % q)
            else:
                self.ddN.addItem('%d' % q)
        self.ddN.setCurrentIndex(S.nUse - 1)
        self.ddN.setEnabled(len(S.tabN) > 1)

    def onPickN(self):
        if self.S is None:
            return
        n = self.ddN.currentIndex() + 1
        if n < 1 or n == self.S.nUse:
            return
        self.lblStatus.setText('Updating...')
        QtWidgets.QApplication.processEvents(
            QtCore.QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents)
        self.S = auto_select_model(self.P, self.S, n)
        self.refreshReport()
        self.drawComponents()
        self.drawResidual()
        if n == self.S.nStar:
            self.lblStatus.setText(
                'In use: %d component(s), the automatic choice.' % n)
        else:
            self.lblStatus.setText(
                'In use: %d component(s), set by hand.' % n)
        self.keepState()

    def pickN(self, n):
        """Choose the number of components as the user would (tests)."""
        self.ddN.setCurrentIndex(int(n) - 1)
        self.safe(self.onPickN, 'Components used')

    def fitName(self):
        if self.S.nUse == self.S.nStar:
            return 'the chosen fit'
        return 'the fit with %d components' % self.S.nUse

    def sameData(self):
        app = self.app
        wl = np.asarray(app.wavelength, float).ravel()
        t = np.asarray(app.delay, float).ravel()
        D = np.asarray(app.deltaA, float)
        return (self.srcWl is not None and wl.shape == self.srcWl.shape
                and np.array_equal(wl, self.srcWl)
                and t.shape == self.srcT.shape
                and np.array_equal(t, self.srcT)
                and D.shape == self.srcD.shape
                and np.array_equal(D, self.srcD, equal_nan=True))

    # ------------------------------------------------------------------
    # use the result
    # ------------------------------------------------------------------
    def doApply(self):
        app, P = self.app, self.P
        if P is None or self.applied:
            return
        if app.aniOn and app.aniView == 'r':
            self.alert('The corrections cannot be applied while the main '
                       'window shows the anisotropy r. Choose the view '
                       '"Parallel" or "Magic angle", run again and apply '
                       'there: both polarizations are then corrected alike.',
                       'Anisotropy view')
            return
        if not P.changed:
            self.inform('The report proposes no correction.',
                        'Nothing to apply')
            return
        if not self.sameData():
            self.alert('The data in the main window changed after the run. '
                       'Press "Run" again.', 'Data changed')
            return
        if not app.confirm(
                'The corrections of the report are written into the '
                'working data and become the new baseline (as with Crop). '
                'The analysis windows close and their stored results are '
                'discarded. "Reset all" in the main window returns to the '
                'loaded data.', 'Apply corrections', yes='Apply',
                no='Cancel', parent=self):
            return
        tK = np.asarray(P.tIn, float)[np.asarray(P.keepT, bool)]
        for key in ('bg', 'chirp', 'crop', 'mask'):
            w = app.childWindows.get(key)
            if w is not None:
                try:
                    w.close()
                except Exception:
                    pass
        # (this window stays open although the data change)
        app.childWindows.pop(self.KEY, None)
        try:
            app.deltaA = np.array(P.Dfull, float)
            app.pendingAuto = app.autoLogEntry(P)
            app.applyCropByRange(float(np.min(P.wl)), float(np.max(P.wl)),
                                 float(tK[0]), float(tK[-1]))
            if P.chirp['shift'] != 0:
                app.applyZeroTimeShift(P.chirp['shift'])
        finally:
            app.childWindows[self.KEY] = self
        self.applied = True
        self.btnApply.setEnabled(False)
        self.srcWl = np.array(app.wavelength, float).ravel()
        self.srcT = np.array(app.delay, float).ravel()
        self.srcD = np.array(app.deltaA, float)
        self.refreshReport()
        self.lblStatus.setText('Corrections applied to the main data.')
        # (the data change dropped every kept result, this one included:
        # keep it again, now for the corrected data)
        self.keepState()
        self.raise_()
        self.activateWindow()

    def doSendGA(self):
        app, P, S = self.app, self.P, self.S
        if S is None:
            return
        if P.changed and not self.applied:
            self.inform('The model belongs to the corrected data. Press '
                        '"Apply corrections to the data" first.',
                        'Apply first')
            return
        if not self.sameData():
            self.alert('The data in the main window changed after the run. '
                       'Press "Run" again.', 'Data changed')
            return
        Mo = S.model
        w = app.childWindows.get('ga')
        if w is not None:
            w.close()
        nT = int(Mo.tau.size)
        app.gaNComp = nT
        app.gaTauInit = [float(v) for v in Mo.tau]
        app.gaTauFixed = [False] * nT
        app.gaBetaInit = [1.0] * nT
        app.gaStretchOn = [False] * nT
        app.gaHasInf = bool(Mo.hasInf)
        app.gaIRFTreat = Mo.irf
        app.gaT0 = float(Mo.t0)
        app.gaFWHM = float(Mo.fwhm)
        if Mo.irf == 'none' or Mo.irfSeed:
            app.gaT0Fixed = True
            app.gaFWHMFixed = True
        else:
            app.gaT0Fixed = False
            app.gaFWHMFixed = bool(Mo.fwFixed)
        app.gaCA = ca_options({'order': int(Mo.caK), 'mode': 'irf',
                               'fwhm': np.nan})
        if np.all(S.cols):
            app.gaFitTMin = None
            app.gaFitTMax = None
        else:
            app.gaFitTMin = float(S.t[0])
            app.gaFitTMax = float(S.t[-1])
        app.clearGaResults()
        app.openGlobalAnalysisWindow()
        self.lblStatus.setText('Model sent: press "Run Fit" there.')

    def doSave(self):
        if not self.L:
            return
        defName = 'TRSpec_auto_report.txt'
        base = self.app.dataStem()
        if base:
            defName = 'TRSpec_%s_auto_report.txt' % base
        full = self.app.putFile('Text (*.txt)', 'Save report', defName,
                                parent=self)
        if not full:
            return
        try:
            fid = open(full, 'w', encoding='utf-8', newline='\n')
        except OSError:
            self.alert('The file could not be opened for writing.',
                       'Save report')
            return
        with fid:
            fid.write('TRSpecAnalyzer, saved %s\n\n'
                      % time.strftime('%Y-%m-%d %H:%M'))
            for ln in self.L:
                fid.write(ln + '\n')
        self.lblStatus.setText('Report saved: %s' % os.path.basename(full))

    # ------------------------------------------------------------------
    # pictures
    # ------------------------------------------------------------------
    def heat(self, cv, ax, cax, wlv, tv, Z):
        ax.__dict__.pop('_tfLeg', None)
        cax.set_visible(True)
        wlv = np.asarray(wlv, float).ravel()
        tv = np.asarray(tv, float).ravel()
        self.maps[ax] = (wlv, tv, np.asarray(Z, float))
        mapOn(self.app, ax, cax, wlv, tv, Z)
        ax.set_xlim(float(np.min(wlv)), float(np.max(wlv)))
        if tv[-1] > tv[0]:
            ax.set_ylim(float(tv[0]), float(tv[-1]))
        cv.relayout()

    def drawCorrections(self):
        app, P = self.app, self.P
        ch = P.chirp
        wlI = np.asarray(P.wlIn, float).ravel()
        tI = np.asarray(P.tIn, float).ravel()
        good = np.asarray(ch['good'], bool)
        vr = np.asarray(P.rowOK, bool) & ~np.asarray(P.maskRows, bool)
        Db = self.runD - np.asarray(P.bg, float)[:, None]
        lo, hi, cRaw = tI[0], tI[-1], np.nan
        if P.hasT0 and np.any(good):
            trg = np.asarray(ch['ridgeT'], float)[good]
            wdt = np.nanmax([P.fwGuess, 2 * ch['dtZ'], 0.0])
            if wdt > 0:
                cRaw = float(np.median(trg))
                lo = max(tI[0], float(np.min(trg)) - 8 * wdt)
                hi = min(tI[-1], float(np.max(trg)) + 12 * wdt)
        cI = (tI >= lo) & (tI <= hi)
        if np.count_nonzero(cI) < 3:
            cI = np.ones(tI.size, bool)
        ax = self.axRaw
        self.heat(self.cvRaw, ax, self.caxRaw, wlI, tI[cI], Db[:, cI])
        hLeg = []
        if np.any(good):
            hLeg += ax.plot(wlI[good], np.asarray(ch['ridgeT'], float)[good],
                            linestyle='none', marker='.', markersize=4,
                            color='k', label='steepest rise')
        if ch['status'] in ('applied', 'found'):
            hLeg += ax.plot(wlI[vr], np.asarray(ch['t0'], float)[vr], '-',
                            color=(0, 0.6, 0), linewidth=1.8,
                            label='fitted time zero')
        if hLeg:
            legend(ax, hLeg)
        setTitle(ax, 'Data before the corrections, near time zero',
                 'Before the corrections, near time zero',
                 'Before the corrections')
        tO = np.asarray(P.t, float).ravel()
        lo2, hi2 = tO[0], tO[-1]
        if np.isfinite(cRaw) and np.isfinite(P.t0Guess):
            lo2 = max(tO[0], P.t0Guess + (lo - cRaw))
            hi2 = min(tO[-1], P.t0Guess + (hi - cRaw))
        cO = (tO >= lo2) & (tO <= hi2)
        if np.count_nonzero(cO) < 3:
            cO = np.ones(tO.size, bool)
        ax = self.axCor
        self.heat(self.cvCor, ax, self.caxCor, P.wl, tO[cO], P.D[:, cO])
        if np.isfinite(P.t0Guess):
            ax.axhline(P.t0Guess, color='k', linestyle=':', linewidth=0.8)
        if P.changed:
            setTitle(ax, 'After the proposed corrections',
                     'After the corrections')
        else:
            setTitle(ax, 'No correction proposed')

        ax = self.axBg
        self._clear(ax)
        if P.bgStatus == 'applied':
            ax.plot(wlI, P.bg, '-', color=lineColor(0), linewidth=1.2)
            setTitle(ax, 'Background: mean of the first %d delays' % P.nBg,
                     'Background (first %d delays)' % P.nBg, 'Background')
        else:
            setTitle(ax, 'No background subtracted')
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel(self.app.sigName())
        ax.set_xlim(float(np.min(wlI)), float(np.max(wlI)))
        sciTicks(ax, 'y')

        ax = self.axNz
        self._clear(ax)
        noise = np.asarray(P.noise, float)
        hLeg = ax.plot(wlI, noise, '-', color=lineColor(0), linewidth=1.0,
                       label='noise')
        with np.errstate(invalid='ignore', divide='ignore'):
            thr = P.opt['maskK'] * noise / np.asarray(P.noiseRatio, float)
        if np.any(np.isfinite(thr)):
            hLeg += ax.plot(wlI, thr, '--', color=(0.5, 0.5, 0.5),
                            label='%.3g x neighbors' % P.opt['maskK'])
        keepW = np.asarray(P.keepW, bool)
        maskRows = np.asarray(P.maskRows, bool)
        mk = maskRows & keepW
        if np.any(mk):
            hLeg += ax.plot(wlI[mk], noise[mk], 'o', color=(0.85, 0.2, 0.2),
                            markerfacecolor='none', label='masked')
        kn = np.asarray(P.noisy, bool) & ~maskRows
        if np.any(kn):
            hLeg += ax.plot(wlI[kn], noise[kn], 's', color=(0.2, 0.55, 0.3),
                            markerfacecolor='none', markersize=5,
                            label='noisy, kept (clear signal)')
        if np.any(~keepW):
            hLeg += ax.plot(wlI[~keepW], noise[~keepW], 'x',
                            color=(0.3, 0.3, 0.3), label='cropped')
        if np.any(noise[np.isfinite(noise)] > 0):
            logY(ax)
        ax.set_xlim(float(np.min(wlI)), float(np.max(wlI)))
        legend(ax, hLeg)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('Noise (same unit as the data)')
        setTitle(ax, 'Noise of every channel')
        for cv in (self.cvRaw, self.cvCor, self.cvBg, self.cvNz):
            cv.draw_idle()

    def drawComponents(self):
        app, P, S = self.app, self.P, self.S
        axs = (self.axSV, self.axRms, self.axDADS, self.axEADS)
        for ax in axs:
            self._clear(ax)
        cvs = (self.cvSV, self.cvRms, self.cvDADS, self.cvEADS)
        if S is None:
            for ax in axs:
                ax.grid(False)
                setTitle(ax, 'No summary')
            self.drawKinetics()
            for cv in cvs:
                cv.draw_idle()
            return
        tuA = app.tUnitAx()
        c0 = lineColor(0)
        ax = self.axSV
        kS = min(15, S.sv.size)
        ax.plot(np.arange(1, kS + 1), S.sv[:kS], 'o-', color=c0,
                markerfacecolor=c0, markersize=5)
        hasFl = bool(np.isfinite(S.svFloor) and S.svFloor > 0)
        hFl = []
        if hasFl:
            hFl = ax.plot([0.5, kS + 0.5], [5 * S.svFloor] * 2, 'k--',
                          linewidth=1.0, label='5 x noise floor')
        logY(ax)
        ax.set_xlim(0.5, kS + 0.5)
        ax.set_xlabel('Component')
        ax.set_ylabel('Singular value')
        setTitle(ax, 'Singular values')
        if hasFl:
            legend(ax, hFl, shorter=['5 x floor'])

        ax = self.axRms
        tabN = np.asarray(S.tabN, float)
        tabRms = np.asarray(S.tabRms, float)
        ax.plot(tabN, tabRms, 'o-', color=c0, markerfacecolor='none',
                markersize=5)
        nmUse = 'in use' if S.nUse != S.nStar else 'chosen'
        hR = ax.plot([S.nUse], [tabRms[S.nUse - 1]], 'o', color=c0,
                     markerfacecolor=c0, markersize=8, label=nmUse)
        short = [nmUse]             # wording for a narrow plot
        if S.nUse != S.nStar:
            short.append('automatic')
            hR += ax.plot([S.nStar], [tabRms[S.nStar - 1]], 's',
                          color=(0.3, 0.3, 0.3), markerfacecolor='none',
                          markersize=12, markeredgewidth=1.2,
                          label='automatic choice')
        if np.isfinite(S.rmsNoise) and S.rmsNoise > 0:
            hR += ax.plot([0.5, np.max(tabN) + 0.5], [S.rmsNoise] * 2, 'k:',
                          linewidth=1.0, label='noise estimate')
            short.append('noise')
        logY(ax)
        ax.set_xticks(tabN)
        ax.set_xlim(0.5, float(np.max(tabN)) + 0.5)
        ax.set_xlabel('Number of components')
        if S.ca.use > 0:
            ax.set_ylabel('RMS of the fit with artifact terms')
        else:
            ax.set_ylabel('RMS of the global fit')
        setTitle(ax, 'RMS of the fits')
        legend(ax, hR, shorter=short)

        Mo = S.model
        wl = np.asarray(P.wl, float).ravel()
        nT = Mo.tau.size
        nc = Mo.DADS.shape[1]
        hD, hE = [], []
        for k in range(nc):
            col = lineColor(k)
            if k < nT:
                nm = '$\\tau_%d$ = %.3g %s' % (k + 1, Mo.tau[k], tuA)
            else:
                nm = 'non-decaying'
            hD += self.axDADS.plot(wl, Mo.DADS[:, k], '-', color=col,
                                   linewidth=1.2, label=nm)
            hE += self.axEADS.plot(wl, Mo.EADS[:, k], '-', color=col,
                                   linewidth=1.2, label='EADS %d' % (k + 1))
        for k in range(Mo.caK):
            hD += self.axDADS.plot(wl, Mo.caSpec[:, k], '--',
                                   color=CA_COLORS[k], linewidth=1.0,
                                   label=CA_NAMES[k])
        for ax, h in ((self.axDADS, hD), (self.axEADS, hE)):
            ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
            ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
            ax.set_xlabel(app.xAxLabel())
            ax.set_ylabel('Amplitude')
            sciTicks(ax, 'y')
            legend(ax, h)
        fn = self.fitName()
        if Mo.caK > 0:
            setTitle(self.axDADS, 'DADS of %s  (dashed: artifact terms)' % fn,
                     'DADS of %s' % fn, 'DADS')
        else:
            setTitle(self.axDADS, 'DADS of %s' % fn, 'DADS')
        setTitle(self.axEADS, 'EADS of %s  (1 → 2 → ...)' % fn,
                 'EADS of %s' % fn, 'EADS')
        self.drawKinetics()
        for cv in cvs:
            cv.draw_idle()

    def kinPick(self, wv, Dk):
        """Channels (indices into wv) whose kinetics are shown."""
        vals = scanNumbers(self.edKin.text().strip())
        isAuto = len(vals) == 0
        if isAuto:
            amp = np.max(np.abs(Dk), axis=1)
            ed = np.linspace(np.min(wv), np.max(wv), 6)
            for b in range(5):
                inB = np.flatnonzero((wv >= ed[b]) & (wv <= ed[b + 1]))
                if inB.size == 0:
                    continue
                vals.append(float(wv[inB[int(np.argmax(amp[inB]))]]))
        span = float(np.max(wv) - np.min(wv))
        idx = []
        for v in vals:
            if not np.isfinite(v):
                continue
            i1 = int(np.argmin(np.abs(wv - v)))
            if abs(wv[i1] - v) > 0.02 * span:
                continue
            if i1 in idx:
                continue
            if isAuto and idx and \
                    np.any(np.abs(wv[idx] - wv[i1]) < 0.03 * span):
                continue
            idx.append(i1)
            if len(idx) >= 6:
                break
        if isAuto:
            idx.sort()
        return idx

    def _layKin(self):
        cv = self.cvKin
        W, H = cv.sizePt()
        # (B and T as the plots beside and below)
        L, B, T, gap = 58.0, 36.0, 24.0, 6.0
        w = W - L - KIN_LEGEND_PT
        h = H - B - T
        if self.kMode == 'split':
            wL = (w - gap) / 2.5
            cv.place(self.axKinL, L, B, wL, h)
            cv.place(self.axKinR, L + wL + gap, B, w - gap - wL, h)
        else:
            cv.place(self.axKinL, L, B, w, h)
            cv.place(self.axKinR, L, B, w, h)

    def drawKinetics(self):
        app, P, S = self.app, self.P, self.S
        axL, axR = self.axKinL, self.axKinR
        for ax in (axL, axR):
            self._clear(ax)
        self.kinIdx = []
        if S is None:
            self.kMode = 'lin'
            axR.set_visible(False)
            axL.set_visible(True)
            axL.grid(False)
            self.lblKin.setTexts(IDLE if P is None else 'No summary')
            self.cvKin.relayout()
            self.cvKin.draw_idle()
            return
        Mo = S.model
        wv = np.asarray(S.wl, float).ravel()
        tk = np.asarray(S.t, float).ravel()
        n = tk.size
        Dk = P.D[np.ix_(S.rows, S.cols)]
        idx = self.kinPick(wv, Dk)
        tz = 0.0
        if Mo.irf == 'numerical' and np.isfinite(Mo.t0):
            tz = Mo.t0
        hit = np.flatnonzero(tk >= tz)
        j0 = int(hit[0]) if hit.size else 0
        jS = None
        if n >= 8 and j0 <= n - 3:
            d = np.diff(tk)
            d0 = np.median(d[j0:min(j0 + 5, n - 1)])
            big = np.flatnonzero(d[j0:] > 3 * d0)
            if big.size:
                jS = int(big[0]) + j0
        kMode = 'lin'
        if jS is not None and tk[jS] > 0 and n - (jS + 1) >= 5 and \
                tk[-1] >= 5 * tk[jS]:
            kMode = 'split'
            if jS + 1 < 5 and tk[0] > 0:
                kMode = 'log'
        if kMode == 'split':
            iL, iR = np.arange(0, jS + 1), np.arange(jS, n)
        elif kMode == 'log':
            iL, iR = np.zeros(0, int), np.arange(n)
        else:
            iL, iR = np.arange(n), np.zeros(0, int)
        self.kMode = kMode
        axL.set_visible(iL.size > 0)
        axR.set_visible(iR.size > 0)
        self.cvKin.relayout()
        axT = axR if iR.size else axL
        if not idx:
            self.lblKin.setTexts('No channel at these positions',
                                 'No channel there')
            self.cvKin.draw_idle()
            return
        self.kinIdx = idx
        hK = []
        yLo, yHi = np.inf, -np.inf
        for q, i1 in enumerate(idx):
            col = lineColor(q)
            nm = '%.4g %s' % (wv[i1], app.xUnitAx())
            h = None
            for ax, ii in ((axL, iL), (axR, iR)):
                if ii.size == 0:
                    continue
                ax.plot(tk[ii], Dk[i1, ii], 'o', color=col,
                        markerfacecolor='none', markersize=3.5)
                h = ax.plot(tk[ii], Mo.fit[i1, ii], '-', color=col,
                            linewidth=1.2, label=nm)[0]
            hK.append(h)
            yLo = min(yLo, float(np.min(Dk[i1, :])),
                      float(np.min(Mo.fit[i1, :])))
            yHi = max(yHi, float(np.max(Dk[i1, :])),
                      float(np.max(Mo.fit[i1, :])))
        pad = 0.05 * (yHi - yLo)
        if not (pad > 0):
            pad = 1.0
        yl = (yLo - pad, yHi + pad)
        tLab = 'Delay time (%s)' % app.tUnitAx()
        if iL.size:
            axL.axhline(0, color='k', linestyle=':', linewidth=0.8)
            axL.set_ylim(*yl)
            if tk[iL[-1]] > tk[iL[0]]:
                axL.set_xlim(tk[iL[0]], tk[iL[-1]])
            axL.set_xlabel(tLab)
            axL.set_ylabel(self.app.sigName())
            sciTicks(axL, 'y')
        if iR.size:
            axR.axhline(0, color='k', linestyle=':', linewidth=0.8)
            axR.set_xscale('log')
            axR.set_ylim(*yl)
            axR.set_xlim(tk[iR[0]], tk[iR[-1]])
            pw = np.arange(int(np.ceil(np.log10(tk[iR[0]]))),
                           int(np.floor(np.log10(tk[iR[-1]]))) + 1)
            if pw.size:
                axR.set_xticks(10.0 ** pw)
                axR.set_xticklabels(['%g' % v for v in 10.0 ** pw])
                axR.minorticks_off()
            if iL.size == 0:
                axR.set_xlabel(tLab)
                axR.set_ylabel(self.app.sigName())
                sciTicks(axR, 'y')
            else:
                axR.tick_params(labelleft=False)
                axR.set_xlabel('log scale')
                axL.set_xlabel(tLab + ', linear')
        legend(axT, hK, loc='upper left', bbox_to_anchor=(1.02, 1.0),
               borderaxespad=0.0)
        byHand = '' if S.nUse == S.nStar else ' with %d' % S.nUse
        self.lblKin.setTexts('dots: data,  lines: ' + self.fitName(),
                             'dots: data,  lines: fit' + byHand,
                             'lines: fit' + byHand)
        self.cvKin.draw_idle()

    def drawResidual(self):
        app, S = self.app, self.S
        for ax in (self.axPow, self.axProf):
            self._clear(ax)
        cvs = (self.cvRes, self.cvPow, self.cvProf)
        if S is None:
            self._mapOff(self.cvRes, self.axRes, self.caxRes)
            for ax in (self.axRes, self.axPow, self.axProf):
                ax.grid(False)
                setTitle(ax, 'No summary')
            for cv in cvs:
                cv.draw_idle()
            return
        O = S.osc
        Rv = S.model.resid
        tv = np.asarray(S.t, float).ravel()
        fn = self.fitName()
        if np.all(np.isfinite(O.win)):
            cR = (tv >= O.win[0]) & (tv <= O.win[1])
            ttl = ('Residual of %s, evenly sampled part' % fn,
                   'Residual of %s' % fn, 'Residual')
        else:
            cR = np.ones(tv.size, bool)
            ttl = ('Residual of %s' % fn, 'Residual')
        if np.count_nonzero(cR) < 3:
            cR = np.ones(tv.size, bool)
        self.heat(self.cvRes, self.axRes, self.caxRes, S.wl, tv[cR],
                  Rv[:, cR])
        setTitle(self.axRes, *ttl)
        if np.size(O.freq) == 0:
            msg = {'short': 'Too few evenly spaced delays',
                   'skipped': 'Oscillations not analyzed'}.get(
                       O.status, 'No oscillation analysis')
            for ax in (self.axPow, self.axProf):
                ax.grid(False)
                setTitle(ax, msg)
            for cv in cvs:
                cv.draw_idle()
            return
        ax = self.axPow
        hP = ax.plot(O.freq, O.power, '-', color=lineColor(0), linewidth=1.0,
                     label='residual')
        hP += ax.plot(O.freq, O.floor, '--', color=(0.5, 0.5, 0.5),
                      linewidth=1.0, label='local floor')
        if np.size(O.pkIdx):
            hP += ax.plot(O.pkF, O.power[O.pkIdx], 'v',
                          color=(0.85, 0.2, 0.2),
                          markerfacecolor=(0.85, 0.2, 0.2), label='detected')
        logY(ax)
        ax.set_xlim(0, O.fMax)
        legend(ax, hP)
        ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        ax.set_ylabel('|FFT|$^2$ summed over the wavelengths '
                      '(noise-weighted)')
        setTitle(ax, 'Spectrum of the residual')
        ax = self.axProf
        if np.size(O.pkIdx) == 0:
            ax.grid(False)
            setTitle(ax, 'No oscillation detected')
            for cv in cvs:
                cv.draw_idle()
            return
        wv = np.asarray(S.wl, float).ravel()
        hQ = []
        for q, kq in enumerate(O.pkIdx):
            pr = O.Pmap[:, kq]
            pm = float(np.max(pr))
            if not (pm > 0):
                pm = 1.0
            hQ += ax.plot(wv, pr / pm, '-', color=lineColor(q),
                          linewidth=1.2,
                          label='%.0f cm$^{-1}$' % O.pkF[q])
        ax.set_xlim(float(np.min(wv)), float(np.max(wv)))
        legend(ax, hQ)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel('|FFT|$^2$ (each scaled to 1)')
        setTitle(ax, 'Where each oscillation is strong')
        for cv in cvs:
            cv.draw_idle()
