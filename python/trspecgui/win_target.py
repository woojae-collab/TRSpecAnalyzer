"""Target Analysis window. Port of openTargetWindow."""
import copy
import os
import re
import traceback

import numpy as np
from matplotlib.patches import FancyBboxPatch
from PySide6 import QtCore, QtWidgets

from trspec._util import Result
from trspec.kinetics import ca_basis, ca_options
from trspec.target import (_str2double, fit_target_analysis,
                           parse_target_scheme, target_scheme_text)

from .diagram import targetDiagramLayout
from .state import ui_str
from .textfit import legend, setTitle
from .widgets import (AX_BG, FONT_PT, BlitManager, DataTable, NumEdit,
                      PlotCanvas, VScroll, button, colorButton, flowbox,
                      label, lineColor, monoFont, panel, runBusy, sciTicks,
                      sidePanel, styleAxes, textArea, tip)
from .win_base import ChildWindow, absMax, fitSize, mapLayout
from .win_kinfit import irfLimitNote, writeTable

RUN_COLOR = (0.60, 0.25, 0.45)

PRESET_KEYS = ['seq3', 'seq3gs', 'seq4', 'par3', 'branch', 'branchLoss',
               'rev', 'exISC', 'exET', 'exSF', 'custom']
PRESET_NAMES = [
    'Sequential:  A -> B -> C',
    'Sequential, all decay:  A -> B -> C -> GS',
    'Sequential, 4 compartments:  A -> B -> C -> D',
    'Parallel:  A, B, C decay independently',
    'Branching:  A -> B,  A -> C',
    'Branching with loss:  A -> B -> GS,  A -> GS',
    'Reversible:  A <-> B -> GS',
    'Example - intersystem crossing:  S1 -> T1, S1 -> S0',
    'Example - energy transfer:  D -> A, D -> GS',
    'Example - singlet fission:  S1 -> TT -> 2 T',
    'Custom']
PRESET_TXT = [
    ['A -> B', 'B -> C'],
    ['A -> B', 'B -> C', 'C -> GS'],
    ['A -> B', 'B -> C', 'C -> D'],
    ['A -> GS', 'B -> GS', 'C -> GS', 'init: A = 1, B = 1, C = 1'],
    ['A -> B', 'A -> C', 'B -> GS', 'C -> GS'],
    ['A -> B', 'B -> GS', 'A -> GS'],
    ['A <-> B : kf, kb', 'B -> GS'],
    ['S1 -> T1 : kisc', 'S1 -> S0 : kS', 'T1 -> S0 : kT'],
    ['D -> A : ket', 'D -> GS : kD', 'A -> GS : kA'],
    ['S1 -> TT', 'TT -> 2 T'],
    []]
SYNTAX_INTRO = ('One statement per line, or several on a line separated by '
                '";".\nNames of compartments and rate parameters are free. '
                'Text after % is a comment.')
SYNTAX_ROWS = [
    ('A -> B', 'first-order step with its own rate parameter'),
    ('A -> B : k1', '... with a named rate parameter'),
    ('B -> C : k1', 'the same name = the same parameter (shared)'),
    ('B -> C : 2*k1', 'rate tied to a parameter (also k1/2)'),
    ('A -> B -> C', 'a chain is split into single steps'),
    ('A <-> B : kf, kb', 'reversible step: two rates'),
    ('B -> 2 C', 'two C per B (stoichiometric factor)'),
    ('A -> GS', 'loss to the ground state; GS, S0 and 0 are the dark sink '
                '(no spectrum)'),
    ('species: A, B, C', 'declare compartments (sets the order; may name '
                         'one without steps)'),
    ('init: A = 0.7, B = 0.3', 'initial populations (default: first '
                               'compartment = 1)'),
    ('zero: B 700 900', 'spectrum of B = 0 from 700 to 900'),
    ('equal: C = B', 'spectrum of C equals that of B'),
    ('equal: C = 0.5*B 600 700', '... proportional, in a range'),
    ('nonneg', 'all spectra >= 0'),
    ('nonneg: A, B', '... only those of A and B')]
SYNTAX_NOTE = ('Ranges are in the units of the spectral axis. Only '
               'first-order steps are supported.')
SYNTAX_HELP = SYNTAX_INTRO + '\n' + ''.join(
    '\n  %s    %s' % r for r in SYNTAX_ROWS) + '\n\n' + SYNTAX_NOTE
IDLE_TEXT = ['Choose or write a scheme, set the starting time',
             'constants, then click "Run fit".']
CA_NAMES = ['Gaussian', '1st derivative', '2nd derivative']


def sigRound(x, n):
    """x with n significant digits (round(x, n, 'significant'))."""
    x = float(x)
    if not np.isfinite(x) or x == 0:
        return x
    return float('%.*g' % (int(n), x))


def irfLine(name, val, sd, isFixed, u):
    if isFixed:
        return '%s = %.4g %s (fixed)' % (name, val, u)
    if np.isinf(sd):
        return '%s = %.4g %s  [not determined]' % (name, val, u)
    if np.isfinite(sd):
        return '%s = %.4g +/- %.2g %s' % (name, val, sd, u)
    return '%s = %.4g %s' % (name, val, u)


class StatusLabel(QtWidgets.QLabel):
    """Two lines of status text. Its height does not depend on its
    width: a label that says so makes Qt give every control of the
    column its preferred height, and the column then needs a scroll
    bar on a low screen instead of using its minimum heights."""

    def __init__(self, text):
        super().__init__(text)
        self.setWordWrap(True)
        self.setFixedHeight(32)
        self.setStyleSheet('color: rgb(51,51,51);')
        self.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft |
                          QtCore.Qt.AlignmentFlag.AlignTop)

    def hasHeightForWidth(self):
        return False

    def heightForWidth(self, w):
        return -1


class SubWindow(QtWidgets.QWidget):
    """A small window that belongs to another one; onClosed is called
    when it goes."""

    def __init__(self, title, w, h, onClosed=None):
        super().__init__(None)
        self.setWindowTitle(title)
        self.resize(*fitSize(w, h))
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._onClosed = onClosed

    def closeEvent(self, ev):
        if self._onClosed is not None:
            try:
                self._onClosed()
            except Exception:
                pass
        ev.accept()


class TargetWindow(ChildWindow):
    KEY = 'target'

    def __init__(self, app):
        super().__init__(app, 'Target Analysis', 1480, 980)
        u = app.tUnitTxt()
        self.sch = None
        self.res = None
        self.svdWin = None
        self.synWin = None
        self.diaWin = None
        self.cvDiaWin = None
        self.hWLspec = None
        self.hWLres = None
        self.dragging = None        # None | 'sas' | 'res'
        self._quiet = False         # the program is filling the controls
        gl = QtWidgets.QHBoxLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setSpacing(10)

        # ================= left: scheme and fit setup ===================
        left, lv = sidePanel('Kinetic scheme and fit setup', 470)
        lv.setSpacing(5)
        gl.addWidget(left)

        def hrow(spacing=6):
            r = QtWidgets.QHBoxLayout()
            r.setSpacing(spacing)
            lv.addLayout(r)
            return r

        row = hrow()
        lb = label('Model:')
        lb.setStyleSheet('font-weight: bold;')
        row.addWidget(lb)
        self.ddPreset = QtWidgets.QComboBox()
        for nm, key in zip(PRESET_NAMES, PRESET_KEYS):
            self.ddPreset.addItem(nm, key)
        self.ddPreset.setMinimumWidth(120)
        row.addWidget(self.ddPreset, 1)

        # (heights between a minimum and the MATLAB height: a low
        # window takes the room from here, not from the result box)
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setMinimumHeight(150)
        self.tabs.setMaximumHeight(218)
        lv.addWidget(self.tabs)
        self.txtScheme = QtWidgets.QPlainTextEdit()
        self.txtScheme.setFont(monoFont())
        self.txtScheme.setLineWrapMode(
            QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        self.txtScheme.setPlainText('\n'.join(PRESET_TXT[0]))
        self.txtScheme.setToolTip(SYNTAX_HELP)
        self.tabs.addTab(self.txtScheme, 'Scheme (text)')
        tabK = QtWidgets.QWidget()
        gK = QtWidgets.QVBoxLayout(tabK)
        gK.setContentsMargins(4, 4, 4, 4)
        gK.setSpacing(4)
        self.tblK = QtWidgets.QTableWidget(0, 0)
        self.tblK.verticalHeader().setDefaultSectionSize(24)
        hk = self.tblK.horizontalHeader()
        hk.setMinimumSectionSize(64)
        hk.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.tblK.setToolTip(tip(
            'Column = compartment a step starts from, row = where it goes. '
            'A cell holds the rate of that step: a parameter name (k1), a '
            'tied rate (2*k1), with "[x2]" appended for two products per '
            'step. Empty = no step. The same name in several cells shares '
            'the parameter. Last row: initial population.'))
        gK.addWidget(self.tblK, 1)
        rk = QtWidgets.QHBoxLayout()
        rk.setSpacing(5)
        gK.addLayout(rk)
        self.edNewSp = QtWidgets.QLineEdit()
        self.edNewSp.setToolTip('Name of a new compartment.')
        rk.addWidget(self.edNewSp, 1)
        self.btnAddSp = button('Add compartment')
        rk.addWidget(self.btnAddSp)
        self.ddRemSp = QtWidgets.QComboBox()
        self.ddRemSp.addItem('-')
        rk.addWidget(self.ddRemSp, 1)
        self.btnRemSp = button('Remove')
        rk.addWidget(self.btnRemSp)
        self.tabs.addTab(tabK, 'K matrix')
        self.cvDia = PlotCanvas()
        self.cvDia.setMinimumSize(120, 80)
        self.axDia = self.cvDia.figure.add_axes([0, 0, 1, 1])
        self.cvDia.setToolTip('Click to enlarge.')
        self.tabs.addTab(self.cvDia, 'Diagram')

        row = hrow()
        self.btnApply = button(
            'Apply scheme',
            tooltip='Read the scheme text and list its rate parameters in '
                    'the table below. Parameters that were already listed '
                    'keep their values. "Run fit" does this as well.')
        row.addWidget(self.btnApply, 1)
        self.btnFromGA = button(
            ui_str('\\tau from GA'),
            tooltip='Copy the time constants of the last Global Analysis '
                    'fit into the table (fastest first, in table order), '
                    'together with its t0 and FWHM.')
        row.addWidget(self.btnFromGA, 1)
        self.btnHelp = button('Syntax...')
        row.addWidget(self.btnHelp, 1)

        self.tblPar = DataTable(
            ui_str(['Rate', '\\tau_{init} (%s)' % u, 'Fixed', '\\tau min',
                    '\\tau max', 'Used in']),
            ['text', 'num', 'bool', 'text', 'text', 'text'],
            editable=[False, True, True, True, True, False],
            widths=[58, 92, 52, 62, 62, 0])
        self.tblPar.verticalHeader().setVisible(False)
        self.tblPar.setMinimumHeight(102)
        self.tblPar.setMaximumHeight(150)
        self.tblPar.setToolTip(tip(ui_str(
            'One row per rate parameter k; the table shows its time '
            'constant \\tau = 1/k. "\\tau min" and "\\tau max" are optional '
            'limits (leave empty for none).')))
        lv.addWidget(self.tblPar)

        t0Def = app.gaResultT0 if app.gaResultT0 is not None else app.gaT0
        fwDef = app.gaResultFWHM if app.gaResultFWHM is not None \
            else app.gaFWHM
        if t0Def is None or not np.isfinite(t0Def):
            t0Def = 0.0
        if fwDef is None or not np.isfinite(fwDef) or fwDef <= 0:
            fwDef = 0.1
        row = hrow()
        lb = label(ui_str('IRF t_0 (%s):' % u))
        lb.setMinimumWidth(110)
        row.addWidget(lb)
        self.edT0 = NumEdit(t0Def)
        row.addWidget(self.edT0, 1)
        self.cbT0fix = QtWidgets.QCheckBox('Fixed')
        self.cbT0fix.setChecked(True)
        row.addWidget(self.cbT0fix)
        row = hrow()
        lb = label('IRF FWHM (%s):' % u)
        lb.setMinimumWidth(110)
        row.addWidget(lb)
        self.edFW = NumEdit(max(fwDef, 1e-6), limits=(1e-6, np.inf))
        row.addWidget(self.edFW, 1)
        self.cbFWfix = QtWidgets.QCheckBox('Fixed')
        self.cbFWfix.setChecked(True)
        row.addWidget(self.cbFWfix)
        row = hrow()
        lb = label('IRF model:')
        lb.setMinimumWidth(110)
        row.addWidget(lb)
        self.ddIRF = QtWidgets.QComboBox()
        self.ddIRF.addItem('Gaussian IRF', 'irf')
        self.ddIRF.addItem('None (IRF-free)', 'none')
        self.ddIRF.setToolTip(tip(
            'Gaussian IRF: every decay mode is the closed-form exponential '
            'convolved with the IRF.\nNone: pure kinetics starting at t0, '
            'no convolution; t0 and FWHM are not fitted. Intended for a '
            'Fit t range that starts after the IRF.'))
        row.addWidget(self.ddIRF, 1)

        caDef = ca_options(None)
        gca = app.gaResultCA
        if isinstance(gca, dict) and gca.get('order', 0) > 0:
            caDef['order'] = int(gca['order'])
            if gca.get('mode') != 'irf':
                caDef['mode'] = 'fixed'
                caDef['fwhm'] = float(gca['fwhm'])
        if not (np.isfinite(caDef['fwhm']) and caDef['fwhm'] > 0):
            caDef['fwhm'] = fwDef
        row = hrow()
        lb = label('Coh. artifact:')
        lb.setMinimumWidth(110)
        row.addWidget(lb)
        self.ddCA = QtWidgets.QComboBox()
        for k, nm in enumerate(['None', '1 term: Gaussian',
                                '2 terms: + 1st derivative',
                                '3 terms: + 2nd derivative']):
            self.ddCA.addItem(nm, k)
        self.ddCA.setCurrentIndex(int(caDef['order']))
        self.ddCA.setToolTip(tip(
            'Signal that exists only while pump and probe overlap. It is '
            'fitted with up to three extra terms at time zero, each with a '
            'free spectrum: a Gaussian and its first and second time '
            'derivative. The spectral constraints of the scheme do not '
            'apply to them.\nWith these terms, steps faster than the IRF '
            'width are less well determined: switch them on only when the '
            'data show an artifact.'))
        row.addWidget(self.ddCA, 1)
        row = hrow()
        self.lblCAW = label('Artifact width:')
        self.lblCAW.setMinimumWidth(110)
        row.addWidget(self.lblCAW)
        self.ddCAW = QtWidgets.QComboBox()
        for nm, key in (('Same as IRF', 'irf'), ('Own, fitted', 'free'),
                        ('Own, fixed', 'fixed')):
            self.ddCAW.addItem(nm, key)
        self.ddCAW.setCurrentIndex(['irf', 'free', 'fixed'].index(
            caDef['mode']))
        self.ddCAW.setToolTip(tip(
            'Width (FWHM) of the artifact terms: the IRF FWHM, or a width '
            'of their own, fitted (within a factor of 5 of the IRF FWHM) or '
            'fixed at the value on the right. With a width of their own '
            'the IRF FWHM is hardly determined by the fit: keep it fixed.'))
        row.addWidget(self.ddCAW, 1)
        self.edCAW = NumEdit(
            max(caDef['fwhm'], 1e-6), limits=(1e-6, np.inf), width=70,
            tooltip='Own width of the artifact terms (FWHM, %s): starting '
                    'value or fixed value.' % u)
        row.addWidget(self.edCAW)
        row = hrow(4)
        row.addWidget(label('Fit t range (%s):' % u))
        self.edFitTmin = NumEdit(
            float(np.min(app.delay)), width=70,
            tooltip='Delay points outside this range are left out of the '
                    'fit (the model is still drawn on the full axis).')
        row.addWidget(self.edFitTmin)
        row.addWidget(label('~'))
        self.edFitTmax = NumEdit(
            float(np.max(app.delay)), width=70,
            tooltip='Upper bound of the fitted delay range.')
        row.addWidget(self.edFitTmax)
        self.btnFull = button('Full', self.setFitTFull,
                              tooltip='Reset to the full delay range.')
        row.addWidget(self.btnFull)
        row.addStretch(1)

        row = hrow()
        self.btnRun = colorButton('Run fit', RUN_COLOR)
        self.btnRun.setStyleSheet(self.btnRun.styleSheet()
                                  + 'QPushButton { font-weight: bold; }')
        row.addWidget(self.btnRun, 12)
        self.btnUse = button(
            ui_str('Fitted \\rightarrow start'),
            tooltip='Copy the fitted time constants, t0 and FWHM into the '
                    'starting values, e.g. to refit with another choice of '
                    'fixed parameters.')
        self.btnUse.setEnabled(False)
        row.addWidget(self.btnUse, 10)
        self.btnReset = button('Reset')
        row.addWidget(self.btnReset, 8)

        row = hrow(5)
        lb = label('Export:')
        lb.setStyleSheet('font-weight: bold;')
        row.addWidget(lb)
        self.btnExpSAS = button('SAS')
        self.btnExpDAS = button('DAS')
        self.btnExpPop = button('Populations')
        self.btnExpFit = button('Fit matrix')
        self.btnExpPar = button('Parameters')
        for b, s in ((self.btnExpSAS, 10), (self.btnExpDAS, 10),
                     (self.btnExpPop, 13), (self.btnExpFit, 12),
                     (self.btnExpPar, 13)):
            b.setEnabled(False)
            row.addWidget(b, s)
        self.lblStatus = StatusLabel('Ready.')
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
        rv = QtWidgets.QVBoxLayout(rin)
        rv.setContentsMargins(4, 0, 4, 2)
        rv.setSpacing(8)
        ro.addWidget(VScroll(rin))
        ctlHost, ctl = flowbox(8)
        rv.addWidget(ctlHost)
        wl = np.asarray(app.wavelength, float)
        self.wlLo, self.wlHi = float(np.min(wl)), float(np.max(wl))
        if not (self.wlHi > self.wlLo):
            self.wlHi = self.wlLo + 1.0
        self.ddSpec = QtWidgets.QComboBox()
        self.ddSpec.addItems(['SAS', 'DAS'])
        self.ddSpec.setMinimumWidth(70)
        self.ddSpec.setToolTip(tip(
            'SAS: spectrum of every compartment of the scheme. DAS: '
            'amplitude spectrum of every decay time of the same fit (what '
            'a global fit with these lifetimes gives).'))
        ctl.addPair(label('Spectra:', right=True), self.ddSpec)
        self.ddScale = QtWidgets.QComboBox()
        self.ddScale.addItems(['Log', 'Linear'])
        self.ddScale.setMinimumWidth(75)
        ctl.addPair(label('Time axis:', right=True), self.ddScale)
        self.edWL = NumEdit(
            min(max(float(app.selWL), self.wlLo), self.wlHi),
            limits=(self.wlLo, self.wlHi), width=70,
            tooltip='Probe position of the kinetics panel. It is marked by '
                    'a line on the spectra panel and on the residual map: '
                    'click or drag on either of them to move it.')
        ctl.addPair(label('Kinetics at %s (%s):' % (
            app.xSymTxt(), app.xUnitTxt()), right=True), self.edWL)
        self.cbNorm = QtWidgets.QCheckBox('Normalize')
        self.cbNorm.setToolTip('Scale every spectrum to a maximum of 1.')
        ctl.addWidget(self.cbNorm)
        self.cbContrib = QtWidgets.QCheckBox('Contributions')
        self.cbContrib.setChecked(True)
        self.cbContrib.setToolTip(tip(
            'Show spectrum x population of every compartment under the '
            'fitted trace.'))
        ctl.addWidget(self.cbContrib)
        self.btnSVD = button(
            'Residual SVD...',
            tooltip='Singular value decomposition of the residual matrix. '
                    'Structure in the first singular vectors is signal the '
                    'model does not describe; pure noise has none.')
        self.btnSVD.setEnabled(False)
        ctl.addWidget(self.btnSVD)

        grid = QtWidgets.QGridLayout()
        grid.setSpacing(8)
        rv.addLayout(grid, 1)
        self.cvSAS = PlotCanvas()
        self.cvPop = PlotCanvas()
        self.cvKin = PlotCanvas()
        self.cvRes = PlotCanvas()
        for cv, (r, c) in ((self.cvSAS, (0, 0)), (self.cvPop, (0, 1)),
                           (self.cvKin, (1, 0)), (self.cvRes, (1, 1))):
            cv.setMinimumHeight(230)
            grid.addWidget(cv, r, c)
        grid.setRowStretch(0, 1)
        grid.setRowStretch(1, 1)
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        self.axSAS = self.cvSAS.figure.add_subplot(111)
        self.axPop = self.cvPop.figure.add_subplot(111)
        self.axKin = self.cvKin.figure.add_subplot(111)
        self.axRes = self.cvRes.figure.add_subplot(111)
        self.caxRes = self.cvRes.figure.add_axes([0.9, 0.1, 0.02, 0.8])

        def lay1(cv, ax):
            # (the three plots without a color bar do not leave the
            # room of the bar of the residual map free: the room right
            # of them is about as wide as the blank left of the label
            # of their y-axis)
            W, H = cv.sizePt()
            cv.place(ax, 58.0, 34.0, W - 58.0 - 22.0, H - 34.0 - 22.0)

        self.cvSAS.layoutFcn = lambda: lay1(self.cvSAS, self.axSAS)
        self.cvPop.layoutFcn = lambda: lay1(self.cvPop, self.axPop)
        self.cvKin.layoutFcn = lambda: lay1(self.cvKin, self.axKin)
        self.cvRes.layoutFcn = self._layRes
        self.cvDia.layoutFcn = self.drawDiagrams
        self.blSAS = BlitManager(self.cvSAS)
        self.blRes = BlitManager(self.cvRes)
        for cv, which in ((self.cvSAS, 'sas'), (self.cvRes, 'res')):
            cv.mpl_connect('button_press_event',
                           lambda ev, w=which: self.onPickDown(w, ev))
            cv.mpl_connect('motion_notify_event',
                           lambda ev, w=which: self.onPickMotion(w, ev))
            cv.mpl_connect('button_release_event',
                           lambda ev: self.onPickUp())
        self.cvDia.mpl_connect(
            'button_press_event',
            lambda ev: self.safe(self.openDiagramWindow, 'Diagram')
            if ev.button == 1 else None)       # (right click: the menu)
        self.cvDia.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)

        # ---- what the controls do --------------------------------------
        self.ddPreset.currentIndexChanged.connect(
            lambda *_: self._user(self.onPreset, 'Model'))
        self.txtScheme.textChanged.connect(self.onSchemeEdited)
        self.tabs.currentChanged.connect(
            lambda *_: self._user(self.onTabChange, 'Scheme'))
        self.tblK.itemChanged.connect(self._onKItem)
        self.btnAddSp.clicked.connect(
            lambda *_: self.safe(self.addCompartment, 'Add compartment'))
        self.btnRemSp.clicked.connect(
            lambda *_: self.safe(self.removeCompartment,
                                 'Remove compartment'))
        self.btnApply.clicked.connect(
            lambda *_: self.safe(self.applyScheme, 'Apply scheme'))
        self.btnFromGA.clicked.connect(
            lambda *_: self.safe(self.fromGA, 'Copy from GA'))
        self.btnHelp.clicked.connect(
            lambda *_: self.safe(self.showSyntaxHelp, 'Syntax'))
        self.ddIRF.currentIndexChanged.connect(
            lambda *_: self.onIRFChange())
        self.ddCA.currentIndexChanged.connect(
            lambda *_: self.syncCAEnable())
        self.ddCAW.currentIndexChanged.connect(
            lambda *_: self.syncCAEnable())
        self.btnRun.clicked.connect(
            lambda *_: self.safe(self.doRun, 'Run fit'))
        self.btnUse.clicked.connect(
            lambda *_: self.safe(self.useFitted, 'Fitted to start'))
        self.btnReset.clicked.connect(
            lambda *_: self.safe(self.doReset, 'Reset'))
        self.btnExpSAS.clicked.connect(
            lambda *_: self.safe(self.exportSAS, 'Export SAS'))
        self.btnExpDAS.clicked.connect(
            lambda *_: self.safe(self.exportDAS, 'Export DAS'))
        self.btnExpPop.clicked.connect(
            lambda *_: self.safe(self.exportPop, 'Export populations'))
        self.btnExpFit.clicked.connect(
            lambda *_: self.safe(self.exportFit, 'Export fit'))
        self.btnExpPar.clicked.connect(
            lambda *_: self.safe(self.exportPar, 'Export parameters'))
        self.btnSVD.clicked.connect(
            lambda *_: self.safe(self.showResidualSVD, 'Residual SVD'))
        self.edWL.valueChanged.connect(
            lambda *_: self.safe(self.onWLChange, 'Kinetics'))
        self.ddScale.currentIndexChanged.connect(
            lambda *_: self._user(self.plotTimePanels, 'Time axis'))
        self.ddSpec.currentIndexChanged.connect(
            lambda *_: self._user(self.plotSpec, 'Spectra'))
        self.cbNorm.toggled.connect(
            lambda *_: self._user(self.plotSpec, 'Spectra'))
        self.cbContrib.toggled.connect(
            lambda *_: self._user(self.plotKin, 'Kinetics'))

        # ---- settings and result of an earlier visit -------------------
        self._quiet = True
        st = app.tgSet
        if isinstance(st, dict) and 'schemeTxt' in st:
            try:
                self._restoreSettings(st)
            except Exception:
                self.lblStatus.setText(
                    'Stored settings could not be restored.')
        self.onIRFChange()
        try:
            self.applyScheme()
        except Exception as err:
            self.lblStatus.setText('Scheme: %s' % err)
        sv = app.tgSaved
        if isinstance(sv, dict) and isinstance(sv.get('res'), dict):
            try:
                r = self._unpackRes(sv['res'])
                if r is not None:
                    self.res = r
                    self.edWL.setValue(min(max(float(sv['curWL']),
                                               self.wlLo), self.wlHi))
                    self.txtResult.setPlainText('\n'.join(sv['lines']))
                    self.lblStatus.setText('%s  [restored]' % sv['status'])
                    self.setResultButtons(True)
            except Exception:
                self.res = None
                self.lblStatus.setText(
                    'Stored result could not be restored -- run again.')
        self._quiet = False
        try:
            self.plotAll()
        except Exception:
            self.lblStatus.setText(
                'Stored result could not be drawn -- run again.')

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def safe(self, fcn, what):
        """Run fcn; an error is shown instead of being lost (safeCb)."""
        try:
            fcn()
        except Exception as err:
            if not isinstance(err, ValueError):     # (not a refused entry)
                traceback.print_exc()
            self.alert('"%s" failed:\n%s' % (what, err), 'Callback error')
            self.lblStatus.setText('Error in %s.' % what)

    def _user(self, fcn, what):
        """A control was changed: by the user (act on it) or by the
        program while it fills the controls (leave it)."""
        if not self._quiet:
            self.safe(fcn, what)

    def _layRes(self):
        mapLayout(self.cvRes, self.axRes, self.caxRes, L=58.0, R=74.0,
                  T=22.0, B=34.0)

    def schemeLines(self):
        return self.txtScheme.toPlainText().split('\n')

    def setSchemeText(self, lines):
        q, self._quiet = self._quiet, True
        try:
            self.txtScheme.setPlainText('\n'.join(lines))
        finally:
            self._quiet = q

    def setPreset(self, key):
        q, self._quiet = self._quiet, True
        try:
            self.ddPreset.setCurrentIndex(PRESET_KEYS.index(key))
        finally:
            self._quiet = q

    def _restoreSettings(self, st):
        app = self.app
        self.setPreset(st['preset'] if st['preset'] in PRESET_KEYS
                       else 'custom')
        txt = st['schemeTxt']
        self.setSchemeText(txt if isinstance(txt, (list, tuple))
                           else str(txt).split('\n'))
        tab = st.get('table')
        if isinstance(tab, (list, tuple)) and \
                all(len(r) == 6 for r in tab):
            self.tblPar.setData([[r[0], float(r[1]), bool(r[2]), r[3], r[4],
                                  r[5]] for r in tab])
        self.edT0.setValue(st['t0'])
        self.edFW.setValue(max(st['fwhm'], 1e-6))
        self.cbT0fix.setChecked(bool(st['t0Fixed']))
        self.cbFWfix.setChecked(bool(st['fwhmFixed']))
        self.ddIRF.setCurrentIndex(0 if st['irf'] == 'irf' else 1)
        if 'ca' in st:
            caSt = ca_options(st['ca'])
            self.ddCA.setCurrentIndex(int(caSt['order']))
            self.ddCAW.setCurrentIndex(['irf', 'free', 'fixed'].index(
                caSt['mode']))
            if np.isfinite(caSt['fwhm']) and caSt['fwhm'] > 0:
                self.edCAW.setValue(caSt['fwhm'])
        lo, hi = self.clampFitRange(st.get('fitTMin'), st.get('fitTMax'))
        self.edFitTmin.setValue(lo)
        self.edFitTmax.setValue(hi)
        self.ddScale.setCurrentText(st['scale'])
        self.cbNorm.setChecked(bool(st['norm']))
        self.cbContrib.setChecked(bool(st['contrib']))
        if st.get('spec') in ('SAS', 'DAS'):
            self.ddSpec.setCurrentText(st['spec'])
        del app

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
    # the scheme
    # ------------------------------------------------------------------
    def onPreset(self):
        k = self.ddPreset.currentIndex()
        if k < 0 or not PRESET_TXT[k]:
            return
        self.setSchemeText(PRESET_TXT[k])
        self.applyScheme()

    def onSchemeEdited(self):
        if not self._quiet:
            self.setPreset('custom')

    def onTabChange(self):
        if self.tabs.currentIndex() in (1, 2):
            self.applyScheme()

    def applyScheme(self):
        app = self.app
        newSch = parse_target_scheme(self.schemeLines())
        old = self.tblPar.data()
        nPar = newSch.n_p
        spanA = float(np.max(app.delay) - np.min(app.delay))
        if not (spanA > 0):
            spanA = 1.0
        tLo = max(2 * self.edFW.value(), 1e-3 * spanA)
        tHi = max(spanA / 3, 2 * tLo)
        defTau = np.logspace(np.log10(tLo), np.log10(tHi), max(nPar, 2))
        rows = []
        for q in range(nPar):
            pn = newSch.params[q]
            use = []
            for stp in newSch.steps:
                if stp['par'] != q:
                    continue
                if stp['mult'] == 1:
                    use.append(stp['label'])
                else:
                    use.append('%s (%s)' % (stp['label'], stp['expr']))
            row = [pn, sigRound(defTau[q], 3), False, '', '',
                   ';  '.join(use)]
            for o in old:
                if len(o) == 6 and o[0] == pn:
                    row[1], row[2], row[3], row[4] = \
                        o[1], bool(o[2]), o[3], o[4]
                    break
            rows.append(row)
        self.tblPar.setData(rows)
        self.sch = newSch
        self.refreshK()
        self.lblStatus.setText(
            'Scheme: %d compartment(s) (%s), %d step(s), %d rate '
            'parameter(s).' % (newSch.n_s, ', '.join(newSch.species),
                               len(newSch.steps), newSch.n_p))
        self.drawDiagrams()

    def kData(self):
        """The K-matrix view as a list of rows of text."""
        t = self.tblK
        return [[t.item(i, j).text() if t.item(i, j) else ''
                 for j in range(t.columnCount())]
                for i in range(t.rowCount())]

    def refreshK(self):
        sch = self.sch
        nC = sch.n_s
        dat = [[''] * nC for _ in range(nC + 2)]
        for sp in sch.steps:
            r = nC if sp['to'] < 0 else sp['to']
            if sp['coef'] == 1:
                dat[r][sp['from']] = sp['expr']
            else:
                dat[r][sp['from']] = '%s [x%g]' % (sp['expr'], sp['coef'])
        for q in range(nC):
            dat[nC + 1][q] = '%g' % sch.c0[q]
        t = self.tblK
        t.blockSignals(True)
        try:
            t.clear()
            t.setColumnCount(nC)
            t.setRowCount(nC + 2)
            t.setHorizontalHeaderLabels(['from ' + n for n in sch.species])
            t.setVerticalHeaderLabels(
                ['to ' + n for n in sch.species]
                + ['to ' + sch.sink_name, 'initial'])
            for i in range(nC + 2):
                for j in range(nC):
                    t.setItem(i, j, QtWidgets.QTableWidgetItem(dat[i][j]))
        finally:
            t.blockSignals(False)
        self.ddRemSp.blockSignals(True)
        self.ddRemSp.clear()
        self.ddRemSp.addItems(list(sch.species))
        self.ddRemSp.setCurrentIndex(nC - 1)
        self.ddRemSp.blockSignals(False)

    def useSchemeStruct(self, s2):
        newTxt = target_scheme_text(s2)
        parse_target_scheme(newTxt)
        self.setSchemeText(newTxt)
        self.setPreset('custom')
        self.applyScheme()

    def _onKItem(self, it):
        r, c, val = it.row(), it.column(), it.text()

        def go():
            try:
                self.safe(lambda: self.onKEdit(r, c, val), 'K matrix')
            except RuntimeError:        # (the window was closed meanwhile)
                pass
        # the edit rebuilds the table: only after the cell editor closed
        QtCore.QTimer.singleShot(0, go)

    def onKEdit(self, r, c, val):
        """A cell of the K-matrix view was edited: row r (compartment
        the step goes to; then the sink, then 'initial'), column c
        (compartment it starts from)."""
        try:
            sch = self.sch
            val = str(val).strip()
            nC = sch.n_s
            s2 = copy.deepcopy(sch)
            s2.c0 = np.array(s2.c0, float)
            if r == nC + 1:
                v0 = 0.0 if val == '' else _str2double(val)
                if not np.isfinite(v0) or v0 < 0:
                    raise ValueError(
                        'The initial population must be a number >= 0.')
                s2.c0[c] = v0
                if not np.any(s2.c0 > 0):
                    raise ValueError(
                        'At least one compartment needs a population > 0.')
            elif r == c:
                raise ValueError('A compartment has no step to itself.')
            else:
                iTo = -1 if r == nC else r
                hit = next((i for i, s in enumerate(s2.steps)
                            if s['from'] == c and s['to'] == iTo), None)
                if val == '':
                    if hit is not None:
                        del s2.steps[hit]
                else:
                    cf, ex = 1.0, val
                    m = re.match(
                        r'^(.*?)\s*\[\s*x\s*([0-9]*\.?[0-9]+)\s*\]$', val)
                    if m:
                        ex = m.group(1).strip()
                        cf = _str2double(m.group(2))
                    if iTo < 0 and cf != 1:
                        raise ValueError('A loss to the ground state takes '
                                         'no "[x..]" factor.')
                    if not np.isfinite(cf) or cf <= 0:
                        raise ValueError(
                            'The factor in "[x..]" must be positive.')
                    stp = {'from': c, 'to': iTo, 'coef': cf, 'par': 0,
                           'mult': 1, 'label': '', 'expr': ex}
                    if hit is None:
                        s2.steps.append(stp)
                    else:
                        s2.steps[hit] = stp
            self.useSchemeStruct(s2)
        except Exception:
            self.refreshK()
            raise

    def addCompartment(self):
        sch = self.sch
        nm = self.edNewSp.text().strip()
        if not nm:
            raise ValueError('Type a name for the new compartment first.')
        if nm in sch.species:
            raise ValueError('A compartment called "%s" already exists.' % nm)
        s2 = copy.deepcopy(sch)
        s2.species = list(s2.species) + [nm]
        s2.c0 = np.append(np.asarray(s2.c0, float), 0.0)
        s2.nonneg = np.append(np.asarray(s2.nonneg, bool),
                              bool(np.all(sch.nonneg)))
        self.useSchemeStruct(s2)
        self.edNewSp.setText('')
        self.lblStatus.setText('Compartment %s added: give it steps in the '
                               'K matrix.' % nm)

    def removeCompartment(self):
        sch = self.sch
        nm = self.ddRemSp.currentText()
        if nm not in sch.species:
            return
        k = list(sch.species).index(nm)
        if sch.n_s < 2:
            raise ValueError('The last compartment cannot be removed.')
        s2 = copy.deepcopy(sch)
        s2.steps = [s for s in s2.steps if s['from'] != k and s['to'] != k]
        for s in s2.steps:
            if s['from'] > k:
                s['from'] -= 1
            if s['to'] > k:
                s['to'] -= 1
        s2.species = [n for i, n in enumerate(s2.species) if i != k]
        s2.c0 = np.delete(np.asarray(s2.c0, float), k)
        s2.nonneg = np.delete(np.asarray(s2.nonneg, bool), k)
        if not np.any(s2.c0 > 0):
            s2.c0[0] = 1.0
        s2.zero = [z for z in s2.zero if z['sp'] != k]
        for z in s2.zero:
            if z['sp'] > k:
                z['sp'] -= 1
        s2.rel = [x for x in s2.rel if x['dep'] != k and x['master'] != k]
        for x in s2.rel:
            if x['dep'] > k:
                x['dep'] -= 1
            if x['master'] > k:
                x['master'] -= 1
        self.useSchemeStruct(s2)
        self.lblStatus.setText('Compartment %s removed, with its steps and '
                               'constraints.' % nm)

    # ------------------------------------------------------------------
    # the table of rate parameters
    # ------------------------------------------------------------------
    def readTable(self):
        dat = self.tblPar.data()
        n = len(dat)
        tauV = np.zeros(n)
        fixV = np.zeros(n, bool)
        loV = np.full(n, np.nan)
        hiV = np.full(n, np.nan)
        for q, row in enumerate(dat):
            v = row[1]
            if isinstance(v, str):
                v = _str2double(v)
            if v is None or not np.isfinite(v) or v <= 0:
                raise ValueError('Parameter %s: the time constant must be a '
                                 'positive number.' % row[0])
            tauV[q] = float(v)
            fixV[q] = bool(row[2])
            loV[q] = self.readLimit(row[3], row[0])
            hiV[q] = self.readLimit(row[4], row[0])
        return tauV, fixV, loV, hiV

    @staticmethod
    def readLimit(c, pname):
        if c is None:
            return np.nan
        if isinstance(c, (int, float, np.floating)):
            if np.isnan(c):
                return np.nan
            v = float(c)
        else:
            c = str(c).strip()
            if not c:
                return np.nan
            v = _str2double(c)
        if not np.isfinite(v) or v <= 0:
            raise ValueError('Parameter %s: a limit must be a positive '
                             'number (or empty for none).' % pname)
        return v

    def fromGA(self):
        app = self.app
        if app.gaResultTau is None or len(app.gaResultTau) == 0:
            self.alert('Run a Global Analysis fit first.', 'No GA result')
            return
        dat = self.tblPar.data()
        tg = np.sort(np.asarray(app.gaResultTau, float).ravel())
        nCopy = min(len(dat), tg.size)
        for q in range(nCopy):
            dat[q][1] = sigRound(tg[q], 4)
        self.tblPar.setData(dat)
        if app.gaResultT0 is not None and np.isfinite(app.gaResultT0):
            self.edT0.setValue(app.gaResultT0)
        if app.gaResultFWHM is not None and app.gaResultFWHM > 0:
            self.edFW.setValue(max(app.gaResultFWHM, 1e-6))
        self.lblStatus.setText('%d time constant(s), t0 and FWHM copied '
                               'from the GA result.' % nCopy)

    def useFitted(self):
        res = self.res
        if res is None:
            return
        dat = self.tblPar.data()
        par = list(res['sch'].params)
        for row in dat:
            if row[0] in par:
                row[1] = sigRound(res['tau'][par.index(row[0])], 5)
        self.tblPar.setData(dat)
        self.edT0.setValue(res['t0'])
        if res['fwhm'] > 0:
            self.edFW.setValue(max(res['fwhm'], 1e-6))
        ca = res.get('ca')
        if isinstance(ca, dict) and ca['mode'] != 'irf' and ca['fwhm'] > 0:
            self.edCAW.setValue(ca['fwhm'])
        self.lblStatus.setText('Fitted values copied to the starting values.')

    def onIRFChange(self):
        none = self.ddIRF.currentData() == 'none'
        self.edFW.setEnabled(not none)
        self.cbFWfix.setEnabled(not none)
        if none:
            self.cbT0fix.setChecked(True)
        self.cbT0fix.setEnabled(not none)
        self.syncCAEnable()

    def syncCAEnable(self):
        hasIRF = self.ddIRF.currentData() != 'none'
        useCA = hasIRF and int(self.ddCA.currentData()) > 0
        self.ddCA.setEnabled(hasIRF)
        self.ddCAW.setEnabled(useCA)
        self.lblCAW.setEnabled(useCA)
        self.edCAW.setEnabled(useCA and self.ddCAW.currentData() != 'irf')

    def setFitTFull(self):
        self.edFitTmin.setValue(float(np.min(self.app.delay)))
        self.edFitTmax.setValue(float(np.max(self.app.delay)))

    def setResultButtons(self, on):
        for b in (self.btnUse, self.btnExpSAS, self.btnExpDAS,
                  self.btnExpPop, self.btnExpFit, self.btnExpPar,
                  self.btnSVD):
            b.setEnabled(bool(on))

    # ------------------------------------------------------------------
    # the fit
    # ------------------------------------------------------------------
    def doRun(self):
        app = self.app
        self.applyScheme()
        tauV, fixV, loV, hiV = self.readTable()
        irfM = self.ddIRF.currentData()
        rLo, rHi = self.fitRangeToStore(self.edFitTmin.value(),
                                        self.edFitTmax.value())
        tR = None if rLo is None else [rLo, rHi]
        caS = {'order': int(self.ddCA.currentData()),
               'mode': self.ddCAW.currentData(),
               'fwhm': float(self.edCAW.value())}
        sch = self.sch
        t0In, fwIn = self.edT0.value(), self.edFW.value()
        t0Fix, fwFix = self.cbT0fix.isChecked(), self.cbFWfix.isChecked()
        self.lblStatus.setText('Fitting...')
        D, wl, t = app.deltaA, app.wavelength, app.delay
        try:
            o = runBusy(app, self, 'Fitting...', lambda: fit_target_analysis(
                D, wl, t, sch, tauV, fixV, t0In, fwIn, t0Fix, fwFix, irfM,
                tR, loV, hiV, caS))
        except Exception as err:
            self.alert('Target fit failed:\n%s' % err, 'Fit error')
            self.lblStatus.setText('Fit failed.')
            return
        r = {
            'sch': sch,
            'tau': np.asarray(o.tau, float),
            'tauFixed': np.asarray(fixV, bool),
            'tauStd': np.asarray(o.tau_std, float),
            'atLimit': np.asarray(o.at_limit, bool),
            't0': float(o.t0), 'fwhm': float(o.fwhm),
            't0Std': float(o.t0_std), 'fwhmStd': float(o.fwhm_std),
            'fwhmLo': float(o.fwhm_lo),
            'fwhmAtLimit': bool(o.fwhm_at_limit),
            'irfRetry': int(o.irf_retry),
            'fwhmHi': float(o.fwhm_hi),
            'fwhmTooWide': bool(o.fwhm_too_wide),
            't0Fixed': bool(t0Fix or irfM == 'none'),
            'fwhmFixed': bool(fwFix or irfM == 'none'),
            'irf': irfM, 'tRange': tR,
            'SAS': np.asarray(o.sas, float),
            'DAS': np.asarray(o.das, float),
            'conc': np.asarray(o.conc, float),
            'fit': np.asarray(o.fit, float),
            'rms': float(o.rms), 'initialRMS': float(o.initial_rms),
            'iters': int(o.iters), 'converged': bool(o.converged),
            'obsTau': np.asarray(o.obs_tau, float),
            'nObs': int(o.n_obs), 'nFree': int(o.n_free),
            'reordered': bool(o.reordered), 'ca': None}
        if o.ca_order > 0:
            r['ca'] = {'order': int(o.ca_order),
                       'spec': np.asarray(o.ca_spec, float),
                       'fwhm': float(o.ca_fwhm),
                       'fwhmStd': float(o.ca_fwhm_std), 'mode': o.ca_mode}
        self.res = r
        self.closeSVD()
        self.txtResult.setPlainText('\n'.join(self.resultLines()))
        if r['converged']:
            self.lblStatus.setText('Done. RMS = %.4g' % r['rms'])
        else:
            self.lblStatus.setText(
                'Stopped at the iteration limit. RMS = %.4g' % r['rms'])
        self.setResultButtons(True)
        self.plotAll()
        self.saveSettings()
        app.tgSaved = {'res': self._packRes(r),
                       'curWL': float(self.edWL.value()),
                       'status': self.lblStatus.text(),
                       'lines': self.txtResult.toPlainText().split('\n')}

    @staticmethod
    def _packRes(r):
        """The result in a form a session file can hold (the scheme as
        its text)."""
        p = {k: v for k, v in r.items() if k != 'sch'}
        p['schemeTxt'] = list(target_scheme_text(r['sch']))
        return p

    def _unpackRes(self, p):
        """Result stored by _packRes; None if it does not belong to the
        data now loaded."""
        app = self.app
        if 'DAS' not in p or 'schemeTxt' not in p:
            return None
        if np.shape(p['SAS'])[0] != app.wavelength.size or \
                np.shape(p['conc'])[0] != app.delay.size:
            return None
        r = {k: v for k, v in p.items() if k != 'schemeTxt'}
        r['sch'] = parse_target_scheme(list(p['schemeTxt']))
        for k in ('tau', 'tauStd', 'obsTau'):
            r[k] = np.atleast_1d(np.asarray(r[k], float))
        for k in ('tauFixed', 'atLimit'):
            r[k] = np.atleast_1d(np.asarray(r[k], bool))
        for k in ('SAS', 'DAS', 'conc', 'fit'):
            r[k] = np.asarray(r[k], float)
        if r.get('tRange') is not None:
            r['tRange'] = [float(v) for v in np.ravel(r['tRange'])]
        if isinstance(r.get('ca'), dict):
            r['ca'] = dict(r['ca'])
            r['ca']['spec'] = np.asarray(r['ca']['spec'], float)
        else:
            r['ca'] = None
        return r

    def resultLines(self):
        app = self.app
        res = self.res
        sch = res['sch']
        u = app.tUnitTxt()
        st1 = sch.steps
        nSt = len(st1)
        tau = res['tau']
        tauStep = [tau[s['par']] / s['mult'] for s in st1]
        L = ['Scheme:']
        wLab = max([len(s['label']) for s in st1] or [1])
        for s in st1:
            L.append('  %s  : %s' % (s['label'].ljust(wLab), s['expr']))
        if nSt == 0:
            L.append('  (no steps)')
        c0s = ['%s = %g' % (sch.species[q], sch.c0[q])
               for q in range(sch.n_s) if sch.c0[q] != 0]
        L.append('Initial population: ' + ', '.join(c0s))
        wl = np.asarray(app.wavelength, float)
        for z in sch.zero:
            nIn = int(np.count_nonzero((wl >= z['lo']) & (wl <= z['hi'])))
            L.append('Spectrum of %s = 0 for %g .. %g %s (%d points)' % (
                sch.species[z['sp']], z['lo'], z['hi'], app.xUnitTxt(), nIn))
        for rl in sch.rel:
            if rl['ratio'] == 1:
                rhs = sch.species[rl['master']]
            else:
                rhs = '%g x %s' % (rl['ratio'], sch.species[rl['master']])
            if np.isfinite(rl['lo']) and np.isfinite(rl['hi']):
                L.append('Spectrum of %s = %s for %g .. %g %s' % (
                    sch.species[rl['dep']], rhs, rl['lo'], rl['hi'],
                    app.xUnitTxt()))
            else:
                L.append('Spectrum of %s = %s' % (sch.species[rl['dep']],
                                                  rhs))
        nonneg = np.asarray(sch.nonneg, bool)
        if nonneg.any():
            L.append('Non-negative spectra: ' + ', '.join(
                n for n, f in zip(sch.species, nonneg) if f))
        L.append('')
        L.append('RMS = %.4g  (start: %.4g), %d iterations' % (
            res['rms'], res['initialRMS'], res['iters']))
        if not res['converged']:
            L.append('** NOT CONVERGED: iteration limit reached.')
            L.append('   Try other starting values or fix more.')
        if res['tRange'] is not None:
            L.append('Fit t range: %.4g .. %.4g %s' % (
                res['tRange'][0], res['tRange'][1], u))
        L.append('')
        L.append('Rate parameters (tau = 1/k):')
        wPar = max([len(p) for p in sch.params] or [1])
        anyUndet = False
        for q in range(sch.n_p):
            pn = sch.params[q].ljust(wPar)
            sd = res['tauStd'][q]
            kTxt = 'k = %.4g /%s' % (1.0 / tau[q], u)
            if res['tauFixed'][q]:
                L.append('  %s  tau = %.4g %s (fixed)   %s' % (
                    pn, tau[q], u, kTxt))
            elif res['atLimit'][q]:
                L.append('  %s  tau = %.4g %s  [at its limit]' % (
                    pn, tau[q], u))
            elif np.isinf(sd):
                anyUndet = True
                L.append('  %s  tau = %.4g %s  [not determined]' % (
                    pn, tau[q], u))
            elif np.isfinite(sd):
                L.append('  %s  tau = %.4g +/- %.2g %s   %s' % (
                    pn, tau[q], sd, u, kTxt))
            else:
                L.append('  %s  tau = %.4g %s   %s' % (pn, tau[q], u, kTxt))
        tied = [q for q in range(nSt) if st1[q]['mult'] != 1]
        if tied:
            L.append('Steps with a tied rate:')
            for q in tied:
                L.append('  %s  tau = %.4g %s  (%s)' % (
                    st1[q]['label'], tauStep[q], u, st1[q]['expr']))
        L.append('')
        if res['irf'] == 'none':
            L.append('No IRF; kinetics start at t0 = %.4g %s' % (
                res['t0'], u))
        else:
            L.append(irfLine('IRF t0  ', res['t0'], res['t0Std'],
                             res['t0Fixed'], u))
            if res.get('fwhmAtLimit') is True:
                L.append('IRF FWHM = %.4g %s  [not resolved]' % (
                    res['fwhm'], u))
            elif res.get('fwhmTooWide') is True:
                L.append('IRF FWHM = %.4g %s  [not determined]' % (
                    res['fwhm'], u))
            else:
                L.append(irfLine('IRF FWHM', res['fwhm'], res['fwhmStd'],
                                 res['fwhmFixed'], u))
            L += irfLimitNote(app, Result(
                fwhm_at_limit=bool(res.get('fwhmAtLimit', False)),
                fwhm_lo=res.get('fwhmLo', np.nan),
                fwhm_too_wide=bool(res.get('fwhmTooWide', False)),
                fwhm_hi=res.get('fwhmHi', np.nan),
                irf_retry=int(res.get('irfRetry', 0))))
        ca = res.get('ca')
        if isinstance(ca, dict) and ca['order'] > 0:
            L.append('')
            L.append('Coherent artifact: %d term(s)' % ca['order'])
            if ca['mode'] == 'irf':
                L.append('  width = IRF FWHM')
            else:
                L.append('  ' + irfLine('width', ca['fwhm'], ca['fwhmStd'],
                                        ca['mode'] == 'fixed', u))
            for q in range(int(ca['order'])):
                col = np.abs(ca['spec'][:, q])
                col = np.where(np.isnan(col), -np.inf, col)
                iMx = int(np.argmax(col))
                L.append('  %-14s largest %.3g at %.4g %s' % (
                    CA_NAMES[q], ca['spec'][iMx, q], wl[iMx],
                    app.xUnitTxt()))
            if ca['mode'] != 'irf' and not res['fwhmFixed'] and \
                    ca['order'] >= 2:
                L += ['  With a width of their own the derivative',
                      '  terms absorb a change of the IRF FWHM: it',
                      '  is then hardly determined. Fix it, or use',
                      '  "Same as IRF".']
        ot = res['obsTau'][np.isfinite(res['obsTau'])]
        if ot.size:
            L.append('')
            L.append('Observable lifetimes (%s):' % u)
            L.append('  ' + ', '.join('%.4g' % x for x in ot))
        first = True
        for q in range(sch.n_s):
            outI = [i for i in range(nSt) if st1[i]['from'] == q]
            if len(outI) < 2:
                continue
            kOut = np.array([1.0 / tauStep[i] for i in outI])
            fr = 100 * kOut / np.sum(kOut)
            parts = []
            for w, i in enumerate(outI):
                dst = 'ground state' if st1[i]['to'] < 0 \
                    else sch.species[st1[i]['to']]
                parts.append('-> %s %.1f %%' % (dst, fr[w]))
            if first:
                L.append('')
                L.append('Branching:')
                first = False
            pIdx = [st1[i]['par'] for i in outI]
            note = ''
            if np.any(np.isinf(res['tauStd'][pIdx])):
                note = '   [not determined]'
            elif np.all(res['tauFixed'][pIdx]) or len(set(pIdx)) == 1:
                note = '   (fixed by the model)'
            L.append('  %s: %s%s' % (sch.species[q], ', '.join(parts), note))
        L.append('')
        L.append('Notes:')
        if anyUndet:
            L.append('- %d free rate parameter(s), %d observable '
                     'lifetime(s):' % (res['nFree'], res['nObs']))
            L += ['  the marked ones only act in combination and',
                  '  cannot be separated by these data. Fix one,',
                  '  tie it to another (k2 = 2*k1), or add a',
                  '  spectral constraint (zero: / equal:).']
        if np.any(res['atLimit']):
            L += ['- A parameter ended on one of its limits: the',
                  '  data would move it further. No error is',
                  '  given for it.']
        if res['reordered']:
            L += ['- The same fit is obtained with time constants',
                  '  interchanged; the order of your starting',
                  '  values was kept.']
        L += ['- The fit tests the decay rates, not the scheme:',
              '  other schemes with the same rates fit equally',
              '  well. The spectra are those implied by the',
              '  scheme chosen here.']
        return L

    # ------------------------------------------------------------------
    # plots
    # ------------------------------------------------------------------
    def plotAll(self):
        self.plotSpec()
        self.plotPop()
        self.plotKin()
        self.plotRes()
        self.drawDiagrams()

    def plotTimePanels(self):
        self.plotPop()
        self.plotKin()
        self.plotRes()

    def applyTimeAxis(self, ax):
        app = self.app
        tA = np.asarray(app.delay, float).ravel()
        isLog = self.ddScale.currentText() == 'Log' and bool(np.any(tA > 0))
        if isLog:
            ax.set_xscale('log')
            xl = [float(np.min(tA[tA > 0])), float(np.max(tA))]
        else:
            ax.set_xscale('linear')
            xl = [float(np.min(tA)), float(np.max(tA))]
        res = self.res
        if res is not None and res['tRange'] is not None:
            fr = list(res['tRange'])
            if isLog:
                fr[0] = max(fr[0], xl[0])
            if fr[0] < fr[1]:
                xl = fr
        if xl[1] > xl[0]:
            ax.set_xlim(xl[0], xl[1])

    def _marker(self, ax, bl, isY):
        wl = self.edWL.value()
        kw = dict(color='k', linewidth=1.0, alpha=0.7)
        h = ax.axhline(wl, **kw) if isY else ax.axvline(wl, **kw)
        bl.setArtists([h])
        return h

    def plotSpec(self):
        app = self.app
        res = self.res
        ax = self.axSAS
        ax.clear()
        styleAxes(ax, grid=True)
        ax.set_xlabel(app.xAxLabel())
        isDAS = self.ddSpec.currentText() == 'DAS'
        norm = self.cbNorm.isChecked()
        if res is None:
            ax.set_ylabel(('Amplitude (%s)' % self.app.sigName()))
            setTitle(ax, 'Spectra (run fit)')
            ax.set_xlim(self.wlLo, self.wlHi)
            self.hWLspec = self._marker(ax, self.blSAS, False)
            self.cvSAS.draw_idle()
            return
        wl = np.asarray(app.wavelength, float)
        if isDAS:
            Y = res['DAS']
            nm = []
            for q in range(Y.shape[1]):
                if np.isfinite(res['obsTau'][q]):
                    nm.append('τ = %.3g %s' % (res['obsTau'][q],
                                               app.tUnitAx()))
                else:
                    nm.append('τ = ∞')
        else:
            Y = res['SAS']
            nm = list(res['sch'].species)

        def scaled(y):
            if norm:
                with np.errstate(invalid='ignore'):
                    m = np.nanmax(np.abs(y)) if np.any(np.isfinite(y)) \
                        else np.nan
                if np.isfinite(m) and m > 0:
                    return y / m
            return y

        for q in range(Y.shape[1]):
            ax.plot(wl, scaled(Y[:, q]), '-', linewidth=1.4,
                    color=lineColor(q), label=nm[q].replace('$', '\\$'))
        ca = res.get('ca')
        if isinstance(ca, dict) and ca['spec'].shape[0] == wl.size:
            caNm = ['artifact: Gaussian', 'artifact: 1st deriv.',
                    'artifact: 2nd deriv.']
            caCol = [(0.30, 0.30, 0.30), (0.50, 0.50, 0.50),
                     (0.68, 0.68, 0.68)]
            for q in range(min(3, ca['spec'].shape[1])):
                ax.plot(wl, scaled(ca['spec'][:, q]), '--', linewidth=1.1,
                        color=caCol[q], label=caNm[q])
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        self.hWLspec = self._marker(ax, self.blSAS, False)
        ax.set_ylabel('Normalized amplitude' if norm else ('Amplitude (%s)' % self.app.sigName()))
        if isDAS:
            setTitle(ax, 'Decay-associated spectra (DAS) of this fit',
                     'DAS of this fit', 'DAS')
        else:
            setTitle(ax, 'Species-associated spectra (SAS)', 'SAS')
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        if not norm:
            sciTicks(ax, 'y')
        legend(ax, frameon=False)
        self.cvSAS.draw_idle()

    def plotPop(self):
        app = self.app
        res = self.res
        ax = self.axPop
        ax.clear()
        styleAxes(ax, grid=True)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel('Population')
        if res is None:
            setTitle(ax, 'Populations (run fit)')
            self.cvPop.draw_idle()
            return
        sch = res['sch']
        for q in range(sch.n_s):
            ax.plot(app.delay, res['conc'][:, q], '-', linewidth=1.4,
                    color=lineColor(q),
                    label=sch.species[q].replace('$', '\\$'))
        setTitle(ax, 'Populations c(t)')
        self.applyTimeAxis(ax)
        legend(ax, frameon=False)
        self.cvPop.draw_idle()

    def plotKin(self):
        app = self.app
        res = self.res
        ax = self.axKin
        ax.clear()
        styleAxes(ax, grid=True)
        wl = np.asarray(app.wavelength, float)
        iw = int(np.argmin(np.abs(wl - self.edWL.value())))
        t = np.asarray(app.delay, float)
        ax.plot(t, app.deltaA[iw, :], 'k.', markersize=5, label='Data')
        if res is not None:
            if self.cbContrib.isChecked():
                sch = res['sch']
                for q in range(sch.n_s):
                    a = res['SAS'][iw, q]
                    if not np.isfinite(a) or a == 0:
                        continue
                    ax.plot(t, a * res['conc'][:, q], '-', linewidth=0.9,
                            color=lineColor(q),
                            label=sch.species[q].replace('$', '\\$'))
                ca = res.get('ca')
                if isinstance(ca, dict) and ca['spec'].shape[0] == wl.size \
                        and np.all(np.isfinite(ca['spec'][iw, :])):
                    caB = ca_basis(t, res['t0'], ca['fwhm'], ca['order'])
                    ax.plot(t, caB @ ca['spec'][iw, :], '--', linewidth=0.9,
                            color=(0.4, 0.4, 0.4), label='artifact')
            ax.plot(t, res['fit'][iw, :], 'r-', linewidth=1.5, label='Fit')
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel(self.app.sigName())
        setTitle(ax, 'Kinetics at %s = %.4g %s' % (
            app.xSym(), wl[iw], app.xUnitAx()))
        self.applyTimeAxis(ax)
        sciTicks(ax, 'y')
        legend(ax, frameon=False)
        self.cvKin.draw_idle()

    def residualInView(self):
        app = self.app
        res = self.res
        tA = np.asarray(app.delay, float).ravel()
        isLog = self.ddScale.currentText() == 'Log' and bool(np.any(tA > 0))
        keep = np.ones(tA.size, bool)
        if isLog:
            keep = tA > 0
        if res['tRange'] is not None:
            keep = keep & (tA >= res['tRange'][0]) & (tA <= res['tRange'][1])
        if np.count_nonzero(keep) < 2:
            keep = np.ones(tA.size, bool)
            isLog = False
        return app.deltaA[:, keep] - res['fit'][:, keep], tA[keep], isLog

    def plotRes(self):
        app = self.app
        res = self.res
        ax, cax = self.axRes, self.caxRes
        ax.clear()
        cax.clear()
        styleAxes(ax)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel(app.xAxLabel())
        if res is None:
            cax.set_visible(False)
            setTitle(ax, 'Residual map (run fit)')
            self.hWLres = None
            self.blRes.setArtists([])
            self._layRes()
            self.cvRes.draw_idle()
            return
        cax.set_visible(True)
        Rk, tK, isLog = self.residualInView()
        wl = np.asarray(app.wavelength, float)
        ax.set_facecolor(AX_BG)
        cl = absMax(Rk)
        if not np.isfinite(cl) or cl == 0:
            cl = 1.0
        h = None
        if Rk.shape[0] > 1 and Rk.shape[1] > 1:
            h = ax.pcolormesh(tK, wl, np.ma.masked_invalid(Rk[:-1, :-1]),
                              shading='flat', cmap=app.getColormap(),
                              vmin=-cl, vmax=cl, rasterized=True)
        ax.set_xscale('log' if isLog else 'linear')
        ax.set_xlim(float(np.min(tK)), float(np.max(tK)))
        ax.set_ylim(float(np.min(wl)), float(np.max(wl)))
        self.hWLres = self._marker(ax, self.blRes, True)
        setTitle(ax, 'Residual (data - fit),  RMS = %.3g' % res['rms'],
                 'Residual,  RMS = %.3g' % res['rms'])
        cax.set_visible(h is not None)
        if h is not None:
            cb = self.cvRes.figure.colorbar(h, cax=cax)
            cb.set_label(self.app.sigName(), fontsize=FONT_PT)
            cb.ax.tick_params(labelsize=FONT_PT, length=3)
            sciTicks(cb.ax, 'y')
        self._layRes()
        self.cvRes.draw_idle()

    # ------------------------------------------------------------------
    # probe position of the kinetics panel
    # ------------------------------------------------------------------
    def onWLChange(self):
        wl = self.edWL.value()
        if self.hWLspec is not None:
            self.hWLspec.set_xdata([wl, wl])
            self.blSAS.update()
        if self.hWLres is not None:
            self.hWLres.set_ydata([wl, wl])
            self.blRes.update()
        self.plotKin()
        if isinstance(self.app.tgSaved, dict) and \
                'curWL' in self.app.tgSaved:
            self.app.tgSaved['curWL'] = float(wl)

    def onPickDown(self, which, ev):
        ax = self.axSAS if which == 'sas' else self.axRes
        if ev.button != 1 or ev.inaxes is not ax:
            return
        if which == 'sas' and self.cvSAS.legendAt(ev) is not None:
            return                               # (moving the legend)
        if which == 'res' and self.res is None:
            return
        self.dragging = which
        self.applyPick(ev.xdata if which == 'sas' else ev.ydata)

    def onPickMotion(self, which, ev):
        if self.dragging != which:
            return
        ax = self.axSAS if which == 'sas' else self.axRes
        if ev.x is None or ev.y is None:
            return
        try:
            x, y = ax.transData.inverted().transform((ev.x, ev.y))
        except Exception:
            return
        self.applyPick(x if which == 'sas' else y)

    def onPickUp(self):
        self.dragging = None

    def applyPick(self, wl):
        if wl is None or not np.isfinite(wl):
            return
        lo, hi = self.edWL.limits()
        self.edWL.setValue(min(max(float(wl), lo), hi))
        self.onWLChange()

    # ------------------------------------------------------------------
    # syntax help, diagram, residual SVD
    # ------------------------------------------------------------------
    def showSyntaxHelp(self):
        if self.synWin is not None:
            self.synWin.raise_()
            self.synWin.activateWindow()
            return
        w = SubWindow('Scheme syntax', 800, 540, self._synClosed)
        g = QtWidgets.QVBoxLayout(w)
        g.setContentsMargins(12, 12, 12, 12)
        g.setSpacing(8)
        lb = label(SYNTAX_INTRO)
        lb.setWordWrap(True)
        g.addWidget(lb)
        tb = QtWidgets.QTableWidget(len(SYNTAX_ROWS), 2)
        tb.setHorizontalHeaderLabels(['Statement', 'Meaning'])
        tb.verticalHeader().setVisible(False)
        tb.verticalHeader().setDefaultSectionSize(24)
        tb.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        tb.setColumnWidth(0, 250)
        tb.horizontalHeader().setStretchLastSection(True)
        for i, (a, b) in enumerate(SYNTAX_ROWS):
            it = QtWidgets.QTableWidgetItem(a)
            it.setFont(monoFont())
            tb.setItem(i, 0, it)
            tb.setItem(i, 1, QtWidgets.QTableWidgetItem(b))
        g.addWidget(tb, 1)
        lb = label(SYNTAX_NOTE, color=(0.3, 0.3, 0.3))
        lb.setWordWrap(True)
        g.addWidget(lb)
        self.synWin = w
        self.synTable = tb
        w.show()

    def _synClosed(self):
        self.synWin = None

    def closeSyntaxHelp(self):
        if self.synWin is not None:
            w, self.synWin = self.synWin, None
            w.close()

    def diagramValues(self):
        """Fitted time constant of every step of the scheme shown, and
        which of them the data do not determine; (None, None) unless
        the scheme shown is the one that was fitted."""
        res, sch = self.res, self.sch
        if res is None or sch is None:
            return None, None
        try:
            rs = res['sch']
            same = list(rs.species) == list(sch.species) and \
                len(rs.steps) == len(sch.steps) and \
                np.array_equal(np.ravel(rs.c0), np.ravel(sch.c0)) and \
                list(rs.zero) == list(sch.zero) and \
                list(rs.rel) == list(sch.rel) and \
                np.array_equal(np.ravel(rs.nonneg), np.ravel(sch.nonneg))
            if same and rs.steps:
                same = all(a['from'] == b['from'] and a['to'] == b['to']
                           and a['coef'] == b['coef']
                           and a['expr'] == b['expr']
                           for a, b in zip(rs.steps, sch.steps))
            if not same:
                return None, None
            tv, uv = [], []
            for s in rs.steps:
                ip = s['par']
                tv.append(res['tau'][ip] / s['mult'])
                uv.append(bool(not res['tauFixed'][ip]
                               and np.isinf(res['tauStd'][ip])))
            return tv, uv
        except Exception:
            return None, None

    def drawDiagram(self, cv, ax, fsMax, hint):
        """The scheme as boxes and arrows on ax, which fills cv."""
        ax.clear()
        ax.set_position([0, 0, 1, 1])
        ax.set_axis_off()
        cv.figure.set_facecolor('white')
        sch = self.sch
        if sch is None:
            return None
        Wd, Hd = float(cv.width()), float(cv.height())
        if Wd < 60 or Hd < 40:
            Wd, Hd = 430.0, 152.0
        tauS, und = self.diagramValues()
        Dg = targetDiagramLayout(sch, Wd, Hd, fsMax, tauS, und,
                                 self.app.tUnitTxt())
        # pixels of the layout -> points of the font: the canvas has
        # 100 pixels per inch; DejaVu Sans is wider than the letters
        # the layout reckons with, hence a little smaller
        pt = 72.0 / 100.0 * 0.92 * min(Wd / Dg['W'], Hd / Dg['H'], 1.0)
        ink = (0.25, 0.25, 0.25)
        for ln in Dg['line']:
            ax.plot(ln['x'], ln['y'], '-', color=ink, linewidth=1.0,
                    solid_capstyle='butt')
        for hd in Dg['head']:
            ax.fill(hd['x'], hd['y'], facecolor=ink, edgecolor='none')
        for rc in Dg['rect']:
            cq = np.array(lineColor(rc['sp']), float)
            x, y, w, h = rc['pos']
            rr = 0.2 * h
            ax.add_patch(FancyBboxPatch(
                (x + rr, y + rr), w - 2 * rr, h - 2 * rr,
                boxstyle='round,pad=%g,rounding_size=%g' % (rr, rr),
                facecolor=1 - 0.18 * (1 - cq), edgecolor=cq,
                linewidth=1 + 1.2 * float(rc['init']),
                mutation_aspect=1.0))
        for tx in Dg['text']:
            tc, tw = (0.40, 0.40, 0.40), 'normal'
            if tx['kind'] == 'name':
                tc, tw = (0.07, 0.07, 0.07), 'bold'
            elif tx['kind'] == 'value':
                tc = RUN_COLOR
            elif tx['kind'] == 'rate':
                tc = (0.20, 0.20, 0.20)
            kw = {}
            if tx['bg']:
                kw['bbox'] = dict(facecolor='white', edgecolor='none',
                                  pad=1.0)
            ax.text(tx['x'], tx['y'], tx['str'].replace('$', '\\$'),
                    fontsize=max(tx['size'] * pt, 3.0), fontweight=tw,
                    color=tc, ha=tx['ha'],
                    va={'middle': 'center_baseline'}.get(tx['va'],
                                                         tx['va']),
                    clip_on=False, **kw)
        if hint:
            ax.text(Dg['W'] - 0.6 * Dg['fs'], Dg['noteY'], hint,
                    fontsize=max(Dg['noteSize'] * pt, 3.0),
                    color=(0.45, 0.45, 0.45), ha='right', va='bottom')
        ax.set_xlim(0, Dg['W'])
        ax.set_ylim(0, Dg['H'])
        return Dg

    def drawDiagrams(self):
        try:
            self.diaLayout = self.drawDiagram(self.cvDia, self.axDia, 11.0,
                                              'click to enlarge')
            self.cvDia.draw_idle()
            if self.diaWin is not None and self.cvDiaWin is not None:
                self.diaWinLayout = self.drawDiagram(
                    self.cvDiaWin, self.axDiaWin, 20.0, '')
                self.cvDiaWin.draw_idle()
        except Exception as err:
            self.lblStatus.setText('Diagram not drawn: %s' % err)

    def openDiagramWindow(self):
        if self.diaWin is not None:
            self.diaWin.raise_()
            self.diaWin.activateWindow()
        else:
            w = SubWindow('Kinetic scheme', 900, 520, self._diaClosed)
            w.setStyleSheet('background: white;')
            g = QtWidgets.QVBoxLayout(w)
            g.setContentsMargins(10, 10, 10, 10)
            self.cvDiaWin = PlotCanvas()
            self.axDiaWin = self.cvDiaWin.figure.add_axes([0, 0, 1, 1])
            self.cvDiaWin.layoutFcn = self.drawDiagrams
            g.addWidget(self.cvDiaWin)
            self.diaWin = w
            w.show()
        self.drawDiagrams()

    def _diaClosed(self):
        self.diaWin = None
        self.cvDiaWin = None
        self.axDiaWin = None

    def closeDiagramWindow(self):
        if self.diaWin is not None:
            w = self.diaWin
            self._diaClosed()
            w.close()

    def closeSVD(self):
        if self.svdWin is not None:
            w, self.svdWin = self.svdWin, None
            w.close()

    def _svdClosed(self):
        self.svdWin = None

    def showResidualSVD(self):
        app = self.app
        res = self.res
        if res is None:
            return
        self.closeSVD()
        tA = np.asarray(app.delay, float).ravel()
        keep = np.ones(tA.size, bool)
        if res['tRange'] is not None:
            keep = (tA >= res['tRange'][0]) & (tA <= res['tRange'][1])
        Rm = np.array(app.deltaA[:, keep] - res['fit'][:, keep], float)
        rowOK = np.any(np.isfinite(Rm), axis=1) & \
            np.all(np.isfinite(res['SAS']), axis=1)
        Rm[~np.isfinite(Rm)] = 0
        Rm[~rowOK, :] = 0
        Us, sv, Vt = np.linalg.svd(Rm, full_matrices=False)
        nShow = min(2, sv.size)
        if nShow < 1:
            return
        tK = tA[keep]
        w = SubWindow('Residual SVD', 1150, 400, self._svdClosed)
        g = QtWidgets.QVBoxLayout(w)
        g.setContentsMargins(10, 10, 10, 10)
        cv = PlotCanvas()
        g.addWidget(cv)
        a1 = cv.figure.add_subplot(131)
        a2 = cv.figure.add_subplot(132)
        a3 = cv.figure.add_subplot(133)

        def lay():
            W, H = cv.sizePt()
            L, R, T, B, gap = 62.0, 16.0, 22.0, 34.0, 60.0
            wu = (W - L - R - 2 * gap) / 2.8
            cv.place(a1, L, B, 0.8 * wu, H - T - B)
            cv.place(a2, L + 0.8 * wu + gap, B, wu, H - T - B)
            cv.place(a3, L + 1.8 * wu + 2 * gap, B, wu, H - T - B)

        cv.layoutFcn = lay
        nSv = min(20, sv.size)
        for a in (a1, a2, a3):
            styleAxes(a, grid=True)
        a1.semilogy(np.arange(1, nSv + 1),
                    np.maximum(sv[:nSv], np.finfo(float).tiny), 'o-',
                    linewidth=1.2, color=(0.2, 0.2, 0.2), markersize=5,
                    markerfacecolor=(0.2, 0.2, 0.2))
        a1.set_xlabel('Index')
        a1.set_ylabel('Singular value')
        setTitle(a1, 'Singular values of the residual', 'Singular values')
        for q in range(nShow):
            a2.plot(app.wavelength, Us[:, q], '-', linewidth=1.3,
                    color=lineColor(q),
                    label='#%d  (s = %.3g)' % (q + 1, sv[q]))
            a3.plot(tK, Vt[q, :], '-', linewidth=1.3, color=lineColor(q),
                    label='#%d' % (q + 1))
        a2.axhline(0, color='k', linestyle=':', linewidth=0.8)
        a3.axhline(0, color='k', linestyle=':', linewidth=0.8)
        a2.set_xlabel(app.xAxLabel())
        a2.set_ylabel('Left singular vector')
        setTitle(a2, 'Spectral structure of the residual',
                 'Spectral structure')
        a2.set_xlim(float(np.min(app.wavelength)),
                    float(np.max(app.wavelength)))
        legend(a2, frameon=False)
        a3.set_xlabel('Delay time (%s)' % app.tUnitAx())
        a3.set_ylabel('Right singular vector')
        setTitle(a3, 'Temporal structure of the residual',
                 'Temporal structure')
        if self.ddScale.currentText() == 'Log' and np.any(tK > 0):
            a3.set_xscale('log')
            a3.set_xlim(float(np.min(tK[tK > 0])), float(np.max(tK)))
        elif np.max(tK) > np.min(tK):
            a3.set_xlim(float(np.min(tK)), float(np.max(tK)))
        legend(a3, frameon=False)
        self.svdWin = w
        self.svdCanvas = cv
        self.svdValues = sv
        w.show()

    # ------------------------------------------------------------------
    # reset, settings, exports
    # ------------------------------------------------------------------
    def doReset(self):
        self.res = None
        self.app.tgSaved = None
        self.closeSVD()
        self.setResultButtons(False)
        self.txtResult.setPlainText('\n'.join(IDLE_TEXT))
        self.plotAll()
        self.lblStatus.setText('Reset.')

    def saveSettings(self):
        sLo, sHi = self.fitRangeToStore(self.edFitTmin.value(),
                                        self.edFitTmax.value())
        self.app.tgSet = {
            'preset': self.ddPreset.currentData(),
            'schemeTxt': self.schemeLines(),
            'table': self.tblPar.data(),
            't0': float(self.edT0.value()),
            'fwhm': float(self.edFW.value()),
            't0Fixed': self.cbT0fix.isChecked(),
            'fwhmFixed': self.cbFWfix.isChecked(),
            'irf': self.ddIRF.currentData(),
            'ca': {'order': int(self.ddCA.currentData()),
                   'mode': self.ddCAW.currentData(),
                   'fwhm': float(self.edCAW.value())},
            'fitTMin': sLo, 'fitTMax': sHi,
            'scale': self.ddScale.currentText(),
            'norm': self.cbNorm.isChecked(),
            'contrib': self.cbContrib.isChecked(),
            'spec': self.ddSpec.currentText()}

    def speciesHeaders(self):
        hdr = []
        for nm in self.res['sch'].species:
            nm = re.sub(r'[^A-Za-z0-9_]', '_', nm)
            if not re.match(r'^[A-Za-z]', nm):
                nm = 'S_' + nm
            hdr.append(nm)
        for q in range(1, len(hdr)):
            while hdr[q] in hdr[:q]:
                hdr[q] += '_'
        return hdr

    def exportSAS(self):
        app, res = self.app, self.res
        if res is None:
            return
        cols = [app.wavelength] + [res['SAS'][:, q]
                                   for q in range(res['SAS'].shape[1])]
        hdr = [app.xHdr()] + self.speciesHeaders()
        ca = res.get('ca')
        if isinstance(ca, dict) and \
                ca['spec'].shape[0] == app.wavelength.size:
            caHdr = ['artifact_Gaussian', 'artifact_1st_derivative',
                     'artifact_2nd_derivative']
            for q in range(min(3, ca['spec'].shape[1])):
                cols.append(ca['spec'][:, q])
                hdr.append(caHdr[q])
        self.doExportTable(cols, hdr, 'target_SAS',
                           'Export species-associated spectra')

    def exportDAS(self):
        app, res = self.app, self.res
        if res is None:
            return
        hdr = []
        for q, v in enumerate(res['obsTau']):
            if np.isfinite(v):
                hdr.append('tau%d=%.4g_%s' % (q + 1, v, app.tUnitHdr()))
            else:
                hdr.append('tau%d=inf' % (q + 1))
        cols = [app.wavelength] + [res['DAS'][:, q]
                                   for q in range(res['DAS'].shape[1])]
        self.doExportTable(cols, [app.xHdr()] + hdr, 'target_DAS',
                           'Export decay-associated spectra')

    def exportPop(self):
        app, res = self.app, self.res
        if res is None:
            return
        cols = [app.delay] + [res['conc'][:, q]
                              for q in range(res['conc'].shape[1])]
        self.doExportTable(cols, ['t_%s' % app.tUnitHdr()]
                           + self.speciesHeaders(), 'target_populations',
                           'Export populations')

    def exportFit(self):
        app, res = self.app, self.res
        if res is None:
            return
        hdr = [app.xHdr()] + ['t=%.4g_%s' % (x, app.tUnitHdr())
                              for x in np.ravel(app.delay)]
        cols = [app.wavelength] + [res['fit'][:, j]
                                   for j in range(res['fit'].shape[1])]
        self.doExportTable(cols, hdr, 'target_fit',
                           'Export fitted data matrix')

    def exportPar(self):
        if self.res is None:
            return
        full = self.app.putFile('Text (*.txt)', 'Export fit parameters',
                                'target_parameters.txt', parent=self)
        if not full:
            return
        try:
            with open(full, 'w', encoding='utf-8', newline='\n') as fh:
                for ln in self.txtResult.toPlainText().split('\n'):
                    fh.write(ln + '\n')
        except OSError:
            self.alert('Cannot write:\n%s' % full, 'Error')
            return
        self.lblStatus.setText('Parameters saved: %s'
                               % os.path.basename(full))

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

    def onClosed(self):
        try:
            self.saveSettings()
        except Exception:
            pass
        self.closeSVD()
        self.closeDiagramWindow()
        self.closeSyntaxHelp()
