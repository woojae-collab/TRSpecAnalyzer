"""Main window of the program. Port of createMainGUI and of the methods
that draw and handle it in TRSpecAnalyzer.m; method names follow the
MATLAB program."""
import json
import os

import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from PySide6 import QtCore, QtWidgets

from . import colormaps
from .state import (AppState, PY_VERSION, VERSION, VERSION_DATE,
                    parse_data_file, ui_str, write_data_file)
from .widgets import (AX_BG, FONT_PT, BlitManager, ElideLabel, NumEdit,
                      PlotCanvas, alert, button, colorButton, flowbox, hbox,
                      inform, label, lineColor, panel, sciTicks, styleAxes,
                      tip)
from . import zoom
from .textfit import legend, setLegendText, setSpan
from .win_base import ANALYSIS_KEYS, fitSize, openChild

Z_LABEL = 'ΔA or Intensity'


class TRSpecAnalyzer(AppState):
    """The program: state plus main window. With headless=True no
    message box is shown (messages are kept in lastMessage)."""

    Version = VERSION
    VersionDate = VERSION_DATE

    def __init__(self, headless=False):
        super().__init__()
        self.headless = headless
        self.lastMessage = None
        self.testPath = None          # stands in for file dialogs (tests)
        self.childWindows = {}
        self.dragging2D = False
        self.drag2DAx = None
        self.ax2D = []
        self.axKin = []
        self.axSpec = None
        self.cbarAx = None
        self.hVLine = []
        self.hHLine = []
        self.hMap = []
        self.createMainGUI()

    # ==================================================================
    # messages
    # ==================================================================
    def alert(self, message, title='TRSpecAnalyzer', parent=None):
        self.lastMessage = (title, message)
        if not self.headless:
            alert(parent or self.fig, message, title)

    def inform(self, message, title='TRSpecAnalyzer', parent=None):
        self.lastMessage = (title, message)
        if not self.headless:
            inform(parent or self.fig, message, title)

    def confirm(self, message, title='Confirm', yes='Yes', no='Cancel',
                parent=None):
        """Ask before a step that cannot be taken back. Headless: the
        answer is testConfirm (yes unless a test sets it)."""
        self.lastMessage = (title, message)
        if self.headless:
            return bool(getattr(self, 'testConfirm', True))
        from .widgets import confirm
        return confirm(parent or self.fig, message, title, yes, no)

    def notYet(self, name):
        self.inform('"%s" is not part of this build of the Python '
                    'version yet.' % name, 'Not available yet')

    # ==================================================================
    # main window
    # ==================================================================
    def createMainGUI(self):
        self.loadPrefs()
        self.fig = QtWidgets.QMainWindow()
        self.fig.setWindowTitle(
            'TR Spectroscopy Analyzer v%s - Python version' % PY_VERSION)
        self.fig.resize(*fitSize(1400, 820))
        self.fig.closeEvent = self._closeEvent
        central = QtWidgets.QWidget()
        self.fig.setCentralWidget(central)
        gl = QtWidgets.QGridLayout(central)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setHorizontalSpacing(8)
        gl.setVerticalSpacing(6)
        gl.setColumnStretch(0, 2)
        gl.setColumnStretch(1, 1)
        gl.setRowStretch(3, 1)
        gl.setRowStretch(4, 1)

        # ---- row 1: data handling -----------------------------------
        # (the three rows of controls continue on a second line when
        # the window is too narrow for them)
        tb, lay = flowbox(4)
        gl.addWidget(tb, 0, 0, 1, 2)
        lay.addWidget(button('Load...', self.loadData, width=80))
        lay.addWidget(button(
            'Reset all', self.resetCorrections, width=90,
            tooltip='Restore the raw data taken at load time. Undoes crop, '
                    'background, chirp, mask, and t=0 shift -- equivalent '
                    'to re-loading the file.'))
        lay.addWidget(button('Crop...', self.openCropWindow, width=80))
        lay.addWidget(button('Background Corr...', self.openBackgroundWindow,
                             width=124))
        lay.addWidget(button('Chirp Corr...', self.openChirpWindow, width=110))
        lay.addWidget(button(
            'Set t=0', self.applyZeroTimeAtCrosshair, width=85,
            tooltip='Make the currently selected crosshair position on the '
                    'main 2D map the new t=0. Click anywhere on the map '
                    'first to move the crosshair to your desired zero-time '
                    'location.'))
        lay.addWidget(button(
            ui_str('Mask \\lambda...'), self.openMaskWindow, width=130,
            tooltip='Mask wavelength regions (e.g., pump scattering). Can '
                    'set them to NaN (excluded from fits) or to 0.'))
        self.ddScale = QtWidgets.QComboBox()
        self.ddScale.addItems(['Linear', 'Log', 'Split'])
        self.ddScale.setCurrentText(self.scaleModeLabel())
        self.ddScale.setMinimumWidth(85)
        self.ddScale.currentTextChanged.connect(self.onScaleModeChange)
        lay.addPair(label('Scale:', right=True), self.ddScale)
        self.lblSplitThreshold = label(
            'Split thr. (%s):' % self.tUnitTxt(), right=True)
        self.edThreshold = NumEdit(self.splitThreshold, limits=(1e-6, np.inf),
                                   width=60)
        self.edThreshold.valueChanged.connect(self.onThresholdChange)
        lay.addPair(self.lblSplitThreshold, self.edThreshold)
        # (the status text takes what is left of the line and is cut
        # with "..." when that is too short: it never adds a line)
        self.lblStatus = ElideLabel('No data loaded', right=True,
                                    color=(0.3, 0.3, 0.3))
        lay.addFill(self.lblStatus, 160)

        # ---- row 2: analysis windows --------------------------------
        tbA, lay = flowbox(6)
        gl.addWidget(tbA, 1, 0, 1, 2)
        lay.addWidget(colorButton(
            'SVD...', (0.20, 0.45, 0.70), width=110,
            tooltip='Singular Value Decomposition: shows singular-value '
                    'spectrum, spectral and temporal basis vectors, and '
                    'reconstruction with N components. Useful for choosing '
                    'the number of exponentials in Global Analysis.'))
        lay.itemAt(lay.count() - 1).widget().clicked.connect(
            lambda *_: self.openSvdWindow())
        b = colorButton('Global Analysis...', (0.55, 0.30, 0.70), width=160)
        b.clicked.connect(lambda *_: self.openGlobalAnalysisWindow())
        lay.addWidget(b)
        b = colorButton(
            'Target Analysis...', (0.60, 0.25, 0.45), width=150,
            tooltip='Target analysis: fit a compartment model of your own '
                    '(any number of compartments; sequential, parallel, '
                    'branching or reversible steps; shared or tied rate '
                    'constants; spectral constraints) and get the '
                    'species-associated spectra and populations it implies.')
        b.clicked.connect(lambda *_: self.openTargetWindow())
        lay.addWidget(b)
        b = colorButton(
            'LDA...', (0.20, 0.55, 0.40), width=110,
            tooltip=ui_str(
                'Lifetime Density Analysis: solve for a continuous '
                'distribution A(\\lambda, \\tau) of amplitudes over a '
                'log-spaced lifetime grid with Tikhonov regularization. '
                'Useful when the kinetics are too complex for a discrete '
                'few-exponentials fit.'))
        b.clicked.connect(lambda *_: self.openLdaWindow())
        lay.addWidget(b)
        b = colorButton(
            'MCR-ALS...', (0.85, 0.50, 0.20), width=110,
            tooltip='Multivariate Curve Resolution by Alternating Least '
                    'Squares: decompose the data as D = C * S^T with '
                    'constraints (non-negativity, unimodality). Most '
                    'model-free soft-modeling approach.')
        b.clicked.connect(lambda *_: self.openMcrWindow())
        lay.addWidget(b)
        b = colorButton(
            'Coherence...', (0.30, 0.65, 0.75), width=130,
            tooltip='Vibrational coherence analysis: FFT of kinetic '
                    'residual at each wavelength gives a 2D oscillation '
                    'map (wavelength x wavenumber). Reveals coherent '
                    'vibrational modes of excited or ground state.')
        b.clicked.connect(lambda *_: self.openCoherenceWindow())
        lay.addWidget(b)
        b = colorButton(
            'Auto Analysis...', (0.20, 0.50, 0.45), width=140,
            tooltip='Automatic analysis: proposes the corrections '
                    '(background, mask, chirp, crop), estimates the number '
                    'of components and their time constants by global '
                    'fits, and looks for oscillations in the residual. The '
                    'result is a report with the plots it rests on; the '
                    'data are changed only when you apply the corrections.')
        b.clicked.connect(lambda *_: self.openAutoWindow())
        lay.addWidget(b)
        lay.addStretch(1)

        # ---- row 3: display range, export, session, units -----------
        tb2, lay = flowbox(5)
        gl.addWidget(tb2, 2, 0, 1, 2)
        self.lblViewTMin = label(
            'Display t_min (%s):' % self.tUnitTxt(), right=True)
        self.edViewTMin = NumEdit(0, width=80)
        self.edViewTMin.valueChanged.connect(
            lambda v: self.onMainViewChange('min', v))
        lay.addPair(self.lblViewTMin, self.edViewTMin)
        self.lblViewTMax = label('t_max (%s):' % self.tUnitTxt(), right=True)
        self.edViewTMax = NumEdit(0, width=80)
        self.edViewTMax.valueChanged.connect(
            lambda v: self.onMainViewChange('max', v))
        lay.addPair(self.lblViewTMax, self.edViewTMax)
        lay.addWidget(button('Reset display', self.resetMainView, width=100))
        b = colorButton('Export Corrected Data...', (0.30, 0.55, 0.80),
                        width=170)
        b.clicked.connect(lambda *_: self.exportCorrected())
        lay.addWidget(b)
        lay.addWidget(button(
            'Save Session...', self.saveSession, width=110,
            tooltip='Save the entire working state (data, corrections, fit '
                    'settings, GA results, view options) to a file so you '
                    'can resume later exactly where you left off.'))
        lay.addWidget(button(
            'Load Session...', self.loadSessionFile, width=110,
            tooltip='Restore a session saved with "Save Session".'))
        self.ddTimeUnit = QtWidgets.QComboBox()
        for txt, dat in (('fs', 'fs'), ('ps', 'ps'), ('ns', 'ns'),
                         ('μs', 'us')):
            self.ddTimeUnit.addItem(txt, dat)
        self.ddTimeUnit.setCurrentIndex(
            max(self.ddTimeUnit.findData(self.timeUnit), 0))
        self.ddTimeUnit.setMinimumWidth(70)
        self.ddTimeUnit.setToolTip(tip(
            'Switch time-axis labels between femtoseconds, picoseconds, '
            'nanoseconds, and microseconds. Numeric values in the loaded '
            'data are NOT converted; only display labels and export '
            'headers change.'))
        self.ddTimeUnit.currentIndexChanged.connect(
            lambda _i: self.onTimeUnitChange(self.ddTimeUnit.currentData()))
        lay.addPair(label('Time unit:', right=True, color=(0.2, 0.2, 0.2)),
                    self.ddTimeUnit)
        self.ddSpecUnit = QtWidgets.QComboBox()
        self.ddSpecUnit.addItem('Wavelength (nm)', 'nm')
        self.ddSpecUnit.addItem('Wavenumber (cm-1)', 'cm-1')
        self.ddSpecUnit.setCurrentIndex(
            max(self.ddSpecUnit.findData(self.specUnit), 0))
        self.ddSpecUnit.setMinimumWidth(165)
        self.ddSpecUnit.setToolTip(tip(
            'Label the spectral axis as wavelength (nm; electronic TA / '
            'TRPL) or wavenumber (cm-1; TR-Raman / TR-IR). Numeric values '
            'in the loaded data are NOT converted; only labels and export '
            'headers change.'))
        self.ddSpecUnit.currentIndexChanged.connect(
            lambda _i: self.onSpecUnitChange(self.ddSpecUnit.currentData()))
        lay.addPair(label('Spectral axis:', right=True,
                          color=(0.2, 0.2, 0.2)), self.ddSpecUnit)
        lay.addWidget(label('full range shown when tmin >= tmax',
                            color=(0.45, 0.45, 0.45), italic=True))
        lay.addStretch(1)

        # ---- 2D map --------------------------------------------------
        self.pnl2D = panel(ui_str('2D Map  -  click to select (\\lambda, t)'))
        gl.addWidget(self.pnl2D, 3, 0, 2, 1)
        v = QtWidgets.QVBoxLayout(self.pnl2D)
        v.setContentsMargins(5, 8, 5, 5)
        v.setSpacing(3)
        ctrl2D, lay = hbox(5)
        v.addWidget(ctrl2D)
        lay.addWidget(label('Colormap:', right=True))
        self.ddColormap = QtWidgets.QComboBox()
        self.ddColormap.addItems(colormaps.NAMES)
        self.ddColormap.setCurrentText(colormaps.normName(self.mapColormap))
        self.ddColormap.setMinimumWidth(110)
        self.ddColormap.currentTextChanged.connect(self.onColormapChange)
        lay.addWidget(self.ddColormap)
        lay.addWidget(label('Z min:', right=True))
        self.edZMin = NumEdit(0, fmt='%.3g', width=90)
        self.edZMin.valueChanged.connect(lambda v_: self.onZChange('min', v_))
        lay.addWidget(self.edZMin)
        lay.addWidget(label('Z max:', right=True))
        self.edZMax = NumEdit(0, fmt='%.3g', width=90)
        self.edZMax.valueChanged.connect(lambda v_: self.onZChange('max', v_))
        lay.addWidget(self.edZMax)
        lay.addWidget(button('Reset Z', self.resetZRange, width=90))
        lay.addStretch(1)
        self.cv2D = PlotCanvas()
        v.addWidget(self.cv2D, 1)
        self.blit2D = BlitManager(self.cv2D)
        self.cv2D.mpl_connect('button_press_event', self.on2DClick)
        self.cv2D.mpl_connect('motion_notify_event', self.on2DMotion)
        self.cv2D.mpl_connect('button_release_event', self.on2DUp)

        # ---- spectrum ------------------------------------------------
        self.pnlSpec = panel('Spectrum  (at selected delay)')
        gl.addWidget(self.pnlSpec, 3, 1)
        v = QtWidgets.QVBoxLayout(self.pnlSpec)
        v.setContentsMargins(5, 8, 5, 5)
        v.setSpacing(3)
        self.cvSpec = PlotCanvas()
        v.addWidget(self.cvSpec, 1)
        ctrlS, lay = hbox(5)
        v.addWidget(ctrlS)
        self.lblSelT = label('Delay (%s):' % self.tUnitTxt())
        lay.addWidget(self.lblSelT)
        self.edSelT = NumEdit(0)
        self.edSelT.valueChanged.connect(lambda v_: self.setSelT(v_))
        lay.addWidget(self.edSelT, 1)
        lay.addWidget(button('Pin', self.pinSpec))
        lay.addWidget(button('Clear pins', self.clearSpecOverlays))

        # ---- kinetics ------------------------------------------------
        self.pnlKin = panel('Kinetics  (at selected wavelength)')
        gl.addWidget(self.pnlKin, 4, 1)
        v = QtWidgets.QVBoxLayout(self.pnlKin)
        v.setContentsMargins(5, 8, 5, 5)
        v.setSpacing(3)
        self.cvKin = PlotCanvas()
        v.addWidget(self.cvKin, 1)
        ctrlK, lay = hbox(5)
        v.addWidget(ctrlK)
        self.lblSelWL = label('%s (%s):' % (self.xSymTxt(), self.xUnitTxt()))
        lay.addWidget(self.lblSelWL)
        self.edSelWL = NumEdit(0)
        self.edSelWL.valueChanged.connect(lambda v_: self.setSelWL(v_))
        lay.addWidget(self.edSelWL, 1)
        lay.addWidget(button('Pin', self.pinKin))
        lay.addWidget(button('Clear pins', self.clearKinOverlays))
        b = colorButton(
            'Fit trace...', (0.55, 0.30, 0.70),
            tooltip='Fit a sum of exponentials (with optional IRF '
                    'convolution) to the kinetic trace at the currently '
                    'selected wavelength.')
        b.clicked.connect(lambda *_: self.openKineticFitWindow())
        lay.addWidget(b)

        self.build2DAxes()
        self.buildSpecAxes()
        self.buildKinAxes()

    def _closeEvent(self, ev):
        self.closeChildWindows()
        ev.accept()

    # ------------------------------------------------------------------
    # axes
    # ------------------------------------------------------------------
    def tLabel(self):
        return 'Delay time (%s)' % self.tUnitAx()

    def build2DAxes(self):
        fig = self.cv2D.figure
        fig.clear()
        self.hVLine = []
        self.hHLine = []
        self.hMap = []
        if self.delayScaleMode == 'split':
            axTop = fig.add_axes([0.10, 0.45, 0.75, 0.45])
            axBot = fig.add_axes([0.10, 0.10, 0.75, 0.25])
            self.ax2D = [axTop, axBot]
        else:
            self.ax2D = [fig.add_axes([0.10, 0.10, 0.75, 0.80])]
        self.cbarAx = fig.add_axes([0.90, 0.10, 0.02, 0.80])
        self.cbarAx.set_visible(False)
        n = len(self.ax2D)
        for k, ax in enumerate(self.ax2D):
            ax.set_facecolor(AX_BG)
            ax.set_ylabel(self.tLabel())
            styleAxes(ax)
            if k == n - 1:
                ax.set_xlabel(self.xAxLabel())
            else:
                # upper panel of the split map: the spectral axis is
                # labeled once, under the lower panel (as when data
                # are shown)
                ax.tick_params(labelbottom=False)
        self.blit2D.setArtists([])
        self.cv2D.layoutFcn = self._layout2D
        self.cv2D.relayout()
        self.cv2D.draw_idle()

    def _layout2D(self):
        cv = self.cv2D
        W, H = cv.sizePt()
        L, R, T, B = 50.0, 74.0, 20.0, 34.0
        w = W - L - R
        if len(self.ax2D) == 2:
            gap = 24.0
            avail = max(H - T - B - gap, 12.0)
            hBot = avail / 3.0
            cv.place(self.ax2D[0], L, B + hBot + gap, w, avail - hBot)
            cv.place(self.ax2D[1], L, B, w, hBot)
        elif self.ax2D:
            cv.place(self.ax2D[0], L, B, w, H - T - B)
        cv.yGroups = [list(self.ax2D)]
        if self.cbarAx is not None:
            cv.place(self.cbarAx, W - R + 8.0, B, 10.0, H - T - B)
        for ax in self.ax2D:              # titles end before the bar
            setSpan(ax, 0.0, W - R + 6.0)

    def buildSpecAxes(self):
        fig = self.cvSpec.figure
        fig.clear()
        self.axSpec = fig.add_axes([0.15, 0.2, 0.8, 0.7])
        self.axSpec.set_xlabel(self.xAxLabel())
        self.axSpec.set_ylabel(Z_LABEL)
        styleAxes(self.axSpec)
        self.cvSpec.layoutFcn = self._layoutSpec
        self.cvSpec.relayout()
        self.cvSpec.draw_idle()

    def _layoutSpec(self):
        W, H = self.cvSpec.sizePt()
        L, R, T, B = 50.0, 16.0, 20.0, 32.0
        self.cvSpec.place(self.axSpec, L, B, W - L - R, H - T - B)

    def buildKinAxes(self):
        fig = self.cvKin.figure
        fig.clear()
        if self.delayScaleMode == 'split':
            axL = fig.add_axes([0.15, 0.2, 0.25, 0.7])
            axR = fig.add_axes([0.45, 0.2, 0.50, 0.7])
            self.axKin = [axL, axR]
        else:
            self.axKin = [fig.add_axes([0.15, 0.2, 0.8, 0.7])]
        for k, ax in enumerate(self.axKin):
            ax.set_xlabel(self.tLabel())
            styleAxes(ax)
            if k == 0:
                ax.set_ylabel(Z_LABEL)
            else:
                # right part of the split plot: it shares the ordinate
                # with the left part (as when data are shown)
                ax.tick_params(labelleft=False)
        self.cvKin.layoutFcn = self._layoutKin
        self.cvKin.relayout()
        self.cvKin.draw_idle()

    def _layoutKin(self):
        cv = self.cvKin
        W, H = cv.sizePt()
        L, R, T, B = 50.0, 16.0, 20.0, 32.0
        h = H - T - B
        if len(self.axKin) == 2:
            gap = 20.0
            w = max(W - L - R - gap, 12.0)
            cv.place(self.axKin[0], L, B, w / 3.0, h)
            cv.place(self.axKin[1], L + w / 3.0 + gap, B, 2.0 * w / 3.0, h)
        elif self.axKin:
            cv.place(self.axKin[0], L, B, W - L - R, h)
        cv.xGroups = [list(self.axKin)]

    # ==================================================================
    # loading
    # ==================================================================
    def loadData(self):
        files, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self.fig, 'Select TA data file(s)', self.startDir(),
            'Data (*.csv *.dat *.txt);;All Files (*)')
        if not files:
            return
        self.rememberDir(os.path.dirname(files[0]))
        if len(files) == 1:
            name = os.path.basename(files[0])
            try:
                wl, t, A = parse_data_file(files[0])
            except Exception as err:
                self.alert('Failed to read file:\n%s' % err, 'Load Error')
                return
            from .win_filecheck import openSingleFileCheck
            wl2, t2, A2, cancel = openSingleFileCheck(self, wl, t, A, name)
            if cancel:
                return
            self.setLoadedData(wl2, t2, A2, 'Loaded: %s' % name)
        else:
            self.openLoadAverageWindow(os.path.dirname(files[0]),
                                       [os.path.basename(f) for f in files])

    def loadFile(self, path, transpose=False):
        """Load one file without any dialog."""
        from .win_filecheck import orient
        wl, t, A = parse_data_file(path)
        wl, t, A = orient(wl, t, A, transpose)
        self.setLoadedData(wl, t, A, 'Loaded: %s' % os.path.basename(path))

    def setLoadedData(self, wl, t, A, sourceDesc):
        self.setLoadedState(wl, t, A, sourceDesc)
        wl, t = self.wavelength, self.delay
        self.edSelWL.setLimits(np.min(wl), np.max(wl))
        self.edSelT.setLimits(np.min(t), np.max(t))
        self.edSelWL.setValue(self.selWL)
        self.edSelT.setValue(self.selT)
        self.edViewTMin.setValue(self.mainViewTMin)
        self.edViewTMax.setValue(self.mainViewTMax)
        fin = np.abs(self.deltaA[np.isfinite(self.deltaA)])
        cl = float(np.max(fin)) if fin.size else 1.0
        if cl == 0:
            cl = 1.0
        self.edZMin.setValue(-cl)
        self.edZMax.setValue(cl)
        self.updateStatusLabel()
        self.dataChanged()
        self.resetPlotViews()
        self.updateAll()

    def resetPlotViews(self):
        """Full view of the three plots (other data: a zoom of the old
        ones means nothing)."""
        for cv in (self.cv2D, self.cvSpec, self.cvKin):
            zoom.reset(cv)

    def updateStatusLabel(self):
        self.lblStatus.setText(self.statusText())

    def dataChanged(self):
        """The working data changed: the analysis windows are closed
        and the results that rest on the old data are dropped."""
        self.closeAnalysisWindows()
        self.invalidateChildResults()

    def closeAnalysisWindows(self):
        for key in ANALYSIS_KEYS:
            w = self.childWindows.get(key)
            if w is not None:
                try:
                    w.close()
                except Exception:
                    pass
                self.childWindows.pop(key, None)

    def closeChildWindows(self):
        for key in list(self.childWindows):
            w = self.childWindows.get(key)
            try:
                w.close()
            except Exception:
                pass
            self.childWindows.pop(key, None)

    def invalidateChildResults(self):
        self.clearChildResults()

    def _syncControlsToData(self, zFromData=None):
        """Limits and values of the edit fields after the axes of the
        data changed (selection clamped, full display range)."""
        wlLo, wlHi = float(np.min(self.wavelength)), float(np.max(self.wavelength))
        tLo, tHi = float(np.min(self.delay)), float(np.max(self.delay))
        self.selWL = float(min(max(self.selWL, wlLo), wlHi))
        self.selT = float(min(max(self.selT, tLo), tHi))
        self.edSelWL.setLimits(wlLo, wlHi)
        self.edSelT.setLimits(tLo, tHi)
        self.edSelWL.setValue(self.selWL)
        self.edSelT.setValue(self.selT)
        self.mainViewTMin, self.mainViewTMax = tLo, tHi
        self.edViewTMin.setValue(tLo)
        self.edViewTMax.setValue(tHi)
        self.mapZMin = None
        self.mapZMax = None
        if zFromData is not None:
            fin = np.abs(zFromData[np.isfinite(zFromData)])
            cl = float(np.max(fin)) if fin.size else 1.0
            if cl == 0:
                cl = 1.0
            self.edZMin.setValue(-cl)
            self.edZMax.setValue(cl)
        self.specOverlays = [v for v in self.specOverlays if tLo <= v <= tHi]
        self.kinOverlays = [v for v in self.kinOverlays if wlLo <= v <= wlHi]

    def _clearCorrectionState(self):
        self.bgApplied = False
        self.bgSpectrum = None
        self.bgN = 0
        self.chirpApplied = False
        self.chirpPts = np.zeros((0, 2))
        self.chirpParams = None
        self.chirpFitRMS = np.nan
        self.ridgeWL = np.zeros(0)
        self.ridgeT = np.zeros(0)
        self.maskedRegions = []
        self.tZeroShift = 0.0

    # ==================================================================
    # corrections that need no window of their own
    # ==================================================================
    def resetCorrections(self):
        """Back to the data as loaded (undoes crop, background, chirp,
        mask and time-zero shift)."""
        if self.originalDeltaA is None:
            return
        self.dataChanged()
        self.wavelength = self.originalWavelength.copy()
        self.delay = self.originalDelay.copy()
        self.deltaA_raw = self.originalDeltaA.copy()
        self._clearCorrectionState()
        self._syncControlsToData(self.originalDeltaA)
        self.updateStatusLabel()
        self.resetPlotViews()
        self.updateAll()

    def applyZeroTimeShift(self, delta):
        """Shift the delay axis so that delta becomes 0."""
        if not np.isfinite(delta) or delta == 0 or self.delay.size == 0:
            return
        delta = float(delta)
        self.dataChanged()
        self.delay = self.delay - delta
        self.tZeroShift += delta
        self.specOverlays = [v - delta for v in self.specOverlays]
        for nm in ('gaFitTMin', 'gaFitTMax', 'kfFitTMin', 'kfFitTMax'):
            if getattr(self, nm) is not None:
                setattr(self, nm, getattr(self, nm) - delta)
        if isinstance(self.tgSet, dict) and \
                self.tgSet.get('fitTMin') is not None:
            self.tgSet['fitTMin'] = self.tgSet['fitTMin'] - delta
            if self.tgSet.get('fitTMax') is not None:
                self.tgSet['fitTMax'] = self.tgSet['fitTMax'] - delta
        self.selT = self.selT - delta
        if self.mainViewTMin is not None:
            self.mainViewTMin -= delta
            self.mainViewTMax -= delta
        tLo, tHi = float(np.min(self.delay)), float(np.max(self.delay))
        self.edSelT.setLimits(tLo, tHi)
        self.selT = float(max(min(self.selT, tHi), tLo))
        self.edSelT.setValue(self.selT)
        if self.mainViewTMin is not None:
            self.edViewTMin.setValue(self.mainViewTMin)
            self.edViewTMax.setValue(self.mainViewTMax)
        self.updateAll()

    def applyZeroTimeAtCrosshair(self):
        if not self.hasData():
            self.alert('Load data first.', 'No data')
            return
        delta = self.selT
        if delta is None or not np.isfinite(delta):
            return
        if delta == 0:
            self.alert('Crosshair is already at t = 0. Click somewhere else '
                       'on the 2D map first to pick a new t_0 location.',
                       'No shift needed')
            return
        self.applyZeroTimeShift(delta)
        self.lblStatus.setText(
            't_0 shift applied: %.4g %s  (total: %.4g %s)' % (
                delta, self.tUnitTxt(), self.tZeroShift, self.tUnitTxt()))

    def applyCropByRange(self, wlMin, wlMax, tMin, tMax):
        """Keep the data inside the two ranges; the result becomes the
        new baseline (corrections made so far are baked in)."""
        if not self.hasData():
            return
        if wlMin > wlMax:
            wlMin, wlMax = wlMax, wlMin
        if tMin > tMax:
            tMin, tMax = tMax, tMin
        wlMask = (self.wavelength >= wlMin) & (self.wavelength <= wlMax)
        tMask = (self.delay >= tMin) & (self.delay <= tMax)
        if not wlMask.any() or not tMask.any():
            self.alert('Empty crop range.', 'Invalid range')
            return
        if np.count_nonzero(wlMask) < 2 or np.count_nonzero(tMask) < 2:
            self.alert('The crop range must keep at least 2 points on each '
                       'axis.', 'Invalid range')
            return
        self.dataChanged()
        cropped = self.deltaA[wlMask, :][:, tMask].copy()
        self.wavelength = self.wavelength[wlMask]
        self.delay = self.delay[tMask]
        self.deltaA = cropped
        self.deltaA_raw = cropped.copy()
        self._clearCorrectionState()
        self._syncControlsToData(self.deltaA_raw)
        self.updateStatusLabel()
        self.updateAll()

    # ==================================================================
    # files: export, session
    # ==================================================================
    def putFile(self, filterSpec, dlgTitle, defName, parent=None):
        """Ask for a file to write; '' if canceled. The folder is
        remembered."""
        if self.headless:
            fn = self.testPath or ''
        else:
            fn, sel = QtWidgets.QFileDialog.getSaveFileName(
                parent or self.fig, dlgTitle,
                os.path.join(self.startDir(), defName), filterSpec)
            if fn and not os.path.splitext(fn)[1]:
                i = sel.rfind('*.')
                if i >= 0:
                    fn += sel[i + 1:].rstrip(') ').split()[0]
        if fn:
            self.rememberDir(os.path.dirname(fn))
        return fn

    def getFile(self, filterSpec, dlgTitle, parent=None):
        if self.headless:
            fn = self.testPath or ''
        else:
            fn, _ = QtWidgets.QFileDialog.getOpenFileName(
                parent or self.fig, dlgTitle, self.startDir(), filterSpec)
        if fn:
            self.rememberDir(os.path.dirname(fn))
        return fn

    def exportCorrected(self):
        if not self.hasData():
            self.alert('Load data first.', 'No data')
            return
        parts = []
        if self.bgApplied:
            parts.append('BG%d' % self.bgN)
        if self.chirpApplied:
            parts.append('chirp')
        if self.isCropped():
            parts.append('cropped')
        suffix = '_'.join(parts) if parts else 'raw'
        nameOnly = self.dataStem() or 'data'
        fullpath = self.putFile(
            'CSV (comma-delimited) (*.csv);;DAT (tab-delimited) (*.dat);;'
            'TXT (tab-delimited) (*.txt)', 'Export corrected data',
            '%s_%s.csv' % (nameOnly, suffix))
        if not fullpath:
            return
        delim = ',' if fullpath.lower().endswith('.csv') else '\t'
        try:
            write_data_file(fullpath, self.wavelength, self.delay,
                            self.deltaA, delim)
            self.inform('Saved %d wavelengths x %d delays.\n\nPipeline: %s'
                        '\n\nFile: %s' % (self.wavelength.size,
                                           self.delay.size, suffix, fullpath),
                        'Export complete')
        except Exception as err:
            self.alert('Export failed:\n%s' % err, 'Export error')

    SESSION_FILTER = 'TRSpecAnalyzer session, Python version (*.trsession)'

    def saveSession(self):
        if not self.hasData():
            self.alert('Load data first.', 'No data')
            return
        base = self.dataStem()
        defName = ('TRSpec_%s_session.trsession' % base if base
                   else 'TRSpec_session.trsession')
        fn = self.putFile(self.SESSION_FILTER, 'Save session', defName)
        if not fn:
            return
        try:
            self.sessionSave(fn)
        except Exception as err:
            self.alert('Save failed:\n%s' % err, 'Session')
            return
        self.lblStatus.setText('Session saved: %s' % os.path.basename(fn))

    def loadSessionFile(self, path=None):
        fn = path or self.getFile(
            self.SESSION_FILTER + ';;All Files (*)', 'Load session')
        if not fn:
            return
        try:
            d = self.sessionLoad(fn)
        except ValueError as err:
            self.alert(str(err), 'Session')
            return
        except Exception as err:
            self.alert('Load failed:\n%s' % err, 'Session')
            return
        self.closeChildWindows()
        self.sessionApply(d)
        self.mapColormap = colormaps.normName(self.mapColormap)
        B = QtCore.QSignalBlocker
        blockers = [B(w) for w in (self.ddScale, self.ddColormap,
                                   self.ddTimeUnit, self.ddSpecUnit)]
        self.ddScale.setCurrentText(self.scaleModeLabel())
        self.ddColormap.setCurrentText(self.mapColormap)
        self.ddTimeUnit.setCurrentIndex(
            max(self.ddTimeUnit.findData(self.timeUnit), 0))
        self.ddSpecUnit.setCurrentIndex(
            max(self.ddSpecUnit.findData(self.specUnit), 0))
        del blockers
        self.edThreshold.setValue(self.splitThreshold)
        self.edViewTMin.setValue(self.mainViewTMin
                                 if self.mainViewTMin is not None else 0)
        self.edViewTMax.setValue(self.mainViewTMax
                                 if self.mainViewTMax is not None else 0)
        self.refreshTimeUnitLabels()
        self.refreshSpecUnitLabels()
        if self.wavelength.size and self.delay.size:
            wlLo, wlHi = float(np.min(self.wavelength)), float(np.max(self.wavelength))
            tLo, tHi = float(np.min(self.delay)), float(np.max(self.delay))
            if self.selWL is None or not np.isfinite(self.selWL):
                self.selWL = float(self.wavelength[
                    max(int(np.floor(self.wavelength.size / 2 + 0.5)) - 1, 0)])
            if self.selT is None or not np.isfinite(self.selT):
                self.selT = float(self.delay[0])
            self.selWL = float(min(max(self.selWL, wlLo), wlHi))
            self.selT = float(min(max(self.selT, tLo), tHi))
            self.edSelWL.setLimits(wlLo, wlHi)
            self.edSelT.setLimits(tLo, tHi)
            self.edSelWL.setValue(self.selWL)
            self.edSelT.setValue(self.selT)
        cMin, cMax = self.computeZLim()
        self.edZMin.setValue(cMin)
        self.edZMax.setValue(cMax)
        self.holdDataChange = True
        try:
            self.build2DAxes()
            self.buildKinAxes()
            self.updateStatusLabel()
            self.resetPlotViews()
            self.updateAll()
        finally:
            self.holdDataChange = False
        self.lblStatus.setText('%s   [session: %s]' % (
            self.lblStatus.text(), os.path.basename(fn)))

    # ==================================================================
    # correction windows
    # ==================================================================
    def _needData(self):
        if not self.hasData():
            self.alert('Load data first.', 'No data')
            return False
        return True

    def openBackgroundWindow(self):
        if self._needData():
            from .win_background import BackgroundWindow
            return openChild(self, BackgroundWindow)

    def openCropWindow(self):
        if self._needData():
            from .win_crop import CropWindow
            return openChild(self, CropWindow)

    def openMaskWindow(self):
        if self._needData():
            from .win_mask import MaskWindow
            return openChild(self, MaskWindow)

    def openChirpWindow(self):
        if self._needData():
            from .win_chirp import ChirpWindow
            return openChild(self, ChirpWindow)

    def openLoadAverageWindow(self, pathBase, fileNames):
        from .win_average import LoadAverageWindow, readScans
        files = readScans(self, pathBase, fileNames)
        if files is None:
            return None
        old = self.childWindows.get('load')
        if old is not None:
            old.close()
        return openChild(self, LoadAverageWindow, pathBase, files)

    # ==================================================================
    # drawing
    # ==================================================================
    def updateAll(self):
        if self.recompute() and not self.holdDataChange:
            self.dataChanged()
        self.drawMap2D()
        self.drawSpectrum()
        self.drawKinetics()

    def getColormap(self):
        return colormaps.getColormap(self.mapColormap)

    def scaleMasks(self, forKin=False):
        """Delay masks and axis scales of the current scale mode (map:
        top/bottom, kinetics: left/right)."""
        t = self.delay
        mode = self.delayScaleMode
        if mode == 'linear':
            return [np.ones(t.shape, bool)], ['linear']
        if mode == 'log':
            pos = t > 0
            if forKin and not pos.any():
                pos = np.ones(t.shape, bool)
            return [pos], ['log']
        thr = self.splitThreshold
        if forKin:
            return [t <= thr, t >= thr], ['linear', 'log']
        return [t >= thr, t <= thr], ['log', 'linear']

    def plotHeatmapOn(self, ax, t, data, wl=None):
        """Color map of data (n_wl x n_t) over (wl, t). Each cell
        takes the color of its lower left corner, as MATLAB's surface
        with flat faces does; NaN cells stay empty."""
        wl = self.wavelength if wl is None else wl
        t = np.asarray(t, float)
        if wl.size < 2 or t.size < 2:
            return None
        C = np.ma.masked_invalid(np.asarray(data, float).T[:-1, :-1])
        h = ax.pcolormesh(wl, t, C, shading='flat')
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        ax.set_ylim(float(np.min(t)), float(np.max(t)))
        return h

    def drawMap2D(self):
        if not self.hasData():
            return
        t = self.delay
        masks, scales = self.scaleMasks()
        cMin, cMax = self.computeZLim()
        cmap = self.getColormap()
        self.hVLine = []
        self.hHLine = []
        self.hMap = []
        n = len(self.ax2D)
        for k, ax in enumerate(self.ax2D):
            ax.clear()
            styleAxes(ax)
            tSub = t[masks[k]]
            if tSub.size == 0:
                continue
            dSub = self.deltaA[:, masks[k]]
            h = self.plotHeatmapOn(ax, tSub, dSub)
            if h is not None:
                h.set_cmap(cmap)
                h.set_clim(cMin, cMax)
                self.hMap.append(h)
            ax.set_yscale(scales[k])
            ylo, yhi = self.clipViewToMask(tSub, scales[k])
            ax.set_ylim(ylo, yhi)
            if k == n - 1:
                ax.set_xlabel(self.xAxLabel())
            else:
                ax.set_xlabel('')
                ax.tick_params(labelbottom=False)
            ax.set_ylabel(self.tLabel())
            ax.set_title(self.mapTitle(scales[k]), fontsize=FONT_PT + 1)
            hv = ax.axvline(self.selWL, color='k', linewidth=1, alpha=0.7)
            hh = ax.axhline(self.selT, color='k', linewidth=1, alpha=0.7)
            hh.set_visible(not (scales[k] == 'log' and self.selT <= 0))
            self.hVLine.append(hv)
            self.hHLine.append(hh)
        self.blit2D.setArtists(self.hVLine + self.hHLine)
        self._drawColorbar(cMin, cMax, cmap)
        self.cv2D.draw_idle()

    def _drawColorbar(self, cMin, cMax, cmap):
        cax = self.cbarAx
        cax.clear()
        cax.set_visible(True)
        sm = ScalarMappable(norm=Normalize(cMin, cMax), cmap=cmap)
        cb = self.cv2D.figure.colorbar(sm, cax=cax)
        cb.set_label(Z_LABEL, fontsize=FONT_PT)
        cb.ax.tick_params(labelsize=FONT_PT, length=3)
        sciTicks(cb.ax, 'y')
        self.cv2D.relayout()

    def applyZLim(self):
        self.drawMap2D()

    def drawSpectrum(self, fast=False):
        """Spectrum at the selected delay, with the pinned ones. With
        fast=True only the current curve and the marker are moved (used
        while the crosshair is dragged)."""
        if not self.hasData():
            return
        if fast and self._fastSpectrum():
            return
        ax = self.axSpec
        ax.clear()
        styleAxes(ax, grid=True)
        legH, legL = [], []
        for k, tv in enumerate(self.specOverlays):
            idx = int(np.argmin(np.abs(self.delay - tv)))
            hO, = ax.plot(self.wavelength, self.deltaA[:, idx], '-',
                          color=lineColor(k), linewidth=1.0)
            legH.append(hO)
            legL.append('t = %.3g %s' % (self.delay[idx], self.tUnitAx()))
        idx = int(np.argmin(np.abs(self.delay - self.selT)))
        hCur, = ax.plot(self.wavelength, self.deltaA[:, idx], 'b-',
                        linewidth=1.6)
        legH.append(hCur)
        legL.append(self._specCurLabel(idx))
        hZ = ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        hV = ax.axvline(self.selWL, color='k', linewidth=1.0, alpha=0.7)
        ax.set_xlabel(self.xAxLabel())
        ax.set_ylabel(Z_LABEL)
        ax.set_xlim(float(np.min(self.wavelength)),
                    float(np.max(self.wavelength)))
        sciTicks(ax, 'y')
        leg = None
        if len(legH) > 1:
            leg = legend(ax, legH, legL, fontsize=8)
        else:
            ax.set_title(legL[-1], fontsize=FONT_PT + 1)
        self._specH = {'ax': ax, 'cur': hCur, 'v': hV, 'zero': hZ,
                       'leg': leg}
        self.cvSpec.draw_idle()

    def _specCurLabel(self, idx):
        return 't = %.3g %s (current)' % (self.delay[idx], self.tUnitAx())

    def _fastSpectrum(self):
        H = getattr(self, '_specH', None)
        if not H or H['ax'] is not self.axSpec or H['cur'].axes is None:
            return False
        ax = self.axSpec
        idx = int(np.argmin(np.abs(self.delay - self.selT)))
        H['cur'].set_ydata(self.deltaA[:, idx])
        H['v'].set_xdata([self.selWL, self.selWL])
        if H['leg'] is not None:
            setLegendText(ax, -1, self._specCurLabel(idx))
        else:
            ax.set_title(self._specCurLabel(idx), fontsize=FONT_PT + 1)
        self._autoY(ax, H['zero'])
        self.cvSpec.draw_idle()
        return True

    @staticmethod
    def _autoY(ax, hZero):
        """y-limits as a fresh plot gets them: from the curves, and
        from the zero line only when zero lies outside them."""
        hZero.set_visible(False)
        ax.relim(visible_only=True)
        ax.autoscale_view(scalex=False, scaley=True)
        lo, hi = ax.get_ybound()
        hZero.set_visible(True)
        if not (lo <= 0 <= hi):
            ax.relim(visible_only=True)
            ax.autoscale_view(scalex=False, scaley=True)

    def _kinCurLabel(self, idx, current=True):
        s = '%s = %.1f %s' % (self.xSym(), self.wavelength[idx],
                              self.xUnitAx())
        return s + ' (current)' if current else s

    def _kinSplitYLim(self, masks):
        """Common y-limits of the two kinetics axes of the split view."""
        yLo, yHi = np.inf, -np.inf
        anyMask = np.zeros(self.delay.shape, bool)
        for m in masks:
            anyMask |= m
        for wv in list(self.kinOverlays) + [self.selWL]:
            idx = int(np.argmin(np.abs(self.wavelength - wv)))
            yj = self.deltaA[idx, anyMask]
            yj = yj[np.isfinite(yj)]
            if yj.size == 0:
                continue
            yLo = min(yLo, float(np.min(yj)))
            yHi = max(yHi, float(np.max(yj)))
        if np.isfinite(yLo) and np.isfinite(yHi) and yLo < yHi:
            pad = 0.05 * (yHi - yLo)
            for ax in self.axKin:
                ax.set_ylim(yLo - pad, yHi + pad)

    def drawKinetics(self, fast=False):
        """Kinetics at the selected wavelength, with the pinned ones
        (fast=True: see drawSpectrum)."""
        if not self.hasData():
            return
        if fast and self._fastKinetics():
            return
        masks, scales = self.scaleMasks(forKin=True)
        nAx = len(self.axKin)
        idxCur = int(np.argmin(np.abs(self.wavelength - self.selWL)))
        kLeg = nAx - 1                 # the legend goes on the widest axes
        handles = []
        leg = None
        for k, ax in enumerate(self.axKin):
            ax.clear()
            styleAxes(ax, grid=True)
            m = masks[k]
            tSub = self.delay[m]
            if tSub.size == 0:
                handles.append(None)
                continue
            legH, legL = [], []
            for j, wv in enumerate(self.kinOverlays):
                idx = int(np.argmin(np.abs(self.wavelength - wv)))
                hO, = ax.plot(tSub, self.deltaA[idx, m], '-',
                              color=lineColor(j), linewidth=1.0)
                legH.append(hO)
                legL.append(self._kinCurLabel(idx, False))
            hCur, = ax.plot(tSub, self.deltaA[idxCur, m], 'r-',
                            linewidth=1.6)
            legH.append(hCur)
            legL.append(self._kinCurLabel(idxCur))
            hZ = ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
            hV = ax.axvline(self.selT, color='k', linewidth=1.0, alpha=0.7)
            hV.set_visible(bool(np.min(tSub) <= self.selT <= np.max(tSub)))
            ax.set_xlabel(self.tLabel())
            if k == 0:
                ax.set_ylabel(Z_LABEL)
            else:
                ax.tick_params(labelleft=False)
            ax.set_xscale(scales[k])
            xlo, xhi = self.clipViewToMask(tSub, scales[k])
            ax.set_xlim(xlo, xhi)
            if scales[k] == 'log' and xlo > 0 and xhi > xlo:
                decLo = int(np.ceil(np.log10(xlo)))
                decHi = int(np.floor(np.log10(xhi)))
                if decHi >= decLo:
                    ax.set_xticks(10.0 ** np.arange(decLo, decHi + 1))
            sciTicks(ax, 'y')
            if k == kLeg and len(legH) > 1:
                leg = legend(ax, legH, legL, fontsize=8)
            if nAx == 1:
                ax.set_title(self._kinCurLabel(idxCur, False),
                             fontsize=FONT_PT + 1)
            handles.append({'ax': ax, 'cur': hCur, 'v': hV, 'zero': hZ,
                            'mask': m,
                            'lo': float(np.min(tSub)),
                            'hi': float(np.max(tSub))})
        if nAx == 2:
            self._kinSplitYLim(masks)
        self._kinH = {'axes': list(self.axKin), 'h': handles, 'leg': leg,
                      'legAx': self.axKin[kLeg], 'masks': masks}
        self.cvKin.draw_idle()

    def _fastKinetics(self):
        H = getattr(self, '_kinH', None)
        if not H or len(H['axes']) != len(self.axKin):
            return False
        for a, b in zip(H['axes'], self.axKin):
            if a is not b:
                return False
        nAx = len(self.axKin)
        idxCur = int(np.argmin(np.abs(self.wavelength - self.selWL)))
        for h in H['h']:
            if h is None:
                continue
            if h['cur'].axes is None:
                return False
            h['cur'].set_ydata(self.deltaA[idxCur, h['mask']])
            h['v'].set_xdata([self.selT, self.selT])
            h['v'].set_visible(bool(h['lo'] <= self.selT <= h['hi']))
            if nAx == 1:
                h['ax'].set_title(self._kinCurLabel(idxCur, False),
                                  fontsize=FONT_PT + 1)
                self._autoY(h['ax'], h['zero'])
        if H['leg'] is not None:
            setLegendText(H['legAx'], -1, self._kinCurLabel(idxCur))
        if nAx == 2:
            self._kinSplitYLim(H['masks'])
        self.cvKin.draw_idle()
        return True

    def updateCrosshairs(self):
        for k, h in enumerate(self.hVLine):
            h.set_xdata([self.selWL, self.selWL])
        for k, h in enumerate(self.hHLine):
            h.set_ydata([self.selT, self.selT])
            isLog = h.axes.get_yscale() == 'log'
            h.set_visible(not (isLog and self.selT <= 0))
        self.blit2D.update()

    # ==================================================================
    # mouse on the 2D map
    # ==================================================================
    def on2DClick(self, ev):
        if not self.hasData() or ev.button != 1:
            return
        if ev.inaxes not in self.ax2D:
            return
        self.drag2DAx = ev.inaxes
        self.dragging2D = True
        self.applyMain2DPoint((ev.xdata, ev.ydata))

    def on2DMotion(self, ev):
        if not self.dragging2D:
            return
        ax = self.drag2DAx
        if getattr(ev, 'inaxes', None) in self.ax2D:
            ax = ev.inaxes          # pointer moved onto the other panel
        if ax is None or ax not in self.ax2D:
            if not self.ax2D:
                return
            ax = self.ax2D[0]
        try:
            x, y = ax.transData.inverted().transform((ev.x, ev.y))
        except Exception:
            return
        if not (np.isfinite(x) and np.isfinite(y)):
            return
        self.applyMain2DPoint((x, y))

    def on2DUp(self, _ev=None):
        self.dragging2D = False

    def applyMain2DPoint(self, pt):
        wl, t = float(pt[0]), float(pt[1])
        wl = min(max(wl, np.min(self.wavelength)), np.max(self.wavelength))
        t = min(max(t, np.min(self.delay)), np.max(self.delay))
        self.setSelWL(wl, False)
        self.setSelT(t, False)
        self.updateCrosshairs()
        self.drawSpectrum(fast=True)
        self.drawKinetics(fast=True)

    def setSelWL(self, val, doUpdate=True):
        if self.wavelength.size == 0:
            return
        val = float(min(max(val, np.min(self.wavelength)),
                        np.max(self.wavelength)))
        self.selWL = val
        self.edSelWL.setValue(val)
        if doUpdate:
            self.updateCrosshairs()
            self.drawSpectrum(fast=True)
            self.drawKinetics(fast=True)

    def setSelT(self, val, doUpdate=True):
        if self.delay.size == 0:
            return
        val = float(min(max(val, np.min(self.delay)), np.max(self.delay)))
        self.selT = val
        self.edSelT.setValue(val)
        if doUpdate:
            self.updateCrosshairs()
            self.drawSpectrum(fast=True)
            self.drawKinetics(fast=True)

    # ==================================================================
    # view controls
    # ==================================================================
    def onScaleModeChange(self, val):
        self.delayScaleMode = str(val).lower()
        self.build2DAxes()
        self.buildKinAxes()
        self.updateAll()
        self.savePrefs()

    def onThresholdChange(self, val):
        if not np.isfinite(val) or val <= 0:
            return
        self.splitThreshold = float(val)
        if self.delayScaleMode == 'split':
            self.updateAll()
        self.savePrefs()

    def onMainViewChange(self, which, val):
        if which == 'min':
            self.mainViewTMin = float(val)
        else:
            self.mainViewTMax = float(val)
        self.drawMap2D()
        self.drawKinetics()

    def resetMainView(self):
        if self.delay.size == 0:
            self.mainViewTMin = None
            self.mainViewTMax = None
            return
        self.mainViewTMin = float(self.delay[0])
        self.mainViewTMax = float(self.delay[-1])
        self.edViewTMin.setValue(self.mainViewTMin)
        self.edViewTMax.setValue(self.mainViewTMax)
        self.drawMap2D()
        self.drawKinetics()

    def onColormapChange(self, val):
        self.mapColormap = val
        self.drawMap2D()
        for w in self.childWindows.values():
            fn = getattr(w, 'onColormapChange', None)
            if fn is not None:
                try:
                    fn()
                except Exception:
                    pass
        self.savePrefs()

    def onZChange(self, which, val):
        if which == 'min':
            newMin, newMax = val, self.edZMax.value()
        else:
            newMin, newMax = self.edZMin.value(), val
        if newMin >= newMax:
            cMin, cMax = self.computeZLim()
            self.edZMin.setValue(cMin)
            self.edZMax.setValue(cMax)
            return
        self.mapZMin, self.mapZMax = float(newMin), float(newMax)
        self.applyZLim()

    def resetZRange(self):
        self.mapZMin = None
        self.mapZMax = None
        cMin, cMax = self.computeZLim()
        self.edZMin.setValue(cMin)
        self.edZMax.setValue(cMax)
        self.applyZLim()

    def onTimeUnitChange(self, newUnit):
        if self.timeUnit == newUnit:
            return
        self.timeUnit = newUnit
        if self.delay.size:
            self.chirpViewTMax = min(self.defaultChirpTMax(),
                                     float(np.max(self.delay)))
        else:
            self.chirpViewTMax = self.defaultChirpTMax()
        self.closeChildWindows()
        self.invalidateChildResults()
        self.refreshTimeUnitLabels()
        self._relabelEmptyAxes()
        self.updateAll()
        self.savePrefs()

    def onSpecUnitChange(self, newUnit):
        if self.specUnit == newUnit:
            return
        self.specUnit = newUnit
        self.closeChildWindows()
        self.invalidateChildResults()
        self.refreshSpecUnitLabels()
        self._relabelEmptyAxes()
        self.updateAll()
        self.savePrefs()

    def _relabelEmptyAxes(self):
        """Axis labels while no data are loaded (nothing is redrawn
        then)."""
        if self.hasData():
            return
        n = len(self.ax2D)
        for k, ax in enumerate(self.ax2D):
            if k == n - 1:          # (the upper panel of the split map
                ax.set_xlabel(self.xAxLabel())      # carries none)
            ax.set_ylabel(self.tLabel())
        self.axSpec.set_xlabel(self.xAxLabel())
        for ax in self.axKin:
            ax.set_xlabel(self.tLabel())
        for cv in (self.cv2D, self.cvSpec, self.cvKin):
            cv.draw_idle()

    def refreshTimeUnitLabels(self):
        u = self.tUnitTxt()
        self.lblSplitThreshold.setText('Split thr. (%s):' % u)
        self.lblViewTMin.setText('Display t_min (%s):' % u)
        self.lblViewTMax.setText('t_max (%s):' % u)
        self.lblSelT.setText('Delay (%s):' % u)

    def refreshSpecUnitLabels(self):
        self.lblSelWL.setText('%s (%s):' % (self.xSymTxt(), self.xUnitTxt()))

    # ==================================================================
    # pins
    # ==================================================================
    def pinSpec(self):
        if not self.hasData():
            return
        idx = int(np.argmin(np.abs(self.delay - self.selT)))
        v = float(self.delay[idx])
        if v not in self.specOverlays:
            self.specOverlays.append(v)
        self.drawSpectrum()

    def pinKin(self):
        if not self.hasData():
            return
        idx = int(np.argmin(np.abs(self.wavelength - self.selWL)))
        v = float(self.wavelength[idx])
        if v not in self.kinOverlays:
            self.kinOverlays.append(v)
        self.drawKinetics()

    def clearSpecOverlays(self):
        self.specOverlays = []
        self.drawSpectrum()

    def clearKinOverlays(self):
        self.kinOverlays = []
        self.drawKinetics()

    # ==================================================================
    # preferences (kept between runs)
    # ==================================================================
    def prefsFile(self):
        d = os.environ.get('TRSPEC_PREFS_DIR')
        if not d:
            base = QtCore.QStandardPaths.writableLocation(
                QtCore.QStandardPaths.StandardLocation.GenericConfigLocation)
            d = os.path.join(base or os.path.expanduser('~'),
                             'TRSpecAnalyzer')
        return os.path.join(d, 'prefs_python.json')

    def loadPrefs(self):
        try:
            with open(self.prefsFile(), 'r', encoding='utf-8') as fh:
                p = json.load(fh)
        except Exception:
            return
        if not isinstance(p, dict):
            return
        v = p.get('lastDir')
        if isinstance(v, str) and v and os.path.isdir(v):
            self.prefLastDir = v
        v = p.get('timeUnit')
        if v in ('fs', 'ps', 'ns', 'us'):
            self.timeUnit = v
            self.chirpViewTMax = self.defaultChirpTMax()
        v = p.get('mapColormap')
        if isinstance(v, str) and (v in colormaps.NAMES or v == 'parula'):
            self.mapColormap = colormaps.normName(v)
        v = p.get('splitThreshold')
        if isinstance(v, (int, float)) and np.isfinite(v) and v > 0:
            self.splitThreshold = float(v)
        v = p.get('delayScaleMode')
        if v in ('linear', 'log', 'split'):
            self.delayScaleMode = v
        v = p.get('specUnit')
        if v in ('nm', 'cm-1'):
            self.specUnit = v

    def savePrefs(self):
        try:
            fn = self.prefsFile()
            os.makedirs(os.path.dirname(fn), exist_ok=True)
            with open(fn, 'w', encoding='utf-8') as fh:
                json.dump({'lastDir': self.prefLastDir,
                           'timeUnit': self.timeUnit,
                           'mapColormap': self.mapColormap,
                           'delayScaleMode': self.delayScaleMode,
                           'splitThreshold': self.splitThreshold,
                           'specUnit': self.specUnit}, fh, indent=1)
        except Exception:
            pass

    def startDir(self):
        d = self.prefLastDir
        if not d or not os.path.isdir(d):
            d = os.getcwd()
        return d

    def rememberDir(self, pth):
        if isinstance(pth, str) and pth and os.path.isdir(pth):
            self.prefLastDir = pth
            self.savePrefs()

    # ==================================================================
    # analysis windows that are not ported yet
    # ==================================================================
    def openKineticFitWindow(self):
        if self._needData():
            from .win_kinfit import KineticFitWindow
            return openChild(self, KineticFitWindow)

    def openSvdWindow(self):
        if self._needData():
            from .win_svd import SvdWindow
            w = openChild(self, SvdWindow)
            if not getattr(w, 'ok', True):
                w.close()
                return None
            return w

    def openGlobalAnalysisWindow(self):
        if self._needData():
            from .win_ga import GlobalAnalysisWindow
            return openChild(self, GlobalAnalysisWindow)

    def openTargetWindow(self):
        if self._needData():
            from .win_target import TargetWindow
            return openChild(self, TargetWindow)

    def openLdaWindow(self):
        if self._needData():
            from .win_lda import LdaWindow
            return openChild(self, LdaWindow)

    def openMcrWindow(self):
        if self._needData():
            from .win_mcr import McrWindow
            return openChild(self, McrWindow)

    def openCoherenceWindow(self):
        if self._needData():
            from .win_coherence import CoherenceWindow
            return openChild(self, CoherenceWindow)

    def openTimeFrequencyWindow(self):
        from .win_coherence import openTimeFrequency
        return openTimeFrequency(self)

    def openAutoWindow(self):
        if self._needData():
            from .win_auto import AutoWindow
            return openChild(self, AutoWindow)
