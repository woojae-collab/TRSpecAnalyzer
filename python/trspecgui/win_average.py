"""Load & Average Datasets window. Port of openLoadAverageWindow."""
import os
import warnings

import numpy as np
from PySide6 import QtCore, QtWidgets

from trspec.dataset import estimate_scan_shifts

from .state import detect_axis_swap, parse_data_file, ui_str, write_data_file
from .textfit import legend
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label, lineColor, sciTicks, styleAxes, tip)
from .win_base import GREEN, ChildWindow, absMax
from .win_filecheck import orient

C_FILE, C_MAX, C_INC, C_SHIFT, C_REF = range(5)


def readScans(app, pathBase, fileNames):
    """Parse the files; returns the list of scans, or None after an
    error message (unreadable file, sizes that differ)."""
    files = []
    for k, name in enumerate(fileNames):
        try:
            wl, t, A = parse_data_file(os.path.join(pathBase, name))
        except Exception as err:
            app.alert('Failed to read %s:\n%s' % (name, err), 'Load Error')
            return None
        files.append({'name': name, 'path': os.path.join(pathBase, name),
                      'wavelength': wl, 'delay': t, 'deltaA': A,
                      'include': True, 'maxAbs': absMaxRaw(A), 'shift': 0.0,
                      'isRef': k == 0})
    if not files:
        app.alert('No files were loaded.', 'Load Error')
        return None
    sz = files[0]['deltaA'].shape
    for f in files[1:]:
        if f['deltaA'].shape != sz:
            app.alert('Shape mismatch: %s is %dx%d, %s is %dx%d.\n'
                      'All files must share the same dimensions.' % (
                          files[0]['name'], sz[0], sz[1], f['name'],
                          f['deltaA'].shape[0], f['deltaA'].shape[1]),
                      'Load Error')
            return None
    return files


def absMaxRaw(A):
    fin = np.abs(A[np.isfinite(A)])
    return float(np.max(fin)) if fin.size else np.nan


def computeAverageOf(files):
    """Mean of the included scans on the axes of the first file; a
    shifted scan is read at t - shift. Returns (wl, t, avg, nUsed);
    avg is None when no file is included."""
    incl = [f for f in files if f['include']]
    if not incl:
        return None, None, None, 0
    wlAv = np.asarray(files[0]['wavelength'], float)
    tAv = np.asarray(files[0]['delay'], float)
    stack = np.full((wlAv.size, tAv.size, len(incl)), np.nan)
    for ii, f in enumerate(incl):
        shi = float(f['shift'])
        sameWL = np.array_equal(f['wavelength'], wlAv)
        if sameWL and shi == 0 and np.array_equal(f['delay'], tAv):
            stack[:, :, ii] = f['deltaA']
            continue
        tQ = tAv - shi
        if sameWL:
            for jj in range(wlAv.size):
                stack[jj, :, ii] = np.interp(tQ, f['delay'], f['deltaA'][jj, :],
                                             left=np.nan, right=np.nan)
        else:
            from scipy.interpolate import RegularGridInterpolator
            rgi = RegularGridInterpolator(
                (np.asarray(f['wavelength'], float),
                 np.asarray(f['delay'], float)), f['deltaA'],
                method='linear', bounds_error=False, fill_value=np.nan)
            WQ, TQ = np.meshgrid(wlAv, tQ, indexing='ij')
            stack[:, :, ii] = rgi(np.stack([WQ, TQ], axis=-1))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)
        avg = np.nanmean(stack, axis=2)
    return wlAv, tAv, avg, len(incl)


class LoadAverageWindow(ChildWindow):
    KEY = 'load'

    def __init__(self, app, pathBase, files):
        super().__init__(app, 'Load & Average Datasets', 1080, 780)
        self.pathBase = pathBase
        self.files = files
        self.filesOrig = [dict(f) for f in files]
        self.szRef = files[0]['deltaA'].shape
        self.applyBusy = False
        self.filling = False
        self.dragS = False
        self.dragK = False
        swapSuggested = detect_axis_swap(files[0]['wavelength'],
                                         files[0]['delay'])
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(10, 10, 10, 10)
        v.setSpacing(8)
        info = label(ui_str(
            'Compare files visually, then check those you want to include '
            'in the average.  Unchecked files are shown dashed / '
            'semi-transparent.  Click + drag on the preview plots to move '
            'Ref \\lambda / Ref t.  If the axes look swapped, toggle '
            '"Transpose axes".  To time-align scans, mark one file as Ref, '
            'enter Shift values manually or click "Auto-align to Ref".'),
            color=(0.15, 0.15, 0.15))
        info.setWordWrap(True)
        v.addWidget(info)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(button('Add more files...', self.onAddMore))
        row.addWidget(button('Remove selected', self.onRemoveSelected))
        row.addWidget(button('Check all', lambda: self.onSetAll(True)))
        row.addWidget(button('Uncheck all', lambda: self.onSetAll(False)))
        self.btnTrans = QtWidgets.QPushButton('Transpose axes')
        self.btnTrans.setCheckable(True)
        self.btnTrans.setChecked(bool(swapSuggested))
        self.btnTrans.setToolTip(tip(
            'Swap wavelength and delay for all files. Use this when the '
            'source files were saved with rows = delay and columns = '
            'wavelength instead of the expected layout.'))
        row.addWidget(self.btnTrans)
        row.addWidget(label('   Signal:', right=True))
        self.ddSig = QtWidgets.QComboBox()
        for k in app.SIG_KINDS:
            self.ddSig.addItem(app.sigShort(k), k)
        self.ddSig.setCurrentIndex(max(self.ddSig.findData(app.sigLoadKind), 0))
        self.ddSig.setFixedWidth(120)
        self.ddSig.setToolTip(tip(
            'What the numbers in the files are: ΔA or ΔT/T (transient '
            'absorption, TR-IR) or an intensity (TR-PL, Raman). Remembered '
            'for the next files.'))
        self.ddSig.currentIndexChanged.connect(lambda *_: self.onSigKind())
        row.addWidget(self.ddSig)
        row.addStretch(1)
        v.addLayout(row)

        self.tbl = QtWidgets.QTableWidget(0, 5)
        self.tbl.setHorizontalHeaderLabels(
            ['File', 'Max |%s|' % app.sigShort(app.sigLoadKind), 'Include',
             'Shift (%s)' % app.tUnitTxt(), 'Ref'])
        hh = self.tbl.horizontalHeader()
        hh.setSectionResizeMode(C_FILE,
                                QtWidgets.QHeaderView.ResizeMode.Stretch)
        for c, w in ((C_MAX, 110), (C_INC, 80), (C_SHIFT, 110), (C_REF, 60)):
            self.tbl.setColumnWidth(c, w)
        self.tbl.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.tbl.setFixedHeight(170)
        self.tbl.itemChanged.connect(self.onEdit)
        v.addWidget(self.tbl)

        wl0, t0 = files[0]['wavelength'], files[0]['delay']
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('Ref %s (%s):' % (app.xSymTxt(), app.xUnitTxt()),
                            right=True))
        self.edRefWL = NumEdit(0, width=100)
        row.addWidget(self.edRefWL)
        row.addWidget(label('Ref t (%s):' % app.tUnitTxt(), right=True))
        self.edRefT = NumEdit(0, width=100)
        row.addWidget(self.edRefT)
        self._setRefDefaults(wl0, t0)
        row.addWidget(label('Kinetics scale:', right=True))
        self.ddKinScale = QtWidgets.QComboBox()
        self.ddKinScale.addItems(['Linear', 'Log'])
        self.ddKinScale.setCurrentText('Log')
        row.addWidget(self.ddKinScale)
        row.addWidget(button(
            'Auto-align to Ref', self.onAutoAlign, tooltip=ui_str(
                'Estimate the time shift of every other included file by '
                'cross-correlating the kinetic trace at "Ref \\lambda" '
                'against the file marked as Ref. Results land in the Shift '
                'column.')))
        row.addWidget(button('Reset shifts', self.onResetShifts,
                             tooltip='Set all shifts back to 0.'))
        self.cbShowLegend = QtWidgets.QCheckBox('Show legend')
        self.cbShowLegend.setToolTip(tip(
            'Toggle the legend on the preview plots. When on, the legend '
            'shows short IDs (Scan 1, Scan 2, ...); full filenames remain '
            'in the table above.'))
        row.addWidget(self.cbShowLegend)
        row.addStretch(1)
        v.addLayout(row)

        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)
        self.cvS = PlotCanvas()
        self.axS = self.cvS.figure.add_axes([0.15, 0.15, 0.8, 0.75])
        self.cvS.layoutFcn = lambda: self._lay(self.cvS, self.axS)
        row.addWidget(self.cvS, 1)
        self.cvK = PlotCanvas()
        self.axK = self.cvK.figure.add_axes([0.15, 0.15, 0.8, 0.75])
        self.cvK.layoutFcn = lambda: self._lay(self.cvK, self.axK)
        row.addWidget(self.cvK, 1)
        v.addLayout(row, 1)
        self.cvS.mpl_connect('button_press_event', self.onPrevSClick)
        self.cvS.mpl_connect('motion_notify_event', self.onPrevSMotion)
        self.cvS.mpl_connect('button_release_event', self.onPrevUp)
        self.cvK.mpl_connect('button_press_event', self.onPrevKClick)
        self.cvK.mpl_connect('motion_notify_event', self.onPrevKMotion)
        self.cvK.mpl_connect('button_release_event', self.onPrevUp)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(button('Save Avg to CSV...', self.onSaveAvg))
        row.addStretch(1)
        row.addWidget(button('Cancel', self.close))
        self.btnLoad = colorButton('Load Average of Checked', GREEN)
        self.btnLoad.clicked.connect(lambda *_: self.onLoadAvg())
        row.addWidget(self.btnLoad)
        v.addLayout(row)

        self.btnTrans.toggled.connect(self.onTransposeToggle)
        self.edRefWL.valueChanged.connect(lambda *_: self.drawPreview())
        self.edRefT.valueChanged.connect(lambda *_: self.drawPreview())
        self.ddKinScale.currentTextChanged.connect(
            lambda *_: self.drawPreview())
        self.cbShowLegend.toggled.connect(lambda *_: self.drawPreview())
        if swapSuggested:
            self.onTransposeToggle(True)
        else:
            self.refreshTable()
            self.drawPreview()

    @staticmethod
    def _lay(cv, ax):
        W, H = cv.sizePt()
        cv.place(ax, 54.0, 34.0, W - 54.0 - 16.0, H - 34.0 - 22.0)

    def _setRefDefaults(self, wl, t):
        self.edRefWL.setLimits(np.min(wl), np.max(wl))
        self.edRefWL.setValue(
            wl[max(int(np.floor(wl.size / 2 + 0.5)) - 1, 0)])
        self.edRefT.setLimits(np.min(t), np.max(t))
        pos = np.flatnonzero(t >= 1)
        i0 = pos[0] if pos.size else max(int(np.floor(t.size / 2 + 0.5)) - 1, 0)
        self.edRefT.setValue(t[i0])

    # ------------------------------------------------------------------
    # table
    # ------------------------------------------------------------------
    def refreshTable(self):
        self.filling = True
        try:
            F = QtCore.Qt.ItemFlag
            self.tbl.setRowCount(len(self.files))
            for i, f in enumerate(self.files):
                it = QtWidgets.QTableWidgetItem(f['name'])
                it.setFlags(F.ItemIsEnabled | F.ItemIsSelectable)
                self.tbl.setItem(i, C_FILE, it)
                it = QtWidgets.QTableWidgetItem('%.3e' % f['maxAbs'])
                it.setFlags(F.ItemIsEnabled | F.ItemIsSelectable)
                self.tbl.setItem(i, C_MAX, it)
                for c, key in ((C_INC, 'include'), (C_REF, 'isRef')):
                    it = QtWidgets.QTableWidgetItem('')
                    it.setFlags(F.ItemIsEnabled | F.ItemIsSelectable |
                                F.ItemIsUserCheckable)
                    it.setCheckState(QtCore.Qt.CheckState.Checked if f[key]
                                     else QtCore.Qt.CheckState.Unchecked)
                    self.tbl.setItem(i, c, it)
                it = QtWidgets.QTableWidgetItem('%.6g' % f['shift'])
                it.setFlags(F.ItemIsEnabled | F.ItemIsSelectable |
                            F.ItemIsEditable)
                it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight |
                                    QtCore.Qt.AlignmentFlag.AlignVCenter)
                self.tbl.setItem(i, C_SHIFT, it)
        finally:
            self.filling = False

    def onEdit(self, item):
        if self.filling:
            return
        row, col = item.row(), item.column()
        if row >= len(self.files):
            return
        f = self.files[row]
        checked = item.checkState() == QtCore.Qt.CheckState.Checked
        if col == C_INC:
            f['include'] = checked
            self.drawPreview()
        elif col == C_SHIFT:
            try:
                val = float(item.text())
            except ValueError:
                val = np.nan
            f['shift'] = float(val) if np.isfinite(val) else 0.0
            self.filling = True
            item.setText('%.6g' % f['shift'])
            self.filling = False
            self.drawPreview()
        elif col == C_REF:
            if checked:
                for jj, g in enumerate(self.files):
                    g['isRef'] = jj == row
            else:
                f['isRef'] = False
            self.refreshTable()

    def setInclude(self, row, val):
        """Tick or clear Include of one row (as a click would)."""
        self.tbl.item(row, C_INC).setCheckState(
            QtCore.Qt.CheckState.Checked if val
            else QtCore.Qt.CheckState.Unchecked)

    def setRef(self, row):
        self.tbl.item(row, C_REF).setCheckState(QtCore.Qt.CheckState.Checked)

    def setShift(self, row, val):
        self.tbl.item(row, C_SHIFT).setText(str(val))

    def onSetAll(self, val):
        for f in self.files:
            f['include'] = bool(val)
        self.refreshTable()
        self.drawPreview()

    def onTransposeToggle(self, doSwap):
        for f, fo in zip(self.files, self.filesOrig):
            wl, t, A = orient(fo['wavelength'], fo['delay'], fo['deltaA'],
                              bool(doSwap))
            f['wavelength'], f['delay'], f['deltaA'] = wl, t, A
            f['maxAbs'] = absMaxRaw(A)
        self._setRefDefaults(self.files[0]['wavelength'],
                             self.files[0]['delay'])
        self.refreshTable()
        self.drawPreview()

    def onAddMore(self):
        names, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, 'Select more TA data file(s)', self.pathBase,
            'Data (*.csv *.dat *.txt);;All Files (*)')
        if not names:
            return
        self.app.rememberDir(os.path.dirname(names[0]))
        self.addFiles(names)

    def addFiles(self, paths):
        swap = self.btnTrans.isChecked()
        for pth in paths:
            name = os.path.basename(pth)
            try:
                wl, t, A = parse_data_file(pth)
            except Exception as err:
                self.alert('Failed to read %s:\n%s' % (name, err),
                           'Load Error')
                continue
            if A.shape != self.szRef:
                self.alert('Shape mismatch with %s; skipping %s.' % (
                    self.files[0]['name'], name), 'Shape Error')
                continue
            base = {'name': name, 'path': pth, 'wavelength': wl, 'delay': t,
                    'deltaA': A, 'include': True, 'maxAbs': absMaxRaw(A),
                    'shift': 0.0, 'isRef': False}
            self.filesOrig.append(dict(base))
            wl2, t2, A2 = orient(wl, t, A, swap)
            cur = dict(base)
            cur.update(wavelength=wl2, delay=t2, deltaA=A2,
                       maxAbs=absMaxRaw(A2))
            self.files.append(cur)
        self.refreshTable()
        self.drawPreview()

    def onRemoveSelected(self):
        rows = sorted({ix.row() for ix in self.tbl.selectedIndexes()
                       if ix.row() < len(self.files)})
        if not rows:
            return
        for r in reversed(rows):
            del self.files[r]
            del self.filesOrig[r]
        if not self.files:
            self.close()
            return
        self.refreshTable()
        self.drawPreview()

    # ------------------------------------------------------------------
    # previews
    # ------------------------------------------------------------------
    def onPrevSClick(self, ev):
        if ev.inaxes is not self.axS or ev.button != 1 or \
                self.cvS.legendAt(ev) is not None:   # (moving the legend)
            return
        self.dragS = True
        self.applyPrevSPoint(ev.xdata)

    def onPrevKClick(self, ev):
        if ev.inaxes is not self.axK or ev.button != 1 or \
                self.cvK.legendAt(ev) is not None:
            return
        self.dragK = True
        self.applyPrevKPoint(ev.xdata)

    def onPrevSMotion(self, ev):
        if self.dragS:
            x, _ = self.axS.transData.inverted().transform((ev.x, ev.y))
            self.applyPrevSPoint(x)

    def onPrevKMotion(self, ev):
        if self.dragK:
            x, _ = self.axK.transData.inverted().transform((ev.x, ev.y))
            self.applyPrevKPoint(x)

    def onPrevUp(self, _ev=None):
        self.dragS = False
        self.dragK = False

    def applyPrevSPoint(self, x):
        if x is None or not np.isfinite(x):
            return
        wl = self.files[0]['wavelength']
        self.edRefWL.setValue(min(max(float(x), np.min(wl)), np.max(wl)))
        self.drawPreview()

    def applyPrevKPoint(self, x):
        if x is None or not np.isfinite(x):
            return
        t = self.files[0]['delay']
        self.edRefT.setValue(min(max(float(x), np.min(t)), np.max(t)))
        self.drawPreview()

    def _legendLabel(self, i, f):
        sid = 'Scan %d' % (i + 1)
        if f['isRef']:
            return '%s [REF]' % sid
        if f['shift'] != 0:
            return '%s (shift=%+.3g)' % (sid, f['shift'])
        return sid

    def onSigKind(self):
        app = self.app
        app.sigLoadKind = self.ddSig.currentData()
        app.savePrefs()
        self.tbl.horizontalHeaderItem(1).setText(
            'Max |%s|' % app.sigShort(app.sigLoadKind))
        self.drawPreview()

    def drawPreview(self):
        app = self.app
        refWL, refT = self.edRefWL.value(), self.edRefT.value()
        wlAx = self.files[0]['wavelength']
        tAx = self.files[0]['delay']
        iW = int(np.argmin(np.abs(wlAx - refWL)))
        showLeg = self.cbShowLegend.isChecked()
        isLog = self.ddKinScale.currentText() == 'Log'

        ax = self.axS
        ax.clear()
        styleAxes(ax, grid=True)
        legH, legL = [], []
        for i, f in enumerate(self.files):
            style, width, alpha = (('-', 1.3, 1.0) if f['include']
                                   else ('--', 0.9, 0.45))
            iT = int(np.argmin(np.abs(f['delay'] - (refT - f['shift']))))
            h, = ax.plot(wlAx, f['deltaA'][:, iT], style,
                         color=lineColor(i) + (alpha,), linewidth=width)
            legH.append(h)
            legL.append(self._legendLabel(i, f))
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        ax.axvline(refWL, color='k', linewidth=1.0, alpha=0.7)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel(self.app.sigNameOf(self.app.sigLoadKind))
        ax.set_xlim(float(np.min(wlAx)), float(np.max(wlAx)))
        ax.set_title('Spectrum at t = %.3g %s' % (refT, app.tUnitAx()),
                     fontsize=FONT_PT + 1)
        sciTicks(ax, 'y')
        if showLeg:
            legend(ax, legH, legL, fontsize=8)
        self.cvS.draw_idle()

        ax = self.axK
        ax.clear()
        styleAxes(ax, grid=True)
        legH, legL = [], []
        for i, f in enumerate(self.files):
            style, width, alpha = (('-', 1.3, 1.0) if f['include']
                                   else ('--', 0.9, 0.45))
            iWi = int(np.argmin(np.abs(f['wavelength'] - refWL)))
            h, = ax.plot(f['delay'] + f['shift'], f['deltaA'][iWi, :], style,
                         color=lineColor(i) + (alpha,), linewidth=width)
            legH.append(h)
            legL.append(self._legendLabel(i, f))
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        if not isLog or refT > 0:
            ax.axvline(refT, color='k', linewidth=1.0, alpha=0.7)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel(self.app.sigNameOf(self.app.sigLoadKind))
        if isLog:
            ax.set_xscale('log')
            pos = tAx > 0
            if pos.any():
                ax.set_xlim(float(np.min(tAx[pos])), float(np.max(tAx)))
        else:
            ax.set_xlim(float(np.min(tAx)), float(np.max(tAx)))
        ax.set_title('Kinetics at %s = %.1f %s' % (
            app.xSym(), wlAx[iW], app.xUnitAx()), fontsize=FONT_PT + 1)
        sciTicks(ax, 'y')
        if showLeg:
            legend(ax, legH, legL, fontsize=8)
        self.cvK.draw_idle()

    # ------------------------------------------------------------------
    # alignment, average
    # ------------------------------------------------------------------
    def onAutoAlign(self):
        files = self.files
        refIdx = next((i for i, f in enumerate(files) if f['isRef']), None)
        if refIdx is None:
            refIdx = 0
            files[0]['isRef'] = True
            self.refreshTable()
        if not files[refIdx]['include']:
            self.alert('The Ref file is unchecked. Check it before running '
                       'auto-align.', 'Auto-align')
            return
        use = [i for i, f in enumerate(files) if f['include']]
        sets = [(files[i]['wavelength'], files[i]['delay'],
                 files[i]['deltaA']) for i in use]
        sh = estimate_scan_shifts(sets, ref=use.index(refIdx),
                                  ref_wl=self.edRefWL.value())
        changed = 0
        for k, i in enumerate(use):
            if i == refIdx:
                continue
            files[i]['shift'] = float(sh[k])
            changed += 1
        self.refreshTable()
        self.drawPreview()
        self.inform('Auto-aligned %d file(s) to "%s".' % (
            changed, files[refIdx]['name']), 'Auto-align')

    def onResetShifts(self):
        for f in self.files:
            f['shift'] = 0.0
        self.refreshTable()
        self.drawPreview()

    def computeAverage(self):
        return computeAverageOf(self.files)

    def onSaveAvg(self):
        wlAv, tAv, avgD, N = self.computeAverage()
        if avgD is None:
            self.alert('Check at least one file.', 'Nothing to save')
            return
        fn = self.app.putFile('CSV (*.csv)', 'Save averaged data',
                              'averaged.csv', parent=self)
        if not fn:
            return
        try:
            write_data_file(fn, wlAv, tAv, avgD, ',')
            self.inform('Saved average of %d files to:\n%s' % (
                N, os.path.basename(fn)), 'Saved')
        except Exception as err:
            self.alert('Save failed:\n%s' % err, 'Save Error')

    def onLoadAvg(self):
        if self.applyBusy:
            return
        wlAv, tAv, avgD, N = self.computeAverage()
        if avgD is None:
            self.alert('Check at least one file.', 'None selected')
            return
        names = [f['name'] for f in self.files if f['include']]
        desc = ('Loaded: %s' % names[0] if N == 1
                else 'Avg of %d files (%s...)' % (N, names[0]))
        self.applyBusy = True
        try:
            self.app.setLoadedData(wlAv, tAv, avgD, desc)
        finally:
            self.applyBusy = False
        self.close()
