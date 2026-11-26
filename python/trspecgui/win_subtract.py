"""'Subtract Data' dialog. Port of openSubtractWindow.

A reference data set measured on the same axes (solvent, substrate) is
subtracted from the data as loaded, before every other correction:
    data - scale * reference(t - shift)
Scale and time shift are typed in or estimated by least squares in a
chosen range."""
import os

import numpy as np
from matplotlib.ticker import MaxNLocator
from PySide6 import QtCore, QtWidgets

from .state import parse_data_file
from .textfit import legend
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label, sciTicks, styleAxes, tip)
from .win_base import GREEN

C_DATA = (0.10, 0.35, 0.75)
C_REF = (0.55, 0.55, 0.55)
C_DIFF = (0.85, 0.20, 0.15)


class SubtractDialog(QtWidgets.QDialog):
    def __init__(self, app):
        super().__init__(app.fig)
        self.app = app
        self.setWindowTitle('Subtract Data')
        from .win_base import fitSize
        self.resize(*fitSize(960, 680))
        self.setModal(True)
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(10, 10, 10, 10)
        v.setSpacing(6)

        # ---- row 1: the file ----
        r = QtWidgets.QHBoxLayout()
        self.btnLoad = button(
            'Load reference...', self.onLoad, width=150,
            tooltip='Data file of the solvent or substrate alone, measured '
                    'at the same wavelengths and delays.')
        r.addWidget(self.btnLoad)
        self.lblFile = label('', color=(0.2, 0.2, 0.2))
        r.addWidget(self.lblFile, 1)
        v.addLayout(r)

        # ---- row 2: factor and shift ----
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('Scale:', right=True))
        self.edScale = NumEdit(app.subScale, fmt='%.6g', width=100)
        self.edScale.setToolTip(tip('Factor of the reference. Below 1 when '
                                    'the sample absorbs part of the pump.'))
        r.addWidget(self.edScale)
        r.addWidget(label('Time shift (%s):' % app.tUnitTxt(), right=True))
        self.edShift = NumEdit(app.subShift, fmt='%.6g', width=100)
        self.edShift.setToolTip(tip('Shift of the reference in time. '
                                    'Positive: the reference is moved to '
                                    'later delays.'))
        r.addWidget(self.edShift)
        self.btnEstS = button(
            'Estimate scale', lambda *_: self.onEstimate(False), width=130,
            tooltip='Least-squares factor in the range given below, with '
                    'the time shift as entered.')
        r.addWidget(self.btnEstS)
        self.btnEstB = button(
            'Estimate scale + shift', lambda *_: self.onEstimate(True),
            width=170,
            tooltip='Least-squares factor and time shift in the range given '
                    'below.')
        r.addWidget(self.btnEstB)
        r.addStretch(1)
        v.addLayout(r)

        # ---- row 3: range used by the estimate ----
        wl, t = app.wavelength, app.delay
        dtm = float(np.median(np.diff(t))) if t.size > 1 else 1.0
        if not np.isfinite(dtm) or dtm <= 0:
            dtm = 1.0
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('Estimate in  λ from:', right=True))
        self.edW1 = NumEdit(float(np.min(wl)), fmt='%.6g', width=80)
        r.addWidget(self.edW1)
        r.addWidget(label('to:', right=True))
        self.edW2 = NumEdit(float(np.max(wl)), fmt='%.6g', width=80)
        r.addWidget(self.edW2)
        r.addWidget(label('t from:', right=True))
        self.edT1 = NumEdit(float(t[0]), fmt='%.6g', width=80)
        r.addWidget(self.edT1)
        r.addWidget(label('to:', right=True))
        self.edT2 = NumEdit(float(t[-1]), fmt='%.6g', width=80)
        r.addWidget(self.edT2)
        r.addWidget(label('Largest shift:', right=True))
        self.edMax = NumEdit(10 * dtm, fmt='%.6g', limits=(0, np.inf),
                             width=80)
        self.edMax.setToolTip(tip('"Estimate scale + shift" looks for the '
                                  'shift between minus and plus this '
                                  'value.'))
        r.addWidget(self.edMax)
        r.addStretch(1)
        v.addLayout(r)

        # ---- row 4: where the two plots look ----
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('Kinetics at λ:', right=True))
        self.edPW = NumEdit(float(app.selWL), fmt='%.6g', width=100)
        r.addWidget(self.edPW)
        r.addWidget(label('Spectrum at t:', right=True))
        self.edPT = NumEdit(float(app.selT), fmt='%.6g', width=100)
        r.addWidget(self.edPT)
        r.addWidget(label('(the kinetics are drawn over the range '
                          '"t from ... to" above)', italic=True,
                          color=(0.45, 0.45, 0.45)))
        r.addStretch(1)
        v.addLayout(r)

        # ---- row 5: plots ----
        self.cv = PlotCanvas()
        v.addWidget(self.cv, 1)
        self.axK = self.cv.figure.add_axes([0.07, 0.12, 0.40, 0.78])
        self.axS = self.cv.figure.add_axes([0.57, 0.12, 0.40, 0.78])
        self.cv.layoutFcn = self._layout

        # ---- row 6: message; row 7: buttons ----
        self.lblMsg = label('', color=(0.25, 0.25, 0.25))
        self.lblMsg.setWordWrap(True)
        self.lblMsg.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        v.addWidget(self.lblMsg)
        r = QtWidgets.QHBoxLayout()
        r.addStretch(1)
        self.btnApply = colorButton('Apply', GREEN, width=110)
        self.btnApply.clicked.connect(lambda *_: self.onApply())
        r.addWidget(self.btnApply)
        self.btnOff = button('Remove subtraction', self.onRemove, width=160)
        r.addWidget(self.btnOff)
        r.addWidget(button('Close', self.close, width=90))
        v.addLayout(r)

        for ed in (self.edScale, self.edShift, self.edW1, self.edW2,
                   self.edT1, self.edT2, self.edPW, self.edPT):
            ed.valueChanged.connect(lambda *_: self.drawPreview())
        self.refresh('')

    # ------------------------------------------------------------------
    def _layout(self):
        W, H = self.cv.sizePt()
        gap, left, right, bottom, top = 62.0, 58.0, 12.0, 36.0, 24.0
        w = (W - left - right - gap) / 2.0
        h = H - bottom - top
        self.cv.place(self.axK, left, bottom, w, h)
        self.cv.place(self.axS, left + w + gap, bottom, w, h)

    def hasRef(self):
        app = self.app
        return (app.refStage is not None and app.deltaA_raw is not None
                and app.refStage.shape == app.deltaA_raw.shape)

    def refresh(self, msg=''):
        app = self.app
        canUse, why = app.subAllowed()
        hasRef = self.hasRef()
        if app.refOrig is None:
            self.lblFile.setText('No reference loaded.')
        else:
            self.lblFile.setText('%s   (%d x %d)' % (
                app.refDesc, app.refOrig.shape[0], app.refOrig.shape[1]))
        self.btnLoad.setEnabled(canUse)
        for b in (self.btnEstS, self.btnEstB, self.btnApply):
            b.setEnabled(canUse and hasRef)
        self.btnOff.setEnabled(bool(app.subApplied and not app.solvData))
        if not canUse:
            msg = why
        elif not msg:
            if app.subApplied:
                msg = ('Subtracted now: data - %.6g x reference, time shift '
                       '%.6g %s.' % (app.subScale, app.subShift,
                                     app.tUnitTxt()))
            elif hasRef:
                msg = ('Nothing is subtracted yet. Set the scale and the '
                       'time shift, or estimate them, then press Apply.')
            elif app.baseRebased:
                msg = ('Load the file of the solvent or substrate alone, as '
                       'measured. The corrections made so far (background, '
                       'chirp, masks, crop) are applied to it before it is '
                       'subtracted.')
            else:
                msg = ('Load the file of the solvent or substrate alone. The '
                       'subtraction is made before background and chirp '
                       'correction.')
        self.lblMsg.setText(msg)
        self.drawPreview()

    def previewIdx(self):
        app = self.app
        iw = int(np.argmin(np.abs(app.wavelength - self.edPW.value())))
        it = int(np.argmin(np.abs(app.delay - self.edPT.value())))
        return iw, it

    def drawPreview(self):
        app = self.app
        axK, axS = self.axK, self.axS
        for ax in (axK, axS):
            ax.clear()
            styleAxes(ax, grid=True)
        iw, it = self.previewIdx()
        tt, ww = app.delay, app.wavelength
        D = app.deltaA_raw
        hk = [axK.plot(tt, D[iw, :], '-', color=C_DATA, linewidth=1.2,
                       label='Data')[0]]
        hs = [axS.plot(ww, D[:, it], '-', color=C_DATA, linewidth=1.2,
                       label='Data')[0]]
        if self.hasRef() and app.subAllowed()[0]:
            Rs = self.edScale.value() * app.shiftRows(app.refStage, tt,
                                                      self.edShift.value())
            hk.append(axK.plot(tt, Rs[iw, :], '-', color=C_REF, linewidth=1.0,
                               label='Reference x scale')[0])
            hk.append(axK.plot(tt, D[iw, :] - Rs[iw, :], '-', color=C_DIFF,
                               linewidth=1.2, label='Difference')[0])
            hs.append(axS.plot(ww, Rs[:, it], '-', color=C_REF, linewidth=1.0,
                               label='Reference x scale')[0])
            hs.append(axS.plot(ww, D[:, it] - Rs[:, it], '-', color=C_DIFF,
                               linewidth=1.2, label='Difference')[0])
        for ax in (axK, axS):
            ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        lo = min(self.edT1.value(), self.edT2.value())
        hi = max(self.edT1.value(), self.edT2.value())
        if hi > lo:
            axK.set_xlim(lo, hi)
        # (a short range, -0.1 ... 0.1 ps: with the ticks matplotlib
        # chooses, -0.100 -0.075 -0.050 ..., the numbers ran into each
        # other at the size the window opens with)
        axK.xaxis.set_major_locator(
            MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))
        axK.set_xlabel('Delay time (%s)' % app.tUnitTxt())
        axK.set_ylabel(app.sigName())
        axS.set_xlabel(app.xAxLabel())
        axS.set_ylabel(app.sigName())
        axK.set_title('Kinetics at %.6g %s' % (ww[iw], app.xUnitAx()),
                      fontsize=FONT_PT + 1)
        axS.set_title('Spectrum at %.6g %s' % (tt[it], app.tUnitTxt()),
                      fontsize=FONT_PT + 1)
        legend(axK, handles=hk, fontsize=FONT_PT)
        legend(axS, handles=hs, fontsize=FONT_PT)
        sciTicks(axK, 'y')
        sciTicks(axS, 'y')
        self.cv.draw_idle()

    # ------------------------------------------------------------------
    def onLoad(self):
        app = self.app
        fn = app.getFile('Data (*.csv *.dat *.txt);;All Files (*)',
                         'Select the reference data file', parent=self)
        if not fn:
            return
        try:
            w2, t2, A2 = parse_data_file(fn)
        except Exception as e:          # noqa: BLE001 (shown to the user)
            self.refresh('Failed to read the file: %s' % e)
            return
        wasOn = app.subApplied
        ok, msg = app.setReference(w2, t2, A2, os.path.basename(fn))
        if not ok:
            self.refresh(msg)
            return
        if wasOn:
            # another reference: nothing is subtracted until Apply
            app.applySubtraction(False)
        self.refresh('')

    def onEstimate(self, withShift):
        app = self.app
        wMask = ((app.wavelength >= min(self.edW1.value(), self.edW2.value()))
                 & (app.wavelength <= max(self.edW1.value(),
                                          self.edW2.value())))
        tMask = ((app.delay >= min(self.edT1.value(), self.edT2.value()))
                 & (app.delay <= max(self.edT1.value(), self.edT2.value())))
        if np.count_nonzero(wMask) < 1 or np.count_nonzero(tMask) < 2:
            self.refresh('The range of the estimate holds too few points.')
            return
        S = app.deltaA_raw[wMask, :]
        R = app.refStage[wMask, :]
        if withShift:
            sc, sh, c0, c1 = app.fitSubtraction(S, R, app.delay, tMask, True,
                                                self.edMax.value())
            self.edShift.setValue(sh)
        else:
            sh = self.edShift.value()
            c1, sc = app.subCost(S, R, app.delay, tMask, sh)
            a = S[:, tMask]
            c0 = float(np.sum(a[np.isfinite(a)] ** 2))
        self.edScale.setValue(sc)
        nPts = max(int(np.count_nonzero(np.isfinite(S[:, tMask]))), 1)
        self.refresh('Estimate: scale %.6g, time shift %.6g %s. RMS in the '
                     'range %.4g before, %.4g after the subtraction. Not '
                     'applied yet.' % (sc, sh, app.tUnitTxt(),
                                       np.sqrt(c0 / nPts), np.sqrt(c1 / nPts)))

    def onApply(self):
        app = self.app
        if app.refStage is None or not app.subAllowed()[0]:
            return
        sc, sh = self.edScale.value(), self.edShift.value()
        if not (np.isfinite(sc) and np.isfinite(sh)):
            self.refresh('Scale and time shift must be numbers.')
            return
        try:
            app.applySubtraction(True, sc, sh)
        except Exception as e:          # noqa: BLE001 (shown to the user)
            self.refresh('Subtraction failed: %s' % e)
            return
        self.refresh('')

    def onRemove(self):
        self.app.applySubtraction(False)
        self.refresh('The subtraction was removed.')
