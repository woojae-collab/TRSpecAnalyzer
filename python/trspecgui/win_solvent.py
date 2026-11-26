"""'Solvent Coherence' dialog of the Coherence window. Port of
openSolventWindow.

The vibrational coherence of the solvent, measured alone, is taken out
of the residual before the FFT:
    residual - A * solvent coherence(t - shift)
A and the time shift are typed in or estimated in the frequency domain,
on the real and imaginary parts of the Fourier transforms in a band that
holds a mode of the solvent (M. Liebel, thesis, Oxford 2014, Fig. 4.6).
The solvent file is the reference data set of "Subtract data...": it is
loaded as measured, the corrections of the data are applied to it.

"Subtract from the data" takes the same coherence out of the data of the
main window (state.solventInData), or, as "Whole solvent signal", hands
the factor and the time shift to the subtraction of "Subtract data..."."""
import os

import numpy as np
from PySide6 import QtCore, QtWidgets

from trspec.coherence import (coh_spectrum, fit_solvent_coherence,
                              solvent_band, solvent_shifted, solvent_vc)

from .state import parse_data_file, ui_str
from .textfit import legend
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label, sciTicks, styleAxes, tip)
from .win_base import GREEN

C_DATA = (0.10, 0.35, 0.75)
C_REF = (0.55, 0.55, 0.55)
C_DIFF = (0.85, 0.20, 0.15)
KCM = 1e12 / 2.99792458e10          # cm^-1 per THz


class SolventDialog(QtWidgets.QDialog):
    def __init__(self, coh):
        app = coh.app
        super().__init__(coh)
        self.coh = coh
        self.app = app
        self.setWindowTitle('Solvent Coherence')
        from .win_base import fitSize
        self.resize(*fitSize(980, 702))
        self.setModal(True)
        self.opts = coh.getOpts()
        o = self.opts
        self.tt = np.asarray(app.delay, float).ravel()
        self.wl = np.asarray(app.wavelength, float).ravel()
        self.inWin = (self.tt >= o['t_min']) & (self.tt <= o['t_max'])
        dtw = float(np.median(np.diff(self.tt[self.inWin])))
        band0 = [float(x) for x in np.ravel(app.solvBand)]
        if np.all(np.isfinite(band0)) and app.solvBandUnit != o['freq_unit']:
            if o['freq_unit'] == 'thz':
                band0 = [x / KCM for x in band0]
            else:
                band0 = [x * KCM for x in band0]
        if not np.all(np.isfinite(band0)):
            band0 = [0.0, 0.0]
        wl0 = [float(x) for x in np.ravel(app.solvWL)]
        if not np.all(np.isfinite(wl0)):
            wl0 = [float(np.min(self.wl)), float(np.max(self.wl))]
        ms0 = app.solvMaxShift
        if not np.isfinite(ms0) or ms0 < 0:
            ms0 = 2 * dtw
        from .win_coherence import freqAxLabel
        self.fuTxt = 'THz' if o['freq_unit'] == 'thz' else ui_str('cm^{-1}')
        self.fLabel = freqAxLabel(o['freq_unit'])
        # ---- state ----
        self.R = None          # residual of the sample (made when first needed)
        self.F = None          # its Fourier transform
        self.freq = None
        self.finfo = None
        self.V = None          # solvent coherence
        self.Vpad = np.nan     # shift it was made for
        self.Rmsg = ''         # why the data are not available

        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(10, 10, 10, 10)
        v.setSpacing(6)

        # ---- row 1: the file ----
        r = QtWidgets.QHBoxLayout()
        self.btnLoad = button(
            'Load solvent data...', self.onLoad, width=170,
            tooltip='Data file of the solvent alone, measured at the same '
                    'wavelengths and delays. It is the reference data set '
                    'of "Subtract data...".')
        r.addWidget(self.btnLoad)
        self.lblFile = label('', color=(0.2, 0.2, 0.2))
        r.addWidget(self.lblFile, 1)
        v.addLayout(r)

        # ---- row 2: where the fit looks ----
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label(ui_str('Fit in the band  \\nu from:'), right=True))
        self.edF1 = NumEdit(band0[0], fmt='%.6g', width=80)
        r.addWidget(self.edF1)
        r.addWidget(label('to:', right=True))
        self.edF2 = NumEdit(band0[1], fmt='%.6g', width=80)
        self.edF2.setToolTip(tip('Frequency band used by the estimate. It '
                                 'should hold a mode of the solvent and no '
                                 'mode of the solute.'))
        r.addWidget(self.edF2)
        r.addWidget(label(self.fuTxt))
        self.btnBand = button(
            'Strongest solvent band', self.onBand, width=170,
            tooltip='Set the band to the strongest mode of the solvent data '
                    '(where its power is above 10 % of the largest value).')
        r.addWidget(self.btnBand)
        r.addWidget(label('%s from:' % app.xSymTxt(), right=True))
        self.edW1 = NumEdit(wl0[0], fmt='%.6g', width=80)
        r.addWidget(self.edW1)
        r.addWidget(label('to:', right=True))
        self.edW2 = NumEdit(wl0[1], fmt='%.6g', width=80)
        self.edW2.setToolTip(tip('Spectral range used by the estimate.'))
        r.addWidget(self.edW2)
        r.addStretch(1)
        v.addLayout(r)

        # ---- row 3: factor and shift ----
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('Factor A:', right=True))
        self.edA = NumEdit(app.solvA, fmt='%.6g', width=90)
        self.edA.setToolTip(tip('Factor of the solvent coherence. Below 1 '
                                'when the solute absorbs part of the pump.'))
        r.addWidget(self.edA)
        r.addWidget(label('Time shift (%s):' % app.tUnitTxt(), right=True))
        self.edShift = NumEdit(app.solvShift, fmt='%.6g', width=90)
        self.edShift.setToolTip(tip('Shift of the solvent coherence in time. '
                                    'Positive: it is moved to later delays.'))
        r.addWidget(self.edShift)
        r.addWidget(label('Largest shift:', right=True))
        self.edMax = NumEdit(ms0, fmt='%.6g', limits=(0, np.inf), width=80)
        self.edMax.setToolTip(tip(
            '"Estimate factor + shift" looks for the shift between minus and '
            'plus this value. Keep it below half a period of the solvent '
            'mode.'))
        r.addWidget(self.edMax)
        self.btnEstA = button(
            'Estimate factor', lambda *_: self.onEstimate(False), width=130,
            tooltip='Least-squares factor in the band, with the time shift '
                    'as entered.')
        r.addWidget(self.btnEstA)
        self.btnEstB = button(
            'Estimate factor + shift', lambda *_: self.onEstimate(True),
            width=170,
            tooltip='Least-squares factor and time shift in the band.')
        r.addWidget(self.btnEstB)
        r.addStretch(1)
        v.addLayout(r)

        # ---- row 4: where the right plot looks ----
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('Traces at %s:' % app.xSymTxt(), right=True))
        lam = coh.curLambda if np.isfinite(coh.curLambda) else float(app.selWL)
        self.edPW = NumEdit(float(lam), fmt='%.6g', width=100)
        r.addWidget(self.edPW)
        r.addWidget(label('(time window %.4g ... %.4g %s and FFT settings: '
                          'those of the Coherence window)'
                          % (o['t_min'], o['t_max'], app.tUnitTxt()),
                          italic=True, color=(0.45, 0.45, 0.45)))
        r.addStretch(1)
        v.addLayout(r)

        # ---- row 5: plots ----
        self.cv = PlotCanvas()
        v.addWidget(self.cv, 1)
        self.axP = self.cv.figure.add_axes([0.07, 0.12, 0.40, 0.78])
        self.axT = self.cv.figure.add_axes([0.57, 0.12, 0.40, 0.78])
        self.cv.layoutFcn = self._layout

        # ---- row 6: message; row 7: this analysis ----
        self.lblMsg = label('', color=(0.25, 0.25, 0.25))
        self.lblMsg.setWordWrap(True)
        self.lblMsg.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        v.addWidget(self.lblMsg)
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addStretch(1)
        r.addWidget(label('From the residual of this analysis:', right=True))
        self.btnApply = colorButton('Apply', GREEN, width=110)
        self.btnApply.clicked.connect(lambda *_: self.onApply())
        r.addWidget(self.btnApply)
        self.btnOff = button('Remove subtraction', self.onRemove, width=160)
        r.addWidget(self.btnOff)
        r.addWidget(button('Close', self.close, width=90))
        v.addLayout(r)

        for ed in (self.edA, self.edShift, self.edMax, self.edF1, self.edF2,
                   self.edW1, self.edW2, self.edPW):
            ed.valueChanged.connect(lambda *_: self.onEdit())
        self.refresh('')

    # ------------------------------------------------------------------
    def _layout(self):
        W, H = self.cv.sizePt()
        gap, left, right, bottom, top = 62.0, 58.0, 12.0, 36.0, 24.0
        w = (W - left - right - gap) / 2.0
        h = H - bottom - top
        self.cv.place(self.axP, left, bottom, w, h)
        self.cv.place(self.axT, left + w + gap, bottom, w, h)

    def needData(self):
        """Residual of the sample, its Fourier transform and the solvent
        coherence: made when first needed. False when not available
        (the reason is in self.Rmsg, if there is one)."""
        app, o = self.app, self.opts
        if app.refStage is None or app.refStage.shape != app.deltaA.shape:
            return False
        if self.R is None:
            self.lblMsg.setText('Computing the residual...')
            QtWidgets.QApplication.processEvents()
            try:
                Rn = self.coh.getResidual()[0]
                F, freq, info = coh_spectrum(Rn, self.tt, **o)
            except Exception as err:        # noqa: BLE001 (shown to the user)
                self.R = None
                self.Rmsg = str(err)
                return False
            self.R, self.F, self.finfo = np.asarray(Rn, float), F, info
            self.freq = np.asarray(freq, float).ravel()
        pad = max(abs(self.edMax.value()), abs(self.edShift.value()))
        if self.V is None or pad != self.Vpad:
            try:
                self.V = solvent_vc(app.solventData(), self.tt, o['t_min'],
                                    o['t_max'], pad)
                self.Vpad = pad
            except Exception as err:        # noqa: BLE001
                self.V = None
                self.Rmsg = str(err)
                return False
        return True

    def fitRows(self):
        lo = min(self.edW1.value(), self.edW2.value())
        hi = max(self.edW1.value(), self.edW2.value())
        return (self.wl >= lo) & (self.wl <= hi)

    def setBand(self):
        """Band of the strongest solvent mode into the two fields."""
        try:
            G0 = coh_spectrum(solvent_shifted(self.V, self.tt, 0.0), self.tt,
                              **self.opts)[0]
            b = solvent_band(G0, self.freq, self.fitRows(), self.finfo.n_fft,
                             self.finfo.n_win)
        except Exception:                   # noqa: BLE001
            b = (np.nan, np.nan)
        if np.all(np.isfinite(b)) and b[1] > b[0]:
            for ed, val in ((self.edF1, b[0]), (self.edF2, b[1])):
                ed.blockSignals(True)
                ed.setValue(val)
                ed.blockSignals(False)

    def refresh(self, msg=''):
        app = self.app
        if app.refOrig is None:
            self.lblFile.setText('No solvent data loaded.')
        else:
            self.lblFile.setText('%s   (%d x %d)' % (
                app.refDesc, app.refOrig.shape[0], app.refOrig.shape[1]))
        self.Rmsg = ''
        ok = self.needData()
        if ok and not (self.edF2.value() > self.edF1.value()):
            self.setBand()             # no band yet: the strongest one
        for b in (self.btnBand, self.btnEstA, self.btnEstB, self.btnApply):
            b.setEnabled(ok)
        self.btnOff.setEnabled(bool(app.solvOn))
        self.btnLoad.setEnabled(not app.subApplied)
        if self.Rmsg:
            msg = self.Rmsg
        elif not msg:
            if app.solvOn:
                msg = ('Subtracted now: residual - %.5g x solvent coherence, '
                       'time shift %.5g %s.' % (app.solvA, app.solvShift,
                                                app.tUnitTxt()))
                if app.solvViewsOK():
                    msg += ('  In the main window, "View" shows the data '
                            'without the solvent: "Solvent coherence '
                            'subtracted" (the oscillating part, in this time '
                            'window) or "Solvent signal subtracted" (the '
                            'whole solvent signal, at all delays).')
            elif ok:
                msg = ('Nothing is subtracted yet. Choose a band that holds '
                       'a mode of the solvent and none of the solute, '
                       'estimate the factor (and the time shift), then press '
                       'Apply. The main window can then show the data '
                       'without the solvent ("View").')
            elif app.refOrig is not None:
                msg = ('The solvent data could not be brought to the state '
                       'of the data shown.')
            else:
                msg = ('Load the file of the solvent alone, as measured. The '
                       'corrections of the data (background, chirp, masks, '
                       'crop) are applied to it.')
        if app.subApplied and not self.Rmsg:
            msg += ('  This file is subtracted in the main window (scale '
                    '%.4g): the factor here is what remains. To use another '
                    'file, remove that subtraction first.' % app.subScale)
        self.lblMsg.setText(msg)
        self.drawPreview()

    def onEdit(self):
        self.needData()
        self.drawPreview()

    def preview(self):
        """(Yd, Gd): the solvent coherence with the shift entered, and
        its Fourier transform; None when there is nothing to show."""
        if self.R is None or self.V is None:
            return None
        A, d = self.edA.value(), self.edShift.value()
        if not (np.isfinite(A) and np.isfinite(d)):
            return None
        try:
            Yd = solvent_shifted(self.V, self.tt, d)
            Gd = coh_spectrum(Yd, self.tt, **self.opts)[0]
        except Exception:                   # noqa: BLE001
            return None
        return Yd, Gd

    def drawPreview(self):
        app = self.app
        axP, axT = self.axP, self.axT
        for ax in (axP, axT):
            ax.clear()
            styleAxes(ax, grid=True)
        axP.set_xlabel(self.fLabel)
        axP.set_ylabel('Power')
        axT.set_xlabel('Delay time (%s)' % app.tUnitAx())
        axT.set_ylabel('Residual')
        pv = self.preview()
        if pv is None:
            self.cv.draw_idle()
            return
        Yd, Gd = pv
        A = self.edA.value()
        rows = self.fitRows()
        if not rows.any():
            rows = np.ones(self.wl.size, bool)
        F = self.F
        P0 = np.sum(np.abs(F[rows, :]) ** 2, axis=0)
        Ps = A ** 2 * np.sum(np.abs(Gd[rows, :]) ** 2, axis=0)
        P1 = np.sum(np.abs(F[rows, :] - A * Gd[rows, :]) ** 2, axis=0)
        hp = [axP.plot(self.freq, P0, '-', color=C_DATA, linewidth=1.2,
                       label='Residual')[0],
              axP.plot(self.freq, Ps, '-', color=C_REF, linewidth=1.0,
                       label='A x solvent')[0],
              axP.plot(self.freq, P1, '-', color=C_DIFF, linewidth=1.2,
                       label='After subtraction')[0]]
        if self.edF2.value() > self.edF1.value():
            axP.axvline(self.edF1.value(), color='k', linestyle=':',
                        linewidth=0.8)
            axP.axvline(self.edF2.value(), color='k', linestyle=':',
                        linewidth=0.8)
        xl = self.coh.getXlimForFreq(self.freq)
        if len(xl) == 2 and xl[1] > xl[0]:
            axP.set_xlim(xl[0], xl[1])
        axP.set_title('Power in the spectral range of the fit',
                      fontsize=FONT_PT + 1)
        legend(axP, handles=hp, fontsize=FONT_PT)
        sciTicks(axP, 'y')
        iw = int(np.argmin(np.abs(self.wl - self.edPW.value())))
        m = self.inWin
        ht = [axT.plot(self.tt[m], self.R[iw, m], '-', color=C_DATA,
                       linewidth=1.2, label='Residual')[0],
              axT.plot(self.tt[m], A * Yd[iw, m], '-', color=C_REF,
                       linewidth=1.0, label='A x solvent coherence')[0],
              axT.plot(self.tt[m], self.R[iw, m] - A * Yd[iw, m], '-',
                       color=C_DIFF, linewidth=1.2,
                       label='After subtraction')[0]]
        axT.axhline(0, color='k', linestyle=':', linewidth=0.8)
        axT.set_title('Traces at %.6g %s' % (self.wl[iw], app.xUnitAx()),
                      fontsize=FONT_PT + 1)
        legend(axT, handles=ht, fontsize=FONT_PT)
        sciTicks(axT, 'y')
        self.cv.draw_idle()

    # ------------------------------------------------------------------
    def hasResult(self):
        return self.coh.res['Pmap'] is not None

    def onLoad(self):
        app = self.app
        fn = app.getFile('Data (*.csv *.dat *.txt);;All Files (*)',
                         'Select the data file of the solvent', parent=self)
        if not fn:
            return
        try:
            w2, t2, A2 = parse_data_file(fn)
        except Exception as e:          # noqa: BLE001 (shown to the user)
            self.refresh('Failed to read the file: %s' % e)
            return
        wasOn = app.solvOn
        ok, msg = app.setReference(w2, t2, A2, os.path.basename(fn))
        if not ok:
            self.refresh(msg)
            return
        self.V = None
        self.edF1.blockSignals(True)
        self.edF2.blockSignals(True)
        self.edF1.setValue(0.0)         # the band of the new file
        self.edF2.setValue(0.0)
        self.edF1.blockSignals(False)
        self.edF2.blockSignals(False)
        em = ''
        if wasOn and self.hasResult():
            # another file: nothing is subtracted until Apply
            em = self.coh.rerunQuiet()
        else:
            self.coh.updateSolvLabel()
        if em:
            self.refresh('The analysis failed: %s' % em)
            return
        self.refresh('')

    def onBand(self):
        if not self.needData():
            self.refresh('')
            return
        self.setBand()
        self.refresh('Band of the strongest solvent mode: %.5g ... %.5g %s.'
                     % (self.edF1.value(), self.edF2.value(), self.fuTxt))

    def onEstimate(self, withShift):
        app = self.app
        if not self.needData():
            self.refresh('')
            return
        if not (self.edF2.value() > self.edF1.value()):
            self.refresh('Give the band of the fit ("from" below "to").')
            return
        rows = self.fitRows()
        if not rows.any():
            self.refresh('The spectral range of the fit holds no point.')
            return
        try:
            q = fit_solvent_coherence(
                self.R, self.V, self.tt,
                [self.edF1.value(), self.edF2.value()], rows,
                self.edMax.value(), withShift, self.edShift.value(),
                **self.opts)
        except Exception as e:          # noqa: BLE001
            self.refresh('Estimate failed: %s' % e)
            return
        self.edA.blockSignals(True)
        self.edShift.blockSignals(True)
        self.edA.setValue(q.A)
        self.edShift.setValue(q.d)
        self.edA.blockSignals(False)
        self.edShift.blockSignals(False)
        msg = ('Estimate: A = %.5g, time shift %.5g %s. Power in the band: '
               '%.4g before, %.4g after (%.1f %% left). Not applied yet.'
               % (q.A, q.d, app.tUnitTxt(), q.P0, q.P1,
                  100 * q.P1 / max(q.P0, np.finfo(float).tiny)))
        mx = self.edMax.value()
        if q.A <= 0:
            msg += ('  The factor is not positive: the band holds nothing '
                    'that matches the solvent.')
        elif withShift and mx > 0 and abs(abs(q.d) - mx) < 1e-6 * mx:
            msg += '  The shift is at the limit "Largest shift": check it.'
        self.refresh(msg)

    def onApply(self):
        app = self.app
        if not self.needData():
            self.refresh('')
            return
        A, d = self.edA.value(), self.edShift.value()
        if not (np.isfinite(A) and np.isfinite(d)):
            self.refresh('Factor and time shift must be numbers.')
            return
        app.solvOn = True
        app.solvA = float(A)
        app.solvShift = float(d)
        app.solvBand = [self.edF1.value(), self.edF2.value()]
        app.solvBandUnit = self.opts['freq_unit']
        app.solvWL = [self.edW1.value(), self.edW2.value()]
        app.solvMaxShift = float(self.edMax.value())
        self.lblMsg.setText('Running the coherence analysis...')
        QtWidgets.QApplication.processEvents()
        em = self.coh.rerunQuiet()
        if em:
            self.refresh('The analysis failed: %s' % em)
            return
        self.refresh('')

    def onRemove(self):
        self.app.solvOn = False
        self.app.clearSolvDef()
        em = ''
        if self.hasResult():
            em = self.coh.rerunQuiet()
        else:
            self.coh.updateSolvLabel()
        if em:
            self.refresh('The analysis failed: %s' % em)
            return
        self.refresh('The subtraction was removed.')
