"""'Anisotropy' dialog. Port of openAnisotropyWindow.

The data of the main window are the ones measured with parallel
polarization. The perpendicular data set, measured on the same axes, is
loaded here; it gets every correction the parallel data got and get.
The main window then shows one of four views (parallel, perpendicular,
magic angle, anisotropy r), chosen under "View"."""
import os

import numpy as np
from PySide6 import QtCore, QtWidgets

from .state import parse_data_file
from .textfit import legend
from .widgets import (FONT_PT, NumEdit, PlotCanvas, button, colorButton,
                      label, sciTicks, styleAxes, tip)
from .win_base import GREEN

C_PAR = (0.10, 0.35, 0.75)
C_PERP = (0.85, 0.20, 0.15)
C_MAGIC = (0.45, 0.45, 0.45)
C_R = (0.20, 0.50, 0.30)
C_FAR = (0.72, 0.72, 0.72)       # a wavelength whose G the noise leaves open
T_GSPEC = 'G of every wavelength in the tail range'


class AnisotropyDialog(QtWidgets.QDialog):
    def __init__(self, app):
        super().__init__(app.fig)
        self.app = app
        self.setWindowTitle('Anisotropy')
        from .win_base import fitSize
        self.resize(*fitSize(1240, 712))
        self.setModal(True)
        v = QtWidgets.QVBoxLayout(self)
        v.setContentsMargins(10, 10, 10, 10)
        v.setSpacing(6)

        v.addWidget(label('Parallel data (the data of the main window):   %s'
                          % app.dataSourceDesc, color=(0.2, 0.2, 0.2)))

        r = QtWidgets.QHBoxLayout()
        r.addWidget(button(
            'Load perpendicular...', self.onLoad, width=170,
            tooltip='Data file measured with perpendicular polarization at '
                    'the same wavelengths and delays. The corrections made '
                    'so far are applied to it.'))
        self.lblFile = label('', color=(0.2, 0.2, 0.2))
        r.addWidget(self.lblFile, 1)
        v.addLayout(r)

        # the range in use while tail matching is on; otherwise a
        # proposal for the delays as they are now
        tl = [float(x) for x in np.ravel(app.aniTailT)]
        if not app.aniTail or not np.all(np.isfinite(tl)):
            tl = app.aniTailProposal()
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('G (factor of the perpendicular data):', right=True))
        self.ddG = QtWidgets.QComboBox()
        self.ddG.addItem('Typed', 'typed')
        self.ddG.addItem('Tail matching', 'tail')
        self.ddG.setCurrentIndex(1 if app.aniTail else 0)
        self.ddG.setFixedWidth(150)
        self.ddG.setToolTip(tip(
            'Typed: the value entered. Tail matching: G is the least-squares '
            'factor that makes G x perpendicular equal to parallel, over all '
            'wavelengths, in the delay range given (the noise of both data '
            'sets is taken into account); use a range in which the '
            'anisotropy has decayed.'))
        r.addWidget(self.ddG)
        r.addWidget(label('G:', right=True))
        self.edG = NumEdit(app.aniG, fmt='%.6g', width=90)
        r.addWidget(self.edG)
        r.addWidget(label('tail t from:', right=True))
        self.edT1 = NumEdit(tl[0], fmt='%.6g', width=90)
        r.addWidget(self.edT1)
        r.addWidget(label('to:', right=True))
        self.edT2 = NumEdit(tl[1], fmt='%.6g', width=90)
        r.addWidget(self.edT2)
        self.lblG = label('', color=(0.25, 0.25, 0.25))
        r.addWidget(self.lblG, 1)
        v.addLayout(r)

        # a time zero that differs between the two measurements
        sr = app.aniShiftProposal()
        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('Time shift of the perpendicular data:',
                          right=True))
        self.edShift = NumEdit(app.aniShift, fmt='%.6g', width=90)
        self.edShift.setToolTip(tip(
            'For two measurements whose time zero differs. The '
            'perpendicular data are moved along the delays by this much '
            'before they are combined with the parallel data: '
            'perp(t - shift); a positive value moves them to later delays. '
            'The values between the delays come from a cubic spline where '
            'the delays are equally spaced and from linear interpolation '
            'elsewhere, as for the time shift of "Subtract data". 0: no '
            'shift.'))
        r.addWidget(self.edShift)
        r.addWidget(label(app.tUnitTxt()))
        self.btnEst = button(
            'Estimate', self.onEstimate, width=100,
            tooltip='The shift for which the perpendicular data follow the '
                    'parallel data best in the delay range given, at every '
                    'wavelength with its own factor and offset. It assumes '
                    'that the two signals have the same shape in the range: '
                    'give a range around the rise of the signal in which the '
                    'anisotropy hardly changes. The value is entered in the '
                    'field; "Apply" uses it.')
        r.addWidget(self.btnEst)
        r.addWidget(label('from the delays', right=True))
        self.edS1 = NumEdit(sr[0], fmt='%.6g', width=90)
        r.addWidget(self.edS1)
        r.addWidget(label('to:', right=True))
        self.edS2 = NumEdit(sr[1], fmt='%.6g', width=90)
        r.addWidget(self.edS2)
        self.cbZoom = QtWidgets.QCheckBox('Show these delays only')
        self.cbZoom.setToolTip(tip(
            'The two plots along the delays show this range only: the rise '
            'of the signal, where a time shift is seen.'))
        r.addWidget(self.cbZoom)
        r.addStretch(1)
        v.addLayout(r)
        self.shiftEst = None            # (the last estimate, for the tests)

        r = QtWidgets.QHBoxLayout()
        r.setSpacing(6)
        r.addWidget(label('Leave r empty where |par + 2 G perp| is below',
                          right=True))
        self.edThr = NumEdit(app.aniThr, fmt='%.6g', limits=(0, np.inf),
                             width=80)
        self.edThr.setToolTip(tip('Where the sum of the two signals is this '
                                  'small, the ratio r is noise only. 0 keeps '
                                  'every value.'))
        r.addWidget(self.edThr)
        r.addWidget(label('x its noise.'))
        r.addWidget(label('      Kinetics at λ:', right=True))
        self.edPW = NumEdit(float(app.selWL), fmt='%.6g', width=100)
        r.addWidget(self.edPW)
        r.addStretch(1)
        v.addLayout(r)

        self.cv = PlotCanvas()
        v.addWidget(self.cv, 1)
        self.axK = self.cv.figure.add_axes([0.05, 0.12, 0.27, 0.78])
        self.axR = self.cv.figure.add_axes([0.38, 0.12, 0.27, 0.78])
        # G of the tail matching along the spectrum
        self.axG = self.cv.figure.add_axes([0.71, 0.12, 0.27, 0.78])
        self.tailSpec = None            # (what the plot shows, for the tests)
        self.cv.layoutFcn = self._layout

        self.lblMsg = label('', color=(0.25, 0.25, 0.25))
        self.lblMsg.setWordWrap(True)
        self.lblMsg.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        v.addWidget(self.lblMsg)
        r = QtWidgets.QHBoxLayout()
        r.addStretch(1)
        self.btnApply = colorButton('Apply', GREEN, width=110)
        self.btnApply.clicked.connect(lambda *_: self.onApply())
        r.addWidget(self.btnApply)
        self.btnOff = button('Remove perpendicular data', self.onRemove,
                             width=210)
        r.addWidget(self.btnOff)
        r.addWidget(button('Close', self.close, width=90))
        v.addLayout(r)

        self.ddG.currentIndexChanged.connect(lambda *_: self.drawPreview())
        self.cbZoom.toggled.connect(lambda *_: self.drawPreview())
        for ed in (self.edG, self.edT1, self.edT2, self.edThr, self.edPW,
                   self.edShift, self.edS1, self.edS2):
            ed.valueChanged.connect(lambda *_: self.drawPreview())
        self.refresh('')

    # ------------------------------------------------------------------
    def _layout(self):
        W, H = self.cv.sizePt()
        # (top: the title of the third plot has two lines)
        gap, left, right, bottom, top = 58.0, 58.0, 12.0, 36.0, 36.0
        w = (W - left - right - 2 * gap) / 3.0
        h = H - bottom - top
        for k, ax in enumerate((self.axK, self.axR, self.axG)):
            self.cv.place(ax, left + k * (w + gap), bottom, w, h)

    def params(self):
        """(G, tailOn, tailT, thr) as entered."""
        return (self.edG.value(), self.ddG.currentData() == 'tail',
                [self.edT1.value(), self.edT2.value()], self.edThr.value())

    def shift(self):
        """The time shift of the perpendicular data as entered."""
        return self.edShift.value()

    def refresh(self, msg=''):
        app = self.app
        if app.perpOrig is None:
            self.lblFile.setText('No perpendicular data loaded.')
        else:
            self.lblFile.setText('%s   (%d x %d)' % (
                app.perpDesc, app.perpOrig.shape[0], app.perpOrig.shape[1]))
        self.btnApply.setEnabled(bool(app.aniOn))
        self.btnOff.setEnabled(bool(app.aniOn))
        self.btnEst.setEnabled(bool(app.aniOn))
        if not msg:
            if app.aniOn:
                msg = ('Both polarizations get the same corrections. The '
                       'main window shows the view chosen under "View" '
                       '(now: %s); that view is what the analysis windows '
                       'work on. Their setup and results are kept for each '
                       'view.' % app.ANI_VIEW_NAMES[app.aniView])
            else:
                msg = ('Load the file measured with perpendicular '
                       'polarization. The anisotropy r = (par - G perp) / '
                       '(par + 2 G perp) is then shown in the main window.')
        self.lblMsg.setText(msg)
        self.drawPreview()

    def drawPreview(self):
        app = self.app
        axK, axR, axG = self.axK, self.axR, self.axG
        for ax in (axK, axR, axG):
            ax.clear()
            styleAxes(ax, grid=True)
        self.tailSpec = None
        typed = self.ddG.currentData() == 'typed'
        self.edG.setEnabled(typed)
        self.edT1.setEnabled(not typed)
        self.edT2.setEnabled(not typed)
        self.lblG.setText('')
        axK.set_xlabel('Delay time (%s)' % app.tUnitTxt())
        axR.set_xlabel('Delay time (%s)' % app.tUnitTxt())
        axK.set_ylabel(app.sigNameOf(app.sigKind))
        axR.set_ylabel('Anisotropy r')
        axG.set_xlabel(app.xAxLabel())
        axG.set_ylabel('G')
        axG.set_title(T_GSPEC, fontsize=FONT_PT + 1)
        if not app.aniOn or app.aniPar is None:
            self.cv.draw_idle()
            return
        G, tailOn, tailT, thr = self.params()
        tt = app.delay
        iw = int(np.argmin(np.abs(app.wavelength - self.edPW.value())))
        # (the perpendicular data with the time shift entered, also
        # before "Apply")
        sh = self.shift()
        perp = app.aniPerpAt(sh if np.isfinite(sh) else app.aniShift)
        _, Gu, _ = app.aniCombine(app.aniPar, perp, tt, 'par', G,
                                  tailOn, tailT, thr)
        rr, _, nB = app.aniCombine(app.aniPar, perp, tt, 'r', G,
                                   tailOn, tailT, thr)
        self.lblG.setText('G = %.5g;  %d values of r left empty' % (Gu, nB))
        hk = [axK.plot(tt, app.aniPar[iw, :], '-', color=C_PAR, linewidth=1.2,
                       label='Parallel')[0],
              axK.plot(tt, Gu * perp[iw, :], '-', color=C_PERP,
                       linewidth=1.2, label='G x perpendicular')[0],
              axK.plot(tt, (app.aniPar[iw, :] + 2 * Gu * perp[iw, :]) / 3,
                       '-', color=C_MAGIC, linewidth=1.0,
                       label='Magic angle')[0]]
        axK.axhline(0, color='k', linestyle=':', linewidth=0.8)
        axR.plot(tt, rr[iw, :], '-', color=C_R, linewidth=1.2)
        axR.axhline(0, color='k', linestyle=':', linewidth=0.8)
        axK.set_title('Kinetics at %.6g %s' % (app.wavelength[iw],
                                               app.xUnitAx()),
                      fontsize=FONT_PT + 1)
        axR.set_title('r at %.6g %s' % (app.wavelength[iw], app.xUnitAx()),
                      fontsize=FONT_PT + 1)
        # the delays of the estimate only: the limits follow what is
        # seen there
        zr = sorted((self.edS1.value(), self.edS2.value()))
        if (self.cbZoom.isChecked() and np.all(np.isfinite(zr))
                and zr[0] < zr[1]):
            mz = (tt >= zr[0]) & (tt <= zr[1])
            if np.count_nonzero(mz) >= 2:
                for ax in (axK, axR):
                    ys = np.concatenate([np.asarray(ln.get_ydata(), float)[mz]
                                         for ln in ax.get_lines()
                                         if np.size(ln.get_xdata()) == tt.size])
                    ys = ys[np.isfinite(ys)]
                    ax.set_xlim(zr[0], zr[1])
                    if ys.size:
                        lo, hi = float(np.min(ys)), float(np.max(ys))
                        pad = 0.06 * max(hi - lo, 1e-12)
                        ax.set_ylim(min(lo, 0.0) - pad, max(hi, 0.0) + pad)
        legend(axK, handles=hk, fontsize=FONT_PT)
        sciTicks(axK, 'y')

        # G of the tail range for every single wavelength, beside the one
        # G in use: a check that one G holds for the whole spectrum
        # (nothing here enters a calculation)
        T = app.aniTailSpectrum(app.aniPar, perp, tt, tailT, Gu)
        self.tailSpec = T
        wlv = np.asarray(app.wavelength, float).ravel()
        hg = []
        far = np.isfinite(T.G) & ~T.ok
        if far.any():
            hg.append(axG.plot(wlv[far], T.G[far], '.', color=C_FAR,
                               markersize=4, label='noise too large')[0])
        if T.ok.any():
            hg.append(axG.errorbar(
                wlv[T.ok], T.G[T.ok], yerr=T.sG[T.ok], fmt='.', color=C_PAR,
                markersize=5, elinewidth=0.8, capsize=0,
                label='one wavelength, +/- its uncertainty'))
        if np.isfinite(Gu):
            hg.append(axG.axhline(Gu, color=C_PERP, linestyle='--',
                                  linewidth=1.2, label='G in use: %.5g' % Gu))
        if wlv.size > 1:
            axG.set_xlim(float(np.min(wlv)), float(np.max(wlv)))
        if T.nOk >= 3:
            yLo = min(float(np.min(T.G[T.ok] - T.sG[T.ok])), Gu)
            yHi = max(float(np.max(T.G[T.ok] + T.sG[T.ok])), Gu)
            pad = 0.08 * max(yHi - yLo, 1e-6 * max(abs(Gu), 1.0))
            axG.set_ylim(yLo - pad, yHi + pad)
        elif np.isfinite(Gu) and Gu != 0:
            axG.set_ylim(*sorted((0.5 * Gu, 1.5 * Gu)))
        if np.isfinite(T.ratio):
            # (two digits below 10, whole numbers above)
            rTxt = ('%.2g' if T.ratio < 9.95 else '%.0f') % T.ratio
            axG.set_title('%s\nscatter: %s x what the noise explains'
                          % (T_GSPEC, rTxt), fontsize=FONT_PT + 1)
        if hg:
            legend(axG, handles=hg, fontsize=FONT_PT)
        self.cv.draw_idle()

    # ------------------------------------------------------------------
    def onLoad(self):
        app = self.app
        fn = app.getFile('Data (*.csv *.dat *.txt);;All Files (*)',
                         'Select the data file of the perpendicular '
                         'polarization', parent=self)
        if not fn:
            return
        try:
            w2, t2, A2 = parse_data_file(fn)
        except Exception as e:          # noqa: BLE001 (shown to the user)
            self.refresh('Failed to read the file: %s' % e)
            return
        ok, msg = app.setPerpendicular(w2, t2, A2, os.path.basename(fn))
        if not ok:
            self.refresh(msg)
            return
        app.aniOn = True
        G, tailOn, tailT, thr = self.params()
        try:
            app.applyAnisotropy(G, tailOn, tailT, thr, self.shift())
        except Exception as e:          # noqa: BLE001
            self.refresh('Failed: %s' % e)
            return
        self.refresh('')

    def onApply(self):
        app = self.app
        if not app.aniOn:
            return
        G, tailOn, tailT, thr = self.params()
        if not (np.isfinite(G) and np.isfinite(thr)) or G == 0:
            self.refresh('G and the threshold must be numbers; G must not '
                         'be 0.')
            return
        if tailOn and not (np.all(np.isfinite(tailT))
                           and tailT[0] != tailT[1]):
            self.refresh('Give the delay range of the tail matching.')
            return
        sh = self.shift()
        if not np.isfinite(sh):
            self.refresh('The time shift of the perpendicular data must be '
                         'a number (0: no shift).')
            return
        app.applyAnisotropy(G, tailOn, tailT, thr, sh)
        self.refresh('Applied: G = %.5g%s%s.' % (
            app.aniGUsed, ' (tail matching)' if tailOn else '',
            '' if app.aniShift == 0 else
            ', time shift of the perpendicular data %.5g %s'
            % (app.aniShift, app.tUnitTxt())))

    def onEstimate(self):
        """Estimate the time shift of the perpendicular data in the
        delay range entered; the value goes to the field (not applied)."""
        app = self.app
        if not app.aniOn or app.aniPar is None or app.aniPerp0 is None:
            return
        R = app.aniShiftEstimate(app.aniPar, app.aniPerp0, app.delay,
                                 [self.edS1.value(), self.edS2.value()])
        self.shiftEst = R
        msg, d = app.aniShiftText(R)
        if np.isfinite(d):
            self.edShift.setValue(d)
        self.refresh(msg)

    def onRemove(self):
        self.app.removeAnisotropy()
        self.refresh('The perpendicular data were removed.')
