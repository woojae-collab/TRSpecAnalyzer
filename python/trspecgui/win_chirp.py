"""Chirp Correction window. Port of openChirpWindow / computeRidge."""
import numpy as np
from matplotlib.patheffects import withStroke
from PySide6 import QtCore, QtWidgets

from trspec.chirp import chirp_model, fit_chirp_params
from trspec.dataset import compute_ridge

from .state import ui_str
from .widgets import (FONT_PT, BlitManager, NumEdit, PlotCanvas, button,
                      colorButton, label, panel, sciTicks, styleAxes, tip)
from .win_base import (BLUE, GREEN, ChildWindow, absMax, mapLayout,
                       mapOn)


class ChirpWindow(ChildWindow):
    KEY = 'chirp'

    def __init__(self, app):
        super().__init__(app, 'Chirp Correction', 1180, 780)
        self.applyBusy = False
        self.candWL = np.nan
        self.candT = np.nan
        self.lastClickWL = None
        self.lastClickT = None
        self.dragging = False
        self.hPts = None
        self.hFitOnMap = None
        self.hRidgeOnMap = None
        self.hCandV = None
        self.hCandH = None
        u = app.tUnitTxt()

        gl = QtWidgets.QGridLayout(self)
        gl.setContentsMargins(10, 10, 10, 10)
        gl.setHorizontalSpacing(8)
        gl.setVerticalSpacing(8)
        gl.setColumnStretch(0, 3)
        gl.setColumnStretch(1, 2)
        gl.setRowStretch(1, 1)
        hdr = label(ui_str(
            'Click + drag on the 2D map to position the (\\lambda, t_0) '
            'crosshair, then press "Save point" to commit it (>= 4 points '
            'needed for fit).  Model:  t_0 = a \\cdot sqrt( (b w^2 - 1)/'
            '(c w^2 - 1) ) + d.'), color=(0.15, 0.15, 0.15))
        hdr.setWordWrap(True)
        gl.addWidget(hdr, 0, 0, 1, 2)

        self.cvMap = PlotCanvas()
        gl.addWidget(self.cvMap, 1, 0)
        self.axMap = self.cvMap.figure.add_axes([0.1, 0.1, 0.75, 0.8])
        self.caxMap = self.cvMap.figure.add_axes([0.9, 0.1, 0.02, 0.8])
        # (upper edge as the fit preview beside it, lower edge as the
        # kinetics plot)
        self.cvMap.layoutFcn = lambda: mapLayout(self.cvMap, self.axMap,
                                                 self.caxMap, T=36.0, B=34.0)
        self.blit = BlitManager(self.cvMap)
        self.cvMap.mpl_connect('button_press_event', self.addPt)
        self.cvMap.mpl_connect('motion_notify_event', self.onChirpDrag)
        self.cvMap.mpl_connect('button_release_event', self.onChirpDragUp)

        right = QtWidgets.QVBoxLayout()
        right.setSpacing(8)
        self.cvFit = PlotCanvas()
        self.axFit = self.cvFit.figure.add_axes([0.15, 0.2, 0.8, 0.65])
        self.cvFit.layoutFcn = lambda: self._lay(self.cvFit, self.axFit, 36.0)
        right.addWidget(self.cvFit, 1)
        self.tbl = QtWidgets.QTableWidget(0, 2)
        self.tbl.setHorizontalHeaderLabels(
            ['%s (%s)' % (app.xSymTxt(), app.xUnitTxt()), 't (%s)' % u])
        self.tbl.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.tbl.horizontalHeader().setStretchLastSection(True)
        self.tbl.setColumnWidth(0, 150)
        self.tbl.setMinimumHeight(90)
        right.addWidget(self.tbl, 1)
        self.cvKin = PlotCanvas()
        self.axKinP = self.cvKin.figure.add_axes([0.15, 0.2, 0.8, 0.65])
        self.cvKin.layoutFcn = lambda: self._lay(self.cvKin, self.axKinP,
                                                 22.0)
        right.addWidget(self.cvKin, 1)
        gl.addLayout(right, 1, 1)

        # ---- ridge helper -------------------------------------------
        ridgePnl = panel(ui_str('Ridge Helper (reference t_0 curve)'))
        gl.addWidget(ridgePnl, 2, 0, 1, 2)
        g = QtWidgets.QGridLayout(ridgePnl)
        g.setContentsMargins(6, 10, 6, 5)
        g.setHorizontalSpacing(5)
        g.setVerticalSpacing(3)
        self.cbShowRidge = QtWidgets.QCheckBox('Show ridge overlay')
        g.addWidget(self.cbShowRidge, 0, 0)
        g.addWidget(label('Method:', right=True), 0, 1)
        self.ddMethod = QtWidgets.QComboBox()
        self.ddMethod.addItems(['max|dA/dt|', 'max|deltaA|'])
        self.ddMethod.setCurrentText(app.ridgeMethod)
        self.ddMethod.setMinimumWidth(130)
        g.addWidget(self.ddMethod, 0, 2)
        g.addWidget(label('Search t (%s):' % u, right=True), 0, 3)
        self.edTmin = NumEdit(app.ridgeTmin, width=70)
        g.addWidget(self.edTmin, 0, 4)
        g.addWidget(label('to'), 0, 5)
        self.edTmax = NumEdit(app.ridgeTmax, width=70)
        g.addWidget(self.edTmax, 0, 6)
        g.addWidget(label(ui_str('\\lambda-smooth:'), right=True), 0, 7)
        self.edSm = NumEdit(app.ridgeSmoothN, limits=(1, 200), width=70,
                            integer=True)
        g.addWidget(self.edSm, 0, 8)
        g.addWidget(button('Recompute', self.recomputeRidge, width=100), 0, 9)
        g.addWidget(label('Auto-suggest N points:', right=True), 1, 1)
        self.edAutoN = NumEdit(8, limits=(4, 100), width=130, integer=True)
        g.addWidget(self.edAutoN, 1, 2)
        g.addWidget(button(
            'Auto-place',
            lambda: self.autoPlace(self.edAutoN.intValue())), 1, 3, 1, 2)
        g.setColumnStretch(10, 1)

        # ---- view range, points, apply ------------------------------
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(6)
        row.addWidget(label('View t_min (%s):' % u, right=True))
        self.edViewTmin = NumEdit(app.chirpViewTMin, width=70)
        row.addWidget(self.edViewTmin)
        row.addWidget(label('t_max (%s):' % u, right=True))
        self.edViewTmax = NumEdit(app.chirpViewTMax, width=70)
        row.addWidget(self.edViewTmax)
        row.addWidget(button('Reset view', self.resetView))
        b = colorButton('Save point', BLUE, tooltip=ui_str(
            'Commit the current crosshair position (\\lambda, t_0) to the '
            'points table. Use after click-and-drag to position the '
            'crosshair on the kinetic rising edge.'))
        b.clicked.connect(lambda *_: self.savePt())
        row.addWidget(b)
        row.addWidget(button('Remove Last', self.undoPt))
        row.addWidget(button('Clear Points', self.clearPts))
        row.addStretch(1)
        row.addWidget(button('Cancel', self.close))
        self.btnApply = colorButton('Fit && Apply', GREEN)
        self.btnApply.clicked.connect(lambda *_: self.applyChirp())
        row.addWidget(self.btnApply)
        gl.addLayout(row, 3, 0, 1, 2)

        self.cbShowRidge.toggled.connect(lambda *_: self.toggleRidge())
        self.ddMethod.currentTextChanged.connect(
            lambda v: self.setAndRecompute('method', v))
        self.edTmin.valueChanged.connect(
            lambda v: self.setAndRecompute('tmin', v))
        self.edTmax.valueChanged.connect(
            lambda v: self.setAndRecompute('tmax', v))
        self.edSm.valueChanged.connect(
            lambda v: self.setAndRecompute('smooth', int(round(v))))
        self.edViewTmin.valueChanged.connect(lambda *_: self.onViewChange())
        self.edViewTmax.valueChanged.connect(lambda *_: self.onViewChange())

        self._initSmallAxes()
        self.drawMapOnce()
        self.redraw()

    # ------------------------------------------------------------------
    @staticmethod
    def _lay(cv, ax, top):
        W, H = cv.sizePt()
        cv.place(ax, 54.0, 34.0, W - 54.0 - 16.0, H - 34.0 - top)

    def _initSmallAxes(self):
        app = self.app
        for ax in (self.axFit, self.axKinP):
            styleAxes(ax, grid=True)
        self.axFit.set_xlabel(app.xAxLabel())
        self.axFit.set_ylabel(self._t0Label())
        self.axFit.set_title('Polynomial fit preview', fontsize=FONT_PT + 1)
        self.axKinP.set_xlabel('Delay time (%s)' % app.tUnitAx())
        self.axKinP.set_ylabel(self.app.sigName())
        self.axKinP.set_title('Kinetics at last-clicked %s' % app.xSym(),
                              fontsize=FONT_PT + 1)

    def _t0Label(self):
        return '$t_0$ (%s)' % self.app.tUnitAx()

    def viewRange(self):
        app = self.app
        tmin, tmax = app.chirpViewTMin, app.chirpViewTMax
        if tmin is None or not np.isfinite(tmin):
            tmin = float(app.delay[0])
        if tmax is None or not np.isfinite(tmax):
            tmax = float(app.delay[-1])
        if tmax <= tmin:
            tmax = tmin + np.spacing(tmin) * 10
            if tmax <= tmin:
                tmax = tmin + 0.1
        return float(tmin), float(tmax)

    def drawMapOnce(self):
        app = self.app
        base = app.getChirpBaseData()
        wl, t = app.wavelength, app.delay
        mapOn(app, self.axMap, self.caxMap, wl, t, base, absMax(base))
        ax = self.axMap
        ax.set_title('Click + drag, then Save point  -  %s base' % (
            'BG-corrected' if app.bgApplied else 'raw'),
            fontsize=FONT_PT + 1)
        ax.set_xlim(float(np.min(wl)), float(np.max(wl)))
        ax.set_ylim(*self.viewRange())
        self.hPts = None
        self.hFitOnMap = None
        self.hRidgeOnMap = None
        self.hCandV = ax.axvline(0, color='k', linewidth=1.0, alpha=0.7,
                                 visible=False)
        self.hCandH = ax.axhline(0, color='k', linewidth=1.0, alpha=0.7,
                                 visible=False)
        self.blit.setArtists([self.hCandV, self.hCandH])
        self.cvMap.relayout()
        self.cvMap.draw_idle()

    # ------------------------------------------------------------------
    # crosshair (candidate point)
    # ------------------------------------------------------------------
    def addPt(self, ev):
        if ev.inaxes is not self.axMap or ev.button != 1:
            return
        self.dragging = True
        self.placeCandidate(ev.xdata, ev.ydata)

    def onChirpDrag(self, ev):
        if not self.dragging:
            return
        try:
            x, y = self.axMap.transData.inverted().transform((ev.x, ev.y))
        except Exception:
            return
        self.placeCandidate(x, y)

    def onChirpDragUp(self, _ev=None):
        self.dragging = False

    def placeCandidate(self, wl, t):
        app = self.app
        tR = self.viewRange()
        self.candWL = float(min(max(wl, np.min(app.wavelength)),
                                np.max(app.wavelength)))
        self.candT = float(min(max(t, tR[0]), tR[1]))
        self.lastClickWL = self.candWL
        self.lastClickT = self.candT
        self.drawCandidate()
        self.updateKinPreview()

    def drawCandidate(self):
        ok = bool(np.isfinite(self.candWL) and np.isfinite(self.candT))
        if self.hCandV is None:
            return
        self.hCandV.set_visible(ok)
        self.hCandH.set_visible(ok)
        if ok:
            self.hCandV.set_xdata([self.candWL, self.candWL])
            self.hCandH.set_ydata([self.candT, self.candT])
        self.blit.update()

    def savePt(self):
        if not (np.isfinite(self.candWL) and np.isfinite(self.candT)):
            self.alert('No candidate point yet -- click on the heatmap '
                       'first.', 'No point')
            return
        app = self.app
        app.chirpPts = np.vstack([np.asarray(app.chirpPts, float).reshape(-1, 2),
                                  [self.candWL, self.candT]])
        self.redraw()

    def undoPt(self):
        app = self.app
        if len(app.chirpPts):
            app.chirpPts = np.asarray(app.chirpPts, float)[:-1, :]
            self.redraw()

    def clearPts(self):
        self.app.chirpPts = np.zeros((0, 2))
        self.redraw()

    # ------------------------------------------------------------------
    def redraw(self):
        """Points, fitted curve and the fit preview."""
        app = self.app
        for h in (self.hPts, self.hFitOnMap):
            if h is not None:
                try:
                    h.remove()
                except Exception:
                    pass
        self.hPts = None
        self.hFitOnMap = None
        pts = np.asarray(app.chirpPts, float).reshape(-1, 2)
        self.tbl.setRowCount(pts.shape[0])
        for i in range(pts.shape[0]):
            for j in range(2):
                it = QtWidgets.QTableWidgetItem('%.5g' % pts[i, j])
                it.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignRight |
                                    QtCore.Qt.AlignmentFlag.AlignVCenter)
                self.tbl.setItem(i, j, it)
        if pts.shape[0]:
            self.hPts, = self.axMap.plot(
                pts[:, 0], pts[:, 1], linestyle='none', marker='+',
                markersize=14, color='w', markeredgewidth=1.8, zorder=5,
                path_effects=[withStroke(linewidth=3.4,
                                         foreground=(0, 0, 0, 0.75))])
        ax = self.axFit
        ax.clear()
        styleAxes(ax, grid=True)
        ax.set_xlabel(app.xAxLabel())
        ax.set_ylabel(self._t0Label())
        nPts = pts.shape[0]
        if nPts >= 4:
            try:
                QtWidgets.QApplication.setOverrideCursor(
                    QtCore.Qt.CursorShape.WaitCursor)
                try:
                    pp, rmsFit = fit_chirp_params(pts)
                finally:
                    QtWidgets.QApplication.restoreOverrideCursor()
                wlFine = np.linspace(np.min(app.wavelength),
                                     np.max(app.wavelength), 400)
                tFine = chirp_model(pp, wlFine)
                self.hFitOnMap, = self.axMap.plot(wlFine, tFine, 'k-',
                                                  linewidth=1.6)
                ax.plot(wlFine, tFine, 'b-', linewidth=1.4)
                ax.plot(pts[:, 0], pts[:, 1], 'o', color='r', markersize=6)
                ax.set_title('N=%d   RMS=%.3g %s\na=%.3g  b=%.3g  c=%.3g  '
                             'd=%.3g' % (nPts, rmsFit, app.tUnitAx(), pp[0],
                                         pp[1], pp[2], pp[3]),
                             fontsize=FONT_PT)
            except Exception as err:
                ax.set_title('Fit failed: %s' % err, fontsize=FONT_PT)
        else:
            ax.set_title('Need %d more point(s)  (>= 4 required)'
                         % max(0, 4 - nPts), fontsize=FONT_PT + 1)
        sciTicks(ax, 'y')
        self.cvFit.draw_idle()
        self.cvMap.draw_idle()

    def updateKinPreview(self):
        if self.lastClickWL is None:
            return
        app = self.app
        base = app.getChirpBaseData()
        iw = int(np.argmin(np.abs(app.wavelength - self.lastClickWL)))
        ax = self.axKinP
        ax.clear()
        styleAxes(ax, grid=True)
        ax.plot(app.delay, base[iw, :], 'k-', linewidth=1.2)
        tr = ax.get_xaxis_transform()
        ax.axvline(self.lastClickT, color='k', linewidth=1.0, alpha=0.7)
        ax.text(self.lastClickT, 0.97, 't = %.3g %s ' % (
            self.lastClickT, app.tUnitAx()), transform=tr, ha='right',
            va='top', fontsize=FONT_PT - 1, clip_on=True)
        if len(app.ridgeWL):
            jr = int(np.argmin(np.abs(np.asarray(app.ridgeWL)
                                      - app.wavelength[iw])))
            tRd = float(np.asarray(app.ridgeT)[jr])
            ax.axvline(tRd, color='g', linestyle='--', linewidth=1.2)
            ax.text(tRd, 0.85, ' $t_{ridge}$ = %.3g %s' % (
                tRd, app.tUnitAx()), transform=tr, ha='left', va='top',
                fontsize=FONT_PT - 1, color=(0, 0.45, 0), clip_on=True)
        ax.axhline(0, color='k', linestyle=':', linewidth=0.8)
        tLo, tHi = app.chirpViewTMin, app.chirpViewTMax
        if tLo is None or not np.isfinite(tLo):
            tLo = float(app.delay[0])
        if tHi is None or not np.isfinite(tHi):
            tHi = float(app.delay[-1])
        if tHi <= tLo:
            tLo, tHi = float(app.delay[0]), float(app.delay[-1])
        ax.set_xlim(tLo, tHi)
        ax.set_xlabel('Delay time (%s)' % app.tUnitAx())
        ax.set_ylabel(self.app.sigName())
        ax.set_title('%s = %.1f %s' % (app.xSym(), app.wavelength[iw],
                                       app.xUnitAx()), fontsize=FONT_PT + 1)
        sciTicks(ax, 'y')
        self.cvKin.draw_idle()

    # ------------------------------------------------------------------
    # ridge helper
    # ------------------------------------------------------------------
    def computeRidge(self):
        app = self.app
        rw, rt = compute_ridge(app.getChirpBaseData(), app.wavelength,
                               app.delay, app.ridgeMethod, app.ridgeTmin,
                               app.ridgeTmax, app.ridgeSmoothN)
        app.ridgeWL, app.ridgeT = rw, rt

    def setAndRecompute(self, field, val):
        app = self.app
        if field == 'method':
            app.ridgeMethod = val
        elif field == 'tmin':
            app.ridgeTmin = float(val)
        elif field == 'tmax':
            app.ridgeTmax = float(val)
        elif field == 'smooth':
            app.ridgeSmoothN = int(val)
        if self.cbShowRidge.isChecked():
            self.recomputeRidge()

    def _dropRidgeLine(self):
        if self.hRidgeOnMap is not None:
            try:
                self.hRidgeOnMap.remove()
            except Exception:
                pass
        self.hRidgeOnMap = None

    def recomputeRidge(self):
        app = self.app
        self.computeRidge()
        self._dropRidgeLine()
        if self.cbShowRidge.isChecked() and len(app.ridgeT):
            self.hRidgeOnMap, = self.axMap.plot(
                app.ridgeWL, app.ridgeT, '-', color=(0, 0.7, 0),
                linewidth=1.2)
        self.cvMap.draw_idle()

    def toggleRidge(self):
        if self.cbShowRidge.isChecked():
            self.recomputeRidge()
        else:
            self._dropRidgeLine()
            self.cvMap.draw_idle()

    def autoPlace(self, N):
        app = self.app
        if not len(app.ridgeT):
            self.computeRidge()
        if not len(app.ridgeT):
            self.alert('Could not compute ridge (check search window).',
                       'No ridge')
            return
        wls = np.linspace(app.wavelength[0], app.wavelength[-1], N + 2)[1:-1]
        newPts = np.zeros((N, 2))
        rw = np.asarray(app.ridgeWL, float)
        rt = np.asarray(app.ridgeT, float)
        for k in range(N):
            j = int(np.argmin(np.abs(rw - wls[k])))
            newPts[k, :] = [rw[j], rt[j]]
        app.chirpPts = newPts
        self.redraw()
        if not self.cbShowRidge.isChecked():
            self.cbShowRidge.setChecked(True)      # draws the ridge

    # ------------------------------------------------------------------
    # view range
    # ------------------------------------------------------------------
    def applyViewRange(self):
        self.axMap.set_ylim(*self.viewRange())
        self.cvMap.draw_idle()

    def onViewChange(self):
        app = self.app
        app.chirpViewTMin = self.edViewTmin.value()
        app.chirpViewTMax = self.edViewTmax.value()
        self.applyViewRange()
        if self.lastClickWL is not None:
            self.updateKinPreview()

    def resetView(self):
        app = self.app
        app.chirpViewTMin = float(app.delay[0])
        app.chirpViewTMax = float(min(app.defaultChirpTMax(), app.delay[-1]))
        self.edViewTmin.setValue(app.chirpViewTMin)
        self.edViewTmax.setValue(app.chirpViewTMax)
        self.applyViewRange()
        if self.lastClickWL is not None:
            self.updateKinPreview()

    # ------------------------------------------------------------------
    def applyChirp(self):
        if self.applyBusy:
            return
        app = self.app
        pts = np.asarray(app.chirpPts, float).reshape(-1, 2)
        if pts.shape[0] < 4:
            self.alert('At least 4 points are required.', 'Not enough points')
            return
        try:
            pp, rmsFit = fit_chirp_params(pts)
        except Exception as err:
            self.alert('Fit failed:\n%s' % err, 'Fit error')
            return
        self.applyBusy = True
        try:
            app.chirpParams = np.asarray(pp, float)
            app.chirpFitRMS = float(rmsFit)
            app.chirpApplied = True
            app.updateAll()
        finally:
            self.applyBusy = False
        self.close()
