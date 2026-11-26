"""'Convert ΔA <-> ΔT/T' dialog. Port of openSignalWindow.

Transient-absorption data are converted by the exact relation
    ΔT/T = 10^(-ΔA) - 1        ΔA = -log10(1 + ΔT/T)        (ΔA in OD)
Like a crop, the conversion takes the data with all corrections made so
far as the new starting point."""
from PySide6 import QtCore, QtWidgets

from .widgets import button, colorButton, label, tip

NOTE = ('The conversion is applied to the data shown, with all corrections '
        'made so far, and the result becomes the new starting point (as '
        'after a crop). Analysis results are discarded. "Reset all" returns '
        'to the data as loaded.')


INT_TXT = ('Intensity data are not converted. If these are '
           'transient-absorption data, choose above what the files contain '
           '(ΔA or ΔT/T).')


class SignalDialog(QtWidgets.QDialog):
    def __init__(self, app):
        super().__init__(app.fig)
        self.app = app
        self.setWindowTitle('Convert ΔA <-> ΔT/T')
        from .win_base import fitSize
        self.resize(*fitSize(540, 360))
        self.setModal(True)
        g = QtWidgets.QGridLayout(self)
        g.setContentsMargins(14, 14, 14, 14)
        g.setVerticalSpacing(9)
        g.setHorizontalSpacing(8)
        g.setColumnMinimumWidth(0, 170)
        g.setColumnStretch(1, 1)

        g.addWidget(label('Loaded files contain:'), 0, 0)
        self.ddLoad = QtWidgets.QComboBox()
        for k in app.SIG_KINDS:
            self.ddLoad.addItem(app.sigShort(k), k)
        self.ddLoad.setToolTip(tip(
            'What the numbers in your data files are (the choice made when '
            'a file is loaded). Changing it here renames the data shown; no '
            'number changes.'))
        g.addWidget(self.ddLoad, 0, 1)

        g.addWidget(label('ΔA is given in:'), 1, 0)
        self.ddUnit = QtWidgets.QComboBox()
        self.ddUnit.addItems(['OD', 'mOD'])
        self.ddUnit.setCurrentText(app.sigAUnit)
        self.ddUnit.setToolTip(tip(
            'Unit of the ΔA values, for the conversion in both directions. '
            'ΔT/T is a fraction (0.01 = 1 %).'))
        g.addWidget(self.ddUnit, 1, 1)

        self.lblF = label('ΔT/T = 10^(-ΔA) - 1          '
                          'ΔA = -log10(1 + ΔT/T)          (ΔA in OD)',
                          color=(0.35, 0.35, 0.35))
        g.addWidget(self.lblF, 2, 0, 1, 2)
        self.lblNow = label('')
        f = self.lblNow.font()
        f.setBold(True)
        self.lblNow.setFont(f)
        g.addWidget(self.lblNow, 3, 0, 1, 2)
        self.btn = colorButton('Convert', (0.30, 0.55, 0.80))
        self.btn.clicked.connect(lambda *_: self.onConvert())
        g.addWidget(self.btn, 4, 0, 1, 2)
        self.lblMsg = label('', color=(0.25, 0.25, 0.25))
        self.lblMsg.setWordWrap(True)
        self.lblMsg.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        self.lblMsg.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop
                                 | QtCore.Qt.AlignmentFlag.AlignLeft)
        g.addWidget(self.lblMsg, 5, 0, 1, 2)
        g.setRowStretch(5, 1)
        row = QtWidgets.QHBoxLayout()
        row.addStretch(1)
        row.addWidget(button('Close', self.close, width=90))
        g.addLayout(row, 6, 0, 1, 2)

        self.refresh(NOTE)
        self.ddLoad.currentIndexChanged.connect(lambda *_: self.onLoadKind())
        self.ddUnit.currentIndexChanged.connect(lambda *_: self.onUnit())

    # ------------------------------------------------------------------
    def _setLoad(self, kind):
        self.ddLoad.blockSignals(True)
        self.ddLoad.setCurrentIndex(max(self.ddLoad.findData(kind), 0))
        self.ddLoad.blockSignals(False)

    def refresh(self, msg=None):
        app = self.app
        if app.deltaA is None:
            self._setLoad(app.sigLoadKind)
            self.lblNow.setText('No data loaded.')
            self.btn.setText('Convert')
            self.btn.setEnabled(False)
        elif app.sigKind == 'int':
            self._setLoad(app.sigOrigKind)
            self.lblNow.setText('The data are now:  %s' % app.sigShort('int'))
            self.btn.setText('Convert')
            self.btn.setEnabled(False)
            if msg == NOTE:
                msg = INT_TXT
        else:
            self._setLoad(app.sigOrigKind)
            other = 'dTT' if app.sigKind == 'dA' else 'dA'
            if app.sigKind == app.sigOrigKind:
                self.lblNow.setText('The data are now:  %s'
                                    % app.sigShort(app.sigKind))
            else:
                self.lblNow.setText(
                    'The data are now:  %s  (converted from %s)'
                    % (app.sigShort(app.sigKind),
                       app.sigShort(app.sigOrigKind)))
            self.btn.setText('Convert to %s' % app.sigShort(other))
            self.btn.setEnabled(True)
        if msg is not None:
            self.lblMsg.setText(msg)

    def onUnit(self):
        self.app.sigAUnit = self.ddUnit.currentText()
        self.app.savePrefs()

    def onLoadKind(self):
        app = self.app
        v = self.ddLoad.currentData()
        if app.deltaA is not None and app.sigKind != app.sigOrigKind:
            # converted data: what they were cannot be renamed
            self.refresh('The data shown were converted. Press "Reset all" '
                         'in the main window first, then change this '
                         'setting.')
            return
        app.sigLoadKind = v
        app.savePrefs()
        if app.deltaA is None:
            app.sigOrigKind = v
            app.sigKind = v
        if app.deltaA is not None and app.sigOrigKind != v:
            # the same numbers under another name: the analyses made so
            # far describe another quantity
            app.dataChanged()
            app.sigOrigKind = v
            app.sigKind = v
            app.updateStatusLabel()
            app.updateAll()
            msg = ('The data are now taken as %s. No number was changed.'
                   % app.sigShort(v))
            if v == 'int':
                msg += '  ' + INT_TXT
            self.refresh(msg)
        else:
            self.refresh(NOTE)

    def onConvert(self):
        app = self.app
        if app.deltaA is None:
            self.refresh('Load data first.')
            return
        if app.sigKind == 'int':
            self.refresh(INT_TXT)
            return
        toKind = 'dTT' if app.sigKind == 'dA' else 'dA'
        import numpy as np
        fin = np.abs(app.deltaA[np.isfinite(app.deltaA)])
        mx = float(np.max(fin)) if fin.size else 0.0
        if toKind == 'dTT' and app.sigAUnit == 'OD' and mx > 2:
            if not app.confirm(
                    'The largest |ΔA| is %.4g. Taken in OD this is a very '
                    'large signal. If the file is in mOD, cancel and choose '
                    'mOD.' % mx, 'Check the unit', yes='Convert as OD',
                    no='Cancel', parent=self):
                return
        elif toKind == 'dA' and mx > 1:
            if not app.confirm(
                    'The largest |ΔT/T| is %.4g. ΔT/T is taken as a fraction '
                    '(0.01 = 1 %%). If the file is in percent or another '
                    'unit, cancel and rescale it first.' % mx,
                    'Check the unit', yes='Convert', no='Cancel',
                    parent=self):
                return
        fromKind = app.sigKind
        try:
            nBad = app.applySignalConversion(toKind)
        except Exception as e:          # noqa: BLE001 (shown to the user)
            self.refresh('Conversion failed: %s' % e)
            return
        msg = 'Converted %s to %s (%s in %s).' % (
            app.sigShort(fromKind), app.sigShort(toKind),
            app.sigShort('dA'), app.sigAUnit)
        if nBad > 0:
            msg += ('  %d values with %s <= -1 have no %s and were left '
                    'empty.' % (nBad, app.sigShort('dTT'), app.sigShort('dA')))
        self.refresh(msg)
