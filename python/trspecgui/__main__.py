"""Start the program:  python -m trspecgui  [data file]"""
import os
import sys
import traceback


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    from PySide6 import QtWidgets
    from .app import TRSpecAnalyzer
    from .widgets import applyLightTheme, prepareEnvironment

    prepareEnvironment()
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication(argv)
    applyLightTheme(qapp)

    def showError(etype, value, tb):
        text = ''.join(traceback.format_exception(etype, value, tb))
        sys.stderr.write(text)
        try:
            box = QtWidgets.QMessageBox()
            box.setIcon(QtWidgets.QMessageBox.Icon.Critical)
            box.setWindowTitle('TRSpecAnalyzer - unexpected error')
            box.setText('An unexpected error occurred. The program keeps '
                        'running; the details below help to find the cause.')
            box.setDetailedText(text)
            box.exec()
        except Exception:
            pass

    sys.excepthook = showError
    app = TRSpecAnalyzer()
    app.fig.show()
    if len(argv) > 1 and os.path.isfile(argv[1]):
        try:
            app.loadFile(argv[1])
        except Exception as err:
            app.alert('Failed to read file:\n%s' % err, 'Load Error')
    return qapp.exec()


if __name__ == '__main__':
    sys.exit(main())
