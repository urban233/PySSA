from src.pyssa.gui.qt import QtCore, pyqtSignal


class WorkerSignals(QtCore.QObject):
  progress = pyqtSignal(object)
  success = pyqtSignal(object)
  error = pyqtSignal(Exception)
  finished = pyqtSignal()
  cancelled = pyqtSignal()
