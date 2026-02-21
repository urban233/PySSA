from src.pyssa.gui.qt import QtCore, pyqtSignal
from concurrent.futures import ProcessPoolExecutor


class ProcessTask(QtCore.QObject):
  success = pyqtSignal(object)
  error = pyqtSignal(Exception)

  def __init__(self, future):
    super().__init__()
    self._future = future
    self._timer = QtCore.QTimer()
    self._timer.setInterval(50)
    self._timer.timeout.connect(self._poll)
    self._timer.start()

  def _poll(self):
    if self._future.done():
      self._timer.stop()
      try:
        result = self._future.result()
        self.success.emit(result)
      except Exception as e:
        self.error.emit(e)

  # Fluent methods
  def on_success(self, fn):
    self.success.connect(fn)
    return self

  def on_error(self, fn):
    self.error.connect(fn)
    return self
