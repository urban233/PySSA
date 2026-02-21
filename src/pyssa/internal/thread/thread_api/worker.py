from src.pyssa.gui.qt import QtCore
from src.pyssa.internal.thread.thread_api import worker_signals


class Worker(QtCore.QRunnable):
  def __init__(self, fn, *args, **kwargs):
    super().__init__()
    self.fn = fn
    self.args = args
    self.kwargs = kwargs
    self.signals = worker_signals.WorkerSignals()
    self._cancelled = False

  def cancel(self):
    self._cancelled = True

  def run(self):
    def progress_callback(value):
      if not self._cancelled:
        self.signals.progress.emit(value)

    try:
      if self._cancelled:
        self.signals.cancelled.emit()
        return

      # Inject progress + cancellation checker
      result = self.fn(
        progress_callback,
        lambda: self._cancelled,
        *self.args,
        **self.kwargs
      )

      if self._cancelled:
        self.signals.cancelled.emit()
      else:
        self.signals.success.emit(result)

    except Exception as e:
      self.signals.error.emit(e)
    finally:
      self.signals.finished.emit()
