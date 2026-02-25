from src.pyssa.gui.qt import QtCore
from src.pyssa.internal.thread.thread_api import worker

class ThreadTask:
  def __init__(self, a_worker: "worker.Worker", a_pool: QtCore.QThreadPool):
    self.worker = a_worker
    self._pool = a_pool
    self._started = False

  def start(self):
    """Starts the task if it hasn't already been started.

    Returns:
        ThreadTask: The current task object.
    """
    if not self._started:
      self._started = True
      self._pool.start(self.worker)
    return self

  def on_success(self, fn):
    self.worker.signals.success.connect(fn)
    return self

  def on_error(self, fn):
    self.worker.signals.error.connect(fn)
    return self

  def on_progress(self, fn):
    self.worker.signals.progress.connect(fn)
    return self
