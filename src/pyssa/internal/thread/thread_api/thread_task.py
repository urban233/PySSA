from src.pyssa.internal.thread.thread_api import worker

class ThreadTask:
  def __init__(self, a_worker: "worker.Worker"):
    self.worker = a_worker

  def on_success(self, fn):
    self.worker.signals.success.connect(fn)
    return self

  def on_error(self, fn):
    self.worker.signals.error.connect(fn)
    return self

  def on_progress(self, fn):
    self.worker.signals.progress.connect(fn)
    return self
