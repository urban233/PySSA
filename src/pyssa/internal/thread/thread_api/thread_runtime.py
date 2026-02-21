from src.pyssa.gui.qt import QtCore
from src.pyssa.internal.thread.thread_api import thread_task
from src.pyssa.internal.thread.thread_api import worker


class ThreadRuntime(QtCore.QObject):
  def __init__(self, thread_count: int = 4):
    super().__init__()
    self._pool = QtCore.QThreadPool().globalInstance()
    # self._pool.setMaxThreadCount(thread_count)

  def run(self, fn, *args, **kwargs) -> "thread_task.ThreadTask":
    tmp_worker = worker.Worker(fn, *args, **kwargs)
    tmp_handle = thread_task.ThreadTask(tmp_worker)
    self._pool.start(tmp_worker)
    return tmp_handle


_singleton_instance: ThreadRuntime | None = None

def get_singleton_thread_runtime() -> ThreadRuntime:
  global _singleton_instance
  if _singleton_instance is None:
    _singleton_instance = ThreadRuntime()
  return _singleton_instance
