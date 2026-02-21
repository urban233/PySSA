from concurrent.futures import ProcessPoolExecutor
from src.pyssa.internal.thread.thread_api import process_task


class ProcessRuntime:
  def __init__(self, max_workers=None):
    self.executor = ProcessPoolExecutor(max_workers=max_workers)

  def run(self, fn, *args):
    future = self.executor.submit(fn, *args)
    return process_task.ProcessTask(future)
