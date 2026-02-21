"""Job scheduler that replaces the legacy ``JobManager``.

The scheduler uses **ThreadRuntime** for orchestration (scheduling, queue
management, status updates on the main thread) and **ProcessRuntime** for
actual job execution — each ``JobType`` queue runs its jobs in a child
process so that different job types achieve true parallelism.

Scheduling rules
----------------
- Jobs of the **same** ``JobType`` run **sequentially** (one after another).
- Jobs of **different** ``JobType``s run **in parallel** (separate processes).
"""
from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any

from src.pyssa.gui.qt import QtCore

from src.pyssa.internal.data_structures.data_classes import (
  job_descriptor as jd_module,
)
from src.pyssa.internal.thread.thread_api.process_runtime import ProcessRuntime
from src.pyssa.internal.thread.thread_api.thread_runtime import (
  get_singleton_thread_runtime,
)
from src.pyssa.model import job_model as jm_module
from src.pyssa.util import enums

logger = logging.getLogger(__name__)


class _TypeQueue:
  """Per-``JobType`` sequential queue managed by the scheduler.

  Attributes:
      pending: FIFO list of ``(row, descriptor)`` tuples waiting to run.
      is_running: ``True`` while a job of this type is executing.
  """

  __slots__ = ("pending", "is_running")

  def __init__(self) -> None:
    self.pending: list[tuple[int, jd_module.JobDescriptor]] = []
    self.is_running: bool = False


class JobScheduler(QtCore.QObject):
  """Orchestrates job execution across multiple ``JobType`` queues.

  Args:
      model: The ``JobModel`` that tracks all jobs.
      process_runtime: The ``ProcessRuntime`` used for executing jobs in
                       child processes.  A default instance is created if
                       ``None``.
      parent: Optional Qt parent.
  """

  def __init__(
      self,
      model: jm_module.JobModel,
      process_runtime: ProcessRuntime | None = None,
      parent: QtCore.QObject | None = None,
  ) -> None:
    super().__init__(parent)
    if model is None:
      raise ValueError("model must not be None")

    self._model = model
    self._process_runtime = process_runtime or ProcessRuntime()
    self._thread_runtime = get_singleton_thread_runtime()
    self._queues: dict[enums.JobType, _TypeQueue] = defaultdict(_TypeQueue)

  # ------------------------------------------------------------------
  # Public API
  # ------------------------------------------------------------------

  def submit(self, descriptor: jd_module.JobDescriptor) -> int:
    """Enqueue a job and start its type-queue if idle.

    Args:
        descriptor: Immutable job metadata describing what to execute.

    Returns:
        The row index assigned to the job in the ``JobModel``.

    Raises:
        ValueError: If *descriptor* is ``None``.
    """
    if descriptor is None:
      raise ValueError("descriptor must not be None")

    row = self._model.add_job(descriptor)
    type_queue = self._queues[descriptor.job_type]
    type_queue.pending.append((row, descriptor))

    logger.info(
      "Job submitted: row=%d, type=%s, project=%s",
      row, descriptor.job_type.value, descriptor.project_name,
    )

    if not type_queue.is_running:
      self._drain_next(descriptor.job_type)

    return row

  def cancel(self, row: int) -> None:
    """Cancel a queued job (basic cancellation).

    If the job is still ``QUEUED`` it is removed from its type queue and
    marked ``CANCELLED`` in the model.  If the job is already ``RUNNING``
    this method is a no-op (advanced cancellation is not implemented).

    Args:
        row: The row index in the ``JobModel``.

    Raises:
        IndexError: If *row* is out of range.
    """
    status = self._model.get_status(row)
    if status != enums.JobStatus.QUEUED:
      logger.debug(
        "Cannot cancel job at row %d — status is %s.", row, status.value,
      )
      return

    descriptor = self._model.get_descriptor(row)
    type_queue = self._queues.get(descriptor.job_type)
    if type_queue is not None:
      type_queue.pending = [
        (r, d) for r, d in type_queue.pending if r != row
      ]

    self._model.update_status(row, enums.JobStatus.CANCELLED)
    logger.info("Job at row %d cancelled.", row)

  def has_running_jobs(self) -> bool:
    """Check whether any queue is currently executing a job.

    Returns:
        ``True`` if at least one type-queue is active.
    """
    return any(q.is_running for q in self._queues.values())

  def shutdown(self) -> None:
    """Gracefully stop all queues.

    Clears every pending queue and marks remaining queued jobs as
    ``CANCELLED``.  Already-running jobs are allowed to finish.
    """
    for job_type, type_queue in self._queues.items():
      for row, _ in type_queue.pending:
        self._model.update_status(row, enums.JobStatus.CANCELLED)
      type_queue.pending.clear()
    logger.info("JobScheduler shutdown: all pending queues cleared.")

  # ------------------------------------------------------------------
  # Internal: draining a type-queue one job at a time
  # ------------------------------------------------------------------

  def _drain_next(self, job_type: enums.JobType) -> None:
    """Pop the next pending job for *job_type* and execute it.

    If the queue is empty the ``is_running`` flag is cleared and the
    method returns.

    Args:
        job_type: The job type whose queue should be drained.
    """
    type_queue = self._queues[job_type]

    if not type_queue.pending:
      type_queue.is_running = False
      logger.debug("Queue for %s is empty — stopping.", job_type.value)
      return

    type_queue.is_running = True
    row, descriptor = type_queue.pending.pop(0)
    self._model.update_status(row, enums.JobStatus.RUNNING)

    logger.info(
      "Executing job: row=%d, type=%s, project=%s",
      row, descriptor.job_type.value, descriptor.project_name,
    )

    process_task = self._process_runtime.run(
      descriptor.run_fn, *descriptor.run_args, **descriptor.run_kwargs,
    )

    process_task.on_success(
      lambda result, _r=row, _d=descriptor, _jt=job_type: self._on_job_success(
        _r, _d, _jt, result,
      ),
    )
    process_task.on_error(
      lambda exc, _r=row, _d=descriptor, _jt=job_type: self._on_job_error(
        _r, _d, _jt, exc,
      ),
    )

  def _on_job_success(
      self,
      row: int,
      descriptor: jd_module.JobDescriptor,
      job_type: enums.JobType,
      result: Any,
  ) -> None:
    """Handle successful completion of a job.

    Updates the model, invokes the hot-project callback if applicable,
    closes the cold handle if applicable, then drains the next job.

    Args:
        row: The row index in the model.
        descriptor: The job's metadata.
        job_type: The type used to look up the queue.
        result: The return value from the worker function.
    """
    self._model.update_status(row, enums.JobStatus.FINISHED, result)
    logger.info("Job at row %d finished successfully.", row)

    if descriptor.is_hot and descriptor.on_result is not None:
      try:
        descriptor.on_result(result)
      except Exception:
        logger.exception(
          "on_result callback failed for job at row %d.", row,
        )

    if descriptor.cold_handle is not None:
      try:
        descriptor.cold_handle.close()
      except Exception:
        logger.exception(
          "Failed to close cold handle for job at row %d.", row,
        )

    self._drain_next(job_type)

  def _on_job_error(
      self,
      row: int,
      descriptor: jd_module.JobDescriptor,
      job_type: enums.JobType,
      exc: Exception,
  ) -> None:
    """Handle a failed job.

    Updates the model, closes the cold handle if applicable, then drains
    the next job.

    Args:
        row: The row index in the model.
        descriptor: The job's metadata.
        job_type: The type used to look up the queue.
        exc: The exception raised by the worker.
    """
    self._model.update_status(row, enums.JobStatus.FAILED, exc)
    logger.error("Job at row %d failed: %s", row, exc)

    if descriptor.cold_handle is not None:
      try:
        descriptor.cold_handle.close()
      except Exception:
        logger.exception(
          "Failed to close cold handle for job at row %d.", row,
        )

    self._drain_next(job_type)
