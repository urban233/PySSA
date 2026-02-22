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
from typing import Any, Callable, TYPE_CHECKING

from src.pyssa.gui.qt import QtCore

from src.pyssa.internal.data_structures.data_classes import job_descriptor
from src.pyssa.internal.thread.thread_api.process_runtime import ProcessRuntime
from src.pyssa.internal.thread.thread_api.thread_runtime import (
  get_singleton_thread_runtime,
)
from src.pyssa.model import job_model
from src.pyssa.util import enums

if TYPE_CHECKING:
  from src.pyssa.io_pyssa.db_pyssa.cold_project_handle import ColdProjectHandle
  from src.pyssa.gui.name_registry import NameRegistry

logger = logging.getLogger(__name__)


class _TypeQueue:
  """Per-``JobType`` sequential queue managed by the scheduler.

  Attributes:
      pending: FIFO list of ``(row, descriptor)`` tuples waiting to run.
      is_running: ``True`` while a job of this type is executing.
  """

  __slots__ = ("pending", "is_running")

  def __init__(self) -> None:
    self.pending: list[tuple[int, job_descriptor.JobDescriptor]] = []
    self.is_running: bool = False


class JobScheduler(QtCore.QObject):
  """Orchestrates job execution across multiple ``JobType`` queues.

  Args:
      model: The ``JobModel`` that tracks all jobs.
      name_registry: The application-wide ``NameRegistry`` used to
                     automatically reserve and release names as jobs
                     are submitted and complete.  Pass ``None`` to
                     disable name-reservation (e.g. in tests).
      process_runtime: The ``ProcessRuntime`` used for executing jobs in
                       child processes.  A default instance is created if
                       ``None``.
      parent: Optional Qt parent.
  """

  def __init__(
      self,
      model: "job_model.JobModel",
      name_registry: "NameRegistry | None" = None,
      process_runtime: ProcessRuntime | None = None,
      parent: QtCore.QObject | None = None,
  ) -> None:
    super().__init__(parent)
    if model is None:
      raise ValueError("model must not be None")

    self._model = model
    self._name_registry = name_registry
    self._process_runtime = process_runtime or ProcessRuntime()
    self._thread_runtime = get_singleton_thread_runtime()
    self._queues: dict[enums.JobType, _TypeQueue] = defaultdict(_TypeQueue)
    self._active_tasks: dict[int, Any] = {}  # row -> ProcessTask, prevents GC

  # ------------------------------------------------------------------
  # Public API
  # ------------------------------------------------------------------

  def submit(self, descriptor: "job_descriptor.JobDescriptor") -> int:
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

    # Reserve the names this job will produce so they cannot be re-used
    # until the job completes, fails, or is cancelled.
    if self._name_registry is not None:
      for category, names in descriptor.reserved_names.items():
        self._name_registry.reserve_many(category, names)

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

    # Release names that were reserved for this cancelled job.
    if self._name_registry is not None:
      for category, names in descriptor.reserved_names.items():
        self._name_registry.release_many(category, names)

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

  def transition_to_cold(
      self,
      project_name: str,
      cold_handle: "ColdProjectHandle",
  ) -> None:
    """Transition all jobs for *project_name* from hot to cold.

    Called by ``AppState.close_project()`` when the user closes a project
    that still has running or queued jobs.  For every matching job:

    - ``descriptor.is_hot`` is set to ``False``.
    - ``descriptor.cold_handle`` is set to the provided handle.
    - ``descriptor.on_result`` is cleared (no UI callback after close).

    When these jobs eventually finish, ``_on_job_success`` will see
    ``is_hot=False`` and skip the ``on_result`` callback.  Instead it
    will close the cold handle, persisting the results in the database.

    Args:
        project_name: The name of the project being closed.
        cold_handle: A ``ColdProjectHandle`` that jobs should use to
                     persist their results.

    Raises:
        ValueError: If *project_name* is empty or *cold_handle* is
                    ``None``.
    """
    if not project_name:
      raise ValueError("project_name must not be empty")
    if cold_handle is None:
      raise ValueError("cold_handle must not be None")

    transitioned = 0

    for type_queue in self._queues.values():
      for row, descriptor in type_queue.pending:
        if descriptor.project_name == project_name and descriptor.is_hot:
          descriptor.is_hot = False
          descriptor.cold_handle = cold_handle
          descriptor.on_result = None
          transitioned += 1

    total_rows = self._model.rowCount()
    for row in range(total_rows):
      status = self._model.get_status(row)
      if status != enums.JobStatus.RUNNING:
        continue
      descriptor = self._model.get_descriptor(row)
      if descriptor.project_name == project_name and descriptor.is_hot:
        descriptor.is_hot = False
        descriptor.cold_handle = cold_handle
        descriptor.on_result = None
        transitioned += 1

    logger.info(
      "Transitioned %d job(s) for project '%s' from hot to cold.",
      transitioned, project_name,
    )

  def transition_to_hot(
      self,
      project_name: str,
      on_result: Callable[[Any], None] | None = None,
  ) -> int:
    """Transition all cold jobs for *project_name* back to hot.

    Called by ``AppState.open_project()`` when the user opens a project
    that still has running or queued cold jobs.  For every matching job:

    - ``descriptor.is_hot`` is set to ``True``.
    - ``descriptor.on_result`` is set to the provided callback.
    - ``descriptor.cold_handle`` is set to ``None`` (the cold DB
      is no longer needed since the hot DB is now open).

    When these jobs eventually finish, ``_on_job_success`` will see
    ``is_hot=True`` and invoke the ``on_result`` callback to update
    the UI model.

    Args:
        project_name: The name of the project being opened.
        on_result: Callback invoked on the **main thread** when a
                   transitioned job finishes.  Receives the worker's
                   return value.  May be ``None`` if the caller only
                   needs the model signals.

    Returns:
        The number of jobs that were transitioned.

    Raises:
        ValueError: If *project_name* is empty.
    """
    if not project_name:
      raise ValueError("project_name must not be empty")

    transitioned = 0

    for type_queue in self._queues.values():
      for row, descriptor in type_queue.pending:
        if descriptor.project_name == project_name and not descriptor.is_hot:
          descriptor.is_hot = True
          descriptor.on_result = on_result
          descriptor.cold_handle = None
          transitioned += 1

    total_rows = self._model.rowCount()
    for row in range(total_rows):
      status = self._model.get_status(row)
      if status != enums.JobStatus.RUNNING:
        continue
      descriptor = self._model.get_descriptor(row)
      if descriptor.project_name == project_name and not descriptor.is_hot:
        descriptor.is_hot = True
        descriptor.on_result = on_result
        descriptor.cold_handle = None
        transitioned += 1

    logger.info(
      "Transitioned %d job(s) for project '%s' from cold to hot.",
      transitioned, project_name,
    )
    return transitioned

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

    # Store reference to prevent garbage collection before signals fire
    self._active_tasks[row] = process_task

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
      descriptor: job_descriptor.JobDescriptor,
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
    # Clean up task reference now that callback has fired
    self._active_tasks.pop(row, None)

    self._model.update_status(row, enums.JobStatus.FINISHED, result)
    logger.info("Job at row %d finished successfully.", row)

    # Release names reserved by this job; they are now part of the project.
    if self._name_registry is not None:
      for category, names in descriptor.reserved_names.items():
        self._name_registry.release_many(category, names)

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
      descriptor: job_descriptor.JobDescriptor,
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
    # Clean up task reference now that callback has fired
    self._active_tasks.pop(row, None)

    self._model.update_status(row, enums.JobStatus.FAILED, exc)
    logger.error("Job at row %d failed: %s", row, exc)

    # Release names reserved by the failed job so they become available again.
    if self._name_registry is not None:
      for category, names in descriptor.reserved_names.items():
        self._name_registry.release_many(category, names)

    if descriptor.cold_handle is not None:
      try:
        descriptor.cold_handle.close()
      except Exception:
        logger.exception(
          "Failed to close cold handle for job at row %d.", row,
        )

    self._drain_next(job_type)
