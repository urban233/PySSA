"""Descriptor that carries all metadata the job scheduler needs to execute one job.

A ``JobDescriptor`` is created by the caller (e.g. a controller) and handed to
``JobScheduler.submit()``.

The descriptor is **mutable** so that the scheduler can transition a running
job from hot to cold when the user closes the project mid-execution (see
``JobScheduler.transition_to_cold``).  Only the scheduler should mutate
a descriptor after submission; callers treat it as effectively read-only.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, TYPE_CHECKING

from src.pyssa.util import enums

if TYPE_CHECKING:
  from src.pyssa.io_pyssa.db_pyssa.cold_project_handle import ColdProjectHandle


@dataclass
class JobDescriptor:
  """Bundle of information required to schedule and execute a job.

  The ``is_hot`` and ``cold_handle`` fields may be mutated by the
  scheduler when a project transitions from hot to cold while jobs
  are still in-flight.

  Args:
      job_type: Determines the queue the job is placed into.
      project_name: Name of the project this job belongs to.
      display_name: Optional human-readable label shown in the job table.
                    May be an empty string when not needed.
      run_fn: The work function submitted to ``ProcessRuntime``.
              Must be picklable (top-level function or staticmethod).
      run_args: Positional arguments forwarded to *run_fn*.
      run_kwargs: Keyword arguments forwarded to *run_fn*.
      is_hot: ``True`` when the job's project is currently open in the UI.
              May be set to ``False`` by the scheduler during a
              hot-to-cold transition.
      cold_handle: Write-only handle for persisting results into a cold
                   project's database.  ``None`` for hot-project jobs.
                   May be set by the scheduler during a hot-to-cold
                   transition.
      on_result: Callback invoked on the **main thread** when a hot-project
                 job finishes successfully.  Receives the worker's return
                 value.  ``None`` for cold-project jobs.
  """

  job_type: enums.JobType
  project_name: str
  display_name: str = ""
  run_fn: Callable[..., Any] = field(repr=False, default=None)
  run_args: tuple = field(default_factory=tuple)
  run_kwargs: dict[str, Any] = field(default_factory=dict)
  is_hot: bool = True
  cold_handle: "ColdProjectHandle | None" = field(default=None, repr=False)
  on_result: "Callable[[Any], None] | None" = field(default=None, repr=False)

