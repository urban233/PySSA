# cold_project_handle.py
"""
ColdProjectHandle: a narrow, write-only interface handed to background jobs
that need to persist results into a project that is not currently open in
the UI.

Design goals
------------
- The job receives only what it needs: the ability to submit writes.
  It does not get a reference to MainWindowController, AppState, or the
  full ProjectDatabase.
- The handle is cheap to pass around (it is just a thin wrapper).
- When the job is finished it calls handle.close(), which drains the write
  queue and releases connections.  The controller is notified via an
  optional on_closed callback so it can remove the entry from _cold_dbs.

Usage (inside a background job worker function)
-----------------------------------------------
    def _my_analysis_worker(progress_callback, is_cancelled, handle, pair_obj):
        # ... do analysis work ...
        handle.submit(WriteOperation(OperationType.INSERT_PROTEIN_PAIR, pair_obj))
        handle.close()
"""
from __future__ import annotations

import logging
from typing import Callable, TYPE_CHECKING

from src.pyssa.io_pyssa.db_pyssa.write_queue import WriteOperation

if TYPE_CHECKING:
  from src.pyssa.io_pyssa.db_pyssa.project_database import ProjectDatabase

logger = logging.getLogger(__name__)


class ColdProjectHandle:
  """Write-only interface to a cold project's database.

  Args:
      project_id:  The stable identifier for this project.
      db:          The underlying ProjectDatabase (kept private).
      on_closed:   Optional callback invoked after close() drains and
                   releases the database.  The controller uses this to
                   remove the entry from _cold_dbs without the job
                   needing to know about the controller.
  """

  def __init__(
          self,
          project_id: str,
          db: "ProjectDatabase",
          on_closed: Callable[[str], None] | None = None,
  ) -> None:
    self._project_id = project_id
    self._db = db
    self._on_closed = on_closed
    self._closed = False

  # ------------------------------------------------------------------
  # Public API exposed to background jobs
  # ------------------------------------------------------------------

  @property
  def project_id(self) -> str:
    return self._project_id

  def submit(self, operation: WriteOperation) -> None:
    """Enqueue a fire-and-forget write. Thread-safe."""
    if self._closed:
      raise RuntimeError(
        f"ColdProjectHandle for '{self._project_id}' is already closed."
      )
    self._db.write_queue.submit(operation)

  def close(self, drain_timeout: float = 10.0) -> None:
    """Drain pending writes, close the database, notify the controller.

    Safe to call multiple times — subsequent calls are no-ops.
    """
    if self._closed:
      return
    self._closed = True
    self._db.close(drain_timeout=drain_timeout)
    if self._on_closed is not None:
      self._on_closed(self._project_id)
    logger.info("ColdProjectHandle for '%s' closed.", self._project_id)
