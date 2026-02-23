# write_queue.py
"""
ProjectWriteQueue: serialises fire-and-forget database writes for one
project using the application's existing ThreadRuntime / Worker API.

Why serialisation matters
--------------------------
QThreadPool can run tasks concurrently.  If we submitted each write
operation as an independent Worker, two writes for the same project could
execute simultaneously, corrupting the database.

This queue solves the problem by using a simple lock + "is a worker
currently running?" flag:

  - When a write is submitted and no worker is running, a new Worker is
    started immediately to drain the queue.
  - When a write is submitted while a worker is already running, it is
    placed in the queue and will be picked up by the running worker on
    its next iteration.
  - When the running worker empties the queue it stops.  The next submit
    will start a new worker.

This means at most one Worker per project is active at any time, giving
us the same serialisation guarantee as the old dedicated thread, but
using the shared thread pool instead.

Usage
-----
    queue = ProjectWriteQueue(db, thread_runtime)

    queue.submit(WriteOperation(OperationType.INSERT_PROTEIN, protein_obj))
    queue.submit(WriteOperation(OperationType.UPDATE_PROTEIN_SESSION, protein_obj))

    queue.drain()   # call before closing the database (blocks until empty)
"""
from __future__ import annotations
from typing import TYPE_CHECKING

import logging
import queue
import threading
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from src.pyssa.internal.thread.thread_api.thread_runtime import ThreadRuntime

if TYPE_CHECKING:
  from src.pyssa.io_pyssa.db_pyssa.project_database import ProjectDatabase

logger = logging.getLogger(__name__)


class OperationType(Enum):
  INSERT_PROTEIN            = "insert_protein"
  DELETE_PROTEIN            = "delete_protein"
  DELETE_CHAIN              = "delete_chain"
  INSERT_PROTEIN_PAIR       = "insert_protein_pair"
  DELETE_PROTEIN_PAIR       = "delete_protein_pair"
  UPDATE_PROTEIN_SESSION    = "update_protein_session"
  UPDATE_PAIR_SESSION       = "update_pair_session"
  UPDATE_PROTEIN_PDB_DATA   = "update_protein_pdb_data"
  INSERT_SEQUENCE           = "insert_sequence"
  DELETE_SEQUENCE           = "delete_sequence"
  UPDATE_SEQUENCE_NAME      = "update_sequence_name"


@dataclass
class WriteOperation:
  """A single unit of deferred database work.

  Args:
      op_type: What to do.
      payload: Data consumed by the handler; type depends on op_type.
               See ProjectWriteQueue._dispatch for expected shapes.
  """
  op_type: OperationType
  payload: Any = field(default=None)


class ProjectWriteQueue:
  """Serialised write queue for one ProjectDatabase, backed by ThreadRuntime.

  Args:
      db:             The ProjectDatabase this queue writes to.
      thread_runtime: The shared ThreadRuntime.
                      Defaults to the application singleton.
      max_queue_size: Maximum number of operations to buffer. When reached,
                      submit() will block until space is available.
  """

  def __init__(
          self,
          db: "ProjectDatabase",
          thread_runtime: ThreadRuntime | None = None,
          max_queue_size: int = 1000,
  ) -> None:
    if db is None:
      raise ValueError("db cannot be None")
    if max_queue_size < 1:
      raise ValueError("max_queue_size must be at least 1")

    self._db = db
    self._runtime = thread_runtime or ThreadRuntime()
    # Use bounded Queue instead of SimpleQueue to prevent memory overflow
    self._queue: queue.Queue[WriteOperation] = queue.Queue(maxsize=max_queue_size)
    self._lock = threading.Lock()
    self._worker_running = False
    self._failed_operations_count = 0
    self._successful_operations_count = 0

    logger.info("ProjectWriteQueue initialized with max_queue_size=%d", max_queue_size)

  # ------------------------------------------------------------------
  # Public API
  # ------------------------------------------------------------------

  def submit(self, operation: WriteOperation, timeout: float | None = 5.0) -> None:
    """Enqueue *operation* for asynchronous execution.

    Thread-safe.  Blocks if the queue is full until space is available
    or timeout expires.

    Args:
        operation: The write operation to enqueue.
        timeout: Maximum seconds to wait if queue is full. None = wait forever.

    Raises:
        queue.Full: If the queue is full and timeout expires.
        ValueError: If operation is None.
    """
    if operation is None:
      raise ValueError("operation cannot be None")

    try:
      self._queue.put(operation, timeout=timeout)
      logger.debug("Queued operation: %s", operation.op_type)
      self._ensure_worker_running()
    except queue.Full:
      logger.error(
        "Write queue is full (max_size=%d). Operation %s was rejected.",
        self._queue.maxsize, operation.op_type
      )
      raise

  def drain(self, timeout: float = 10.0) -> bool:
    """Block until the queue is empty or *timeout* seconds have passed.

    Call this before closing the database to guarantee all pending
    writes have been committed.

    Returns:
        True if the queue was successfully drained, False if timed out.
    """
    done = threading.Event()
    success = [False]  # Use list to allow mutation in nested function

    def _check() -> None:
      # Poll until the queue is empty and no worker is active.
      import time
      deadline = time.monotonic() + timeout
      while time.monotonic() < deadline:
        with self._lock:
          if self._queue.empty() and not self._worker_running:
            success[0] = True
            done.set()
            logger.info(
              "Queue drained successfully. Stats: %d successful, %d failed",
              self._successful_operations_count, self._failed_operations_count
            )
            return
        time.sleep(0.05)
      logger.warning(
        "ProjectWriteQueue.drain() timed out after %.1fs. "
        "Queue size: %d, Worker running: %s",
        timeout, self._queue.qsize(), self._worker_running
      )
      done.set()

    t = threading.Thread(target=_check, daemon=True)
    t.start()
    done.wait()
    return success[0]

  # ------------------------------------------------------------------
  # Internal: worker lifecycle
  # ------------------------------------------------------------------

  def _ensure_worker_running(self) -> None:
    with self._lock:
      if self._worker_running:
        return          # existing worker will pick up the new item
      self._worker_running = True

    # Create stats dict to share with worker
    stats = {
      "successful": 0,
      "failed": 0,
    }

    (
      self._runtime
      .run(_drain_worker, self._db, self._queue, stats)
      .on_success(lambda _: self._on_worker_finished(stats))
      .on_error(self._on_worker_error)
    )

  def _on_worker_finished(self, stats: dict[str, int]) -> None:
    # Called on the main thread after the worker returns.
    # Update global stats
    self._successful_operations_count += stats["successful"]
    self._failed_operations_count += stats["failed"]

    # Check if new items arrived while the worker was finishing.
    with self._lock:
      if self._queue.empty():
        self._worker_running = False
        return
      # Items arrived in the narrow window — reuse the "running" flag
      # and start a new worker immediately.

    # Create new stats dict for new worker
    new_stats = {
      "successful": 0,
      "failed": 0,
    }

    (
      self._runtime
      .run(_drain_worker, self._db, self._queue, new_stats)
      .on_success(lambda _: self._on_worker_finished(new_stats))
      .on_error(self._on_worker_error)
    )

  def _on_worker_error(self, exc: Exception) -> None:
    logger.exception("ProjectWriteQueue worker raised an error.", exc_info=exc)
    with self._lock:
      self._worker_running = False
      self._failed_operations_count += 1
    # Re-start if items are still waiting — the failed operation is lost
    # (already logged) but subsequent operations should still execute.
    if not self._queue.empty():
      logger.info("Restarting worker after error. Queue size: %d", self._queue.qsize())
      self._ensure_worker_running()

  # ------------------------------------------------------------------
  # Statistics
  # ------------------------------------------------------------------

  def get_queue_size(self) -> int:
    """Return the current number of pending operations."""
    return self._queue.qsize()

  def get_stats(self) -> dict[str, int]:
    """Return statistics about queue operations.

    Returns:
        Dictionary with keys: 'successful', 'failed', 'pending', 'max_size'
    """
    return {
      "successful": self._successful_operations_count,
      "failed": self._failed_operations_count,
      "pending": self._queue.qsize(),
      "max_size": self._queue.maxsize,
    }


# ---------------------------------------------------------------------------
# Worker function — runs on a QThreadPool thread via ThreadRuntime.
#
# Worker injects (progress_callback, is_cancelled) as the first two positional
# args before our own args.
# ---------------------------------------------------------------------------

def _drain_worker(
        progress_callback,                          # injected by Worker
        is_cancelled,                               # injected by Worker
        db: "ProjectDatabase",
        work_queue: queue.Queue[WriteOperation],
        stats: dict[str, int],  # Shared stats dictionary
) -> None:
  """Drain as many operations as possible from *work_queue*.

  Stops when the queue is empty or cancellation is requested.
  Errors in individual operations are logged and skipped so that one
  bad write does not block subsequent ones.

  Args:
      progress_callback: Injected by Worker (unused).
      is_cancelled: Injected by Worker to check cancellation.
      db: The ProjectDatabase to write to.
      work_queue: The queue of pending operations.
      stats: Dictionary to track successful/failed operation counts.
  """
  operations_processed = 0
  while not is_cancelled():
    try:
      op = work_queue.get_nowait()
    except queue.Empty:
      logger.debug("Worker processed %d operations before queue empty", operations_processed)
      return          # queue exhausted — signal back via on_success

    try:
      _dispatch(db, op)
      stats["successful"] += 1
      operations_processed += 1
      logger.debug("Successfully executed operation: %s", op.op_type)
    except Exception:
      stats["failed"] += 1
      logger.exception(
        "Write operation %s failed and was skipped.", op.op_type
      )


# ---------------------------------------------------------------------------
# Dispatch — maps OperationType to ProjectDatabase calls
# ---------------------------------------------------------------------------

def _dispatch(db: "ProjectDatabase", op: WriteOperation) -> None:
  match op.op_type:
    case OperationType.INSERT_PROTEIN:
      db.insert_protein_full(op.payload)

    case OperationType.DELETE_PROTEIN:
      db.delete_protein_full(op.payload)          # payload: protein_id

    case OperationType.DELETE_CHAIN:
      protein_id, chain_id = op.payload
      db.delete_chain(protein_id, chain_id)

    case OperationType.INSERT_PROTEIN_PAIR:
      db.insert_protein_pair_full(op.payload)

    case OperationType.DELETE_PROTEIN_PAIR:
      db.delete_protein_pair_full(op.payload)     # payload: pair_id

    case OperationType.UPDATE_PROTEIN_SESSION:
      # payload: Protein domain object
      db.update_protein_session(
        op.payload.get_id(), op.payload.pymol_session
      )

    case OperationType.UPDATE_PROTEIN_PDB_DATA:
      # payload: Protein domain object
      db.replace_pdb_atoms(
          op.payload.get_id(),
          op.payload.get_pdb_data()
      )

    case OperationType.UPDATE_PAIR_SESSION:
      # payload: ProteinPair domain object
      db.update_protein_pair_session(
        op.payload.get_id(), op.payload.pymol_session
      )

    case OperationType.INSERT_SEQUENCE:
      # payload: SeqRecord
      db.insert_sequence(
        str(op.payload.id),
        str(op.payload.seq),
        op.payload.name,
        project_id=1,       # single-project databases always use id=1
      )

    case OperationType.DELETE_SEQUENCE:
      db.delete_sequence(op.payload.name)         # payload: SeqRecord

    case OperationType.UPDATE_SEQUENCE_NAME:
      new_name, old_name, seq = op.payload
      db.update_sequence_name(new_name, old_name, seq)

    case _:
      logger.warning("Unknown OperationType: %s", op.op_type)
