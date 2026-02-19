"""Dedicated PyMOL worker process with a clean, thread-safe public API.

Design goals:
    Maximum performance: PyMOL sessions are expensive to load, so each worker
        caches the last loaded session path and only reloads when it changes.
    Thread safety: A threading.Lock serialises all pipe I/O so multiple threads
        in the host process can share one worker without races.
    Parallel workers: Workers are independent OS processes; any number can run
        concurrently.
    Clean public API: Callers never touch raw pipe connections or WorkerCommand
        objects directly.

Example:
    Typical usage with an explicit lifecycle::

        worker = PmlWorker()
        worker.start()
        worker.set_session_path("/path/to/session.pse")

        # Fire-and-forget (async)
        worker.do("zoom")

        # Synchronous call that waits for completion
        result = worker.do("select", args=("sele", "resi 50"), sync=True)

        worker.stop()

    For one-off tasks prefer one_shot_do() or the context manager::

        with PmlWorker.session("/path/to/session.pse") as worker:
            worker.do("zoom", sync=True)
"""
from __future__ import annotations

import logging
import pathlib
import threading
from contextlib import contextmanager
from multiprocessing import Process
from multiprocessing.connection import Connection, Pipe
from typing import Any, Generator

import pymol2

from pyssa.io_pyssa import binary_data
from src.pyssa.gui.user_pymol import UserPyMOL
from src.pyssa.internal.pymol.worker_command import CommandType, WorkerCommand
from src.pyssa.util import constants

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


class PmlWorker:
  """A handle to a dedicated PyMOL worker sub-process.

  All public methods are thread-safe. Concurrent callers from different
  threads are serialised through an internal lock and will not corrupt each
  other's communication. Independent PmlWorker instances are fully parallel —
  they each own their own OS process and pipe.
  """

  def __init__(self) -> None:
    self._process: Process | None = None
    self._conn: Connection | None = None  # host-side pipe end
    self._session_path: str | None = None
    # Serialises all pipe send/recv calls so multiple threads can share
    # one worker instance safely.
    self._lock = threading.Lock()

  # ------------------------------------------------------------------
  # Lifecycle
  # ------------------------------------------------------------------

  def start(self) -> None:
    """Spawns the background PyMOL process. Must be called before do()."""
    if self._process is not None and self._process.is_alive():
      raise RuntimeError("Worker is already running.")
    host_conn, worker_conn = Pipe()
    self._conn = host_conn
    self._process = Process(
      target=_worker_loop,
      args=(worker_conn,),
      daemon=True,  # Automatically reaped if the host process dies.
      name="PmlWorker",
    )
    self._process.start()
    # The worker-side connection object is no longer needed in *this*
    # process; closing it here prevents fd leaks.
    worker_conn.close()
    logger.debug("PmlWorker started (pid=%d).", self._process.pid)

  def stop(self, timeout: float = 10.0) -> None:
    """Gracefully shuts down the worker process.

    Sends a shutdown command and waits up to ``timeout`` seconds for the
    process to terminate. If it does not exit in time it is forcibly killed.

    Args:
        timeout: Seconds to wait for a clean shutdown before killing the process.
    """
    if self._process is None or not self._process.is_alive():
      return
    with self._lock:
      try:
        self._conn.send(WorkerCommand.shutdown())
      except OSError:
        pass  # Pipe already broken; process likely already dead.
      finally:
        # Closing the host-side pipe end is what actually unblocks the
        # worker's recv() via EOFError. Do it inside the lock so no
        # other thread can attempt a send after this point.
        self._conn.close()
    self._process.join(timeout=timeout)
    if self._process.is_alive():
      logger.warning("Worker did not exit cleanly; killing it.")
      self._process.kill()
      self._process.join()
    logger.debug("PmlWorker stopped.")

  # ------------------------------------------------------------------
  # Session management
  # ------------------------------------------------------------------

  def set_session_path(self, session_path: str | pathlib.Path) -> None:
    """Sets (or changes) the PyMOL session file the worker should use.

    Args:
        session_path: Path to the ``.pse`` session file.
    """
    self._session_path = str(session_path)

  @property
  def session_path(self) -> str | None:
    """The currently configured session path, or None if unset."""
    return self._session_path

  # ------------------------------------------------------------------
  # Command dispatch
  # ------------------------------------------------------------------

  def do(
          self,
          command_name: str,
          args: tuple[Any, ...] = (),
          *,
          sync: bool = False,
  ) -> str | None:
    """Sends a PyMOL command to the worker process.

    Args:
        command_name: A valid PyMOL command, e.g. ``"zoom"`` or ``"select"``.
        args: Positional arguments appended to the command string.
        sync: When True, blocks until the worker sends back an
            acknowledgement and returns it as a string. When False
            (default), returns immediately with None.

    Returns:
        The worker's acknowledgement string when ``sync`` is True,
        otherwise None.

    Raises:
        RuntimeError: If no session path has been set, or the worker has
            not been started.
    """
    self._assert_ready()
    cmd = WorkerCommand.do(
      session_path=self._session_path,
      pymol_command=command_name,
      args=args,
      sync=sync,
    )
    return self._send(cmd)

  # ------------------------------------------------------------------
  # Convenience helpers
  # ------------------------------------------------------------------
  @staticmethod
  def cache_session(session_data: str, session_id: str) -> str:
    """Saves the PyMOL session in the base64 string format to the application cache directory.

    Args:
        session_data: The base64 string of the session.
        session_id: A unique identifier used as the filename stem.

    Returns:
        Absolute path of the saved ``.pse`` file.
    """
    tmp_cache_dir = pathlib.Path(constants.CACHE_PYMOL_SESSION_DIR)
    tmp_cache_dir.mkdir(parents=True, exist_ok=True)
    tmp_session_path = tmp_cache_dir / f"{session_id}.pse"
    binary_data.write_binary_file_from_base64_string(
      tmp_session_path, session_data
    )
    return str(tmp_session_path)

  @staticmethod
  def cache_user_session(user_pymol: UserPyMOL, session_id: str) -> str:
    """Saves the current PyMOL session to the application cache directory.

    Args:
        user_pymol: The live PyMOL session object to read from.
        session_id: A unique identifier used as the filename stem.

    Returns:
        Absolute path of the saved ``.pse`` file.
    """
    tmp_cache_dir = pathlib.Path(constants.CACHE_PYMOL_SESSION_DIR)
    tmp_cache_dir.mkdir(parents=True, exist_ok=True)
    tmp_session_path = tmp_cache_dir / f"{session_id}.pse"
    user_pymol.get_cmd_module().save(str(tmp_session_path))
    return str(tmp_session_path)

  @classmethod
  @contextmanager
  def session(
          cls, session_path: str | pathlib.Path
  ) -> Generator[PmlWorker, None, None]:
    """Context manager that starts a worker, binds a session, and stops it on exit.

    Args:
        session_path: Path to the ``.pse`` session file to bind.

    Yields:
        A started PmlWorker instance with the session path already set.

    Example::

        with PmlWorker.session("/path/to/file.pse") as worker:
            worker.do("zoom", sync=True)
    """
    worker = cls()
    worker.start()
    worker.set_session_path(session_path)
    try:
      yield worker
    finally:
      worker.stop()

  # ------------------------------------------------------------------
  # Internal helpers
  # ------------------------------------------------------------------

  def _assert_ready(self) -> None:
    if self._process is None or not self._process.is_alive():
      raise RuntimeError(
        "Worker is not running.  Call start() before do()."
      )
    if self._session_path is None:
      raise RuntimeError(
        "No session path set.  Call set_session_path() before do()."
      )

  def _send(self, cmd: WorkerCommand) -> str | None:
    """Sends a command over the pipe, optionally blocking for a reply.

    Args:
        cmd: The WorkerCommand to send.

    Returns:
        The worker's reply string if ``cmd.sync`` is True, otherwise None.
    """
    with self._lock:
      self._conn.send(cmd)
      if cmd.sync:
        return self._conn.recv()
    return None


# ---------------------------------------------------------------------------
# Module-level convenience functions
# ---------------------------------------------------------------------------


def one_shot_do(
        session_path: str | pathlib.Path,
        command_name: str,
        args: tuple[Any, ...] = (),
        *,
        sync: bool = False,
) -> str | None:
  """Starts a worker, executes a single command, then stops the worker.

  This is convenient for infrequent, isolated operations. For repeated
  commands against the same session prefer PmlWorker directly so the session
  is only loaded once.

  Args:
      session_path: Path to the ``.pse`` file.
      command_name: PyMOL command to execute.
      args: Arguments forwarded to the command.
      sync: Whether to wait for and return the worker's acknowledgement.

  Returns:
      The worker acknowledgement string when ``sync`` is True, otherwise None.
  """
  with PmlWorker.session(session_path) as worker:
    return worker.do(command_name, args=args, sync=sync)


# ---------------------------------------------------------------------------
# Worker-side event loop  (runs inside the subprocess)
# ---------------------------------------------------------------------------


def _worker_loop(conn: Connection) -> None:
  """Entry point executed inside the worker sub-process.

  Keeps one PyMOL instance alive for the lifetime of this process, caches
  the last loaded session path and only reloads when it changes (the key
  performance optimisation), dispatches ``do`` commands via ``cmd.do()``,
  replies on the pipe when ``sync=True``, and cleanly stops PyMOL on
  shutdown.

  Args:
      conn: The worker-side end of the multiprocessing Pipe used to receive
          commands from and send replies to the host process.
  """
  pymol_instance = pymol2.PyMOL()
  pymol_instance.start()

  loaded_session_path: str | None = None

  try:
    while True:
      try:
        cmd: WorkerCommand = conn.recv()
      except (EOFError, OSError):
        # The host closed its pipe end — this is the normal teardown
        # path when using the context manager. Exit cleanly.
        break

      if cmd.command_type is CommandType.SHUTDOWN:
        logger.debug("Worker: received shutdown command.")
        break

      if cmd.command_type is CommandType.LOAD_SESSION:
        _load_if_needed(pymol_instance, cmd.session_path, loaded_session_path)
        loaded_session_path = cmd.session_path
        if cmd.sync:
          conn.send(f"Session loaded: {cmd.session_path}")
        continue

      if cmd.command_type is CommandType.DO:
        # Load session only when it has changed – this is the key
        # performance optimisation for repeated commands on the same
        # session.
        if cmd.session_path != loaded_session_path:
          _load_if_needed(
            pymol_instance, cmd.session_path, loaded_session_path
          )
          loaded_session_path = cmd.session_path

        pymol_cmd_str = cmd.build_pymol_string()

        # TODO: Integrate this more smoothly
        match cmd.pymol_command.value:
          case "get_scene_list":
            result = pymol_instance.cmd.get_scene_list()
          case "get_model":
            result = pymol_instance.cmd.get_model(cmd.build_pymol_args())
          case _:
            pymol_instance.cmd.do(pymol_cmd_str)
            result = f"OK: {pymol_cmd_str}"

        logger.debug("Worker executed: %s", pymol_cmd_str)

        if cmd.sync:
          conn.send(result)
        continue

      logger.warning("Worker: unknown command type %s; skipping.", cmd.command_type)

  finally:
    pymol_instance.stop()
    conn.close()
    logger.debug("Worker: PyMOL stopped, pipe closed.")


def _load_if_needed(
        pymol_instance: pymol2.PyMOL,
        new_path: str | None,
        current_path: str | None,
) -> None:
  """Loads a session into PyMOL only when the path has changed.

  Args:
      pymol_instance: The active PyMOL instance to load the session into.
      new_path: Path to the ``.pse`` file to load, or None to skip.
      current_path: The path of the session currently loaded in
          ``pymol_instance``, or None if nothing has been loaded yet.
  """
  if new_path is None:
    return
  if new_path == current_path:
    return
  logger.debug("Worker: loading session %s", new_path)
  pymol_instance.cmd.load(new_path)
