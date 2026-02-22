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
import os
import numpy as np
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

  @classmethod
  def one_shot_do(
      cls,
      command_name: Any,
      args: tuple[Any, ...] = (),
      *,
      sync: bool = True,
  ) -> Any:
      """Starts a worker, executes a single command, then stops the worker.
      
      Args:
          command_name: The command or macro to execute.
          args: Positional arguments for the command.
          sync: Whether to block until completion.
          
      Returns:
          The worker's reply.
      """
      worker = cls()
      worker.start()
      try:
          cmd = WorkerCommand.do(
              session_path=None,
              pymol_command=command_name,
              args=args,
              sync=sync,
          )
          return worker._send(cmd)
      finally:
          worker.stop()


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
          case "consolidate_molecule":
            pymol_instance.cmd.reinitialize()
            loaded_session_path = None
            result = _consolidate_molecule(pymol_instance, *cmd.args)
          case "clean_protein_update_structure":
            # Pass args explicitly
            result = _clean_protein(pymol_instance, *cmd.args)
            loaded_session_path = None
          case "create_new_session":
            pymol_instance.cmd.reinitialize()
            loaded_session_path = None
            result = _create_new_session(pymol_instance, *cmd.args)
          case "get_chains":
            pymol_instance.cmd.reinitialize()
            loaded_session_path = None
            result = _get_chains(pymol_instance, *cmd.args)
          case "distance_analysis":
            pymol_instance.cmd.reinitialize()
            loaded_session_path = None
            result = _do_distance_analysis(pymol_instance, *cmd.args)
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


def _consolidate_molecule(pymol_instance: pymol2.PyMOL, a_pdb_filepath: str) -> str:
  if not a_pdb_filepath or not os.path.exists(a_pdb_filepath):
    return ""
  try:
    pymol_instance.cmd.load(filename=str(a_pdb_filepath))
    tmp_protein_name = pathlib.Path(a_pdb_filepath).name.replace(".pdb", "").replace(" ", "_")
  except Exception as e:
    logger.error(f"Failed to load PDB: {e}")
    return ""

  tmp_states = pymol_instance.cmd.count_states(tmp_protein_name)
  if tmp_states > 1:
    try:
      pymol_instance.cmd.create("new_object", tmp_protein_name, 1, 1)
      pymol_instance.cmd.delete(tmp_protein_name)
      pymol_instance.cmd.set_name("new_object", tmp_protein_name)
      tmp_pdb_cache_filepath = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein_name}.pdb")
      pymol_instance.cmd.save(str(tmp_pdb_cache_filepath))
      return str(tmp_pdb_cache_filepath)
    except Exception as e:
      logger.error(f"Protein states consolidation failed: {e}")
      return ""
  return a_pdb_filepath


def _clean_protein(pymol_instance: pymol2.PyMOL, a_pymol_session: str, a_protein_name: str) -> tuple[str, str]:
  if not a_pymol_session or not a_protein_name:
    return "", ""
  tmp_session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/temp_session.pse")
  binary_data.write_binary_file_from_base64_string(tmp_session_filepath, a_pymol_session)

  try:
    pymol_instance.cmd.load(str(tmp_session_filepath))
    tmp_all_object_names = pymol_instance.cmd.get_names()
    if len(tmp_all_object_names) == 0 or a_protein_name not in tmp_all_object_names:
      return "", ""
    pymol_instance.cmd.remove("solvent")
    pymol_instance.cmd.remove("organic")
    tmp_export_session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/export_temp_session.pse")
    tmp_export_pdb_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/export_temp_clean.pdb")
    pymol_instance.cmd.save(str(tmp_export_session_filepath))
    pymol_instance.cmd.save(str(tmp_export_pdb_filepath))
    base64_string = binary_data.create_base64_string_from_file(tmp_export_session_filepath)
  except Exception as e:
    logger.error(f"Clean protein failed: {e}")
    return "", ""
  finally:
    if os.path.exists(str(tmp_session_filepath)):
      os.remove(str(tmp_session_filepath))
  return base64_string, str(tmp_export_pdb_filepath)


def _create_new_session(pymol_instance: pymol2.PyMOL, a_pdb_filepath: str) -> str:
  if not a_pdb_filepath or not os.path.exists(a_pdb_filepath):
    return ""
  try:
    tmp_protein_name = pathlib.Path(a_pdb_filepath).name.replace(".pdb", "")
    pymol_instance.cmd.load(filename=str(a_pdb_filepath), object=tmp_protein_name)
    pymol_instance.cmd.bg_color("black")
    
    pymol_instance.cmd.set("valence", 0)
    pymol_instance.cmd.set("scene_buttons", 0)
    pymol_instance.cmd.set("ray_trace_mode", 1)
    pymol_instance.cmd.set("antialias", 1)
    pymol_instance.cmd.set("ambient", 0.5)
    pymol_instance.cmd.set("cartoon_fancy_helices", 1)
    pymol_instance.cmd.set("cartoon_discrete_colors", 1)
    pymol_instance.cmd.set("cartoon_sampling", 14)
    pymol_instance.cmd.set("spec_power", 350)
    pymol_instance.cmd.set("spec_reflect", 0.2)
    pymol_instance.cmd.set("ray_transparency_contrast", 0.4)
    pymol_instance.cmd.set("ray_transparency_oblique", 1.0)
    pymol_instance.cmd.set("ray_transparency_oblique_power", 4.0)
    pymol_instance.cmd.set("ray_trace_color", "black")
    pymol_instance.cmd.unset("depth_cue")
    pymol_instance.cmd.color("green", tmp_protein_name)
    pymol_instance.cmd.reset()
    pymol_instance.cmd.scene("base", action="store")
    
    session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/{tmp_protein_name}_session.pse")
    pymol_instance.cmd.save(str(session_filepath))
    base64_string = binary_data.create_base64_string_from_file(session_filepath)
    if os.path.exists(str(session_filepath)):
      os.remove(str(session_filepath))
    return base64_string
  except Exception as e:
    logger.error(f"Create new session failed: {e}")
    return ""


def _get_chains(pymol_instance: pymol2.PyMOL, a_pdb_filepath: str) -> list[tuple]:
  if not a_pdb_filepath or not os.path.exists(a_pdb_filepath):
    return []
  try:
    pymol_instance.cmd.load(filename=str(a_pdb_filepath))
    tmp_protein_name = pathlib.Path(a_pdb_filepath).name.replace(".pdb", "")
    tmp_chains = pymol_instance.cmd.get_chains()
    chains_of_protein = []
    
    protein_resns = {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'}
    
    for tmp_chain in tmp_chains:
      sequence_of_chain = pymol_instance.cmd.get_model(f"chain {tmp_chain}")
      if not sequence_of_chain or not sequence_of_chain.atom:
        continue
      is_protein_chain = sequence_of_chain.atom[0].resn in protein_resns
      
      fasta_sequence_of_chain = pymol_instance.cmd.get_fastastr(f"chain {tmp_chain}")
      fasta_sequence_of_chain_without_header = fasta_sequence_of_chain[fasta_sequence_of_chain.find("\\n"):]
      complete_sequence_of_chain = (tmp_protein_name, fasta_sequence_of_chain_without_header.replace("\\n", ""))
      
      chains_of_protein.append((tmp_chain, complete_sequence_of_chain, "PROTEIN" if is_protein_chain else "NON_PROTEIN"))
    return chains_of_protein
  except Exception as e:
    logger.error(f"Get chains failed: {e}")
    return []


def _do_distance_analysis(
  pymol_instance: pymol2.PyMOL,
  the_protein_pair_name: str,
  a_protein_1_pdb_cache_filepath: str,
  a_protein_2_pdb_cache_filepath: str,
  a_protein_1_pymol_selection_string: str,
  a_protein_2_pymol_selection_string: str,
  a_cutoff: float,
  the_cycles: int,
) -> tuple[tuple, str]:
  tmp_protein_1_name = pathlib.Path(a_protein_1_pdb_cache_filepath).name.replace(".pdb", "")
  tmp_protein_2_name = pathlib.Path(a_protein_2_pdb_cache_filepath).name.replace(".pdb", "")

  pymol_instance.cmd.load(filename=str(a_protein_1_pdb_cache_filepath), object=tmp_protein_1_name)
  pymol_instance.cmd.load(filename=str(a_protein_2_pdb_cache_filepath), object=tmp_protein_2_name)

  pymol_instance.cmd.color("green", tmp_protein_1_name)
  pymol_instance.cmd.color("cyan", tmp_protein_2_name)

  pymol_instance.cmd.zoom("all")
  pymol_instance.cmd.scene("base", action="store")

  pymol_instance.cmd.bg_color("white")
  pymol_instance.cmd.set("valence", 0)
  pymol_instance.cmd.set("scene_buttons", 0)
  pymol_instance.cmd.set("ray_trace_mode", 1)
  pymol_instance.cmd.set("antialias", 1)
  pymol_instance.cmd.set("ambient", 0.5)
  pymol_instance.cmd.set("cartoon_fancy_helices", 1)
  pymol_instance.cmd.set("cartoon_discrete_colors", 1)
  pymol_instance.cmd.set("cartoon_sampling", 14)
  pymol_instance.cmd.set("spec_power", 350)
  pymol_instance.cmd.set("spec_reflect", 0.2)
  pymol_instance.cmd.set("ray_transparency_contrast", 0.4)
  pymol_instance.cmd.set("ray_transparency_oblique", 1.0)
  pymol_instance.cmd.set("ray_transparency_oblique_power", 4.0)

  seq_len_1_atoms = pymol_instance.cmd.get_model(f"{tmp_protein_1_name} and n. CA").atom
  seq_len_protein_1 = len(seq_len_1_atoms) if seq_len_1_atoms else 0
  seq_len_2_atoms = pymol_instance.cmd.get_model(f"{tmp_protein_2_name} and n. CA").atom
  seq_len_protein_2 = len(seq_len_2_atoms) if seq_len_2_atoms else 0

  if seq_len_protein_1 > seq_len_protein_2:
    tmp_align_results = pymol_instance.cmd.align(
      f"{tmp_protein_2_name} and {a_protein_2_pymol_selection_string}",
      f"{tmp_protein_1_name} and {a_protein_1_pymol_selection_string}",
      cycles=the_cycles,
      cutoff=a_cutoff,
      object="aln",
    )
  else:
    tmp_align_results = pymol_instance.cmd.align(
      f"{tmp_protein_1_name} and {a_protein_1_pymol_selection_string}",
      f"{tmp_protein_2_name} and {a_protein_2_pymol_selection_string}",
      cycles=the_cycles,
      cutoff=a_cutoff,
      object="aln",
    )
      
  rmsd_dict = {
    "rmsd": str(round(tmp_align_results[0], 2)),
    "aligned_residues": f"{str(tmp_align_results[1])} / {seq_len_protein_1}",
  }

  index_list, ref_chain_list, ref_pos_list, ref_resi_list = [], [], [], []
  model_chain_list, model_pos_list, model_resi_list, distance_list = [], [], [], []

  idx2resi = []
  pymol_instance.cmd.iterate("aln", "idx2resi.append((model, chain, resi, resn))", space={"idx2resi": idx2resi})
  
  prot_1_indices = []
  prot_2_indices = []
  for tmp_prot_atom in idx2resi:
    if tmp_prot_atom[0] == tmp_protein_1_name:
      try: tmp_residue_number = int(tmp_prot_atom[2])
      except ValueError: tmp_residue_number = int(tmp_prot_atom[2][:-1])
      prot_1_indices.append((tmp_prot_atom[1], tmp_residue_number, tmp_prot_atom[3]))
    if tmp_prot_atom[0] == tmp_protein_2_name:
      try: tmp_residue_number = int(tmp_prot_atom[2])
      except ValueError: tmp_residue_number = int(tmp_prot_atom[2][:-1])
      prot_2_indices.append((tmp_prot_atom[1], tmp_residue_number, tmp_prot_atom[3]))

  for resi_no in range(len(prot_1_indices)):
    atom1 = f"/{tmp_protein_1_name}//{prot_1_indices[resi_no][0]}/{prot_1_indices[resi_no][1]}/CA"
    atom2 = f"/{tmp_protein_2_name}//{prot_2_indices[resi_no][0]}/{prot_2_indices[resi_no][1]}/CA"
    distance = round(pymol_instance.cmd.get_distance(atom1, atom2, state=1), 2)

    ref_chain_list.append(prot_1_indices[resi_no][0])
    ref_pos_list.append(int(prot_1_indices[resi_no][1]))
    ref_resi_list.append(prot_1_indices[resi_no][2])
    model_chain_list.append(prot_2_indices[resi_no][0])
    model_pos_list.append(int(prot_2_indices[resi_no][1]))
    model_resi_list.append(prot_2_indices[resi_no][2])
    distance_list.append(distance)
    index_list.append(resi_no)

  pymol_instance.cmd.show("cartoon")
  pymol_instance.cmd.hide("cgo", "all")
  pymol_instance.cmd.reset()
  pymol_instance.cmd.zoom("all")
  pymol_instance.cmd.scene(key=f"{tmp_protein_1_name}-{tmp_protein_2_name}", action="store")

  j, i = 0, 0
  for distance_value in distance_list:
    if float(distance_value) > float(a_cutoff):
      ref_pos = ref_pos_list[i]
      ref_chain = ref_chain_list[i]
      model_pos = model_pos_list[i]
      model_chain = model_chain_list[i]

      measurement_obj = f"measure{j}"
      atom1 = f"/{tmp_protein_1_name}//{ref_chain}/{ref_pos}/CA"
      atom2 = f"/{tmp_protein_2_name}//{model_chain}/{model_pos}/CA"
        
      pymol_instance.cmd.zoom(f"/{tmp_protein_1_name}//{ref_chain}/{ref_pos}", 10)
      pymol_instance.cmd.distance(measurement_obj, atom1, atom2)
      pymol_instance.cmd.label(atom1, r"'%s-%s' % (resn, resi)")
      pymol_instance.cmd.label(atom2, r"'%s-%s' % (resn, resi)")
      pymol_instance.cmd.set("label_position", (0, 0, 10))
      pymol_instance.cmd.scene(key=f"{ref_pos}-{model_pos}", action="store")
      pymol_instance.cmd.hide("labels", atom1)
      pymol_instance.cmd.hide("labels", atom2)
      pymol_instance.cmd.hide("labels", measurement_obj)
      pymol_instance.cmd.hide("dashes", measurement_obj)
      j += 1
    i += 1

  session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/{the_protein_pair_name}_session.pse")
  pymol_instance.cmd.save(str(session_filepath))
  base64_string = binary_data.create_base64_string_from_file(session_filepath)
  if os.path.exists(str(session_filepath)):
    os.remove(str(session_filepath))

  result_lists_hashtable = {
    "index": index_list,
    "ref_chain": ref_chain_list,
    "ref_pos": ref_pos_list,
    "ref_resi": ref_resi_list,
    "model_chain": model_chain_list,
    "model_pos": model_pos_list,
    "model_resi": model_resi_list,
    "distance": distance_list,
  }
  
  distance_analysis_results_object_values = (
    result_lists_hashtable,
    base64_string,
    float(rmsd_dict["rmsd"]),
    rmsd_dict["aligned_residues"],
  )
  return distance_analysis_results_object_values, base64_string
