#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module for all asynchronous functions that are related to the project object."""
import copy
import logging
import pathlib
import platform
import time
import shutil
from typing import Optional

from pyssa.controller import database_manager, watcher, interface_manager, pymol_session_manager
from pyssa.internal.data_structures import project, protein
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.thread import database_thread
from pyssa.internal.thread.async_pyssa import custom_signals
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def create_new_project(
    the_project_name: str,
    the_workspace_path: pathlib.Path,
    the_watcher: "watcher.Watcher",
    the_interface_manager: "interface_manager.InterfaceManager",
) -> tuple[
    str,
    Optional["project.Project"],
    Optional["watcher.Watcher"],
    Optional["interface_manager.InterfaceManager"],
]:
  """Creates a new project object.

  Args:
      the_project_name (str): The name of the new project.
      the_workspace_path (pathlib.Path): The path of the workspace where the project will be created.
      the_watcher (watcher.Watcher): An instance of the "watcher.Watcher" class used for setting up blacklists.
      the_interface_manager (interface_manager.InterfaceManager): The instance of the "interface_manager.InterfaceManager" class used in the application.

  Returns:
      A tuple containing:
      - A string denoting the result of creating the new project.
      - An optional instance of the "project.Project" class representing the newly created project.
      - An optional instance of the "watcher.Watcher" class.
      - An optional instance of the "interface_manager.InterfaceManager" class.
  """
  # <editor-fold desc="Checks">
  if the_project_name is None or the_project_name == "":
    logger.error("the_project_name is either None or an empty string.")
    return "", None, None, None
  if the_workspace_path is None:
    logger.error("the_workspace_path is None.")
    return "", None, None, None
  if the_watcher is None:
    logger.error("the_watcher is None.")
    return "", None, None, None
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None, None, None

  # </editor-fold>

  try:
    tmp_project = project.Project(the_project_name, the_workspace_path)

    tmp_database_filepath = str(
        pathlib.Path(f"{the_workspace_path}/{the_project_name}.db")
    )
    with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
      tmp_project.set_id(
          db_manager.insert_new_project(
              tmp_project.get_project_name(), platform.system()
          )
      )
    constants.PYSSA_LOGGER.info("Create empty project finished.")
    the_watcher.setup_blacklists(
        tmp_project,
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(
            enums.JobType.DISTANCE_ANALYSIS
        ),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
  except Exception as e:
    logger.error(e)
    return "", None, None, None
  else:
    return "result", tmp_project, the_watcher, the_interface_manager


def create_use_project(
    the_project_name: str,
    the_workspace_path: pathlib.Path,
    the_proteins_to_add: list,
    the_watcher: "watcher.Watcher",
    the_interface_manager: "interface_manager.InterfaceManager",
) -> tuple[
    str,
    Optional["project.Project"],
    Optional["watcher.Watcher"],
    Optional["interface_manager.InterfaceManager"],
]:
  """Creates a new project object.

  Args:
      the_project_name (str): The name of the project.
      the_workspace_path (pathlib.Path): The path to the workspace.
      the_proteins_to_add (list): The list of proteins to add.
      the_watcher (watcher.Watcher): The watcher object.
      the_interface_manager (interface_manager.InterfaceManager): The interface manager object.

  Returns:
      A tuple containing the result, the project object, the watcher object, and the interface manager object.
  """
  # <editor-fold desc="Checks">
  if the_project_name is None or the_project_name == "":
    logger.error("the_project_name is either None or an empty string.")
    return "", None, None, None
  if the_workspace_path is None:
    logger.error("the_workspace_path is None.")
    return "", None, None, None
  if the_proteins_to_add is None:
    logger.error("the_proteins_to_add is None.")
    return "", None, None, None
  if the_watcher is None:
    logger.error("the_watcher is None.")
    return "", None, None, None
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None, None, None

  # </editor-fold>

  try:
    tmp_project = project.Project(the_project_name, the_workspace_path)

    tmp_database_filepath = str(
        pathlib.Path(f"{the_workspace_path}/{the_project_name}.db")
    )
    with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
      tmp_project.set_id(
          db_manager.insert_new_project(
              tmp_project.get_project_name(), platform.system()
          )
      )
      for tmp_protein in the_proteins_to_add:
        tmp_protein_copy = copy.deepcopy(tmp_protein)
        tmp_protein_copy.db_project_id = tmp_project.get_id()
        tmp_project.add_existing_protein(tmp_protein_copy)
        tmp_protein_copy.set_id(db_manager.insert_new_protein(tmp_protein_copy))

    constants.PYSSA_LOGGER.info("Use project finished.")
    the_watcher.setup_blacklists(
        tmp_project,
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(
            enums.JobType.DISTANCE_ANALYSIS
        ),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
  except Exception as e:
    logger.error(e)
    return "", None, None, None
  else:
    return "result", tmp_project, the_watcher, the_interface_manager


def import_project(
    a_project_name: str,
    a_filepath: str,
    the_interface_manager: "interface_manager.InterfaceManager",
) -> tuple[str, Optional[str]]:
  """Imports a project from a file.

  Args:
      a_project_name (str): The name of the project being imported.
      a_filepath (str): The filepath of the project database file being imported.
      the_interface_manager (interface_manager.InterfaceManager): The instance of the InterfaceManager class.

  Returns:
      A tuple containing a result string and the temporary project database file path.
      The result string represents the success or failure of the import process.
      The temporary project database file path is the path of the imported project database file, or None if an error occurred.
  """
  # <editor-fold desc="Checks">
  if a_project_name is None or a_project_name == "":
    logger.error("a_project_name is either None or an empty string.")
    return "", None
  if a_filepath is None or a_filepath == "":
    logger.error("a_filepath is either None or an empty string.")
    return "", None
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None

  # </editor-fold>

  try:
    tmp_project_database_filepath = str(
        pathlib.Path(
            f"{the_interface_manager.get_application_settings().workspace_path}/{a_project_name}.db",
        ),
    )
    shutil.copyfile(a_filepath, tmp_project_database_filepath)
    with database_manager.DatabaseManager(
        tmp_project_database_filepath
    ) as db_manager:
      db_manager.update_project_name(a_project_name)
      tmp_project = db_manager.get_project_as_object(
          a_project_name,
          the_interface_manager.get_application_settings().workspace_path,
          the_interface_manager.get_application_settings(),
      )
    the_interface_manager.set_new_project(tmp_project)
    the_interface_manager.refresh_workspace_model()
    the_interface_manager.pymol_session_manager.reinitialize_session()
  except Exception as e:
    logger.error(e)
    return "", None
  return "result", tmp_project_database_filepath


def open_project(
    tmp_project_name: str,
    tmp_project_database_filepath: str,
    the_interface_manager: "interface_manager.InterfaceManager",
    the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
    the_custom_progress_signal: "custom_signals.ProgressSignal",
    the_watcher: "watcher.Watcher",
) -> tuple[
    str,
    Optional["project.Project"],
    Optional["interface_manager.InterfaceManager"],
    Optional["watcher.Watcher"],
]:
  """Opens a project.

  Args:
      tmp_project_name (str): The name of the temporary project.
      tmp_project_database_filepath (str): The filepath of the temporary project database.
      the_interface_manager (interface_manager.InterfaceManager): The interface manager object.
      the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PyMOL session manager object.
      the_custom_progress_signal (custom_signals.ProgressSignal): The custom progress signal object.
      the_watcher (watcher.Watcher): The watcher object.

  Returns:
      A tuple containing the result as a string, the temporary project object, the interface manager object, and the watcher object.
      If an exception occurs during the execution of the method, an empty string and None values will be returned.
  """
  # <editor-fold desc="Checks">
  if tmp_project_name is None or tmp_project_name == "":
    logger.error("tmp_project_name is either None or an empty string.")
    return "", None, None, None
  if (
      tmp_project_database_filepath is None
      or tmp_project_database_filepath == ""
  ):
    logger.error(
        "tmp_project_database_filepath is either None or an empty string."
    )
    return "", None, None, None
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None, None, None
  if the_pymol_session_manager is None:
    logger.error("the_pymol_session_manager is None.")
    return "", None, None, None
  if the_custom_progress_signal is None:
    logger.error("the_custom_progress_signal is None.")
    return "", None, None, None
  if the_watcher is None:
    logger.error("the_watcher is None.")
    return "", None, None, None

  # </editor-fold>

  try:
    the_custom_progress_signal.emit_signal("Opening database ...", 10)
    with database_manager.DatabaseManager(
        tmp_project_database_filepath
    ) as db_manager:
      the_custom_progress_signal.emit_signal("Setting up project ...", 30)
      tmp_project = db_manager.get_project_as_object(
          tmp_project_name,
          the_interface_manager.get_application_settings().workspace_path,
          the_interface_manager.get_application_settings(),
          the_custom_progress_signal,
      )
    the_interface_manager.set_new_project(tmp_project)
    the_custom_progress_signal.emit_signal(
        "Reinitializing PyMOL session ...", 96
    )
    the_pymol_session_manager.reinitialize_session()
    the_watcher.setup_blacklists(
        tmp_project,
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(
            enums.JobType.DISTANCE_ANALYSIS
        ),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
  except Exception as e:
    logger.error(e)
    return "", None, None, None
  else:
    return "result", tmp_project, the_interface_manager, the_watcher


def close_project(
    the_database_thread: "database_thread.DatabaseThread",
    the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
) -> tuple[str, str]:
  """Closes a project by putting a close project operation into the database thread's queue and waiting for it to finish.

  Args:
      the_database_thread (database_thread.DatabaseThread): The database thread which will execute the close project operation.
      the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The Pymol session manager used to reinitialize the session.

  Returns:
      A tuple containing the result and status message of the operation.

  Raises:
      exception.IllegalArgumentError: If any of the arguments are None.
  """
  # <editor-fold desc="Checks">
  if the_database_thread is None:
    logger.error("the_database_thread is None.")
    raise exception.IllegalArgumentError("the_database_thread is None.")
  if the_pymol_session_manager is None:
    logger.error("the_pymol_session_manager is None.")
    raise exception.IllegalArgumentError("the_pymol_session_manager is None.")

  # </editor-fold>

  try:
    the_database_thread.put_database_operation_into_queue(
        database_operation.DatabaseOperation(
            enums.SQLQueryType.CLOSE_PROJECT, (0, "")
        ),
    )
    print(the_database_thread.queue_is_running())
    while the_database_thread.queue_is_running():
      time.sleep(0.5)
    the_pymol_session_manager.reinitialize_session()
  except Exception as e:
    logger.error(
        f"Unknown error occurred while waiting for the database thread to finish: {e}."
    )
    return "", "Waiting for database thread queue failed!"
  else:
    logger.info("Waiting for database thread queue finished.")
    return "result", "Waiting for database thread queue finished."


def add_protein_from_pdb_to_project(
    tmp_protein_name: str,
    the_interface_manager: "interface_manager.InterfaceManager",
) -> tuple[str, Optional["protein.Protein"]]:
  """Adds a protein from a pdb file to the protein model.

  Args:
      tmp_protein_name (str): The name of the protein to be added from the PDB.
      the_interface_manager (interface_manager.InterfaceManager): An instance of the "interface_manager.InterfaceManager" class that manages the user interface.

  Returns:
      A tuple of a string and an optional instance of the "protein.Protein" class. The string represents the result of the method call,
      and the optional protein object represents the added protein if successful, or None if not successful.
  """
  # <editor-fold desc="Checks">
  if tmp_protein_name is None or tmp_protein_name == "":
    logger.error("tmp_protein_name is either None or an empty string.")
    return "", None
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None

  # </editor-fold>

  try:
    the_main_socket, the_general_purpose_socket = (
        the_interface_manager.job_manager.get_general_purpose_socket_pair()
    )
    with database_manager.DatabaseManager(
        str(the_interface_manager.get_current_project().get_database_filepath())
    ) as db_manager:
      tmp_ref_protein = protein.Protein(tmp_protein_name.upper())
      tmp_ref_protein.set_id(db_manager.get_next_id_of_protein_table())
      tmp_ref_protein.db_project_id = (
          the_interface_manager.get_current_project().get_id()
      )
      tmp_ref_protein.add_protein_structure_data_from_pdb_db(
          tmp_protein_name.upper(), the_main_socket, the_general_purpose_socket
      )
      tmp_ref_protein.add_id_to_all_chains(
          db_manager.get_next_id_of_chain_table()
      )
    tmp_ref_protein.create_new_pymol_session(
        the_main_socket, the_general_purpose_socket
    )
    the_interface_manager.add_protein_to_proteins_model(tmp_ref_protein)
  except Exception as e:
    logger.error(e)
    return "", None
  else:
    return "result", tmp_ref_protein


def add_protein_from_local_filesystem_to_project(
    tmp_protein_name: str,
    the_interface_manager: "interface_manager.InterfaceManager",
) -> tuple[str, Optional["protein.Protein"]]:
  """Adds a protein from the local filesystem to the protein model.

  Args:
      tmp_protein_name (str): The name of the temporary protein file.
      the_interface_manager (interface_manager.InterfaceManager): The interface manager object.

  Returns:
      A tuple containing the result message and the protein object if it was successfully added, or None if there was an error.
  """
  # <editor-fold desc="Checks">
  if tmp_protein_name is None or tmp_protein_name == "":
    logger.error("tmp_protein_name is either None or an empty string.")
    return "", None
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None

  # </editor-fold>

  try:
    pdb_filepath = pathlib.Path(tmp_protein_name)
    tmp_ref_protein = protein.Protein(
        pdb_filepath.name.replace(".pdb", ""),
    )
    with database_manager.DatabaseManager(
        str(the_interface_manager.get_current_project().get_database_filepath())
    ) as db_manager:
      tmp_ref_protein.set_id(db_manager.get_next_id_of_protein_table())
      tmp_ref_protein.db_project_id = (
          the_interface_manager.get_current_project().get_id()
      )
      the_main_socket, the_general_purpose_socket = (
          the_interface_manager.job_manager.get_general_purpose_socket_pair()
      )
      tmp_ref_protein.add_protein_structure_data_from_local_pdb_file(
          pdb_filepath, the_main_socket, the_general_purpose_socket
      )
      tmp_ref_protein.add_id_to_all_chains(
          db_manager.get_next_id_of_chain_table()
      )
    tmp_ref_protein.create_new_pymol_session(
        the_main_socket, the_general_purpose_socket
    )
    the_interface_manager.add_protein_to_proteins_model(tmp_ref_protein)
  except Exception as e:
    logger.error(e)
    return "", None
  else:
    return "result", tmp_ref_protein
