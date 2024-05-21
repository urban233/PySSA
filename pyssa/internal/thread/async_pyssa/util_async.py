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
"""Module for util functions that are used by async functions."""
import logging
import os
import pathlib
import shutil
import subprocess
from io import BytesIO
from typing import Optional

import requests
import zmq
import pygetwindow

from pyssa.controller import database_manager, interface_manager, pymol_session_manager
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.thread import database_thread
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, enums, exception
from pyssa.util import globals

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def open_documentation_on_certain_page(
    a_page_name: str,
    the_docs_window,
) -> tuple:
  """Opens the pyssa documentation on the given page name.

  Args:
      a_page_name (str): The name of the page to open in the documentation.
      the_docs_window: The window of the documentation.

  Returns:
      A tuple containing the result, response, and temporary documentation window.
  """
  # <editor-fold desc="Checks">
  if a_page_name is None or a_page_name == "":
    logger.error("a_page_name is either None or an empty string.")
    return "", None, None
  if the_docs_window is None:
    logger.error("the_docs_window is None.")
    return "", None, None

  # </editor-fold>

  try:
    if the_docs_window is not None:
      tmp_is_docs_open_flag = True
    else:
      tmp_is_docs_open_flag = False

    tmp_docs_window = None
    if tmp_is_docs_open_flag:
      # Docs window is already open
      logger.info("Update currently open documentation.")
      tmp_docs_window = the_docs_window
    elif globals.g_server_status == enums.DocsServerStatus.PENDING:
      # Docs are getting build
      logger.warning(
          "Docs need to be build. The build process was started from another process!"
      )
      flag = False
      while flag is not True:
        if (
            len(
                pygetwindow.getWindowsWithTitle(
                    constants.WINDOW_TITLE_OF_HELP_CENTER
                )
            )
            == 1
            or globals.g_server_status == enums.DocsServerStatus.ACTIVE
        ):
          flag = True
          globals.g_server_status = enums.DocsServerStatus.ACTIVE
          tmp_docs_window = pygetwindow.getWindowsWithTitle(
              constants.WINDOW_TITLE_OF_HELP_CENTER
          )[0]
          tmp_docs_window.minimize()
    else:
      # Docs are need to be built in this process
      logger.info("Trying to run the mkdocs serve command ...")
      try:
        os.chdir(constants.PLUGIN_DOCS_PATH)
        subprocess.Popen(
            [
                r"C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Scripts\mkdocs.exe",
                "serve",
                "-a",
                "127.0.0.1:7018",
            ],
            creationflags=subprocess.CREATE_NO_WINDOW,
        )
        subprocess.Popen(
            [r"C:\ProgramData\pyssa\extra_tools\browser.exe"],
            creationflags=subprocess.CREATE_NO_WINDOW,
        )
      except subprocess.CalledProcessError as e:
        logger.error(f"Error starting mkdocs serve: {e}")
        tmp_docs_window = None
      else:
        logger.info("Running the mkdocs serve command finished without errors.")
        globals.g_server_status = enums.DocsServerStatus.PENDING
        flag = False
        while flag is not True:
          if (
              len(
                  pygetwindow.getWindowsWithTitle(
                      constants.WINDOW_TITLE_OF_HELP_CENTER
                  )
              )
              == 1
          ):
            flag = True
            globals.g_server_status = enums.DocsServerStatus.ACTIVE
            tmp_docs_window = pygetwindow.getWindowsWithTitle(
                constants.WINDOW_TITLE_OF_HELP_CENTER
            )[0]
            tmp_docs_window.minimize()

    # After the docs window is ready, it can be accessed by ZeroMQ
    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect("tcp://127.0.0.1:7017")
    socket.send_string(f"update_url http://127.0.0.1:7018/{a_page_name}")
    response = socket.recv_string()
    logger.debug(f"Response from server: {response}")
    socket.close()
  except Exception as e:
    logger.error(e)
    return "", None, None
  else:
    return "result", response, tmp_docs_window


def start_documentation_server(
    placeholder_1: int,
    placeholder_2: int,
) -> tuple:
  """Opens the pyssa documentation on the given page name.

  Args:
      placeholder_1: just a placeholder to be able to use LegacyTask.
      placeholder_2: just a placeholder to be able to use LegacyTask.

  Returns:
      A tuple with ("result", response or None).
  """
  try:
    logger.info("Trying to run the mkdocs serve command ...")
    try:
      os.chdir(constants.PLUGIN_DOCS_PATH)
      subprocess.Popen(
          [
              r"C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Scripts\mkdocs.exe",
              "serve",
              "-a",
              "127.0.0.1:7018",
          ],
          creationflags=subprocess.CREATE_NO_WINDOW,
      )
      subprocess.Popen(
          [r"C:\ProgramData\pyssa\extra_tools\browser.exe"],
          creationflags=subprocess.CREATE_NO_WINDOW,
      )
    except subprocess.CalledProcessError as e:
      logger.error(f"Error starting mkdocs serve: {e}")
    else:
      logger.info("Running the mkdocs serve command finished without errors.")
      globals.g_server_status = enums.DocsServerStatus.PENDING

    flag = False
    while flag is not True:
      if (
          len(
              pygetwindow.getWindowsWithTitle(
                  constants.WINDOW_TITLE_OF_HELP_CENTER
              )
          )
          == 1
      ):
        flag = True
        globals.g_server_status = enums.DocsServerStatus.ACTIVE
    tmp_docs_window = pygetwindow.getWindowsWithTitle(
        constants.WINDOW_TITLE_OF_HELP_CENTER
    )[0]
    tmp_docs_window.minimize()

    context = zmq.Context()
    socket = context.socket(zmq.REQ)
    socket.connect("tcp://127.0.0.1:7017")
    socket.send_string("update_url http://127.0.0.1:7018/help/")
    response = socket.recv_string()
    logger.debug(f"Response from server: {response}")
    socket.close()
  except Exception as e:
    logger.error(e)
    return "", None
  else:
    return "result", tmp_docs_window


def add_proteins_to_project_and_model(
    the_interface_manager: "interface_manager.InterfaceManager",
    tmp_protein_names: list,
) -> tuple[str, Optional["interface_manager.InterfaceManager"]]:
  """Adds proteins to the current project and model.

  Args:
      the_interface_manager (interface_manager.InterfaceManager): The interface manager object.
      tmp_protein_names (list): A list of temporary protein names to add.

  Returns:
      A tuple containing the result status and the updated interface manager object, or None if an exception occurred.
  """
  # <editor-fold desc="Checks">
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None
  if tmp_protein_names is None:
    logger.error("tmp_protein_names is None.")
    return "", None

  # </editor-fold>

  try:
    the_interface_manager.project_lock.lock()
    with database_manager.DatabaseManager(
        str(the_interface_manager.get_current_project().get_database_filepath())
    ) as db_manager:
      for tmp_protein_name in tmp_protein_names:
        tmp_protein = db_manager.get_protein_as_object(tmp_protein_name)
        the_interface_manager.add_protein_to_current_project(tmp_protein)
        the_interface_manager.add_protein_to_proteins_model(tmp_protein)
    the_interface_manager.project_lock.unlock()
    the_interface_manager.watcher.setup_blacklists(
        the_interface_manager.get_current_project(),
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(
            enums.JobType.DISTANCE_ANALYSIS
        ),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
  except Exception as e:
    logger.error(e)
    the_interface_manager.project_lock.unlock()
    return "", None
  else:
    return "result", the_interface_manager


def add_protein_pairs_to_project_and_model(
    the_interface_manager: "interface_manager.InterfaceManager",
    tmp_protein_pair_names: list,
) -> tuple[str, Optional["interface_manager.InterfaceManager"]]:
  """Adds protein pairs to the current project and model.

  Args:
      the_interface_manager (interface_manager.InterfaceManager): The instance of InterfaceManager.
      tmp_protein_pair_names (list): The list of protein pair names.

  Returns:
      A tuple containing the result message and the updated instance of InterfaceManager.
  """
  # <editor-fold desc="Checks">
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None
  if tmp_protein_pair_names is None:
    logger.error("tmp_protein_pair_names is None.")
    return "", None

  # </editor-fold>

  try:
    the_interface_manager.project_lock.lock()
    with database_manager.DatabaseManager(
        str(the_interface_manager.get_current_project().get_database_filepath())
    ) as db_manager:
      for tmp_protein_pair_name in tmp_protein_pair_names:
        print(tmp_protein_pair_name)
        tmp_protein_pair = db_manager.get_protein_pair_as_object(
            tmp_protein_pair_name,
            the_interface_manager.get_current_project(),
            the_interface_manager.get_settings_manager().settings,
        )
        print(tmp_protein_pair)
        the_interface_manager.add_protein_pair_to_current_project(
            tmp_protein_pair
        )
        the_interface_manager.add_protein_pair_to_protein_pairs_model(
            tmp_protein_pair
        )
    the_interface_manager.project_lock.unlock()
    the_interface_manager.watcher.setup_blacklists(
        the_interface_manager.get_current_project(),
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(
            enums.JobType.DISTANCE_ANALYSIS
        ),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
  except Exception as e:
    logger.error(e)
    the_interface_manager.project_lock.unlock()
    return "", None
  else:
    return "result", the_interface_manager


def add_proteins_and_protein_pairs_to_project_and_model(
    the_interface_manager: "interface_manager.InterfaceManager",
    tmp_protein_names: list,
    tmp_protein_pair_names: list,
) -> tuple[str, Optional["interface_manager.InterfaceManager"]]:
  """Adds proteins and protein pairs to a project and model.

  Args:
      the_interface_manager (interface_manager.InterfaceManager): The instance of InterfaceManager.
      tmp_protein_names (list): A list of temporary protein names to add.
      tmp_protein_pair_names (list): The list of protein pair names.

  Returns:
      A tuple containing a result string and an optional instance of the interface_manager.InterfaceManager class.
      The result string indicates the outcome of the method execution.
      The optional InterfaceManager instance is returned only if the method execution is successful.
  """
  # <editor-fold desc="Checks">
  if the_interface_manager is None:
    logger.error("the_interface_manager is None.")
    return "", None
  if tmp_protein_names is None:
    logger.error("tmp_protein_names is None.")
    return "", None
  if tmp_protein_pair_names is None:
    logger.error("tmp_protein_pair_names is None.")
    return "", None

  # </editor-fold>

  try:
    add_proteins_to_project_and_model(the_interface_manager, tmp_protein_names)
    add_protein_pairs_to_project_and_model(
        the_interface_manager, tmp_protein_pair_names
    )
  except Exception as e:
    logger.error(e)
    return "", None
  return "result", the_interface_manager


def close_project_automatically(
    a_project_is_open: bool,
    the_database_thread: "database_thread.DatabaseThread",
    the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
    a_project_name: str,
) -> tuple[str, int]:
  """Closes a project automatically.

  Args:
      a_project_is_open (bool): Whether a project is currently open.
      the_database_thread (database_thread.DatabaseThread): The database thread responsible for executing database operations.
      the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PyMOL session manager.
      a_project_name (str): The name of the project.

  Returns:
      A tuple containing the project name and a status code indicating whether the project was closed successfully (0) or encountered an error (non-zero).
  """
  # <editor-fold desc="Checks">
  if a_project_is_open is None:
    logger.error("a_project_is_open is None.")
    return a_project_name, 1
  if the_database_thread is None:
    logger.error("the_database_thread is None.")
    return a_project_name, 1
  if the_pymol_session_manager is None:
    logger.error("the_pymol_session_manager is None.")
    return a_project_name, 1
  if a_project_name is None or a_project_name == "":
    logger.error("a_project_name is either None or an empty string.")
    return "", 1

  # </editor-fold>

  if a_project_is_open:
    try:
      the_database_thread.put_database_operation_into_queue(
          database_operation.DatabaseOperation(
              enums.SQLQueryType.CLOSE_PROJECT, (0, "")
          ),
      )
      the_pymol_session_manager.reinitialize_session()
    except Exception as e:
      logger.error(
          f"Unknown error occurred while waiting for the database thread to finish: {e}."
      )
      return a_project_name, 0
    else:
      logger.info("Waiting for database thread queue finished.")
      return a_project_name, 0
  else:
    return a_project_name, 0


def download_demo_projects(
    the_workspace_path: str, a_placeholder: int
) -> tuple[bool]:
  """Downloads demo projects to the current workspace.

  Args:
      the_workspace_path: The path of the workspace where the demo projects will be downloaded and extracted.
      a_placeholder: A placeholder parameter that serves as an example.

  Returns:
      A tuple containing a boolean value indicating the success of the download process.
  """
  if the_workspace_path is None or the_workspace_path == "":
    logger.error("the_workspace_path is either None or an empty string.")
    return (False,)
  if not os.path.exists(the_workspace_path):
    return (False,)
  try:
    import zipfile

    download_dest = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip")
    if os.path.exists(download_dest):
      os.remove(download_dest)

    # download demo projects
    url = "https://w-hs.sciebo.de/s/ZHJa6XB9SKWtqGi/download"
    try:
      response = requests.get(url)
      response.raise_for_status()  # Check for errors
      zipfile = zipfile.ZipFile(BytesIO(response.content))
      zipfile.extractall(
          pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
      )
    except requests.exceptions.HTTPError as errh:
      constants.PYSSA_LOGGER.error(f"HTTP Error: {errh}")
      return (False,)
    except requests.exceptions.ConnectionError as errc:
      constants.PYSSA_LOGGER.error(f"Error Connecting: {errc}")
      return (False,)
    except requests.exceptions.Timeout as errt:
      constants.PYSSA_LOGGER.error(f"Timeout Error: {errt}")
      return (False,)
    except requests.exceptions.RequestException as err:
      constants.PYSSA_LOGGER.error(f"Error: {err}")
      return (False,)
    else:
      constants.PYSSA_LOGGER.info(
          "Demo projects downloaded and extracted successfully."
      )
    try:
      path_of_demo_projects = pathlib.Path(
          f"{constants.SETTINGS_DIR}/demo-projects"
      )
      for tmp_filename in os.listdir(path_of_demo_projects):
        # Copy db file into new workspace
        tmp_project_database_filepath = str(
            pathlib.Path(
                f"{the_workspace_path}/{tmp_filename}",
            ),
        )
        tmp_src_filepath = str(
            pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}")
        )
        shutil.copyfile(tmp_src_filepath, tmp_project_database_filepath)
      constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
    except Exception as e:
      constants.PYSSA_LOGGER.error(
          f"Import process of demo projects finished with the error: {e}."
      )
      return (False,)
    else:
      return (True,)
  except Exception as e:
    logger.error(f"An error occurred: {e}")
    return (False,)
