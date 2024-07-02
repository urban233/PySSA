#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
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
"""Module for the PyMOL session manager class."""
import collections
import logging
import os.path
import pathlib
from typing import Optional, Callable

from PyQt5 import QtCore

from src.application_process import application_process_manager
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.data_structures import protein, protein_pair
from src.pyssa.internal.data_structures.data_classes import residue_color_config
from src.pyssa.io_pyssa import binary_data
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants, exception, protein_pair_util
from src.pyssa_pymol import user_pymol_connector
from src.pyssa_pymol import pymol_enums
from src.tea.thread import action, task_result_factory, task_manager, task_scheduler
from src.tea.thread import task_result

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class PymolSessionManager:
  """Holds information and manages the active PyMOL session."""

  # <editor-fold desc="Class attributes">
  session_name: str
  """The name of the active PyMOL session."""

  session_object_type: str
  """The object type (protein or protein_pair) of the active PyMOL session."""

  session_objects: list
  """A list of all objects in the active PyMOL session."""

  current_scene_name: str
  """The name of the current scene."""

  all_scenes: list[str]
  """A list of all scenes in the active PyMOL session."""

  # </editor-fold>

  def __init__(
      self,
      the_app_process_manager: "application_process_manager.ApplicationProcessManager",
  ) -> None:
    """Constructor.

    Args:
        the_app_process_manager (application_process_manager.ApplicationProcessManager): An instance of the ApplicationProcessManager class that manages the application processes.

    Raises:
        exception.IllegalArgumentError: If `the_app_process_manager` is None.
    """
    # <editor-fold desc="Checks">
    if the_app_process_manager is None:
      logger.error("the_app_process_manager is None.")
      raise exception.IllegalArgumentError("the_app_process_manager is None.")

    # </editor-fold>

    self.session_name = ""
    self.session_object_type = ""
    self.session_objects = []
    self.current_scene_name: str = ""
    self.all_scenes: list[str] = []
    self._app_process_manager = the_app_process_manager
    self.user_pymol_connector: "user_pymol_connector.UserPyMOLConnector" = (
        user_pymol_connector.UserPyMOLConnector(
            self._app_process_manager,
        )
    )
    self.lock_user_pymol_connector: QtCore.QMutex = QtCore.QMutex()

  # <editor-fold desc="Private methods">
  def _check_session_integrity(self, a_protein_name: str) -> bool:
    """Checks the session integrity.

    Args:
        a_protein_name (str): The name of a protein.

    Returns:
        A boolean value indicating the integrity of the session.

    Raises:
        exception.IllegalArgumentError: If `a_protein_name` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )

    # </editor-fold>

    tmp_result = self.user_pymol_connector.get_all_object_names()
    if tmp_result == {}:
      logger.warning("get_all_object_names returned an empty dict.")
      return False
    if tmp_result["success"]:
      if a_protein_name in tmp_result["data"]:
        return True
      else:
        return False
    else:
      logger.warning(tmp_result["message"])
      return False

  def _load_pymol_session(self, a_pymol_session: str) -> None:
    """Loads a PyMOL session from a base64 string.

    Args:
        a_pymol_session (str): Base64 encoded string representing the PyMOL session.

    Raises:
        exception.IllegalArgumentError: If `a_pymol_session` is either None or an empty string.
        FileIsEmptyError: If the PyMOL session file is empty.
    """
    # <editor-fold desc="Checks">
    if a_pymol_session is None or a_pymol_session == "":
      logger.error("a_pymol_session is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_pymol_session is either None or an empty string."
      )

    # </editor-fold>

    if self.session_name != "":
      tmp_session_path = pathlib.Path(
          f"{constants.CACHE_PYMOL_SESSION_DIR}/session_of_{self.session_name}.pse",
      )
    else:
      tmp_session_path = pathlib.Path(
          f"{constants.CACHE_PYMOL_SESSION_DIR}/temp_session_file.pse",
      )
    binary_data.write_binary_file_from_base64_string(
        tmp_session_path, a_pymol_session
    )

    if self._app_process_manager.pymol_crashed() is False:
      try:
        tmp_result = self.user_pymol_connector.load_pymol_session(
            str(tmp_session_path)
        )
      except exception.FileIsEmptyError:
        logger.warning(
            "load_pymol_session raised an 'exception.FileIsEmptyError'."
        )
        return
      if tmp_result == {}:
        logger.warning("load_pymol_session returned an empty dict.")
        return
      if tmp_result["success"]:
        self.show_sequence_view()
        logger.info("Session loaded.")
      else:
        logger.warning(tmp_result["message"])
        logger.error("Session loaded failed!")

  def _convert_pymol_session_to_base64_string(
      self, pymol_molecule_object: str
  ) -> str:
    """Converts a pymol session file into a base64 string.

    Args:
        pymol_molecule_object (str): PyMOL molecule object to be converted.

    Raises:
        exception.IllegalArgumentError: If `pymol_molecule_object` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if pymol_molecule_object is None or pymol_molecule_object == "":
      logger.error("pymol_molecule_object is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "pymol_molecule_object is either None or an empty string."
      )

    # </editor-fold>

    session_filepath = pathlib.Path(
        f"{constants.SCRATCH_DIR}/{pymol_molecule_object}_session.pse"
    )
    tmp_result = self.user_pymol_connector.save_pymol_session(
        str(session_filepath)
    )
    if tmp_result == {}:
      logger.warning("save_pymol_session returned an empty dict.")
      return ""
    if tmp_result["success"]:
      base64_string = binary_data.create_base64_string_from_file(
          session_filepath
      )
      os.remove(session_filepath)
      return base64_string
    else:
      logger.warning(tmp_result["message"])
      return ""

  # </editor-fold>

  # <editor-fold desc="Public methods">
  # <editor-fold desc="Non-cmd methods">
  def is_the_current_pymol_scene_base(self) -> bool:
    """Checks if the current scene in PyMOL is the base scene.

    Returns:
        True if the current scene is base, False otherwise.
    """
    if self.current_scene_name == "base":
      return True
    return False

  def is_the_current_session_empty(self) -> bool:
    """Checks if the manager is in an empty session state."""
    if (
        self.session_name == ""
        and self.session_object_type == ""
        and self.session_objects == []
    ):
      return True
    else:
      return False

  def is_the_current_protein_in_session(
      self, the_name_of_the_selected_protein: str
  ) -> bool:
    """Determines if the current protein is in the session.

    Args:
        the_name_of_the_selected_protein (str): The name of the selected protein.

    Returns:
        True if the current protein is in the session, False otherwise.

    Raises:
        exception.IllegalArgumentError: If `the_name_of_the_selected_protein` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if (
        the_name_of_the_selected_protein is None
        or the_name_of_the_selected_protein == ""
    ):
      logger.error(
          "the_name_of_the_selected_protein is either None or an empty string."
      )
      raise exception.IllegalArgumentError(
          "the_name_of_the_selected_protein is either None or an empty string."
      )

    # </editor-fold>

    if (
        self.session_object_type == "protein"
        and self.session_name == the_name_of_the_selected_protein
    ):
      return True
    else:
      return False

  def is_the_current_protein_pair_in_session(
      self, the_name_of_the_selected_protein_pair
  ) -> bool:
    """Determines if the current protein pair is in the session.

    Args:
        the_name_of_the_selected_protein_pair (str): The name of the selected protein pair.

    Returns:
        True if the current protein pair is in the session, False otherwise.

    Raises:
        exception.IllegalArgumentError: If `the_name_of_the_selected_protein_pair` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if (
        the_name_of_the_selected_protein_pair is None
        or the_name_of_the_selected_protein_pair == ""
    ):
      logger.error(
          "the_name_of_the_selected_protein_pair is either None or an empty string."
      )
      raise exception.IllegalArgumentError(
          "the_name_of_the_selected_protein_pair is either None or an empty string."
      )

    # </editor-fold>

    if (
        self.session_object_type == "protein_pair"
        and self.session_name == the_name_of_the_selected_protein_pair
    ):
      return True
    else:
      return False

  # </editor-fold>
  
  def cmd(self, a_command_name: "pymol_enums.CommandEnum", the_args: tuple = ()) -> dict:
    """Runs a PyMOL cmd command.
    
    Args:
      a_command_name (pymol_enums.CommandEnum): The name of the command to be executed.
      the_args (tuple): (Optional) The arguments to be passed to the command. Defaults to ().

    Returns:
      The result of the command execution as a dictionary.
    
    Raises:
      exception.PyMOLCommandFailedError: If the pymol command failed.
    """
    with QtCore.QMutexLocker(self.lock_user_pymol_connector):
      tmp_success_flag, tmp_result = self.user_pymol_connector.run_command(a_command_name, the_args)
    if tmp_success_flag:
      return tmp_result
    raise exception.PyMOLCommandFailedError(f"The pymol command {a_command_name.value} with the args {the_args} failed.")
  
  def async_cmd(self,
                the_task_manager: "task_manager.TaskManager",
                the_task_scheduler: "task_scheduler.TaskScheduler",
                a_command_name: "pymol_enums.CommandEnum",
                the_args: tuple = (),
                an_await_function: Optional[Callable] = None) -> None:
    """Runs a PyMOL cmd command in an asynchronous manner.
    
    The PyMOL command gets executed in the User PyMOL.
    The results will be given to the await function as argument.
    Therefore, the await function needs an argument to use the results.
    
    Args:
      the_task_manager (task_manager.TaskManager): The task manager instance of the interface manager.
      the_task_scheduler (task_scheduler.TaskScheduler): The task scheduler instance of the interface manager.
      a_command_name (pymol_enums.CommandEnum): The name of the command to be executed.
      the_args (tuple): (Optional) The arguments to be passed to the command. Defaults to ().
      an_await_function (Callable): An await function to call after the async function finished. Defaults to None.
      
    Raises:
      exception.IllegalArgumentError: If any of the arguments is None except an_await_function.
    """
    # <editor-fold desc="Checks">
    if the_task_manager is None:
      logger.error("the_task_manager is None.")
      raise exception.IllegalArgumentError("the_task_manager is None.")
    if the_task_scheduler is None:
      logger.error("the_task_scheduler is None.")
      raise exception.IllegalArgumentError("the_task_scheduler is None.")
    if a_command_name is None:
      logger.error("a_command_name is None.")
      raise exception.IllegalArgumentError("a_command_name is None.")
    if the_args is None:
      logger.error("the_args is None.")
      raise exception.IllegalArgumentError("the_args is None.")
    
    # </editor-fold>
    
    the_task_manager.append_task_result(
      task_result_factory.TaskResultFactory.run_task_result(
        a_task_result=task_result.TaskResult.from_action(
          an_action=action.Action(
            a_target=self.cmd,
            args=(
              a_command_name,
              the_args,
            ),
          ),
          an_await_function=an_await_function,
        ),
        a_task_scheduler=the_task_scheduler,
      )
    )

  def async_cmds(self,
                 the_task_manager: "task_manager.TaskManager",
                 the_task_scheduler: "task_scheduler.TaskScheduler",
                 the_command_names: tuple["pymol_enums.CommandEnum"],
                 the_args: tuple[tuple],
                 an_await_function: Optional[Callable] = None) -> None:
    """Runs multiple PyMOL cmd command sequentially in an asynchronous manner.

    The PyMOL command gets executed in the User PyMOL.
    The results will be given to the await function as argument.
    Therefore, the await function needs an argument to use the results.

    Args:
      the_task_manager (task_manager.TaskManager): The task manager instance of the interface manager.
      the_task_scheduler (task_scheduler.TaskScheduler): The task scheduler instance of the interface manager.
      the_command_names (tuple[pymol_enums.CommandEnum]): The name of the command to be executed.
      the_args tuple[tuple]: The arguments to be passed to the command.
      an_await_function (Callable): An await function to call after the async function finished. Defaults to None.
    
    Raises:
      exception.IllegalArgumentError: If any of the arguments is None except `an_await_function` or `the_command_names` or `the_args` are an empty list.
    """
    # <editor-fold desc="Checks">
    if the_task_manager is None:
      logger.error("the_task_manager is None.")
      raise exception.IllegalArgumentError("the_task_manager is None.")
    if the_task_scheduler is None:
      logger.error("the_task_scheduler is None.")
      raise exception.IllegalArgumentError("the_task_scheduler is None.")
    if the_command_names is None or len(the_command_names) == 0:
      logger.error("the_command_names is either None or an empty list.")
      raise exception.IllegalArgumentError("the_command_names is either None or an empty list.")
    if the_args is None or len(the_args) == 0:
      logger.error("the_args is either None or an empty list.")
      raise exception.IllegalArgumentError("the_args is either None or an empty list.")
    
    # </editor-fold>
    
    tmp_actions: collections.deque = collections.deque()
    i = 0
    for tmp_command_name in the_command_names:
      tmp_actions.append(
        action.Action(
          a_target=self.cmd,
          args=(tmp_command_name, the_args[i])
        )
      )
      i += 1
    
    the_task_manager.append_task_result(
      task_result_factory.TaskResultFactory.run_task_result(
        a_task_result=task_result.TaskResult.from_actions(
          the_actions=tuple(tmp_actions),
          an_await_function=an_await_function,
        ),
        a_task_scheduler=the_task_scheduler,
      )
    )
    
  # <editor-fold desc="Cmd-depended methods">
  def reinitialize_session(self) -> None:
    """Reinitialize the pymol session and class attributes."""
    # reset class attributes
    self.session_name = ""
    self.session_object_type = ""
    self.session_objects: list = []
    # reset actual pymol session
    self.user_pymol_connector.reinitialize_session()

  def load_protein_session(self, a_protein: "protein.Protein") -> None:
    """Loads a protein session into the current session.

    Args:
        a_protein (protein.Protein): The protein session to load.

    Raises:
        exception.IllegalArgumentError: If `a_protein` is None.
        ProteinNotFoundInPyMOLSession: If loading the PyMOL session fails because the protein cannot be found in the PyMOL object list.
    """
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    self.session_name = a_protein.get_molecule_object()
    self.session_object_type = "protein"
    self.session_objects = [a_protein]
    self._load_pymol_session(a_protein.pymol_session)

    # <editor-fold desc="Integrity check">
    if not self._check_session_integrity(a_protein.get_molecule_object()):
      raise exception.ProteinNotFoundInPyMOLSession(
          f"Loading the PyMOL session failed, because the protein {self.session_name} can not be found in the PyMOL object list."
      )

    # </editor-fold>

  def load_protein_pair_session(
      self, a_protein_pair: "protein_pair.ProteinPair"
  ) -> None:
    """Loads a protein session into the current session.

    Args:
        a_protein_pair (protein_pair.ProteinPair): An instance of the ProteinPair class representing the protein pair to be loaded into the session.

    Raises:
        exception.IllegalArgumentError: If `a_protein_pair` is None.
        exception.ProteinNotFoundInPyMOLSession: If any of the proteins in the protein pair cannot be found in the PyMOL object list.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    # </editor-fold>

    self.session_name = a_protein_pair.name
    self.session_object_type = "protein_pair"
    self.session_objects = [a_protein_pair]
    self._load_pymol_session(a_protein_pair.pymol_session)
    self.user_pymol_connector.scene("base", an_action="recall")
    self.user_pymol_connector.reset()
    self.user_pymol_connector.scene("base", an_action="update")

    # <editor-fold desc="Integrity check">
    if (
        a_protein_pair.protein_1.get_molecule_object()
        == a_protein_pair.protein_2.get_molecule_object()
    ):
      if not self._check_session_integrity(
          f"{a_protein_pair.protein_1.get_molecule_object()}_1"
      ):
        raise exception.ProteinNotFoundInPyMOLSession(
            f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()}_1 can not be found in the PyMOL object list."
        )
      if not self._check_session_integrity(
          f"{a_protein_pair.protein_2.get_molecule_object()}_2"
      ):
        raise exception.ProteinNotFoundInPyMOLSession(
            f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_2.get_molecule_object()}_2 can not be found in the PyMOL object list."
        )
    else:
      if not self._check_session_integrity(
          a_protein_pair.protein_1.get_molecule_object()
      ):
        raise exception.ProteinNotFoundInPyMOLSession(
            f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()} can not be found in the PyMOL object list."
        )
      if not self._check_session_integrity(
          a_protein_pair.protein_2.get_molecule_object()
      ):
        raise exception.ProteinNotFoundInPyMOLSession(
            f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_2.get_molecule_object()} can not be found in the PyMOL object list."
        )

    # </editor-fold>

  def load_scene(self, a_scene_name: str) -> bool:
    """Loads the scene with the given scene name.

    Args:
        a_scene_name: A string representing the name of the scene to be loaded.

    Returns:
        A boolean value indicating whether the scene was successfully loaded.

    Raises:
        exception.IllegalArgumentError: If `a_scene_name` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_scene_name is None or a_scene_name == "":
      logger.error("a_scene_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_scene_name is either None or an empty string."
      )

    # </editor-fold>

    tmp_result = self.user_pymol_connector.load_scene(a_scene_name)
    if tmp_result == {}:
      logger.warning("save_pymol_session returned an empty dict.")
      return False
    return tmp_result["success"]

  def load_current_scene(self) -> bool:
    """Loads the current scene.

    Returns:
        True if the current scene is loaded successfully, False otherwise.
    """
    tmp_result = self.user_pymol_connector.load_scene(self.current_scene_name)
    if tmp_result == {}:
      logger.warning("save_pymol_session returned an empty dict.")
      return False
    return tmp_result["success"]

  def save_current_pymol_session_as_pse_cache_file(self) -> pathlib.Path:
    """Save the current Pymol session as a PSE cache file.

    Returns:
        The path to the saved PSE cache file.
    """
    tmp_session_path = pathlib.Path(
        f"{constants.CACHE_PYMOL_SESSION_DIR}/{self.session_name}.pse",
    )
    tmp_session_base64 = self._convert_pymol_session_to_base64_string(
        self.session_name
    )
    binary_data.write_binary_file_from_base64_string(
        tmp_session_path, tmp_session_base64
    )
    return tmp_session_path

  def save_current_session_as_base64(self) -> str:
    """Saves the current session as a Base64-encoded string.

    Returns:
        The Base64-encoded string representation of the saved session.
    """
    tmp_session_filepath = pathlib.Path(
        f"{constants.CACHE_PYMOL_SESSION_DIR}/session_of_{self.session_name}.pse"
    )
    tmp_session_base64 = self._convert_pymol_session_to_base64_string(
        self.session_name
    )
    binary_data.write_binary_file_from_base64_string(
        tmp_session_filepath, tmp_session_base64
    )
    return binary_data.create_base64_string_from_file(tmp_session_filepath)

  def get_all_scenes_in_current_session(self) -> None:
    """Clears the list of all scenes in the current session and retrieves the list of scenes from the user_pymol_connector."""
    self.all_scenes.clear()
    tmp_result = self.user_pymol_connector.get_scene_list()
    if tmp_result == {}:
      logger.warning("save_pymol_session returned an empty dict.")
      self.all_scenes = []
      return
    if tmp_result["success"]:
      self.all_scenes = tmp_result["data"]
    else:
      self.all_scenes = []

  def set_all_scenes_for_current_session(self, all_scenes: list) -> None:
    """Sets all scenes for the current session.

    Args:
        all_scenes (list): A list of scenes to set for the current session.

    Raises:
        exception.IllegalArgumentError: If `all_scenes` is None.
    """
    # <editor-fold desc="Checks">
    if all_scenes is None:
      logger.error("all_scenes is None.")
      raise exception.IllegalArgumentError("all_scenes is None.")

    # </editor-fold>

    self.all_scenes.clear()
    self.all_scenes = all_scenes

  def show_sequence_view(self) -> None:
    """Sets the custom setting "seq_view" to 1.

    This method is used to show the sequence view in PyMOL for the current session.
    """
    self.user_pymol_connector.set_custom_setting("seq_view", 1)

  def hide_sequence_view(self) -> None:
    """Sets the custom setting "seq_view" to 0.

    This method is used to hides the sequence view in PyMOL for the current session.
    """
    self.user_pymol_connector.set_custom_setting("seq_view", 0)

  # </editor-fold>

  def show_specific_representation(
      self, a_representation: str, a_selection_string: str
  ) -> None:
    """Shows a specific representation of the selection.

    Args:
        a_representation: (str) The specific representation to be shown in PyMOL.
        a_selection_string: (str) The selection string specifying the atoms or residues to apply the representation to.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_representation is None or a_representation == "":
      logger.error("a_representation is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_representation is either None or an empty string."
      )
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    self.user_pymol_connector.show_custom_representation(
        a_representation, a_selection_string
    )

  def hide_specific_representation(
      self, a_representation: str, a_selection_string: str
  ) -> None:
    """Hides a specific representation of the selection.

    Args:
        a_representation: (str) The specific representation to be hidden in PyMOL.
        a_selection_string: (str) The selection string specifying the atoms or residues to apply the representation to.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_representation is None or a_representation == "":
      logger.error("a_representation is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_representation is either None or an empty string."
      )
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    self.user_pymol_connector.hide_custom_representation(
        a_representation, a_selection_string
    )

  def get_residue_colors(self, a_selection_string: str) -> Optional[dict]:
    """Gets all the colors of a residue in the selection.

    Args:
        a_selection_string (str): The selection string used to specify the residues.

    Returns:
        A dictionary containing the residue colors. Returns None if an error occurred.

    Raises:
        exception.IllegalArgumentError: If `a_selection_string` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    tmp_result = self.user_pymol_connector.get_residue_colors(
        a_selection_string
    )
    if tmp_result["success"]:
      return tmp_result["data"]
    return None

  def get_chain_color(
      self, a_selection_string: str, chain_letter: str
  ) -> Optional[dict]:
    """Gets all the colors of a chain in the selection.

    Args:
        a_selection_string (str): The selection string to specify the chains to be considered.
        chain_letter (str): The chain letter to specify the chain color to retrieve.

    Returns:
        A dictionary containing the chain color information if the operation is successful, otherwise returns None.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )
    if chain_letter is None or chain_letter == "":
      logger.error("chain_letter is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "chain_letter is either None or an empty string."
      )

    # </editor-fold>

    tmp_result = self.user_pymol_connector.get_chain_color(
        a_selection_string, chain_letter
    )
    print(tmp_result)
    if tmp_result["success"]:
      return tmp_result["data"]
    return None

  def get_residue_color_config_of_a_given_selection(
      self, a_protein_name: str, chain_letter: str
  ) -> "residue_color_config.ResidueColorConfig":
    """Gets the residue color configuration for a specific chain.

    Args:
        a_protein_name (str): The name of the protein.
        chain_letter (str): The letter representing the chain of the protein.

    Returns:
        The color configuration of the specified protein selection.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )
    if chain_letter is None or chain_letter == "":
      logger.error("chain_letter is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "chain_letter is either None or an empty string."
      )

    # </editor-fold>

    tmp_result = self.user_pymol_connector.get_residue_color_config(
        a_protein_name, chain_letter
    )
    if tmp_result["success"]:
      tmp_color_config: list = tmp_result["data"]
      return residue_color_config.ResidueColorConfig(
          tmp_color_config[0], tmp_color_config[1], tmp_color_config[2]
      )
    return residue_color_config.ResidueColorConfig("", "", "")

  def get_chain_repr_state(
      self, a_selection_string: str, chain_letter: str
  ) -> Optional[dict]:
    """Gets the representation state for a certain chain.

    Args:
        a_selection_string (str): A selection in the PyMOL session.
        chain_letter (str): The chain letter of the protein in the PyMOL session.

    Returns:
        A dictionary representing the state of the chain representation in the PyMOL session.
        Returns None if the chain representation state is not available or an error occurred.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )
    if chain_letter is None or chain_letter == "":
      logger.error("chain_letter is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "chain_letter is either None or an empty string."
      )

    # </editor-fold>

    tmp_result = self.user_pymol_connector.get_chain_repr_state(
        a_selection_string, chain_letter
    )
    if tmp_result == {}:
      logger.warning("save_pymol_session returned an empty dict.")
      return None
    if tmp_result["success"]:
      return tmp_result["data"]
    return None

  def show_protein_selection_as_balls_and_sticks(self, selection: str) -> None:
    """Displays the specified protein selection as balls and sticks representation.

    Args:
        selection (str): The selection of atoms and/or residues in the protein to be displayed as balls and sticks.

    Raises:
        exception.IllegalArgumentError: If `selection` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if selection is None or selection == "":
      logger.error("selection is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "selection is either None or an empty string."
      )

    # </editor-fold>

    self.user_pymol_connector.show_custom_representation("sticks", selection)

  def hide_protein_selection_as_balls_and_sticks(self, selection: str) -> None:
    """Hides the balls and sticks representation for the specified protein selection.

    Args:
        selection (str): The selection of atoms and/or residues in the protein to hide the balls and sticks from.

    Raises:
        exception.IllegalArgumentError: If `selection` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if selection is None or selection == "":
      logger.error("selection is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "selection is either None or an empty string."
      )

    # </editor-fold>

    self.user_pymol_connector.hide_custom_representation("sticks", selection)

  def zoom_to_residue_in_protein_position(self, selection: str) -> None:
    """Zooms to the specified selection.

    Args:
        selection (str): A PyMOL selection string.

    Raises:
        exception.IllegalArgumentError: If `selection` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if selection is None or selection == "":
      logger.error("selection is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "selection is either None or an empty string."
      )

    # </editor-fold>

    self.user_pymol_connector.zoom_with_custom_parameters(selection)

  def color_protein(self, pymol_color: str, a_selection_string: str) -> None:
    """Colors a specific protein selection with a given PyMOL color.

    Args:
        pymol_color (str): A color which is available in PyMOL.
        a_selection_string (str): A PyMOL conform selection string.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or an empty string.
        ValueError: If `pymol_color` is not part of the constants.PYMOL_COLORS.
    """
    # <editor-fold desc="Checks">
    if pymol_color is None or pymol_color == "":
      logger.error("pymol_color is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "pymol_color is either None or an empty string."
      )
    if pymol_color not in constants.PYMOL_COLORS:
      raise ValueError(f"An illegal color argument. {pymol_color}")
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    self.user_pymol_connector.color_selection(pymol_color, a_selection_string)

  def color_protein_pair_by_rmsd(
      self, a_protein_pair: "protein_pair.ProteinPair"
  ) -> None:
    """Colors a specific protein pair based on their rmsd value.

    Args:
        a_protein_pair (protein_pair.ProteinPair): The protein pair to color.

    Raises:
        exception.IllegalArgumentError: If `a_protein_pair` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    # </editor-fold>

    cutoff_1 = 0.5
    cutoff_2 = 1.0
    cutoff_3 = 2
    cutoff_4 = 4
    cutoff_5 = 6

    color_1 = "br0"
    color_2 = "br2"
    color_3 = "br4"
    color_4 = "br6"
    color_5 = "br8"
    color_6 = "red"

    self.user_pymol_connector.color_selection(
        "hydrogen", a_protein_pair.protein_2.get_molecule_object()
    )

    i: int = 0
    for (
        distance_value
    ) in a_protein_pair.distance_analysis.analysis_results.distance_data.get(
        "distance"
    ):
      if distance_value <= cutoff_1:
        atom_info = protein_pair_util.get_chain_and_position(
            a_protein_pair.distance_analysis.analysis_results.distance_data,
            i,
        )
        # create two atoms for the get_distance command
        atom1: str = (
            f"/{a_protein_pair.protein_1.get_molecule_object()}//"
            f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
        )
        atom2: str = (
            f"/{a_protein_pair.protein_2.get_molecule_object()}//"
            f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
        )
        # coloring
        self.user_pymol_connector.color_selection(color_1, atom1)
        self.user_pymol_connector.color_selection(color_1, atom2)
        i += 1

      elif distance_value <= cutoff_2:
        atom_info = protein_pair_util.get_chain_and_position(
            a_protein_pair.distance_analysis.analysis_results.distance_data,
            i,
        )
        # create two atoms for the get_distance command
        atom1: str = (
            f"/{a_protein_pair.protein_1.get_molecule_object()}//"
            f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
        )
        atom2: str = (
            f"/{a_protein_pair.protein_2.get_molecule_object()}//"
            f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
        )
        # coloring
        self.user_pymol_connector.color_selection(color_2, atom1)
        self.user_pymol_connector.color_selection(color_2, atom2)
        i += 1

      elif distance_value <= cutoff_3:
        atom_info = protein_pair_util.get_chain_and_position(
            a_protein_pair.distance_analysis.analysis_results.distance_data,
            i,
        )
        # create two atoms for the get_distance command
        atom1: str = (
            f"/{a_protein_pair.protein_1.get_molecule_object()}//"
            f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}/CA"
        )
        atom2: str = (
            f"/{a_protein_pair.protein_2.get_molecule_object()}//"
            f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}/CA"
        )
        # coloring
        self.user_pymol_connector.color_selection(color_3, atom1)
        self.user_pymol_connector.color_selection(color_3, atom2)
        i += 1

      elif distance_value <= cutoff_4:
        atom_info = protein_pair_util.get_chain_and_position(
            a_protein_pair.distance_analysis.analysis_results.distance_data,
            i,
        )
        # create two atoms for the get_distance command
        atom1: str = (
            f"/{a_protein_pair.protein_1.get_molecule_object()}//"
            f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
        )
        atom2: str = (
            f"/{a_protein_pair.protein_2.get_molecule_object()}//"
            f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
        )
        # coloring
        self.user_pymol_connector.color_selection(color_4, atom1)
        self.user_pymol_connector.color_selection(color_4, atom2)
        i += 1

      elif distance_value <= cutoff_5:
        atom_info = protein_pair_util.get_chain_and_position(
            a_protein_pair.distance_analysis.analysis_results.distance_data,
            i,
        )
        # create two atoms for the get_distance command
        atom1: str = (
            f"/{a_protein_pair.protein_1.get_molecule_object()}//"
            f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
        )
        atom2: str = (
            f"/{a_protein_pair.protein_2.get_molecule_object()}//"
            f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
        )
        # coloring
        self.user_pymol_connector.color_selection(color_5, atom1)
        self.user_pymol_connector.color_selection(color_5, atom2)
        i += 1

      elif distance_value > cutoff_5:
        atom_info = protein_pair_util.get_chain_and_position(
            a_protein_pair.distance_analysis.analysis_results.distance_data,
            i,
        )
        # create two atoms for the get_distance command
        atom1: str = (
            f"/{a_protein_pair.protein_1.get_molecule_object()}//"
            f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
        )
        atom2: str = (
            f"/{a_protein_pair.protein_2.get_molecule_object()}//"
            f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
        )
        # coloring
        self.user_pymol_connector.color_selection(color_6, f"({atom1})")
        self.user_pymol_connector.color_selection(color_6, f"({atom2})")
        i += 1

  def setup_default_session_graphic_settings(self) -> None:
    """This functions modifies the pymol session to look fancy."""
    self.user_pymol_connector.set_background_color(
        constants.PYMOL_DEFAULT_BACKGROUND_COLOR
    )
    self.user_pymol_connector.set_default_graphic_settings()

  def setup_default_image_graphic_settings(
      self, ray_shadows: bool, opaque_background: int = 0
  ) -> None:
    """Sets up the default image graphic settings for PyMOL.

    Args:
        ray_shadows (bool): False if no shadows, True if shadows should be displayed.
        opaque_background (int, optional): 0 for a transparent background and 1 for a white background.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if ray_shadows is None:
      logger.error("ray_shadows is None.")
      raise exception.IllegalArgumentError("ray_shadows is None.")
    if opaque_background is None:
      logger.error("opaque_background is None.")
      raise exception.IllegalArgumentError("opaque_background is None.")

    # </editor-fold>

    if not ray_shadows:
      opt_ray_shadows: str = "off"
    else:
      opt_ray_shadows: str = "on"
    self.user_pymol_connector.set_background_color(
        constants.PYMOL_DEFAULT_BACKGROUND_COLOR
    )
    self.user_pymol_connector.set_custom_setting(
        "ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE
    )
    self.user_pymol_connector.set_custom_setting(
        "antialias", constants.PYMOL_DEFAULT_ANTIALIAS
    )
    self.user_pymol_connector.set_custom_setting("ray_shadows", opt_ray_shadows)
    self.user_pymol_connector.set_custom_setting(
        "ray_opaque_background", opaque_background
    )

  def setup_default_graphic_settings_for_interesting_regions(self) -> None:
    """Sets up the default graphic settings for interesting regions."""
    self.user_pymol_connector.set_background_color(
        constants.PYMOL_DEFAULT_BACKGROUND_COLOR
    )
    self.user_pymol_connector.set_custom_setting("label_size", str(14))
    self.user_pymol_connector.set_custom_setting("label_font_id", str(13))
    self.user_pymol_connector.set_custom_setting("label_color", "hotpink")
    self.user_pymol_connector.set_custom_setting("depth_cue", str(0))
    # interacts directly with molecule objects in the session
    self.user_pymol_connector.hide_custom_representation("cartoon", "all")
    self.user_pymol_connector.show_custom_representation("ribbon", "all")

  def check_if_sele_is_empty(self) -> bool:
    """Check if the sele object in PyMOL is empty.

    Returns:
        True if the sele object is empty or there is an error retrieving the sele object. False if the sele object is not empty.
    """
    tmp_result = self.user_pymol_connector.get_model("sele")
    if tmp_result["success"] and tmp_result["data"] is None:
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "Please select at least one residue from the sequence view.",
          "PyMOL Selection",
          custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
      )
      tmp_dialog.exec_()
      return True
    elif tmp_result["success"] and tmp_result["data"] is not None:
      return False
    else:
      return True
    # else:
    #     # gets thrown if no sele object exists in pymol
    #     tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         "Please select at least one residue from the sequence view.",
    #         "PyMOL Selection",
    #         custom_message_box.CustomMessageBoxIcons.INFORMATION.value
    #     )
    #     tmp_dialog.exec_()
    #     return True
    # try:
    #     tmp_selection.atom[0].resi
    # except IndexError:
    #     # gets thrown if sele object is empty
    #     tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         "Please select at least one residue from the sequence view.",
    #         "PyMOL Selection",
    #         custom_message_box.CustomMessageBoxIcons.INFORMATION.value
    #     )
    #     tmp_dialog.exec_()
    #     return True
    # return False

  def check_if_specific_selection_is_empty(
      self, a_selection_string: str
  ) -> bool:
    """Checks if a specific selection is empty.

    Args:
        a_selection_string (str): The string representation of the selection.

    Returns:
        True if the selection is empty, False otherwise.

    Raises:
        exception.IllegalArgumentError: If `a_selection_string` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    tmp_result = self.user_pymol_connector.get_model(a_selection_string)
    if tmp_result["success"]:
      tmp_selection = tmp_result["data"]
    else:
      # gets thrown if no sele object exists in pymol
      return True
    try:
      tmp_selection.atom[0].resi
    except IndexError:
      # gets thrown if sele object is empty
      return True
    return False

  # </editor-fold>
