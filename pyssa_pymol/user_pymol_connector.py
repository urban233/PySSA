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
"""Module that functions as client for the User PyMOL interface."""
import logging
import os
from typing import Union, Callable

import pygetwindow
import zmq

from application_process import application_process_manager
from pyssa.internal.data_structures.data_classes import pymol_command
from pyssa.logging_pyssa import log_handlers
from pyssa.util import pyssa_exception, exception, constants
from pyssa_pymol import pymol_enums


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class UserPyMOLConnector:
  """Functions as connector to connect PySSA to PyMOL.

  Notes:
      This class is used within PySSA to communicate PyMOL commands.
  """

  def __init__(
      self,
      an_app_process_manager: "application_process_manager.ApplicationProcessManager",
  ) -> None:
    """Initializes an instance of the class.

    Args:
        an_app_process_manager: An instance of the ApplicationProcessManager class.

    Raises:
        ValueError: If an_app_process_manager is None.
    """
    # <editor-fold desc="Checks">
    if an_app_process_manager is None:
      raise ValueError("an_app_process_manager is None.")

    # </editor-fold>

    self._app_process_manager: (
        "application_process_manager.ApplicationProcessManager"
    ) = an_app_process_manager
    
    self.commands: dict["pymol_enums.CommandEnum", Callable] = {
      pymol_enums.CommandEnum.REINITIALIZE_SESSION: self.reinitialize_session,
      pymol_enums.CommandEnum.LOAD_PYMOL_SESSION: self.load_pymol_session,
      pymol_enums.CommandEnum.SAVE_PYMOL_SESSION: self.save_pymol_session,
      pymol_enums.CommandEnum.GET_ALL_OBJECT_NAMES: self.get_all_object_names,
      pymol_enums.CommandEnum.GET_MODEL: self.get_model,
      pymol_enums.CommandEnum.SELECT: self.select,
      pymol_enums.CommandEnum.SCENE: self.scene,
      pymol_enums.CommandEnum.LOAD_SCENE: self.load_scene,
      pymol_enums.CommandEnum.GET_SCENE_LIST: self.get_scene_list,
      pymol_enums.CommandEnum.SET_CUSTOM_SETTINGS: self.set_custom_setting,
      pymol_enums.CommandEnum.GET_RESIDUE_COLORS: self.get_residue_colors,
      pymol_enums.CommandEnum.GET_CHAIN_COLOR: self.get_chain_color,
      pymol_enums.CommandEnum.GET_RESIDUE_COLOR_CONFIG: self.get_residue_color_config,
      pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE: self.get_chain_repr_state,
      pymol_enums.CommandEnum.SHOW_CUSTOM_REPRESENTATION: self.show_custom_representation,
      pymol_enums.CommandEnum.HIDE_CUSTOM_REPRESENTATION: self.hide_custom_representation,
      pymol_enums.CommandEnum.ZOOM_WITH_CUSTOM_PARAMETERS: self.zoom_with_custom_parameters,
      pymol_enums.CommandEnum.COLOR_SELECTION: self.color_selection,
      pymol_enums.CommandEnum.SET_BACKGROUND_COLOR: self.set_background_color,
      pymol_enums.CommandEnum.SET_DEFAULT_GRAPHIC_SETTINGS: self.set_default_graphic_settings,
      pymol_enums.CommandEnum.RAY: self.ray,
      pymol_enums.CommandEnum.DRAW: self.draw,
      pymol_enums.CommandEnum.PNG: self.png,
    }
    
    context = zmq.Context()
    self._main_socket = context.socket(zmq.REQ)
    self._main_socket.setsockopt(zmq.REQ_RELAXED, 1)
    self._main_socket.connect("tcp://127.0.0.1:9070")
    # self._poller = zmq.Poller()
    # self._poller.register(self._main_socket, zmq.POLLIN)

    self._sender_socket = context.socket(zmq.PUSH)
    self._sender_socket.connect("tcp://127.0.0.1:9071")
    self._recv_socket = context.socket(zmq.PULL)
    self._recv_socket.connect("tcp://127.0.0.1:9072")

    self._poller = zmq.Poller()
    self._poller.register(self._recv_socket, zmq.POLLIN)
  
  def run_command(self, a_command_name: "pymol_enums.CommandEnum", the_args: tuple = ()) -> tuple[bool, dict]:
    """Runs a PyMOL command in the User PyMOL.
    
    Args:
      a_command_name (pymol_enums.CommandEnum): The name of the PyMOL command.
      the_args (tuple[bool, dict]): A tuple of the command arguments.
    
    Returns:
      A tuple with two elements:
        - A boolean indicating success or failure.
        - A dict containing the reply from the PyMOL command execution.
    
    Raises:
      exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_command_name is None:
      logger.error("a_command_name is None.")
      raise exception.IllegalArgumentError("a_command_name is None.")
    if the_args is None:
      logger.error("the_args is None.")
      raise exception.IllegalArgumentError("the_args is None.")
    
    # </editor-fold>
    
    try:
      tmp_command_function: Callable = self.commands[a_command_name]
      if len(the_args) == 0:
        tmp_result: dict = tmp_command_function()
      else:
        tmp_result: dict = tmp_command_function(*the_args)
    except Exception as e:
      logger.error(e)
      return False, {}
    else:
      return True, tmp_result
    
  def reinitialize_session(self) -> dict:
    """Reinitializes the PyMOL session.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.
    """
    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.REINITIALIZE_SESSION,
        (0, 0),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def load_pymol_session(self, a_session_filepath: str) -> dict:
    """Loads a PyMOL session from a specified file.

    Args:
        a_session_filepath (str): A string specifying the filepath of the PyMOL session to load.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_session_filepath is None or an empty string.
        exception.IllegalArgumentError: If a_session_filepath could not be found.
        exception.FileIsEmptyError: If a_session_filepath links to an empty file.
    """
    # <editor-fold desc="Checks">
    if a_session_filepath is None or a_session_filepath == "":
      logger.error("a_session_filepath is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_session_filepath is either None or an empty string."
      )
    if not os.path.exists(a_session_filepath):
      logger.error("a_session_filepath could not be found.")
      raise exception.IllegalArgumentError(
          "a_session_filepath could not be found."
      )
    if os.path.getsize(a_session_filepath) == 0:
      logger.error("a_session_filepath links to an empty file.")
      raise exception.FileIsEmptyError(
          "a_session_filepath links to an empty file."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.LOAD_PYMOL_SESSION,
        (0, str(a_session_filepath)),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def save_pymol_session(self, a_session_filepath: str) -> dict:
    """Saves the PyMOL session to the specified file path.

    Args:
        a_session_filepath (str): A string representing the file path to save the PyMOL session.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_session_filepath is None or an empty string.
        exception.IllegalArgumentError: If a_session_filepath could not be found.

    """
    # <editor-fold desc="Checks">
    if a_session_filepath is None or a_session_filepath == "":
      logger.error("a_session_filepath is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_session_filepath is either None or an empty string."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.SAVE_PYMOL_SESSION,
        (0, str(a_session_filepath)),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def get_all_object_names(self) -> dict:
    """Gets all object names of the current PyMOL session.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.
    """
    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.GET_ALL_OBJECT_NAMES,
        (0, 0),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def get_model(self, a_selection_string: str) -> dict:
    """Gets the model of the given selection string.

    Args:
        a_selection_string (str): A PyMOL conform selection string.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If the a_selection_string is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.GET_MODEL,
        (0, a_selection_string),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def select(self, a_name: str, a_selection_string: str) -> dict:
    """Wrapper for the select command of PyMOL.

    Args:
        a_name (str): The name of the selection to be created or modified.
        a_selection_string (str): A PyMOL conform selection string.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_name is either None or an empty string.
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None or a_name == "":
      logger.error("a_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_name is either None or an empty string."
      )
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.SELECT,
        (a_name, a_selection_string),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def scene(self, a_key: str, an_action: str) -> dict:
    """Wrapper for the scene command of PyMOL.

    Args:
        a_key (str): The key/name of the scene.
        an_action (str): The action to perform on the scene.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_key is either None or an empty string.
        exception.IllegalArgumentError: If an_action is either None or an empty string.
        exception.IllegalArgumentError: If an_action cannot be renamed because a new scene name cannot be defined.
        exception.IllegalArgumentError: If an_action is an invalid option.
    """
    # <editor-fold desc="Checks">
    if a_key is None or a_key == "":
      logger.error("a_key is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_key is either None or an empty string."
      )
    if an_action is None or an_action == "":
      logger.error("an_action is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "an_action is either None or an empty string."
      )
    if an_action == "rename":
      logger.error(
          "an_action cannot be rename because a new scene name cannot be defined."
      )
      raise exception.IllegalArgumentError(
          "an_action cannot be renamed because a new scene name cannot be defined."
      )
    if an_action not in [
        "store",
        "recall",
        "insert_after",
        "insert_before",
        "next",
        "previous",
        "update",
        "rename",
        "clear",
        "append",
    ]:
      logger.error("an_action is an invalid option.")
      raise exception.IllegalArgumentError("an_action is an invalid option.")

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.SCENE,
        (a_key, an_action),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def load_scene(self, a_scene_name: str) -> dict:
    """Loads a scene by name.

    Args:
        a_scene_name (str): The name of the scene to be loaded.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_scene_name is either None or an empty string.

    Example usage:
        success, error_message = load_scene("Scene1")
        if success:
            print("Scene loaded successfully.")
        else:
            print(f"Failed to load scene: {error_message}")
    """
    # <editor-fold desc="Checks">
    if a_scene_name is None or a_scene_name == "":
      logger.error("a_scene_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_scene_name is either None or an empty string."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.LOAD_SCENE,
        (0, str(a_scene_name)),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def get_scene_list(self) -> dict:
    """Retrieves a list of scene names of the current PyMOL session.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.
    """
    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.GET_SCENE_LIST,
        (0, 0),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def set_custom_setting(self, a_setting_name: str, a_value: str) -> dict:
    """Sets a custom setting with the given name and value.

    Args:
        a_setting_name (str): The name of the custom setting.
        a_value (str): The value to set for the custom setting.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_setting_name is either None or an empty string.
        exception.IllegalArgumentError: If a_value is None.
    """
    # <editor-fold desc="Checks">
    if a_setting_name is None or a_setting_name == "":
      logger.error("a_setting_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_setting_name is either None or an empty string."
      )
    if a_value is None:
      logger.error("a_value is None.")
      raise exception.IllegalArgumentError("a_value is None.")

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.SET_CUSTOM_SETTINGS,
        (a_setting_name, a_value),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def get_residue_colors(self, a_selection_string: str) -> dict:
    """Gets a dict of colors of the given selection.

    Args:
        a_selection_string (str): A PyMOL conform selection string.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string"
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.GET_RESIDUE_COLORS,
        (0, a_selection_string),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def get_chain_color(
      self, a_selection_string: str, a_chain_letter: str
  ) -> dict:
    """Gets a tuple of colors of the given selection.

    Args:
        a_selection_string (str): A PyMOL conform selection string.
        a_chain_letter (str): A letter representing the specific chain within the protein.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
        exception.IllegalArgumentError: If a_chain_letter is either None or an empty string.
        exception.IllegalArgumentError: If a_chain_letter is not part of the chain_dict.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string"
      )
    if a_chain_letter is None or a_chain_letter == "":
      logger.error("a_chain_letter is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_chain_letter is either None or an empty string."
      )
    if a_chain_letter not in constants.chain_dict.values():
      logger.error("a_chain_letter is not part of the chain_dict.")
      raise exception.IllegalArgumentError(
          "a_chain_letter is not part of the chain_dict."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.GET_CHAIN_COLOR,
        (a_selection_string, a_chain_letter),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def get_residue_color_config(
      self, a_protein_name: str, a_chain_letter: str
  ) -> dict:
    """Gets the colors of C-, N-, and O-atoms for the first residue of the given selection.

    Args:
        a_protein_name (str): The name of the protein.
        a_chain_letter (str): The letter representing the chain of the protein.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_protein_name is either None or an empty string.
        exception.IllegalArgumentError: If a_chain_letter is either None or an empty string.
        exception.IllegalArgumentError: If a_chain_letter is not part of the chain_dict.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )
    if a_chain_letter is None or a_chain_letter == "":
      logger.error("a_chain_letter is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_chain_letter is either None or an empty string."
      )
    if a_chain_letter not in constants.chain_dict.values():
      logger.error(f"a_chain_letter (value: {a_chain_letter}) is not part of the chain_dict.")
      raise exception.IllegalArgumentError(
          "a_chain_letter is not part of the chain_dict."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.GET_RESIDUE_COLOR_CONFIG,
        (a_protein_name, a_chain_letter),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def get_chain_repr_state(
      self, a_selection_string: str, a_chain_letter: str
  ) -> dict:
    """Returns the representation state of a specific chain in PyMOL.

    Args:
        a_selection_string (str): A PyMOL conform selection string.
        a_chain_letter (str): The chain letter to search for in the representation state.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
        exception.IllegalArgumentError: If a_chain_letter is either None or an empty string.
        exception.IllegalArgumentError: If a_chain_letter is not part of the chain_dict.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string"
      )
    if a_chain_letter is None or a_chain_letter == "":
      logger.error("a_chain_letter is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_chain_letter is either None or an empty string"
      )
    if a_chain_letter not in constants.chain_dict.values():
      logger.error(f"a_chain_letter (value: {a_chain_letter}) is not part of the chain_dict.")
      raise exception.IllegalArgumentError(
          "a_chain_letter is not part of the chain_dict"
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE,
        (a_selection_string, a_chain_letter),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def show_custom_representation(
      self, a_representation: str, a_selection_string: str
  ) -> dict:
    """Wrapper for the show command with the options `representation` and `selection`.

    Args:
        a_representation (str): The custom representation to show.
        a_selection_string (str): A PyMOL conform selection string.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_representation is either None or an empty string.
        exception.IllegalArgumentError: If a_representation is not found in the PYMOL_REPS_WITH_INDICES dict.
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_representation is None or a_representation == "":
      logger.error("a_representation is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_representation is either None or an empty string."
      )
    if a_representation not in constants.PYMOL_REPS_WITH_INDICES.values():
      logger.error(
          "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict."
      )
      raise exception.IllegalArgumentError(
          "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict."
      )
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.SHOW_CUSTOM_REPRESENTATION,
        (a_representation, a_selection_string),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def hide_custom_representation(
      self, a_representation: str, a_selection_string: str
  ) -> dict:
    """Wrapper for the hide command with the options `representation` and `selection`.

    Args:
        a_representation (str): The custom representation to show.
        a_selection_string (str): A PyMOL conform selection string.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_representation is either None or an empty string.
        exception.IllegalArgumentError: If a_representation is not found in the PYMOL_REPS_WITH_INDICES dict.
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_representation is None or a_representation == "":
      logger.error("a_representation is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_representation is either None or an empty string."
      )
    if a_representation not in constants.PYMOL_REPS_WITH_INDICES.values():
      logger.error(
          "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict."
      )
      raise exception.IllegalArgumentError(
          "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict."
      )
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.HIDE_CUSTOM_REPRESENTATION,
        (a_representation, a_selection_string),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def zoom_with_custom_parameters(
      self,
      a_selection_string: str,
      a_buffer_size: float = 8.0,
      a_state: int = 0,
      a_complete_flag: int = 0,
  ) -> dict:
    """Zooms with custom parameters.

    Args:
        a_selection_string (str): A PyMOL conform selection string.
        a_buffer_size (float): The buffer size around the selection to include in the zoom.
        a_state (int): The state number to apply the zoom to. Default is 0.
        a_complete_flag (int): Flag indicating whether to complete the zoom. Default is 0.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
        exception.IllegalArgumentError: If a_buffer_size is either None or a value less than 0.
        exception.IllegalArgumentError: If a_state is either None or a value less than -1.
        exception.IllegalArgumentError: If a_complete_flag is either None or invalid (it can only take 0 or 1).
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )
    if a_buffer_size is None or a_buffer_size < 0.0:
      logger.error("a_buffer_size is either None or a value less than 0.")
      raise exception.IllegalArgumentError(
          "a_buffer_size is either None or a value less than 0."
      )
    if a_state is None or a_state < -1:
      logger.error("a_state is either None or a value less than -1.")
      raise exception.IllegalArgumentError(
          "a_state is either None or a value less than -1."
      )
    if a_complete_flag is None or a_complete_flag < 0 or a_complete_flag > 1:
      logger.error(
          "a_complete_flag is either None or invalid (it can only take 0 or 1)."
      )
      raise exception.IllegalArgumentError(
          "a_complete_flag is either None or invalid (it can only take 0 or 1)."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.ZOOM_WITH_CUSTOM_PARAMETERS,
        (a_selection_string, a_buffer_size, a_state, a_complete_flag),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def color_selection(
      self, a_pymol_color: str, a_selection_string: str
  ) -> dict:
    """Color the specified PyMOL selection.

    Args:
        a_pymol_color (str): The color to apply to the selection. Must be a valid PyMOL color.
        a_selection_string (str): A PyMOL conform selection string.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_pymol_color is either None or an empty string.
        exception.IllegalArgumentError: If a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict.
        exception.IllegalArgumentError: If a_selection_string is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_pymol_color is None or a_pymol_color == "":
      logger.error("a_pymol_color is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_pymol_color is either None or an empty string."
      )
    if (
        a_pymol_color not in constants.PYMOL_COLORS_WITH_INDICES.values()
        and a_pymol_color != "atomic"
    ):
      logger.error(
          "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict."
      )
      raise exception.IllegalArgumentError(
          "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict."
      )
    if a_selection_string is None or a_selection_string == "":
      logger.error("a_selection_string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_selection_string is either None or an empty string."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.COLOR_SELECTION,
        (a_pymol_color, a_selection_string),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def set_background_color(self, a_pymol_color: str) -> dict:
    """Sets the background color in PyMOL.

    Args:
        a_pymol_color (str): The color to set as the background. Must be a valid PyMOL color.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_pymol_color is either None or an empty string.
        exception.IllegalArgumentError: If a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict.
    """
    # <editor-fold desc="Checks">
    if a_pymol_color is None or a_pymol_color == "":
      logger.error("a_pymol_color is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_pymol_color is either None or an empty string."
      )
    if a_pymol_color not in constants.PYMOL_COLORS_WITH_INDICES.values():
      logger.error(
          "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict."
      )
      raise exception.IllegalArgumentError(
          "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.SET_BACKGROUND_COLOR,
        (0, a_pymol_color),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def set_default_graphic_settings(self) -> dict:
    """Sets the default graphics settings for PyMOL.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.
    """
    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.SET_DEFAULT_GRAPHIC_SETTINGS,
        (0, 0),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def ray(self, a_width: int, a_height: int, a_renderer: int) -> dict:
    """Wrapper for the ray command with the options `width`, `height` and `renderer`.

    Args:
        a_width (int): The width parameter.
        a_height (int): The height parameter.
        a_renderer (int): The renderer parameter.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_width is either None or has a value less than 0.
        exception.IllegalArgumentError: If a_height is either None or has a value less than 0.
        exception.IllegalArgumentError: If a_renderer is either None or is invalid (only 0 - 3 is valid)
    """
    # <editor-fold desc="Checks">
    if a_width is None or a_width < 0:
      logger.error("a_width is either None or has a value less than 0.")
      raise exception.IllegalArgumentError(
          "a_width is either None or has a value less than 0."
      )
    if a_height is None or a_height < 0:
      logger.error("a_height is either None or has a value less than 0.")
      raise exception.IllegalArgumentError(
          "a_height is either None or has a value less than 0."
      )
    if a_renderer is None or a_renderer < 0 or a_renderer > 3:
      logger.error(
          "a_renderer is either None or is invalid (only 0 - 3 is valid)."
      )
      raise exception.IllegalArgumentError(
          "a_renderer is either None or is invalid (only 0 - 3 is valid)"
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.RAY,
        (a_width, a_height, a_renderer),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def draw(self, a_width: int, a_height: int, an_antialias_value: int) -> dict:
    """Wrapper for the draw command with the options `width`, `height` and `antialias`.

    Args:
        a_width (int): The width of the drawing.
        a_height (int): The height of the drawing.
        an_antialias_value (int): The level of antialiasing to be applied.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If a_width is either None or has a value less than 0.
        exception.IllegalArgumentError: If a_height is either None or has a value less than 0.
        exception.IllegalArgumentError: If an_antialias_value is either None or has a value less than 0.
    """
    # <editor-fold desc="Checks">
    if a_width is None or a_width < 0:
      logger.error("a_width is either None or has a value less than 0.")
      raise exception.IllegalArgumentError(
          "a_width is either None or has a value less than 0."
      )
    if a_height is None or a_height < 0:
      logger.error("a_height is either None or has a value less than 0.")
      raise exception.IllegalArgumentError(
          "a_height is either None or has a value less than 0."
      )
    if an_antialias_value is None or an_antialias_value < 0:
      logger.error(
          "an_antialias_value is either None or has a value less than 0."
      )
      raise exception.IllegalArgumentError(
          "an_antialias_value is either None or has a value less than 0."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.DRAW,
        (a_width, a_height, an_antialias_value),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def png(self, an_image_filepath: str, a_dpi_value: int) -> dict:
    """Wrapper for the png command with the options `filename` and `dpi`.

    Args:
        an_image_filepath (str): The filepath of the image file.
        a_dpi_value (int): The DPI (dots per inch) value to set for the PNG file.

    Returns:
        A dictionary containing the reply from PyMOL or an empty dict if PyMOL crashed.

    Raises:
        exception.IllegalArgumentError: If an_image_filepath is either None or an empty string.
        exception.IllegalArgumentError: If a_dpi_value is either None or has a value less than 0.

    Notes:
        The a_dpi_value argument must be greater than 0. Therefore, the default value of the PyMOL cmd command
        cannot be used.
    """
    # <editor-fold desc="Checks">
    if an_image_filepath is None or an_image_filepath == "":
      logger.error("an_image_filepath is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "an_image_filepath is either None or an empty string."
      )
    if a_dpi_value is None or a_dpi_value < 0:
      logger.error("a_dpi_value is either None or has a value less than 0.")
      raise exception.IllegalArgumentError(
          "a_dpi_value is either None or has a value less than 0."
      )

    # </editor-fold>

    tmp_pymol_command = pymol_command.PyMOLCommand(
        pymol_enums.CommandEnum.PNG,
        (an_image_filepath, a_dpi_value),
    )
    try:
      tmp_reply = self.send_command_to_pymol(
          tmp_pymol_command, self._poller, self._app_process_manager
      )
    except pyssa_exception.PyMOLNotRespondingError as e:
      logger.error(e)
      return {}
    else:
      return tmp_reply

  def send_command_to_pymol(
      self,
      a_pymol_command: "pymol_command.PyMOLCommand",
      the_poller: zmq.Poller,
      the_app_process_manager: "application_process_manager.ApplicationProcessManager",
      a_timeout: int = 1000,
      a_timeout_cycle_number: int = 10,
  ) -> dict:
    """Sends a command to PyMOL and retrieves the response.

    Args:
        the_main_socket: The main socket to communicate with PyMOL.
        a_pymol_command: An instance of PyMOLCommand containing the command to send.
        the_poller: The polling object to check for events.
        the_app_process_manager: The application process manager for PyMOL.
        a_timeout: The timeout value for polling PyMOL response (default is 1000).
        a_timeout_cycle_number: The number of cycles to wait for PyMOL response (default is 10).

    Returns:
        A dictionary containing the response from PyMOL.

    Raises:
        PyMOLNotRespondingError: If PyMOL is crashed and not responding.

    Notes:
        The total timeout is a combination of a_timeout and a_timeout_cycle_number.
        They will be multiplied to generate the total timeout value.
    """
    # the_main_socket.send_json(a_pymol_command.get_command())
    # tmp_reply_is_ready = False
    # i = 0
    # logger.debug(a_pymol_command.command)
    # while tmp_reply_is_ready is False:
    #     logger.debug(f"Waiting for an event to poll (with timeout of {a_timeout})  ...")
    #     tmp_events = the_poller.poll(a_timeout)
    #     print(f"{i} - {tmp_events} - {a_pymol_command.command}")
    #     if tmp_events:
    #         logger.debug("An event is received.")
    #         tmp_reply_is_ready = True
    #     if the_app_process_manager.pymol_crashed():
    #         logger.warning("PyMOL crashed and an exception will now be raised.")
    #         raise pyssa_exception.PyMOLNotRespondingError()
    #     if i == 20:
    #         logger.warning("PyMOL cannot handle the request anymore and an exception will now be raised.")
    #         raise pyssa_exception.PyMOLNotRespondingError()
    #     i += 1
    #
    # logger.debug("Exit while-loop.")
    # if tmp_reply_is_ready:
    #     logger.debug("Return response from PyMOL.")
    #     return the_main_socket.recv_json()
    # self._main_socket.send_json(a_pymol_command.get_command())
    self._sender_socket.send_json(a_pymol_command.get_command())
    tmp_reply_is_ready = False
    i = 0
    logger.debug(a_pymol_command.command)
    while tmp_reply_is_ready is False:
      logger.debug(
          f"Waiting for an event to poll (with timeout of {a_timeout})  ..."
      )
      tmp_events = self._recv_socket.poll(a_timeout)
      logger.debug(f"{i} - {tmp_events} - {a_pymol_command.command}")
      if tmp_events:
        logger.debug("An event is received.")
        tmp_reply_is_ready = True
      if the_app_process_manager.pymol_crashed():
        logger.warning("PyMOL crashed and an exception will now be raised.")
        raise pyssa_exception.PyMOLNotRespondingError()
      if i == a_timeout_cycle_number:
        logger.warning(
            "PyMOL cannot handle anymore requests. An exception will now be raised and PyMOL will be restarted.."
        )
        pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[
            0
        ].close()
        raise pyssa_exception.PyMOLNotRespondingError()
      i += 1

    logger.debug("Exit while-loop.")
    if tmp_reply_is_ready:
      logger.debug("Return response from PyMOL.")
      tmp_reply = self._recv_socket.recv_json()
      logger.debug(tmp_reply)
      return tmp_reply
    return {}

  def reset_connection(self) -> None:
    """Reset the connection by closing the existing sockets and re-initializing them."""
    # Close the sockets
    self._main_socket.close()
    # Re-initialize the sockets
    context = zmq.Context()
    self._main_socket = context.socket(zmq.REQ)
    self._main_socket.setsockopt(zmq.REQ_RELAXED, 1)
    self._main_socket.connect("tcp://127.0.0.1:9070")
    self._poller = zmq.Poller()
    self._poller.register(self._main_socket, zmq.POLLIN)
