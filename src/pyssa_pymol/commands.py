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
"""Module with commands for PyMOL."""
import os
from typing import Optional, Union

import pymol
from pymol import cmd
from pymol import CmdException

from src.pyssa.util import constants
from src.pyssa_pymol import local_logging

logger = local_logging.setup_logger(__file__)
__docformat__ = "google"


def reinitialize_session() -> tuple[bool, str]:
  """Reinitializes PyMOL session.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if
          the command was successfully executed.
  """
  logger.info("Executing command.")
  try:
    cmd.reinitialize()
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def load_pymol_session(a_filepath: str) -> tuple[bool, str]:
  """Loads a PyMOL session from a given filepath.

  Args:
      a_filepath (str): The filepath of the PyMOL session to load.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if
          the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_filepath is None or a_filepath == "":
    logger.error("a_filepath is either None or an empty string.")
    return False, "a_filepath is either None or an empty string."
  if not os.path.exists(a_filepath):
    logger.error("a_filepath could not be found.")
    return False, "a_filepath could not be found."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.load(a_filepath)
    # Hide all hydrogens
    cmd.hide("everything", "h.")
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def save_pymol_session(a_filepath: str) -> tuple[bool, str]:
  """Saves the current PyMOL session to the specified file.

  Args:
      a_filepath (str): The file path to save the PyMOL session to.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if
          the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_filepath is None or a_filepath == "":
    logger.error("a_filepath is either None or an empty string.")
    return False, "a_filepath is either None or an empty string."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.save(a_filepath)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def get_all_object_names() -> tuple[bool, str, Optional[list]]:
  """Retrieves all object names in the current loaded PyMOL session.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if
          the command was successfully executed.
          - A list of all object names in the current loaded PyMOL session or None if an error occurred.
  """
  logger.info("Executing command.")
  try:
    tmp_names = cmd.get_names()
  except CmdException as e:
    return False, e.message, None
  else:
    logger.info("Command executed successfully.")
    return True, "", tmp_names


def get_model(a_selection_string: str) -> tuple[bool, str, Optional[int]]:
  """Gets a PyMOL model from the given selection string.

  Args:
      a_selection_string (str): A PyMOL conform selection string.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
          - An int of value 0 if the model is in the session or None if an error occurred.
  """
  # <editor-fold desc="Checks">
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string.", None

  # </editor-fold>

  logger.info("Executing command.")
  try:
    tmp_model = cmd.get_model(a_selection_string)
    tmp_resi = tmp_model.atom[0].resi
  except CmdException as e:
    return True, e.message, None
  except IndexError as e:
    return True, str(e), None
  else:
    logger.info("Command executed successfully.")
    return True, "", 0


def select(a_name: str, a_selection_string: str) -> tuple[bool, str]:
  """Wrapper for the select command of PyMOL.

  Args:
      a_name (str): The name of the selection to be created or modified.
      a_selection_string (str): A PyMOL conform selection string.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_name is None or a_name == "":
    logger.error("a_name is either None or an empty string.")
    return False, "a_name is either None or an empty string."
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string"

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.select(name=a_name, selection=a_selection_string)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def scene(a_key: str, an_action: str) -> tuple[bool, str]:
  """Wrapper for the scene command of PyMOL.

  Args:
      a_key (str): The key/name of the scene.
      an_action (str): The action to perform on the scene.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_key is None or a_key == "":
    logger.error("a_key is either None or an empty string.")
    return False, "a_key is either None or an empty string."
  if an_action is None or an_action == "":
    logger.error("an_action is either None or an empty string.")
    return False, "an_action is either None or an empty string."
  if an_action == "rename":
    logger.error(
        "an_action cannot be rename because a new scene name cannot be defined."
    )
    return (
        False,
        "an_action cannot be rename because a new scene name cannot be defined.",
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
    return False, "an_action is an invalid option."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.refresh()  # This could make a difference but this is not sure.
    cmd.scene(key=a_key, action=an_action)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def load_scene(a_scene_name: str) -> tuple[bool, str]:
  """Loads a scene by name.

  Args:
      a_scene_name (str): The name of the scene to be loaded.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.

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
    return False, "a_scene_name is either None or an empty string."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.scene(a_scene_name, "recall")
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def get_scene_list() -> tuple[bool, str, Optional[list]]:
  """Retrieves a list of scene names of the current PyMOL session.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
          - The list of scene names, if the command executed successfully. If an error occurred None.
  """
  logger.info("Executing command.")
  try:
    tmp_scene_names = cmd.get_scene_list()
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message, None
  else:
    logger.info("Command executed successfully.")
    return True, "", tmp_scene_names


def set_custom_setting(
    a_setting_name: str, a_value: Union[int, float, str]
) -> tuple[bool, str]:
  """Sets a custom setting with the given name and value.

  Args:
      a_setting_name (str): The name of the custom setting.
      a_value (str): The value to set for the custom setting.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_setting_name is None or a_setting_name == "":
    logger.error("a_setting_name is either None or an empty string.")
    return False, "a_setting_name is either None or an empty string."
  if a_value is None:
    logger.error("a_value is None.")
    return False, "a_value is None."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.set(a_setting_name, a_value)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def get_residue_colors(
    a_selection_string: str,
) -> tuple[bool, str, Optional[dict]]:
  """Gets a dict of colors of the given selection.

  Args:
      a_selection_string (str): A PyMOL conform selection string.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
          - dict: A dictionary where the keys are tuples representing a residue (chain, resi, name) and the values are colors for the residue.
  """
  # <editor-fold desc="Checks">
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string", None

  # </editor-fold>

  logger.info("Executing command.")
  try:
    pymol.stored.colors = []
    cmd.iterate(
        a_selection_string, "stored.colors.append((chain, resi, name, color))"
    )
    res_colors: dict = {}
    for chain, resi, name, color_index in pymol.stored.colors:
      if name == "CA":  # c-alpha atom
        res_colors[(chain, resi, name)] = constants.PYMOL_COLORS_WITH_INDICES[
            color_index
        ]
  except Exception as e:
    logger.error(f"Command failed with error: {e}")
    return False, "", None
  else:
    logger.info("Command executed successfully.")
    return True, "", res_colors


def get_chain_color(
    a_selection_string: str, a_chain_letter: str
) -> tuple[bool, str, Optional[tuple]]:
  """Gets a tuple of colors of the given selection.

  Args:
      a_selection_string (str): A PyMOL conform selection string.
      a_chain_letter (str): A letter representing the specific chain within the protein.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
          - A tuple containing the chain color, a list of colors for CA atoms, and a boolean flag indicating if coloring is done by elements.
  """
  # <editor-fold desc="Checks">
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string", None
  if a_chain_letter is None or a_chain_letter == "":
    logger.error("a_chain_letter is either None or an empty string.")
    return False, "a_chain_letter is either None or an empty string", None
  if a_chain_letter not in constants.chain_dict.values():
    logger.error("a_chain_letter is not part of the chain_dict.")
    return False, "a_chain_letter is not part of the chain_dict", None

  # </editor-fold>

  logger.info("Executing command.")
  try:
    pymol.stored.colors = []
    cmd.iterate(
        a_selection_string, "stored.colors.append((chain, resi, name, color))"
    )
    tmp_ca_atom_colors = []
    tmp_residue_atom_colors = []
    tmp_is_colored_by_elements = False
    for chain, resi, name, color_index in pymol.stored.colors:
      if chain == a_chain_letter and name == "CA":  # c-alpha atom
        try:
          tmp_ca_atom_colors.append(
              constants.PYMOL_COLORS_WITH_INDICES[color_index]
          )
        except KeyError:
          tmp_ca_atom_colors.append("By RMSD")
      else:
        try:
          tmp_ca_atom_colors.append(
              constants.PYMOL_COLORS_WITH_INDICES[color_index]
          )
        except KeyError:
          tmp_is_colored_by_elements = True
    tmp_chain_color = list(set(tmp_ca_atom_colors))[0]
  except Exception as e:
    logger.error(f"Command failed with error: {e}")
    return False, "", None
  else:
    logger.info("Command executed successfully.")
    return (
        True,
        "",
        (tmp_chain_color, tmp_ca_atom_colors, tmp_is_colored_by_elements),
    )


def get_residue_color_config(
    a_protein_name: str, a_chain_letter: str
) -> tuple[bool, str, Optional[list]]:
  """Gets the colors of C-, N-, and O-atoms for the first residue of the given selection.

  Args:
      a_protein_name (str): The name of the protein.
      a_chain_letter (str): The letter representing the chain of the protein.

  Returns:
      A tuple containing:
      - A boolean value indicating if the command was executed successfully.
      - An empty string (not used in this method).
      - A list of three elements representing the colors of the residue elements (C, N, O) or None if an error occurred.
  """
  # <editor-fold desc="Checks">
  if a_protein_name is None or a_protein_name == "":
    logger.error("a_protein_name is either None or an empty string.")
    return False, "a_protein_name is either None or an empty string.", None
  if a_chain_letter is None or a_chain_letter == "":
    logger.error("a_chain_letter is either None or an empty string.")
    return False, "a_chain_letter is either None or an empty string", None
  if a_chain_letter not in constants.chain_dict.values():
    logger.error("a_chain_letter is not part of the chain_dict.")
    return False, "a_chain_letter is not part of the chain_dict", None

  # </editor-fold>

  logger.info("Executing command.")
  try:
    pymol.stored.colors = []
    tmp_selections = [
        f"first e. N and chain {a_chain_letter} and {a_protein_name}",
        f"first e. O and chain {a_chain_letter} and {a_protein_name}",
        f"first e. C and chain {a_chain_letter} and {a_protein_name}",
    ]  # example
    tmp_residue_element_colors = ["", "", ""]
    for tmp_selection in tmp_selections:
      cmd.iterate(
          tmp_selection,
          "stored.colors.append((chain, resi, name, color, elem))",
      )
      for chain, resi, name, color_index, element in pymol.stored.colors:
        if chain == a_chain_letter and element == "C":  # C atom
          tmp_residue_element_colors[0] = constants.PYMOL_COLORS_WITH_INDICES[
              color_index
          ]
        elif chain == a_chain_letter and element == "N":  # N atom
          tmp_residue_element_colors[1] = constants.PYMOL_COLORS_WITH_INDICES[
              color_index
          ]
        elif chain == a_chain_letter and element == "O":  # O atom
          tmp_residue_element_colors[2] = constants.PYMOL_COLORS_WITH_INDICES[
              color_index
          ]
  except Exception as e:
    logger.error(f"Command failed with error: {e}")
    return False, "", None
  else:
    logger.info("Command executed successfully.")
    return True, "", tmp_residue_element_colors


def get_chain_repr_state(
    a_selection_string: str, a_chain_letter: str
) -> tuple[bool, str, Optional[dict]]:
  """Returns the representation state of a specific chain in PyMOL.

  Args:
      a_selection_string (str): A PyMOL conform selection string.
      a_chain_letter (str): The chain letter to search for in the representation state.

  Returns:
      A tuple with three elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or
            an empty string if the command was successfully executed.
          - A dictionary representing the representation state for the chain. The dictionary
            is structured as: {representation_index: representation_state}. If no matching
            chain is found, it returns {-1: ""}.
  """
  # <editor-fold desc="Checks">
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string", None
  if a_chain_letter is None or a_chain_letter == "":
    logger.error("a_chain_letter is either None or an empty string.")
    return False, "a_chain_letter is either None or an empty string", None
  if a_chain_letter not in constants.chain_dict.values():
    logger.error("a_chain_letter is not part of the chain_dict.")
    return False, "a_chain_letter is not part of the chain_dict", None

  # </editor-fold>

  logger.info("Executing command.")
  try:
    pymol.stored.reps = []
    cmd.iterate(
        a_selection_string, "stored.reps.append((chain, resi, name, reps))"
    )
    for chain, resi, name, reps_index in pymol.stored.reps:
      if chain == a_chain_letter:  # c-alpha atom
        logger.debug("Representation index: " + str(reps_index))
        logger.info("Command executed successfully.")
        return True, "", constants.PYMOL_REPR_STATES_WITH_INDICES[reps_index]
    logger.info("Command executed successfully.")
  except Exception as e:
    logger.error(f"Command failed with error: {e}")
    return False, "", None
  else:
    return True, "", {-1: ""}


def show_custom_representation(
    a_representation: str, a_selection_string: str
) -> tuple[bool, str]:
  """Wrapper for the show command with the options `representation` and `selection`.

  Args:
      a_representation (str): The custom representation to show.
      a_selection_string (str): A PyMOL conform selection string.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_representation is None or a_representation == "":
    logger.error("a_representation is either None or an empty string.")
    return False, "a_representation is either None or an empty string."
  if a_representation not in constants.PYMOL_REPS_WITH_INDICES.values():
    logger.error(
        "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict."
    )
    return (
        False,
        "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict.",
    )
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.show(representation=a_representation, selection=a_selection_string)
    # Hide all hydrogens
    cmd.hide(a_representation, "h.")
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def hide_custom_representation(
    a_representation: str, a_selection_string: str
) -> tuple[bool, str]:
  """Wrapper for the hide command with the options `representation` and `selection`.

  Args:
      a_representation (str): The custom representation to show.
      a_selection_string (str): A PyMOL conform selection string.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_representation is None or a_representation == "":
    logger.error("a_representation is either None or an empty string.")
    return False, "a_representation is either None or an empty string."
  if a_representation not in constants.PYMOL_REPS_WITH_INDICES.values():
    logger.error(
        "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict."
    )
    return (
        False,
        "a_representation is not found in the PYMOL_REPS_WITH_INDICES dict.",
    )
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.hide(representation=a_representation, selection=a_selection_string)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def zoom_with_custom_parameters(
    a_selection_string: str,
    a_buffer_size: float = 8.0,
    a_state: int = 0,
    a_complete_flag: int = 0,
) -> tuple[bool, str]:
  """Zooms with custom parameters.

  Args:
      a_selection_string (str): A PyMOL conform selection string.
      a_buffer_size (float): The buffer size around the selection to include in the zoom.
      a_state (int): The state number to apply the zoom to. Default is 0.
      a_complete_flag (int): Flag indicating whether to complete the zoom. Default is 0.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string."
  if a_buffer_size is None or a_buffer_size < 0.0:
    logger.error("a_buffer_size is either None or a value less than 0.")
    return False, "a_buffer_size is either None or a value less than 0."
  if a_state is None or a_state < -1:
    logger.error("a_state is either None or a value less than -1.")
    return False, "a_state is either None or a value less than -1."
  if a_complete_flag is None or a_complete_flag < 0 or a_complete_flag > 1:
    logger.error(
        "a_complete_flag is either None or invalid (it can only take 0 or 1)."
    )
    return (
        False,
        "a_complete_flag is either None or invalid (it can only take 0 or 1).",
    )

  # </editor-fold>

  try:
    cmd.zoom(
        selection=a_selection_string,
        buffer=a_buffer_size,
        state=a_state,
        complete=a_complete_flag,
    )
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def color_selection(
    a_pymol_color: str, a_selection_string: str
) -> tuple[bool, str]:
  """Color the specified PyMOL selection.

  Args:
      a_pymol_color (str): The color to apply to the selection. Must be a valid PyMOL color.
      a_selection_string (str): A PyMOL conform selection string.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_pymol_color is None or a_pymol_color == "":
    logger.error("a_pymol_color is either None or an empty string.")
    return False, "a_pymol_color is either None or an empty string."
  if (
      a_pymol_color not in constants.PYMOL_COLORS_WITH_INDICES.values()
      and a_pymol_color != "atomic"
  ):
    logger.error(
        "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict."
    )
    return (
        False,
        "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict.",
    )
  if a_selection_string is None or a_selection_string == "":
    logger.error("a_selection_string is either None or an empty string.")
    return False, "a_selection_string is either None or an empty string."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.color(a_pymol_color, a_selection_string)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def set_background_color(a_pymol_color: str) -> tuple[bool, str]:
  """Sets the background color in PyMOL.

  Args:
      a_pymol_color (str): The color to set as the background. Must be a valid PyMOL color.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_pymol_color is None or a_pymol_color == "":
    logger.error("a_pymol_color is either None or an empty string.")
    return False, "a_pymol_color is either None or an empty string."
  if a_pymol_color not in constants.PYMOL_COLORS_WITH_INDICES.values():
    logger.error(
        "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict."
    )
    return (
        False,
        "a_pymol_color could not be found in the PYMOL_COLORS_WITH_INDICES dict.",
    )

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.bg_color(a_pymol_color)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def set_default_graphics_settings() -> tuple[bool, str]:
  """Sets the default graphics settings for PyMOL.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  logger.info("Executing command.")
  try:
    cmd.set("valence", 0)
    cmd.set("scene_buttons", 0)
    cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
    cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
    cmd.set("ambient", constants.PYMOL_DEFAULT_AMBIENT)
    cmd.set("cartoon_fancy_helices", constants.PYMOL_DEFAULT_FANCY_HELICES)
    cmd.set(
        "cartoon_discrete_colors",
        constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS,
    )
    cmd.set("cartoon_sampling", constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
    cmd.set("spec_power", constants.PYMOL_DEFAULT_SPEC_POWER)
    cmd.set("spec_reflect", constants.PYMOL_DEFAULT_SPEC_REFLECT)
    cmd.set(
        "ray_transparency_contrast",
        constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST,
    )
    cmd.set(
        "ray_transparency_oblique",
        constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE,
    )  # noqa: E501
    cmd.set(
        "ray_transparency_oblique_power",
        constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER,
    )
    cmd.set("ray_trace_color", constants.PYMOL_DEFAULT_RAY_COLOR)
    cmd.unset("depth_cue")
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def ray(a_width: int, a_height: int, a_renderer: int) -> tuple[bool, str]:
  """Wrapper for the ray command with the options `width`, `height` and `renderer`.

  Args:
      a_width (int): The width parameter.
      a_height (int): The height parameter.
      a_renderer (int): The renderer parameter.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_width is None or a_width < 0:
    logger.error("a_width is either None or has a value less than 0.")
    return False, "a_width is either None or has a value less than 0."
  if a_height is None or a_height < 0:
    logger.error("a_height is either None or has a value less than 0.")
    return False, "a_height is either None or has a value less than 0."
  if a_renderer is None or a_renderer < 0 or a_renderer > 3:
    logger.error(
        "a_renderer is either None or is invalid (only 0 - 3 is valid)."
    )
    return (
        False,
        "a_renderer is either None or is invalid (only 0 - 3 is valid)",
    )

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.ray(width=a_width, height=a_height, renderer=a_renderer)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def draw(
    a_width: int, a_height: int, an_antialias_value: int
) -> tuple[bool, str]:
  """Wrapper for the draw command with the options `width`, `height` and `antialias`.

  Args:
      a_width (int): The width of the drawing.
      a_height (int): The height of the drawing.
      an_antialias_value (int): The level of antialiasing to be applied.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.
  """
  # <editor-fold desc="Checks">
  if a_width is None or a_width < 0:
    logger.error("a_width is either None or has a value less than 0.")
    return False, "a_width is either None or has a value less than 0."
  if a_height is None or a_height < 0:
    logger.error("a_height is either None or has a value less than 0.")
    return False, "a_height is either None or has a value less than 0."
  if an_antialias_value is None or an_antialias_value < 0:
    logger.error(
        "an_antialias_value is either None or has a value less than 0."
    )
    return (
        False,
        "an_antialias_value is either None or has a value less than 0.",
    )

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.draw(width=a_width, height=a_height, antialias=an_antialias_value)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""


def png(an_image_filepath: str, a_dpi_value: int) -> tuple[bool, str]:
  """Wrapper for the png command with the options `filename` and `dpi`.

  Args:
      an_image_filepath (str): The filepath of the image file.
      a_dpi_value (int): The DPI (dots per inch) value to set for the PNG file.

  Returns:
      A tuple with two elements:
          - A boolean indicating whether the command was executed successfully (True) or not (False).
          - A string providing additional error information in case the command failed or an empty string if the command was successfully executed.

  Notes:
      The a_dpi_value argument must be greater than 0. Therefore, the default value of the PyMOL cmd command
      cannot be used.
  """
  # <editor-fold desc="Checks">
  if an_image_filepath is None or an_image_filepath == "":
    logger.error("an_image_filepath is either None or an empty string.")
    return False, "an_image_filepath is either None or an empty string."
  if a_dpi_value is None or a_dpi_value < 0:
    logger.error("a_dpi_value is either None or has a value less than 0.")
    return False, "a_dpi_value is either None or has a value less than 0."

  # </editor-fold>

  logger.info("Executing command.")
  try:
    cmd.png(filename=an_image_filepath, dpi=a_dpi_value)
  except CmdException as e:
    logger.error(f"Command failed with error: {e}")
    return False, e.message
  else:
    logger.info("Command executed successfully.")
    return True, ""
