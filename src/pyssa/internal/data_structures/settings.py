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
"""Module for handling everything related to the settings.xml."""
import json
import logging
import os
import pathlib

from src.pyssa.io_pyssa import safeguard
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Settings:
  """Holds the setting from the settings.json in memory."""

  def __init__(self, dir_settings: str, filename: str) -> None:
    """Constructor.

    Args:
        dir_settings (str): The directory where the settings.json is stored.
        filename (str): The name of the settings.json.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if dir_settings is None:
      logger.error("dir_settings is None.")
      raise exception.IllegalArgumentError("dir_settings is None.")
    if filename is None:
      logger.error("filename is None.")
      raise exception.IllegalArgumentError("filename is None.")

    # </editor-fold>

    if not os.path.exists(constants.SETTINGS_DIR):
      os.mkdir(constants.SETTINGS_DIR)
    if not os.path.exists(constants.DEFAULT_WORKSPACE_PATH):
      os.mkdir(constants.DEFAULT_WORKSPACE_PATH)

    self.workspace_path = constants.DEFAULT_WORKSPACE_PATH
    self.cycles: int = 0
    self.cutoff: float = 1.0
    self.app_launch = 1
    self.dir_settings: str = dir_settings
    self.filename: str = filename
    self.local_colabfold: int = 1
    self.wsl_install: int = 1
    self.color_vision_mode: str = "normal"
    self.ask_save_pymol_session: int = 0
    self.image_background_color: str = "black"
    self.image_renderer: str = "-1"  # or "0"
    self.image_ray_trace_mode: int = 1  # ranges from 0 to 3
    self.image_ray_texture: int = 0  # ranges from 0 to 5
    self.proteins_tab_use_toggle: int = 1  # or "0"
    self.proteins_tab_use_combobox_for_colors: int = 0  # or "1"
    self.protein_pairs_tab_use_toggle: int = 1  # or "0"
    self.protein_pairs_tab_use_combobox_for_colors: int = 0  # or "1"

  def serialize_settings(self) -> None:
    """Serializes the protein object."""
    if not os.path.exists(constants.SETTINGS_DIR):
      os.mkdir(constants.SETTINGS_DIR)
    tmp_settings_filepath = (
        f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}"
    )
    if os.path.exists(tmp_settings_filepath):
      os.remove(tmp_settings_filepath)
    settings_dict = self.__dict__
    print(settings_dict)
    update = {
        "workspace_path": str(self.workspace_path),
        "dir_settings": str(self.dir_settings),
    }
    settings_dict.update(update)
    with open(
        f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}",
        "w",
        encoding="utf-8",
    ) as settings_file:
      json.dump(settings_dict, settings_file, indent=4)

  @staticmethod
  def deserialize_settings() -> "Settings":
    """Constructs the settings object from the json file.

    Returns:
        A complete settings object deserialized from a json file.
    """
    if not os.path.exists(
        f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}"
    ):
      raise FileNotFoundError()  # TODO: more precise error is needed
    with open(
        pathlib.Path(
            f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}"
        ),
        "r",
        encoding="utf-8",
    ) as settings_obj_file:
      settings_dict = json.load(settings_obj_file)
    try:
      tmp_settings: Settings = Settings(
          settings_dict.get("dir_settings"), settings_dict.get("filename")
      )

      tmp_settings.workspace_path = Settings._check_integrity_of_workspace_path(
          settings_dict.get("workspace_path"),
      )
      tmp_settings.cycles = Settings._check_integrity_of_cycles(
          int(settings_dict.get("cycles"))
      )
      tmp_settings.cutoff = Settings._check_integrity_of_cutoff(
          float(settings_dict.get("cutoff"))
      )
      tmp_settings.app_launch = Settings._check_integrity_of_app_start_value(
          settings_dict.get("app_launch")
      )
      tmp_settings.wsl_install = Settings._check_integrity_of_wsl_install_flag(
          int(settings_dict.get("wsl_install")),
      )
      tmp_settings.local_colabfold = (
          Settings._check_integrity_of_colabfold_install_flag(
              int(settings_dict.get("local_colabfold")),
          )
      )
      tmp_settings.color_vision_mode = (
          Settings._check_integrity_of_color_blindness_value(
              settings_dict.get("color_vision_mode"),
          )
      )
      tmp_settings.ask_save_pymol_session = (
          Settings._check_integrity_of_save_pymol_session_flag(
              int(settings_dict.get("ask_save_pymol_session")),
          )
      )
      tmp_settings.image_background_color = (
          Settings._check_integrity_of_bg_color_value(
              settings_dict.get("image_background_color"),
          )
      )
      tmp_settings.image_renderer = Settings._check_integrity_of_renderer_value(
          settings_dict.get("image_renderer"),
      )
      tmp_settings.image_ray_trace_mode = (
          Settings._check_integrity_of_ray_trace_mode_value(
              settings_dict.get("image_ray_trace_mode"),
          )
      )
      tmp_settings.image_ray_texture = (
          Settings._check_integrity_of_ray_texture_value(
              settings_dict.get("image_ray_texture"),
          )
      )
      tmp_settings.proteins_tab_use_toggle = (
          Settings._check_integrity_of_proteins_tab_use_toggle_flag(
              int(settings_dict.get("proteins_tab_use_toggle")),
          )
      )
      tmp_settings.proteins_tab_use_combobox_for_colors = Settings._check_integrity_of_proteins_tab_use_combobox_for_colors_flag(
          int(settings_dict.get("proteins_tab_use_combobox_for_colors")),
      )
      # integrity checks are also valid for protein pairs
      # TODO: create more abstract integrity checks
      tmp_settings.protein_pairs_tab_use_toggle = (
          Settings._check_integrity_of_proteins_tab_use_toggle_flag(
              int(settings_dict.get("protein_pairs_tab_use_toggle")),
          )
      )
      tmp_settings.protein_pairs_tab_use_combobox_for_colors = Settings._check_integrity_of_proteins_tab_use_combobox_for_colors_flag(
          int(settings_dict.get("protein_pairs_tab_use_combobox_for_colors")),
      )

    except ValueError as e:
      raise AttributeError(
          f"An error occurred during the deserialization of the settings: {e}"
      )
    return tmp_settings

  # <editor-fold desc="Integrity checks">
  @staticmethod
  def _check_integrity_of_workspace_path(a_workspace_path: str) -> str:
    """Checks integrity of the workspace path.

    Args:
        a_workspace_path: A string representing the path to the workspace directory.

    Returns:
        A string representing the valid workspace path if it passes the integrity check.

    Raises:
        exception.IllegalArgumentError: If a_workspace_path is None.
        ValueError: If a_workspace_path is not a valid path.
    """
    # <editor-fold desc="Checks">
    if a_workspace_path is None:
      logger.error("a_workspace_path is None.")
      raise exception.IllegalArgumentError("a_workspace_path is None.")

    # </editor-fold>

    if safeguard.Safeguard.check_filepath(pathlib.Path(a_workspace_path)):
      return a_workspace_path
    raise ValueError("a_workspace_path is not a valid path.")

  @staticmethod
  def _check_integrity_of_cycles(a_cycles_value: int) -> int:
    """Checks the integrity of the cycles value.

    Args:
        a_cycles_value: An integer representing the value of cycles.

    Returns:
        An integer representing the value of cycles if it is a positive number.

    Raises:
        exception.IllegalArgumentError: If a_cycles_value is None.
        ValueError: If the value of cycles is not a positive number.
    """
    # <editor-fold desc="Checks">
    if a_cycles_value is None:
      logger.error("a_cycles_value is None.")
      raise exception.IllegalArgumentError("a_cycles_value is None.")

    # </editor-fold>

    if safeguard.Safeguard.check_if_number_is_positive(a_cycles_value):
      return a_cycles_value
    raise ValueError("a_cycles_value is a negative value.")

  @staticmethod
  def _check_integrity_of_cutoff(a_cutoff_value: float) -> float:
    """Checks the integrity of the cutoff value.

    Args:
        a_cutoff_value (float): A float representing the cutoff value to be checked for integrity.

    Returns:
        The cutoff value if it passes the integrity check.

    Raises:
        exception.IllegalArgumentError: If a_cutoff_value is None.
        ValueError: If the cutoff value is not a positive number.
    """
    # <editor-fold desc="Checks">
    if a_cutoff_value is None:
      logger.error("a_cutoff_value is None.")
      raise exception.IllegalArgumentError("a_cutoff_value is None.")

    # </editor-fold>

    if safeguard.Safeguard.check_if_number_is_positive(a_cutoff_value):
      return a_cutoff_value
    raise ValueError("a_cutoff_value is a negative value.")

  @staticmethod
  def _check_integrity_of_app_start_value(an_app_start_value: int) -> int:
    """Checks the integrity of the application start value.

    Args:
        an_app_start_value (int): The application start value to check.

    Returns:
        The validated application start value.

    Raises:
        exception.IllegalArgumentError: If an_app_start_value is None.
        ValueError: If the application start value is not 0 or 1.
    """
    # <editor-fold desc="Checks">
    if an_app_start_value is None:
      logger.error("an_app_start_value is None.")
      raise exception.IllegalArgumentError("an_app_start_value is None.")

    # </editor-fold>

    if an_app_start_value == 0 or an_app_start_value == 1:
      return an_app_start_value
    raise ValueError("an_app_start_value is not 0 or 1.")

  @staticmethod
  def _check_integrity_of_wsl_install_flag(a_wsl_install_flag: int) -> int:
    """Checks the integrity of the Windows Subsystem for Linux (WSL) install flag.

    Args:
        a_wsl_install_flag (int): The WSL install flag to check.

    Returns:
        The verified WSL install flag.

    Raises:
        exception.IllegalArgumentError: If a_wsl_install_flag is None.
        ValueError: If the provided WSL install flag is not 0 or 1.
    """
    # <editor-fold desc="Checks">
    if a_wsl_install_flag is None:
      logger.error("a_wsl_install_flag is None.")
      raise exception.IllegalArgumentError("a_wsl_install_flag is None.")

    # </editor-fold>

    if a_wsl_install_flag == 0 or a_wsl_install_flag == 1:
      return a_wsl_install_flag
    raise ValueError("a_wsl_install_flag is not 0 or 1.")

  @staticmethod
  def _check_integrity_of_colabfold_install_flag(
      a_colabfold_install_flag: int,
  ) -> int:
    """Checks the integrity of the ColabFold installation flag.

    Args:
        a_colabfold_install_flag (int): The ColabFold installation flag to be checked.

    Returns:
        The validated ColabFold installation flag.

    Raises:
        exception.IllegalArgumentError: If a_colabfold_install_flag is None.
        ValueError: If the ColabFold installation flag is neither 0 nor 1.
    """
    # <editor-fold desc="Checks">
    if a_colabfold_install_flag is None:
      logger.error("a_colabfold_install_flag is None.")
      raise exception.IllegalArgumentError("a_colabfold_install_flag is None.")

    # </editor-fold>

    if a_colabfold_install_flag == 0 or a_colabfold_install_flag == 1:
      return a_colabfold_install_flag
    raise ValueError("a_colabfold_install_flag is neither 0 nor 1.")

  @staticmethod
  def _check_integrity_of_color_blindness_value(
      a_color_blindness_value: str,
  ) -> str:
    """Checks the integrity of a color blindness value.

    Args:
        a_color_blindness_value (str): The color blindness value to check.

    Returns:
        The color blindness value.

    Raises:
        exception.IllegalArgumentError: If a_color_blindness_value is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_color_blindness_value is None or a_color_blindness_value == "":
      logger.error("a_color_blindness_value is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_color_blindness_value is either None or an empty string."
      )

    # </editor-fold>

    return a_color_blindness_value

  @staticmethod
  def _check_integrity_of_save_pymol_session_flag(
      save_pymol_session_flag: int,
  ) -> int:
    """Checks the integrity of `save_pymol_session_flag`.

    Args:
        save_pymol_session_flag (int): An integer representing the save flag for a Pymol session.

    Returns:
        The save flag value passed as the input parameter.

    Raises:
        exception.IllegalArgumentError: If save_pymol_session_flag is None.
        ValueError: If the save_pymol_session_flag is not 0 or 1.
    """
    # <editor-fold desc="Checks">
    if save_pymol_session_flag is None:
      logger.error("save_pymol_session_flag is None.")
      raise exception.IllegalArgumentError("save_pymol_session_flag is None.")

    # </editor-fold>

    if save_pymol_session_flag == 0 or save_pymol_session_flag == 1:
      return save_pymol_session_flag
    raise ValueError("save_pymol_session_flag is not 0 or 1")

  @staticmethod
  def _check_integrity_of_bg_color_value(a_bg_color_value: str) -> str:
    """Checks the integrity of the background color value.

    Args:
        a_bg_color_value (str): The background color value to be checked.

    Returns:
        The validated background color value.

    Raises:
        exception.IllegalArgumentError: If a_bg_color_value is None.
        ValueError: If the background color value is not "black" or "white".
    """
    # <editor-fold desc="Checks">
    if a_bg_color_value is None:
      logger.error("a_bg_color_value is None.")
      raise exception.IllegalArgumentError("a_bg_color_value is None.")

    # </editor-fold>

    if a_bg_color_value == "white" or a_bg_color_value == "black":
      return a_bg_color_value
    raise ValueError(
        f"Invalid background color! The value is not black or white: {a_bg_color_value}"
    )

  @staticmethod
  def _check_integrity_of_renderer_value(a_renderer_value: str) -> str:
    """Checks the integrity of the given renderer value.

    Args:
        a_renderer_value (str): A string representing the renderer value.

    Returns:
        The renderer value if it is valid.

    Raises:
        exception.IllegalArgumentError: If a_renderer_value is None.
        ValueError: If the renderer value is not "0" or "-1".
    """
    # <editor-fold desc="Checks">
    if a_renderer_value is None:
      logger.error("a_renderer_value is None.")
      raise exception.IllegalArgumentError("a_renderer_value is None.")
    # </editor-fold>

    if a_renderer_value == "0" or a_renderer_value == "-1":
      return a_renderer_value
    raise ValueError(
        f"Invalid renderer! The value is not 0 or -1: {a_renderer_value}"
    )

  @staticmethod
  def _check_integrity_of_ray_trace_mode_value(
      a_ray_trace_mode_value: int,
  ) -> int:
    """Checks the integrity of `a_ray_trace_mode_value`.

    Args:
        a_ray_trace_mode_value (int): An integer representing the ray trace mode value.

    Returns:
        An integer representing of the ray trace mode value if it is valid.

    Raises:
        exception.IllegalArgumentError: If a_ray_trace_mode_value is None.
        ValueError: If the ray trace mode value is not between 0 and 3.
    """
    # <editor-fold desc="Checks">
    if a_ray_trace_mode_value is None:
      logger.error("a_ray_trace_mode_value is None.")
      raise exception.IllegalArgumentError("a_ray_trace_mode_value is None.")

    # </editor-fold>

    tmp_possible_values: list = list(range(3))
    if a_ray_trace_mode_value in tmp_possible_values:
      return a_ray_trace_mode_value
    raise ValueError(
        f"Invalid ray trace mode! The value is not between 0 and 3: {a_ray_trace_mode_value}"
    )

  @staticmethod
  def _check_integrity_of_ray_texture_value(a_ray_texture_value: int) -> int:
    """Checks the integrity of the ray texture value.

    Args:
        a_ray_texture_value (int): Integer value representing the ray texture mode.

    Returns:
        The ray texture value passed as an argument, if it is valid.

    Raises:
        exception.IllegalArgumentError: If a_ray_texture_value is None.
        ValueError: If the ray texture value is not between 0 and 5.
    """
    # <editor-fold desc="Checks">
    if a_ray_texture_value is None:
      logger.error("a_ray_texture_value is None.")
      raise exception.IllegalArgumentError("a_ray_texture_value is None.")

    # </editor-fold>

    tmp_possible_values: list = list(range(6))
    if a_ray_texture_value in tmp_possible_values:
      return a_ray_texture_value
    raise ValueError(
        f"Invalid ray trace mode! The value is not between 0 and 5: {a_ray_texture_value}"
    )

  @staticmethod
  def _check_integrity_of_proteins_tab_use_toggle_flag(
      proteins_tab_use_toggle_flag: int,
  ) -> int:
    """Checks the integrity of the `proteins_tab_use_toggle_flag` parameter.

    Args:
        proteins_tab_use_toggle_flag (int): The flag indicating whether to use the proteins tab.

    Returns:
        The `proteins_tab_use_toggle_flag` parameter if it is either 0 or 1.

    Raises:
        exception.IllegalArgumentError: If proteins_tab_use_toggle_flag is None.
        ValueError: If the `proteins_tab_use_toggle_flag` parameter is not 0 or 1.
    """
    # <editor-fold desc="Checks">
    if proteins_tab_use_toggle_flag is None:
      logger.error("proteins_tab_use_toggle_flag is None.")
      raise exception.IllegalArgumentError(
          "proteins_tab_use_toggle_flag is None."
      )

    # </editor-fold>

    if proteins_tab_use_toggle_flag == 0 or proteins_tab_use_toggle_flag == 1:
      return proteins_tab_use_toggle_flag
    raise ValueError("proteins_tab_use_toggle_flag parameter is not 0 or 1.")

  @staticmethod
  def _check_integrity_of_proteins_tab_use_combobox_for_colors_flag(
      proteins_tab_use_combobox_for_colors_flag: int,
  ) -> int:
    """Checks the integrity of proteins_tab_use_combobox_for_colors_flag.

    Args:
        proteins_tab_use_combobox_for_colors_flag (int): The flag indicating whether to use a combobox for colors in the proteins tab.

    Returns:
        The input flag if it is either 0 or 1, otherwise raises a ValueError.

    Raises:
        exception.IllegalArgumentError: If proteins_tab_use_combobox_for_colors_flag is None.
        ValueError: If proteins_tab_use_combobox_for_colors_flag is not 0 or 1.
    """
    # <editor-fold desc="Checks">
    if proteins_tab_use_combobox_for_colors_flag is None:
      logger.error("proteins_tab_use_combobox_for_colors_flag is None.")
      raise exception.IllegalArgumentError(
          "proteins_tab_use_combobox_for_colors_flag is None."
      )

    # </editor-fold>

    if (
        proteins_tab_use_combobox_for_colors_flag == 0
        or proteins_tab_use_combobox_for_colors_flag == 1
    ):
      return proteins_tab_use_combobox_for_colors_flag
    raise ValueError("proteins_tab_use_combobox_for_colors_flag is not 0 or 1.")

  # </editor-fold>

  def get_workspace_path(self) -> pathlib.Path:
    """Gets the value of the `workspace_path` attribute.

    Returns:
        The workspace path.
    """
    return self.workspace_path

  def set_workspace_path(self, value: pathlib.Path) -> None:
    """Gets the value of the `workspace_path` attribute.

    Args:
        value (pathlib.Path): The value of the new workspace path.

    Raises:
        exception.IllegalArgumentError: If `value` is None.
        exception.DirectoryNotFoundError: If `value` does not exist.
    """
    # <editor-fold desc="Checks">
    if value is None:
      logger.error("value is None.")
      raise exception.IllegalArgumentError("value is None.")
    if not os.path.exists(value):
      logger.error("value does not exist.")
      raise exception.DirectoryNotFoundError("value does not exist.")

    # </editor-fold>

    self.workspace_path = value

  def get_cycles(self) -> int:
    """Gets the value of the `cycles` attribute.

    Returns:
        The `cycles` attribute.
    """
    return int(self.cycles)

  def set_cycles(self, value: int) -> None:
    """Sets the value for the cycles attribute.

    Args:
        value (int): The value of the new cycles attribute.

    Raises:
        exception.IllegalArgumentError: If value is None.
    """
    # <editor-fold desc="Checks">
    if value is None:
      logger.error("value is None.")
      raise exception.IllegalArgumentError("value is None.")
    # </editor-fold>

    self.cycles = value

  def get_cutoff(self) -> float:
    """Gets the value of the cutoff attribute.

    Returns:
        The `cutoff` attribute.
    """
    return float(self.cutoff)

  def set_cutoff(self, value: float) -> None:
    """Sets the value for the cutoff attribute.

    Args:
        value (int): The value of the new cutoff attribute.

    Raises:
        exception.IllegalArgumentError: If value is None.
    """
    # <editor-fold desc="Checks">
    if value is None:
      logger.error("value is None.")
      raise exception.IllegalArgumentError("value is None.")
    # </editor-fold>

    self.cutoff = value

  # def get_app_launch(self) -> int:
  #     """This function gets the value of the app_launch variable.
  #
  #     Returns (int):
  #        app_launch
  #     """
  #     return int(self.app_launch)
  #
  # def set_app_launch(self, value: int) -> None:
  #     """This function gets the value of the app_launch variable."""
  #     self.app_launch = value

  def restore_settings(self, dir_settings: str, filename: str) -> None:
    """Resets the settings to the default values.

    Args:
        dir_settings (str): The directory where the settings.json should be stored.
        filename (str): The name of the settings.json.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if dir_settings is None:
      logger.error("dir_settings is None.")
      raise exception.IllegalArgumentError("dir_settings is None.")
    if filename is None:
      logger.error("filename is None.")
      raise exception.IllegalArgumentError("filename is None.")

    # </editor-fold>

    self.workspace_path = constants.DEFAULT_WORKSPACE_PATH
    self.cycles: int = 0
    self.cutoff: float = 1.0
    self.app_launch = 1
    self.dir_settings: str = dir_settings
    self.filename: str = filename
    self.wsl_install: int = 0
    self.local_colabfold: int = 0
    self.ask_save_pymol_session: int = 0
    self.image_background_color: str = "black"
    self.image_renderer: str = "-1"  # or "0"
    self.image_ray_trace_mode: int = 1  # ranges from 0 to 3
    self.image_ray_texture: int = 0  # ranges from 0 to 5
    self.proteins_tab_use_toggle: int = 1  # or "0"
    self.proteins_tab_use_combobox_for_colors: int = 0  # or "1"
    self.protein_pairs_tab_use_toggle: int = 1  # or "0"
    self.protein_pairs_tab_use_combobox_for_colors: int = 0  # or "1"

    self.serialize_settings()
