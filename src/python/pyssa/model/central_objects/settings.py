"""Module for handling everything related to the settings.json."""
import json
import logging
import pathlib
from typing import Union
from pyssa.model.preference import model_definitions
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class Settings:
  """Class for storing all application settings."""

  def __init__(self) -> None:
    """Constructor."""
    # <editor-fold desc="Instance attributes">
    self._workspace_path: pathlib.Path = model_definitions.ModelDefinitions.DEFAULT_WORKSPACE_PATH
    """The path to the current workspace"""
    self.cycles: int = 0
    """The number of cycles used in the distance analysis process"""
    self.cutoff: float = 1.0
    """The cutoff value used in the distance analysis process"""
    self.app_launch: int = 1
    """An integer flag indicating whether the application is launched for the first time or not"""
    self.wsl_install: int = 1
    """An integer flag indicating whether WSL2 is installed or not"""
    self.local_colabfold: int = 1
    """An integer flag indicating whether ColabFold is installed or not"""
    # </editor-fold>

  # <editor-fold desc="Integrity check methods">
  def check_if_path_exists(self, a_path: pathlib.Path) -> pathlib.Path:
    """Checks integrity of a given path.

    Args:
      a_path: The path to check

    Raises:
      exception.NoneValueError: If `a_path` is None.
      exception.IllegalArgumentError: If `a_path` could not be found.
    """
    # <editor-fold desc="Checks">
    if a_path is None:
      default_logging.append_to_log_file(logger, "a_path is None.", logging.ERROR)
      raise exception.NoneValueError("a_path is None.")
    if not a_path.exists():
      default_logging.append_to_log_file(logger, "a_path could not be found.", logging.ERROR)
      raise exception.IllegalArgumentError("a_path could not be found.")
    # </editor-fold>
    return a_path

  def check_if_number_is_greater_zero(self, a_number: Union[int, float]) -> Union[int, float]:
    """Checks integrity of a given number.

    Args:
      a_number: The number to check

    Raises:
      exception.NoneValueError: If `a_number` is None.
      exception.IllegalArgumentError: If `a_number` has a value below zero.
    """
    # <editor-fold desc="Checks">
    if a_number is None:
      default_logging.append_to_log_file(logger, "a_number is None.", logging.ERROR)
      raise exception.NoneValueError("a_number is None.")
    if a_number < 0:
      default_logging.append_to_log_file(logger, "a_number has a value below zero.", logging.ERROR)
      raise exception.IllegalArgumentError("a_number has a value below zero.")
    # </editor-fold>
    return a_number

  def check_if_number_is_between_one_and_zero(self, a_number: Union[int, float]) -> Union[int, float]:
    """Checks if a given number is between 1 and 0.

    Args:
      a_number: The number to check

    Raises:
      exception.NoneValueError: If `a_number` is None.
      exception.IllegalArgumentError: If `a_number` has a value below zero.
    """
    # <editor-fold desc="Checks">
    if a_number is None:
      default_logging.append_to_log_file(logger, "a_number is None.", logging.ERROR)
      raise exception.NoneValueError("a_number is None.")
    if a_number < 0 or a_number > 1:
      default_logging.append_to_log_file(logger, "a_number has either a value below zero or above one.", logging.ERROR)
      raise exception.IllegalArgumentError("a_number has either a value below zero or above one")
    # </editor-fold>
    return a_number
  # </editor-fold>

  # <editor-fold desc="Getter & setter">
  def get_workspace_path(self) -> pathlib.Path:
    """Gets the value of the `workspace_path` attribute.

    Returns:
        The workspace path.
    """
    return self._workspace_path

  def set_workspace_path(self, a_new_workspace_path: pathlib.Path) -> None:
    """Gets the value of the `workspace_path` instance attribute.

    Args:
      a_new_workspace_path (pathlib.Path): The new workspace path

    Raises:
      exception.NoneValueError: If `a_new_workspace_path` is None.
      exception.DirectoryNotFoundError: If a_new_workspace_path does not exist.
    """
    # <editor-fold desc="Checks">
    if a_new_workspace_path is None:
      logger.error("a_new_workspace_path is None.")
      raise exception.NoneValueError("a_new_workspace_path is None.")
    if not a_new_workspace_path.exists():
      logger.error("a_new_workspace_path does not exist.")
      raise exception.DirectoryNotFoundError("a_new_workspace_path does not exist.")
    # </editor-fold>
    self._workspace_path = a_new_workspace_path

  def get_cycles(self) -> int:
    """Gets the value of the `cycles` instance attribute.

    Returns:
      The `cycles` attribute.
    """
    return int(self.cycles)

  def set_cycles(self, a_value: int) -> None:
    """Sets the value for the `cycles` instance attribute.

    Args:
      a_value (int): The value of the new `cycles` instance attribute.

    Raises:
      exception.NoneValueError: If `a_value` is None.
    """
    # <editor-fold desc="Checks">
    if a_value is None:
      logger.error("a_value is None.")
      raise exception.NoneValueError("a_value is None.")
    # </editor-fold>
    self.cycles = a_value

  def get_cutoff(self) -> float:
    """Gets the value of the `cutoff` instance attribute.

    Returns:
      The `cutoff` attribute.
    """
    return float(self.cutoff)

  def set_cutoff(self, a_value: float) -> None:
    """Sets the value for the `cutoff` instance attribute.

    Args:
      a_value (int): The value of the new `cutoff` instance attribute.

    Raises:
      exception.NoneValueError: If `a_value` is None.
    """
    # <editor-fold desc="Checks">
    if a_value is None:
      logger.error("a_value is None.")
      raise exception.NoneValueError("a_value is None.")
    # </editor-fold>
    self.cutoff = a_value

  def get_app_launch(self) -> int:
    """Gets the value of the `app_launch` instance attribute.

    Returns:
     The `app_launch` attribute.
    """
    return int(self.app_launch)

  def set_app_launch(self, a_value: int) -> None:
    """Gets the value of the app_launch instance attribute.

    Args:
      a_value (int): The value of the new `app_launch` attribute.

    Raises:
      exception.NoneValueError: If `a_value` is None.
    """
    # <editor-fold desc="Checks">
    if a_value is None:
      logger.error("a_value is None.")
      raise exception.NoneValueError("a_value is None.")
    # </editor-fold>
    self.app_launch = a_value
  # </editor-fold>

  def serialize_settings(self) -> None:
    """Serializes the settings object to a json file."""
    model_definitions.ModelDefinitions.DEFAULT_SETTINGS_PATH.mkdir(parents=True, exist_ok=True)
    model_definitions.ModelDefinitions.DEFAULT_SETTINGS_FILEPATH.unlink(missing_ok=True)
    settings_dict = self.__dict__
    settings_dict.update({
      "workspace_path": str(self._workspace_path)
    })
    with open(
            model_definitions.ModelDefinitions.DEFAULT_SETTINGS_FILEPATH,
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
    if not model_definitions.ModelDefinitions.DEFAULT_SETTINGS_FILEPATH.exists():
      raise FileNotFoundError(f"File: {model_definitions.ModelDefinitions.DEFAULT_SETTINGS_FILEPATH} can not be found!")
    with open(
            model_definitions.ModelDefinitions.DEFAULT_SETTINGS_FILEPATH,
            "r",
            encoding="utf-8",
    ) as settings_obj_file:
      settings_dict = json.load(settings_obj_file)
    try:
      tmp_settings: Settings = Settings()
      tmp_settings._workspace_path = tmp_settings.check_if_path_exists(
        pathlib.Path(settings_dict.get("workspace_path"))
      )
      tmp_settings.cycles = tmp_settings.check_if_number_is_greater_zero(
        int(settings_dict.get("cycles"))
      )
      tmp_settings.cutoff = tmp_settings.check_if_number_is_greater_zero(
        float(settings_dict.get("cutoff"))
      )
      tmp_settings.app_launch = tmp_settings.check_if_number_is_between_one_and_zero(
        int(settings_dict.get("app_launch"))
      )
      tmp_settings.wsl_install = tmp_settings.check_if_number_is_between_one_and_zero(
        int(settings_dict.get("wsl_install"))
      )
      tmp_settings.local_colabfold = tmp_settings.check_if_number_is_between_one_and_zero(
        int(settings_dict.get("local_colabfold"))
      )
    except exception.IllegalArgumentError as e:
      raise AttributeError(
        f"An error occurred during the deserialization of the settings: {e}"
      )
    return tmp_settings

  def restore_settings(self) -> None:
    """Resets the settings to the default values."""
    self._workspace_path: pathlib.Path = model_definitions.ModelDefinitions.DEFAULT_WORKSPACE_PATH
    self.cycles: int = 0
    self.cutoff: float = 1.0
    self.app_launch: int = 1
    self.wsl_install: int = 0
    self.local_colabfold: int = 0
    self.serialize_settings()
