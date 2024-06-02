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
"""Module for the application process manager class."""
import os
import subprocess
import time
from typing import Callable, Optional

from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.util import constants, exception

__docformat__ = "google"


class ApplicationProcessManager:
  """Manages the application's processes, including the initiation, checking, and closing of the PyMOL and PySSA processes.

  Instances of this class provide an interface for managing these processes.
  """

  # <editor-fold desc="Class attributes">
  _reset_pymol_session_func: Optional[Callable] = None
  """A function that can reset the pymol session."""

  pymol_process: Optional[subprocess.Popen] = None
  """The pymol process object."""

  _should_exit: bool = False
  """A boolean flag that determines if the PyMOL process should exit."""

  _is_crashed: bool = False
  """A boolean flag that determines if the PyMOL process crashed."""

  # </editor-fold>

  def __init__(self, the_reset_pymol_session_func: Callable) -> None:
    """Constructor.

    Args:
        the_reset_pymol_session_func (Callable): The function used to reset the PyMOL session.

    Raises:
        exception.IllegalArgumentError: If the_reset_pymol_session_func is None.
    """
    # <editor-fold desc="Checks">
    if the_reset_pymol_session_func is None:
      raise exception.IllegalArgumentError(
          "the_reset_pymol_session_func is None."
      )

    # </editor-fold>

    self._reset_pymol_session_func = the_reset_pymol_session_func
    # self.pymol_process = None
    # self._should_exit = False
    # self._is_crashed = False

  def arrange_windows(self) -> None:
    """Runs a script to arrange the PySSA and PyMOL window on the screen.

    If the script for arranging the windows does not exist, a custom message box is displayed
    with an error message. If the script exists, it is executed using subprocess.Popen(), which
    runs the script in a separate process without displaying a console window.
    """
    if not os.path.exists(constants.ARRANGE_WINDOWS_EXE_FILEPATH):
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "The script for arranging the windows could not be found!",
          "Arrange Windows",
          custom_message_box.CustomMessageBoxIcons.ERROR.value,
      )
      tmp_dialog.exec_()
    else:
      subprocess.Popen(
          [constants.ARRANGE_WINDOWS_EXE_FILEPATH],
          creationflags=subprocess.CREATE_NO_WINDOW,
      )

  def start_pymol(self) -> None:
    """Starts PyMOL application."""
    self.pymol_process = subprocess.Popen(
        [f"{constants.PLUGIN_PATH}\\scripts\\batch\\start_pymol.bat"],
        creationflags=subprocess.CREATE_NO_WINDOW,
    )
    if self.pymol_process.poll() is None:
      print("PyMOL from ApplicationProcessManager class started correctly.")
      self._is_crashed = False
    else:
      print("PyMOL failed to start.")
      self._is_crashed = True

  def close_manager(self) -> None:
    """Closes the manager.

    This method sets the internal flag `_should_exit` to `True`, indicating that the manager should exit.
    """
    self._should_exit = True

  def pymol_closed(self) -> bool:
    """Gets the `_should_exit` flag, indicating that the manager should exit.

    Returns:
        A boolean that is True if the application should be closed, and False otherwise.
    """
    return self._should_exit

  def pymol_crashed(self) -> bool:
    """Returns the class attribute `_is_crashed`.

    Returns:
        A boolean that is True if the PyMOL process has crashed, False otherwise.
    """
    return self._is_crashed

  def check_process(
      self, placeholder_1: int, placeholder_2: int
  ) -> tuple[str, str]:
    """Checks the status of a PyMOL process, looped until the process exists or crashes.

    The method also initiates the start of the PyMOL.

    Args:
        placeholder_1 (int): The Placeholder that is needed for use with LegacyTask.
        placeholder_2 (int): The Placeholder that is needed for use with LegacyTask.

    Returns:
        A tuple containing two empty strings needed for the task class.
    """
    if self.pymol_process is None:
      self.start_pymol()
    while self._should_exit is False and self._is_crashed is False:
      if self.pymol_process.poll() is not None:
        print("PyMOL crashed!")  # TODO: change to logger message
        self._is_crashed = True
      else:
        time.sleep(2)
    print("Closing check_process method.")  # TODO: change to logger message
    return "", ""  # These two empty strings are needed for the task class
