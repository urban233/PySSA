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
"""Module for the add scene view controller."""
import logging

from PyQt5 import QtCore

from src.pyssa.controller import interface_manager
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class AddSceneViewController(QtCore.QObject):
  """Class for the AddSceneViewController."""

  user_input = QtCore.pyqtSignal(tuple)
  """Singal used to transfer data back to the previous window."""

  def __init__(
      self, the_interface_manager: "interface_manager.InterfaceManager"
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): The InterfaceManager object.

    Raises:
        exception.IllegalArgumentError: If `the_interface_manager` is None.
    """
    # <editor-fold desc="Checks">
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")

    # </editor-fold>

    super().__init__()
    self._interface_manager = the_interface_manager
    self._view = the_interface_manager.get_add_scene_view()
    self._interface_manager.pymol_session_manager.get_all_scenes_in_current_session()
    self._all_current_scenes = (
        self._interface_manager.pymol_session_manager.all_scenes
    )
    self._view.lbl_status.setStyleSheet("color: #ba1a1a; font-size: 11px;")
    self._connect_all_ui_elements_to_slot_functions()

  def restore_ui(self) -> None:
    """Restores the UI."""
    self._view.line_edit_scene_name.clear()
    self._view.line_edit_scene_name.setStyleSheet(
        """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
    )
    self._view.lbl_status.setText("")

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.line_edit_scene_name.textChanged.connect(
        self._validate_scene_name
    )
    self._view.btn_add_scene.clicked.connect(self._add_scene)

  def _add_scene(self) -> None:
    """Adds a new scene based on the user input."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked."
    )
    self._view.close()
    self.user_input.emit((self._view.line_edit_scene_name.text(), True))

  def _validate_scene_name(self, text: str) -> None:
    """Validates the scene name entered by the user.

    Args:
        text (str): The input text to be validated.

    Raises:
        exception.IllegalArgumentError: If `text` is None.
    """
    # <editor-fold desc="Checks">
    if text is None:
      logger.error("text is None.")
      raise exception.IllegalArgumentError("text is None.")

    # </editor-fold>

    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
    new_text = "".join(char for char in text)
    self._view.line_edit_scene_name.setText(new_text)
    if new_text in self._all_current_scenes:
      self._view.btn_add_scene.setEnabled(False)
      self._view.lbl_status.setText(
          "This scene name already exists. Please enter another name."
      )
      self._view.line_edit_scene_name.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
    else:
      self._view.btn_add_scene.setEnabled(True)
      self._view.lbl_status.setText("")
      self._view.line_edit_scene_name.setStyleSheet(
          """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
      )
