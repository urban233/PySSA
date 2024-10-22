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
"""Module for the open project view controller."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from src.pyssa.controller import interface_manager
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import util_async
from src.pyssa.util import constants, ui_util, exception
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class OpenProjectViewController(QtCore.QObject):
  """Class for the OpenProjectViewController."""

  return_value = QtCore.pyqtSignal(str)
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
    self._view = the_interface_manager.get_open_view()
    self._fill_projects_list_view()
    self._project_names = self._convert_model_into_set()
    self._connect_all_ui_elements_to_slot_functions()
    self.restore_default_view()

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    self._interface_manager.help_manager.open_open_project_page()

  def restore_default_view(self) -> None:
    """Restores the default UI."""
    self._view.ui.label_28.hide()
    self._view.ui.txt_open_search.setPlaceholderText("Search")
    self._view.ui.txt_open_search.clear()
    self._view.ui.txt_open_selected_project.clear()
    self._view.ui.btn_open_project.setEnabled(False)

  def _fill_projects_list_view(self) -> None:
    """Lists all projects."""
    self._view.ui.projects_list_view.setModel(
        self._interface_manager.get_workspace_model()
    )

  def _convert_model_into_set(self) -> set:
    """Converts the model into a set of project names.

    Returns:
        A set containing project names.
    """
    tmp_project_names = []
    for tmp_row in range(self._view.ui.projects_list_view.model().rowCount()):
      tmp_project_names.append(
          self._view.ui.projects_list_view.model()
          .index(tmp_row, 0)
          .data(Qt.DisplayRole),
      )
    return set(tmp_project_names)

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.txt_open_search.textChanged.connect(
        self._validate_open_search
    )
    self._view.ui.projects_list_view.clicked.connect(
        self._select_project_from_open_list
    )
    self._view.ui.txt_open_selected_project.textChanged.connect(
        self._activate_open_button
    )
    self._view.ui.btn_open_project.clicked.connect(self._open_selected_project)
    self._view.ui.projects_list_view.doubleClicked.connect(
        self._open_selected_project
    )
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

  def _validate_open_search(self, the_entered_text: str) -> None:
    """Validates the input of the project name in real-time.

    Args:
        the_entered_text (str): The text entered by the user for the open search.

    Raises:
        exception.IllegalArgumentError: If `the_entered_text` is None.
    """
    # <editor-fold desc="Checks">
    if the_entered_text is None:
      logger.error("the_entered_text is None.")
      raise exception.IllegalArgumentError("the_entered_text is None.")

    # </editor-fold>

    ui_util.select_matching_string_in_q_list_view(
        self._view.ui.txt_open_search.text(),
        self._view.ui.projects_list_view,
        self._view.ui.txt_open_selected_project,
    )

  def _select_project_from_open_list(self) -> None:
    """Sets the selected project name in the text box."""
    tmp_project_name = self._view.ui.projects_list_view.model().data(
        self._view.ui.projects_list_view.currentIndex(),
        Qt.DisplayRole,
    )
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        f"The project '{tmp_project_name}' from the list was clicked.",
    )
    try:
      self._view.ui.txt_open_selected_project.setText(tmp_project_name)
    except AttributeError:
      self._view.ui.txt_open_selected_project.setText("")

  def _activate_open_button(self) -> None:
    """Activates the open button."""
    if self._view.ui.txt_open_selected_project.text() == "":
      self._view.ui.btn_open_project.setEnabled(False)
      # styles.color_button_not_ready(self._view.ui.btn_open_project)
    else:
      self._view.ui.btn_open_project.setEnabled(True)
      # styles.color_button_ready(self._view.ui.btn_open_project)

  def _open_selected_project(self) -> None:
    """Opens the selected project by sending the `return_value` signal and closing the dialog."""
    tmp_project_name = self._view.ui.projects_list_view.model().data(
        self._view.ui.projects_list_view.currentIndex(),
        Qt.DisplayRole,
    )
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        f"The project '{tmp_project_name}' from the list was double-clicked or the 'Open' button was clicked.",
    )
    self._view.close()
    self.return_value.emit(self._view.ui.txt_open_selected_project.text())
