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
"""Module for the delete project view controller."""
import logging
import os

from src.pyssa.gui import app_state
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import Qt
from src.pyssa.gui.ui.views import delete_project_view

from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.util import constants, enums, ui_util, exception
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class DeleteProjectViewController(QtCore.QObject):
  """Class for the DeleteProjectViewController."""

  def __init__(
      self, the_app_state: "app_state.AppState"
  ) -> None:
    """Constructor.

    Args:
        the_app_state (app_state.AppState): The AppState object.

    Raises:
        exception.IllegalArgumentError: If `the_app_state` is None.
    """
    # <editor-fold desc="Checks">
    if the_app_state is None:
      logger.error("the_app_state is None.")
      raise exception.IllegalArgumentError("the_app_state is None.")

    # </editor-fold>
    super().__init__()
    self._app_state = the_app_state
    self._view = delete_project_view.DeleteProjectView()
    self._fill_projects_list_view()
    self._connect_all_ui_elements_to_slot_functions()
    self.restore_default_view()

  def get_view(self):
    return self._view

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )

  def restore_default_view(self) -> None:
    """Restores the default UI."""
    self._view.ui.label_31.hide()
    self._view.ui.txt_delete_search.setPlaceholderText("Search")
    self._view.ui.txt_delete_search.clear()
    self._view.ui.txt_delete_selected_projects.clear()
    self._view.ui.btn_delete_delete_project.setEnabled(False)

  def _fill_projects_list_view(self) -> None:
    """Lists all projects."""
    self._view.ui.list_delete_projects_view.setModel(
        self._app_state.workspace.get_model()
    )

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.txt_delete_search.textChanged.connect(
        self.validate_delete_search
    )
    self._view.ui.list_delete_projects_view.clicked.connect(
        self.select_project_from_delete_list
    )
    self._view.ui.txt_delete_selected_projects.textChanged.connect(
        self.activate_delete_button
    )
    self._view.ui.btn_delete_delete_project.clicked.connect(self.delete_project)
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

  def validate_delete_search(self) -> None:
    """Validates the input of the project name in real-time."""
    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
    projects_list_view = self._view.ui.list_delete_projects_view
    ui_util.select_matching_string_in_q_list_view(
        self._view.ui.txt_delete_search.text(),
        projects_list_view,
        self._view.ui.txt_delete_selected_projects,
    )

  def select_project_from_delete_list(self) -> None:
    """Selects a project from the project list on the delete page."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A project from the list of existing projects was clicked.",
    )
    try:
      if (
          len(
              self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()
          )
          == 1
      ):
        self._view.ui.txt_delete_selected_projects.setText(
            self._view.ui.list_delete_projects_view.model().data(
                self._view.ui.list_delete_projects_view.currentIndex(),
                Qt.DisplayRole,
            )
        )
      elif (
          len(
              self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()
          )
          > 1
      ):
        self._view.ui.txt_delete_selected_projects.setText(
            "Multiple projects selected."
        )
      else:
        self._view.ui.txt_delete_selected_projects.setText("")
    except AttributeError:
      self._view.ui.txt_delete_selected_projects.setText("")

  def activate_delete_button(self) -> None:
    """Activates the delete button."""
    if self._view.ui.txt_delete_selected_projects.text() == "":
      self._view.ui.btn_delete_delete_project.setEnabled(False)
    else:
      self._view.ui.btn_delete_delete_project.setEnabled(True)

  def delete_project(self) -> None:
    """Deletes an existing project."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete' button was clicked."
    )

    selected_indexes = self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()
    if not selected_indexes:
      return

    if len(selected_indexes) > 1:
      tmp_dialog = custom_message_box.CustomMessageBoxDelete(
          "Are you sure you want to delete these projects?",
          "Delete Projects",
          custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec()
      if not tmp_dialog.response:
        return
      
      names_to_delete = []
      for idx in selected_indexes:
        names_to_delete.append(self._view.ui.list_delete_projects_view.model().data(idx, Qt.DisplayRole))
        
      for name in names_to_delete:
        self._app_state.workspace.delete_project(name)

    elif len(selected_indexes) == 1:
      tmp_dialog = custom_message_box.CustomMessageBoxDelete(
          "Are you sure you want to delete this project?",
          "Delete Project",
          custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec()
      if not tmp_dialog.response:
        return
        
      project_name = self._view.ui.list_delete_projects_view.model().data(selected_indexes[0], Qt.DisplayRole)
      self._app_state.workspace.delete_project(project_name)

    self.restore_default_view()
