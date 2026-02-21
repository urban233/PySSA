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
"""Module for the welcome screen view controller."""
import logging

from src.pyssa.controller import create_project_view_controller
from src.pyssa.gui import app_state, main_window
from src.pyssa.gui.qt import QtCore, QtWidgets
from src.pyssa.gui.qt import Qt
from src.pyssa.gui.ui.views import welcome_screen_view
from src.pyssa.model import psa_objects_model
from src.pyssa.util import ui_util, exception
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.internal.thread.thread_api import thread_runtime
from src.pyssa.io_pyssa.db_pyssa import ProjectDatabase

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class WelcomeScreenViewController(QtCore.QObject):
  """Controller for the Open Project dialog."""

  def __init__(
          self,
          the_main_window: "main_window.MainWindow",
          the_app_state: "app_state.AppState"
  ) -> None:
    """Constructor.

    Args:
        the_app_state: The AppState instance.
    """
    super().__init__()
    self._app_state = the_app_state
    self._view = welcome_screen_view.WelcomeScreenView(the_main_window)
    self._fill_projects_list_view()
    self._project_names = self._convert_model_into_set()
    self._connect_all_ui_elements_to_slot_functions()
    self.restore_default_view()

  def get_view(self):
    return self._view

  # ------------------------------------------------------------------
  # Setup
  # ------------------------------------------------------------------

  def restore_default_view(self) -> None:
    """Restores the default UI."""
    self._set_ui_loading(False)
    self._view.ui.btn_open_project.setEnabled(False)

  def _fill_projects_list_view(self) -> None:
    """Lists all projects from the workspace model."""
    self._view.ui.projects_list_view.setModel(
      self._app_state.workspace.get_model()
    )

  def _convert_model_into_set(self) -> set:
    """Converts the list model into a set of project name strings."""
    model = self._view.ui.projects_list_view.model()
    return {
      model.index(row, 0).data(Qt.DisplayRole)
      for row in range(model.rowCount())
    }

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their slot functions."""
    self._view.ui.projects_list_view.doubleClicked.connect(self._open_selected_project)
    self._view.ui.projects_list_view.clicked.connect(self._activate_open_button)
    self._view.ui.btn_new_project.clicked.connect(self._create_new_project)
    self._view.ui.btn_open_project.clicked.connect(self._open_selected_project)
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

  # ------------------------------------------------------------------
  # Slot functions
  # ------------------------------------------------------------------

  def _open_help_for_dialog(self) -> None:
    """Opens the help page for this dialog."""
    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked.")

  def _activate_open_button(self) -> None:
    """Enables the Open button when a project name is present."""
    has_selection = bool(self._view.ui.projects_list_view.selectionModel().hasSelection())
    self._view.ui.btn_open_project.setEnabled(has_selection)

  def _open_selected_project(self) -> None:
    """Starts the async project load and closes the dialog on success."""
    project_name = self._view.ui.projects_list_view.currentIndex().data(Qt.DisplayRole)
    db_path = str(self._app_state.workspace.construct_project_db_path(project_name))
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      f"Opening project '{project_name}'.",
    )
    self._set_ui_loading(True)

    def load_project(progress_callback, is_cancelled):
      tmp_db = ProjectDatabase(db_path=db_path, project_id=project_name)
      if is_cancelled():
        tmp_db.close()
        raise InterruptedError("Cancelled before loading project data.")
      tmp_project = tmp_db.service.load_project(
        project_name=project_name,
        workspace_path=db_path,
        app_settings=self._app_state.get_settings(),
        progress_signal=_ProgressCallable(progress_callback),
      )
      if is_cancelled():
        tmp_db.close()
        raise InterruptedError("Cancelled after loading project data.")
      tmp_pyssa_objects_model = psa_objects_model.PSAObjectsModel()
      tmp_pyssa_objects_model.build_model(tmp_project)
      return tmp_project, tmp_db, tmp_pyssa_objects_model

    def on_success(result):
      tmp_project, tmp_db, tmp_pyssa_objects_model = result
      self._app_state.pyssa_objects_model = tmp_pyssa_objects_model
      self._app_state.open_project(tmp_project, tmp_db)
      self._app_state.status_bar_manager.show_permanent_message("", False)
      self._app_state.status_bar_manager.show_temporary_message("Project loaded.")

    def on_error(exc):
      logger.exception("Failed to open project.", exc_info=exc)
      self._set_ui_loading(False)
      from src.pyssa.gui.qt import QtWidgets
      QtWidgets.QMessageBox.critical(
        self._view,
        "Failed to open project",
        f"Could not open the project:\n{exc}",
      )
      self._app_state.status_bar_manager.show_error_message(
        "Failed to open project", True
      )

    (
      thread_runtime.get_singleton_thread_runtime()
      .run(load_project)
      .on_success(on_success)
      .on_error(on_error)
    )
    self._view.close()
    self._app_state.status_bar_manager.show_permanent_message(
      "Loading project ...", True
    )

  def _create_new_project(self):
    self._create_project_controller = create_project_view_controller.CreateProjectViewController(
      self._app_state, self._view
    )
    self._create_project_controller.restore_default_view()
    self._create_project_controller.get_view().exec()
    # TODO: Fix the issue with the cancel of the create project dialog
    # also closes the welcome screen
    self._view.close()

  # ------------------------------------------------------------------
  # UI helpers
  # ------------------------------------------------------------------

  def _set_ui_loading(self, loading: bool) -> None:
    self._view.ui.btn_open_project.setEnabled(not loading)
    self._view.ui.projects_list_view.setEnabled(not loading)


class _ProgressCallable:
  """Minimal adapter so Worker's progress_callback satisfies the
  emit_signal(msg, pct) interface expected by ProjectService.load_project.

  Once ProjectService is updated to accept a plain callable, this class
  can be deleted and progress_callback passed directly.
  """

  __slots__ = ("_callback",)

  def __init__(self, callback) -> None:
    self._callback = callback

  def emit_signal(self, msg: str, pct: int) -> None:
    self._callback(pct)
