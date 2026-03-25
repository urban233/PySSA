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
"""Module for the rename protein view controller."""
import logging
from src.pyssa.gui.qt import QtCore

from typing import TYPE_CHECKING
if TYPE_CHECKING:
  from src.pyssa.gui import app_state
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class RenameProteinViewController(QtCore.QObject):
  """Class for the RenameProteinViewController."""

  def __init__(
      self, 
      the_app_state: "app_state.AppState", 
      a_protein: "protein.Protein",
      a_parent=None
  ) -> None:
    """Constructor.

    Args:
        the_app_state (app_state.AppState): The AppState object.
        a_parent: Parent widget to pass to the view.

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
    self._protein = a_protein
    from src.pyssa.gui.ui.views import rename_protein_view
    self._view = rename_protein_view.RenameProteinView(a_parent)
    self._connect_all_ui_elements_to_slot_functions()

  def get_view(self):
    return self._view

  def restore_ui(self) -> None:
    """Restores the UI."""
    self._view.ui.le_name.clear()
    self._view.ui.lbl_status.setText("")

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_rename.clicked.connect(self._rename_protein)

  def _rename_protein(self) -> None:
    """Renames a protein and updates the state."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Rename' button was clicked."
    )
    new_name = self._view.ui.le_name.text()
    if not new_name:
      return

    self._view.close()

    from src.pyssa.internal.thread.thread_api import thread_runtime

    def rename_protein(progress_callback, is_cancelled):
      # Update the database
      old_name = self._protein.get_molecule_object()
      self._app_state.hot_db.update_protein_name(new_name, old_name, self._protein.get_id())
      
      # Update PyMOL Session mapping (trigger update)
      self._app_state.hot_db.update_protein_session(self._protein.get_id(), "")
      
      # Update the protein in memory
      self._protein.set_molecule_object(new_name)
      return

    def on_success(result):
      self._app_state._on_state_changed()

    def on_error(exc):
      logger.error(f"Error during protein rename: {exc}")
      from src.pyssa.gui.qt import QtWidgets
      QtWidgets.QMessageBox.critical(
        self._view,
        "Rename Failed",
        f"Could not rename protein:\n{exc}",
      )

    (
      thread_runtime.get_singleton_thread_runtime()
      .run(rename_protein)
      .on_success(on_success)
      .on_error(on_error)
      .start()
    )
