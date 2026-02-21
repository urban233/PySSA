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
"""Module for the rename sequence view controller."""
import logging
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import Qt

from typing import TYPE_CHECKING
if TYPE_CHECKING:
  from src.pyssa.gui import app_state
from src.pyssa.util import input_validator, exception
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class RenameSequenceViewController(QtCore.QObject):
  """Class for the RenameSequenceViewController."""

  def __init__(
      self, 
      the_app_state: "app_state.AppState", 
      a_sequence: "sequence.Sequence",
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
    self._sequence = a_sequence
    from src.pyssa.gui.ui.views import rename_sequence_view
    self._view = rename_sequence_view.RenameSequenceView(a_parent)
    self._sequence_names = self._convert_sequence_model_into_set()
    self._connect_all_ui_elements_to_slot_functions()

  def get_view(self):
    return self._view

  def restore_ui(self) -> None:
    """Restores the UI."""
    self._view.ui.le_name.clear()
    self._view.ui.lbl_status.setText("")
    self._view.ui.le_name.setStyleSheet(
        """QTextEdit {color: #000000; border-color: #DCDBE3;}""",
    )

  def _convert_sequence_model_into_set(self) -> set:
    """Converts the sequence model into a set of sequence names.

    Returns:
        A set of sequence names.
    """
    tmp_sequence_names = []
    if self._app_state.has_open_project():
      for seq in self._app_state.project.sequences:
        tmp_sequence_names.append(seq.name)
    return set(tmp_sequence_names)

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_rename.clicked.connect(self._rename_sequence)
    self._view.ui.le_name.textChanged.connect(self._validate_protein_name)

  def _validate_protein_name(self, the_entered_text: str) -> None:
    """Validates the input of the protein name in real-time.

    Args:
        the_entered_text (str): The input text to be validated.

    Raises:
        exception.IllegalArgumentError: If `the_entered_text` is None.
    """
    # <editor-fold desc="Checks">
    if the_entered_text is None:
      logger.error("the_entered_text is None.")
      raise exception.IllegalArgumentError("the_entered_text is None.")

    # </editor-fold>

    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
    tmp_input_validator = input_validator.InputValidator(self._view.ui.le_name)
    tmp_validate_flag, tmp_message = (
        tmp_input_validator.validate_input_for_sequence_name(
            the_entered_text,
            self._sequence_names,
        )
    )
    if tmp_validate_flag:
      self._view.ui.lbl_status.setText("")
      self._view.ui.btn_rename.setEnabled(True)
    else:
      self._view.ui.lbl_status.setText(tmp_message)
      self._view.ui.btn_rename.setEnabled(False)

  def _rename_sequence(self) -> None:
    """Renames the sequence and updates the state."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Rename' button was clicked."
    )
    new_name = self._view.ui.le_name.text()
    self._view.close()

    from src.pyssa.internal.thread.thread_api import thread_runtime
    import copy

    def rename_sequence(progress_callback, is_cancelled):
      # Update the database
      self._app_state.hot_db.update_sequence_name(new_name, self._sequence.name, self._sequence.seq)
      
      # Update the sequence in memory
      tmp_project = copy.deepcopy(self._app_state.project)
      for seq in tmp_project.sequences:
        if seq.name == self._sequence.name:
          seq.name = new_name
      self._sequence.name = new_name
      return tmp_project

    def on_success(tmp_project):
      self._app_state.project = tmp_project
      # Trigger UI refresh
      self._app_state._on_state_changed()

    def on_error(exc):
      logger.error(f"Error during sequence rename: {exc}")
      from src.pyssa.gui.qt import QtWidgets
      QtWidgets.QMessageBox.critical(
        self._view,
        "Rename Failed",
        f"Could not rename sequence:\n{exc}",
      )

    (
      thread_runtime.get_singleton_thread_runtime()
      .run(rename_sequence)
      .on_success(on_success)
      .on_error(on_error)
    )
