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
"""Module for the add sequence view controller."""
import logging
import os
import pathlib
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui.ui.styles import styles
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import util_async
from src.pyssa.io_pyssa import safeguard
from src.pyssa.util import input_validator, constants, exception
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class AddSequenceViewController(QtCore.QObject):
  """Class for the AddSequenceViewController."""

  return_value = QtCore.pyqtSignal(tuple)
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
    self._view = the_interface_manager.get_add_sequence_view()

    self._view.ui.lbl_status.setStyleSheet("color: #ba1a1a; font-size: 11px;")
    self._view.ui.lbl_status_seq.setStyleSheet(
        "color: #ba1a1a; font-size: 11px;"
    )
    self._sequence_names = self._convert_sequence_model_into_set()
    self._connect_all_ui_elements_to_slot_functions()
    self.restore_default_view()

  # <editor-fold desc="Util methods">
  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    self._interface_manager.help_manager.open_sequence_add_page()

  # </editor-fold>

  def _convert_sequence_model_into_set(self) -> set:
    """Converts the sequence model into a set of sequence names.

    Returns:
        A set of sequence names.
    """
    tmp_sequence_names = []
    for tmp_row in range(
        self._interface_manager.get_main_view()
        .ui.seqs_list_view.model()
        .rowCount()
    ):
      tmp_sequence_names.append(
          self._interface_manager.get_main_view()
          .ui.seqs_list_view.model()
          .index(tmp_row, 0)
          .data(Qt.DisplayRole),
      )
    return set(tmp_sequence_names)

  def restore_default_view(self) -> None:
    """Restores the default UI."""
    self._view.ui.le_seq_name.clear()
    self._view.ui.le_protein_seq.clear()
    self._view.ui.lbl_status.setText("")
    self._view.ui.lbl_status_seq.setText("")
    self._switch_ui_to_sequence_name_input()
    self._view.ui.btn_next.setEnabled(False)
    self._view.ui.btn_add.setEnabled(False)
    self._view.ui.le_protein_seq.setStyleSheet(
        """QTextEdit {color: #000000; border-color: #DCDBE3;}""",
    )
    self._view.ui.le_seq_name.setStyleSheet(
        """QTextEdit {color: #000000; border-color: #DCDBE3;}""",
    )

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.le_seq_name.textChanged.connect(self._validate_protein_name)
    self._view.ui.le_protein_seq.textChanged.connect(
        self._validate_protein_sequence
    )
    self._view.ui.btn_next.clicked.connect(self._switch_ui_to_sequence_input)
    self._view.ui.btn_back.clicked.connect(
        self._switch_ui_to_sequence_name_input
    )
    self._view.ui.btn_add.clicked.connect(self._add_sequence)
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

  def _validate_protein_name(self, the_entered_text: str) -> None:
    """Validates the input of the protein name in real-time.

    Args:
        the_entered_text (str): A string representing the text entered by the user.

    Raises:
        exception.IllegalArgumentError: If `the_entered_text` is None.
    """
    # <editor-fold desc="Checks">
    if the_entered_text is None:
      logger.error("the_entered_text is None.")
      raise exception.IllegalArgumentError("the_entered_text is None.")

    # </editor-fold>

    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
    tmp_input_validator = input_validator.InputValidator(
        self._view.ui.le_seq_name
    )
    tmp_validate_flag, tmp_message = (
        tmp_input_validator.validate_input_for_sequence_name(
            the_entered_text,
            self._sequence_names,
        )
    )
    if tmp_validate_flag:
      self._view.ui.lbl_status.setText("")
      self._view.ui.btn_next.setEnabled(True)
    else:
      self._view.ui.lbl_status.setText(tmp_message)
      self._view.ui.btn_next.setEnabled(False)

  def _validate_protein_sequence(self) -> None:
    """Validates the input of the protein sequence in real-time."""
    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
    tmp_input_validator = input_validator.InputValidator(
        self._view.ui.le_protein_seq
    )
    the_entered_text = self._view.ui.le_protein_seq.toPlainText()
    tmp_validate_flag, tmp_message = (
        tmp_input_validator.validate_input_for_protein_sequence(
            the_entered_text,
        )
    )
    if tmp_validate_flag:
      self._view.ui.lbl_status_seq.setText("")
      self._view.ui.btn_add.setEnabled(True)
    else:
      self._view.ui.lbl_status_seq.setText(tmp_message)
      self._view.ui.btn_add.setEnabled(False)

  def _switch_ui_to_sequence_input(self) -> None:
    """Switches the user interface to the sequence input mask."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Next' button was clicked."
    )
    self._view.ui.le_seq_name.setEnabled(False)
    self._view.ui.btn_next.setEnabled(False)
    self._view.ui.lbl_protein_seq.show()
    self._view.ui.le_protein_seq.show()
    self._view.ui.lbl_status_seq.show()
    self._view.ui.btn_back.show()
    self._activate_add_button()

  def _switch_ui_to_sequence_name_input(self) -> None:
    """Switches the user interface to the sequence name input mask."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Back' button was clicked."
    )
    self._view.ui.le_seq_name.setEnabled(True)
    self._view.ui.btn_next.setEnabled(True)
    self._view.ui.lbl_protein_seq.hide()
    self._view.ui.le_protein_seq.hide()
    self._view.ui.lbl_status_seq.hide()
    self._view.ui.btn_back.hide()
    self._view.ui.btn_add.setEnabled(False)

  def _activate_add_button(self) -> None:
    """Activates the open button."""
    if self._view.ui.le_protein_seq.toPlainText() == "":
      self._view.ui.le_protein_seq.setStyleSheet(
          """QTextEdit {color: #000000; border-color: #DCDBE3;}""",
      )
      self._view.ui.lbl_status_seq.setText("")
      self._view.ui.btn_add.setEnabled(False)
    else:
      self._validate_protein_sequence()

  def _add_sequence(self) -> None:
    """Adds a protein to the global variable and closes the dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked."
    )
    self.return_value.emit(
        (
            self._view.ui.le_seq_name.text(),
            self._view.ui.le_protein_seq.toPlainText(),
        )
    )
    self._view.close()
