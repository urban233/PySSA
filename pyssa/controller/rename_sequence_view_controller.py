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
"""Module for the rename sequence view controller."""
import glob
import logging
import os
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.util import input_validator, constants
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class RenameSequenceViewController(QtCore.QObject):
    """Class for the RenameSequenceViewController class"""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_rename_sequence_view()
        self._sequence_names = self._convert_sequence_model_into_set()
        self._connect_all_ui_elements_to_slot_functions()

    def restore_ui(self):
        self._view.ui.le_name.clear()
        self._view.ui.lbl_status.setText("")
        self._view.ui.le_name.setStyleSheet(
            """QTextEdit {color: #000000; border-color: #DCDBE3;}"""
        )

    def _convert_sequence_model_into_set(self) -> set:
        tmp_sequence_names = []
        for tmp_row in range(self._interface_manager.get_main_view().ui.seqs_list_view.model().rowCount()):
            tmp_sequence_names.append(
                self._interface_manager.get_main_view().ui.seqs_list_view.model().index(tmp_row, 0).data(Qt.DisplayRole)
            )
        return set(tmp_sequence_names)

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_rename.clicked.connect(self._rename_sequence)
        self._view.ui.le_name.textChanged.connect(self._validate_protein_name)

    def _validate_protein_name(self, the_entered_text: str) -> None:
        """Validates the input of the protein name in real-time."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
        tmp_input_validator = input_validator.InputValidator(self._view.ui.le_name)
        tmp_validate_flag, tmp_message = tmp_input_validator.validate_input_for_sequence_name(
            the_entered_text, self._sequence_names
        )
        if tmp_validate_flag:
            self._view.ui.lbl_status.setText("")
            self._view.ui.btn_rename.setEnabled(True)
        else:
            self._view.ui.lbl_status.setText(tmp_message)
            self._view.ui.btn_rename.setEnabled(False)

    def _rename_sequence(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Rename' button was clicked.")
        self._view.close()
        self.user_input.emit((self._view.ui.le_name.text(), True))
