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
"""Module for the fasta file import preview view controller."""
import logging

from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.internal.data_structures.data_classes import basic_seq_info
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class FastaFileImportPreviewViewController(QtCore.QObject):
    """Class for the FastaFileImportPreviewViewController."""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self,
                 the_interface_manager: "interface_manager.InterfaceManager",
                 the_parsed_sequences: list[basic_seq_info.BasicSeqInfo]):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_fasta_file_import_preview_view()
        self._parsed_sequences: list[basic_seq_info.BasicSeqInfo] = the_parsed_sequences
        self._connect_all_ui_elements_to_slot_functions()

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_save.clicked.connect(self.__slot_save_sequence_changes)
        self._view.ui.table.itemSelectionChanged.connect(self.__slot_enable_table_edit_buttons)
        self._view.ui.table.cellChanged.connect(self._check_if_sequences_can_be_saved)
        self._view.ui.btn_delete.clicked.connect(self.__slot_delete_sequence_from_table)
        self._view.ui.btn_duplicate.clicked.connect(self.__slot_duplicate_selected_sequence)
        self._view.ui.btn_restore.clicked.connect(self._restore_default_table_values)

    def restore_ui(self):
        self._view.ui.lbl_description.setText("Sequences Overview")
        self._view.ui.lbl_description_2.setText("The chain letter will be lost, if you do the import! (The chain letter is important to store the sequences in a specific order.)")
        self._view.ui.btn_delete.setEnabled(False)
        self._view.ui.btn_duplicate.setEnabled(False)
        self._view.ui.btn_add.setEnabled(True)
        self._view.ui.btn_add.hide()
        self._view.ui.btn_insert.setEnabled(False)
        self._view.ui.btn_insert.hide()
        self._view.ui.btn_save.setEnabled(True)
        styles.color_bottom_frame_button(self._view.ui.btn_save)
        # restore table
        self._restore_table()

    def _restore_table(self):
        self._view.ui.table.clear()
        self._view.ui.table.setColumnCount(3)
        self._view.ui.table.setHorizontalHeaderLabels(["Name", "Chain", "Sequence"])

    def _restore_default_table_values(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Restore' button was clicked.")
        self._restore_table()
        self.fill_sequence_table()

    def fill_sequence_table(self):
        self._view.ui.table.setRowCount(len(self._parsed_sequences))
        i = 0
        for tmp_seq_info in self._parsed_sequences:
            tmp_name_item = QtWidgets.QTableWidgetItem(tmp_seq_info.name)
            tmp_chain_item = QtWidgets.QTableWidgetItem(tmp_seq_info.chain)
            tmp_sequence_item = QtWidgets.QTableWidgetItem(tmp_seq_info.seq)
            #tmp_sequence_item.setFlags(tmp_sequence_item.flags() & ~Qt.ItemIsEditable)
            self._view.ui.table.setItem(i, 0, tmp_name_item)
            self._view.ui.table.setItem(i, 1, tmp_chain_item)
            self._view.ui.table.setItem(i, 2, tmp_sequence_item)
            i += 1
        self._view.ui.table.resizeColumnsToContents()

    def _check_if_sequences_can_be_saved(self):
        if self._view.ui.table.rowCount() == 0:
            self._view.ui.btn_save.setEnabled(False)
            return
        else:
            self._view.ui.btn_save.setEnabled(True)
        
        try:
            for tmp_row_no in range(self._view.ui.table.rowCount()):
                if self._view.ui.table.item(tmp_row_no, 1).text() == "Enter a chain letter!" or self._view.ui.table.item(tmp_row_no, 1).text() == "Invalid input!":
                    self._view.ui.btn_save.setEnabled(False)
                    return
                elif self._view.ui.table.item(tmp_row_no, 0).text() == "Invalid input!":
                    self._view.ui.btn_save.setEnabled(False)
                    return
        except AttributeError:
            # Is needed because the signal cellChanged triggers this function even if the item is None
            return
        else:
            self._view.ui.btn_save.setEnabled(True)

    # @SLOT
    def __slot_save_sequence_changes(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Save' button was clicked.")
        self._parsed_sequences = []
        for tmp_row_no in range(self._view.ui.table.rowCount()):
            tmp_name = self._view.ui.table.item(tmp_row_no, 0).text()
            tmp_chain = self._view.ui.table.item(tmp_row_no, 1).text()
            tmp_sequence = self._view.ui.table.item(tmp_row_no, 2).text()
            tmp_seq_info = basic_seq_info.BasicSeqInfo(tmp_name, tmp_chain, tmp_sequence)
            self._parsed_sequences.append(tmp_seq_info)

        self._view.close()
        self.user_input.emit((0, self._parsed_sequences))

    def __slot_delete_sequence_from_table(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete' button was clicked.")
        selected_rows = sorted(set(index.row() for index in self._view.ui.table.selectedIndexes()))
        # Reverse the selected rows list to ensure correct removal when removing rows
        for row in reversed(selected_rows):
            self._view.ui.table.removeRow(row)
        self._view.ui.btn_delete.setEnabled(False)
        self._check_if_sequences_can_be_saved()

    def __slot_enable_table_edit_buttons(self) -> None:
        if len(self._view.ui.table.selectedItems()) == 0:
            self._view.ui.btn_delete.setEnabled(False)
            self._view.ui.btn_duplicate.setEnabled(False)
        else:
            self._view.ui.btn_delete.setEnabled(True)
            self._view.ui.btn_duplicate.setEnabled(True)
        self._check_if_sequences_can_be_saved()

    def __slot_duplicate_selected_sequence(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Duplicate' button was clicked.")
        tmp_name_item = QtWidgets.QTableWidgetItem(self._view.ui.table.item(self._view.ui.table.currentRow(), 0).text())
        tmp_chain_item = QtWidgets.QTableWidgetItem("Enter a valid chain letter!")
        tmp_seq_item = QtWidgets.QTableWidgetItem(self._view.ui.table.item(self._view.ui.table.currentRow(), 2).text())
        #tmp_seq_item.setFlags(tmp_seq_item.flags() & ~Qt.ItemIsEditable)
        self._view.ui.table.insertRow(self._view.ui.table.currentRow() + 1)
        self._view.ui.table.setItem(self._view.ui.table.currentRow() + 1, 0, tmp_name_item)
        self._view.ui.table.setItem(self._view.ui.table.currentRow() + 1, 1, tmp_chain_item)
        self._view.ui.table.setItem(self._view.ui.table.currentRow() + 1, 2, tmp_seq_item)

        self._view.ui.btn_duplicate.setEnabled(False)
        self._view.ui.btn_delete.setEnabled(False)
        self._view.ui.btn_save.setEnabled(False)
        self._view.ui.table.resizeColumnsToContents()