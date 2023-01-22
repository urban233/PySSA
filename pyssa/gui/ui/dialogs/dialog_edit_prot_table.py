#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
from PyQt5.QtWidgets import QTableWidgetItem
from pymol import Qt
from util import input_validator
from pyssa.gui.ui.forms.auto_generated.auto_dialog_edit_prot_table import Ui_Dialog


class DialogEditProtTable(Qt.QtWidgets.QDialog):

    def __init__(self, parent=None):
        """Constructor

        Args:
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.validator = input_validator.InputValidator()

        self.protein_name = None
        self.chain = None
        self.sequence = None
        self.table = None
        self.current_row = None

        self.ui.lbl_edit_prot_table_status_protein_name.setText("")
        self.ui.lbl_edit_prot_table_status_chain.setText("")
        self.ui.lbl_edit_prot_table_status_seq.setText("")

        self.ui.lbl_edit_prot_table_protein_name.hide()
        self.ui.lbl_edit_prot_table_status_protein_name.hide()
        self.ui.txt_edit_prot_table_protein_name.hide()

        self.ui.lbl_edit_prot_table_chain.hide()
        self.ui.lbl_edit_prot_table_status_chain.hide()
        self.ui.txt_edit_prot_table_chain.hide()

        self.ui.btn_edit_prot_table_cancel.clicked.connect(self.close_dialog)
        self.ui.btn_edit_prot_table_save.clicked.connect(self.set_parameters_in_table)
        self.ui.txt_edit_prot_table_seq.textChanged.connect(self.check_sequence)

        self.setMaximumSize(500, 125)
        self.setWindowTitle("Edit protein")

    # @SLOT
    def close_dialog(self):
        self.close()

    def set_all_parameters(self, table, current_row):
        self.table = table
        self.current_row = current_row
        self.ui.txt_edit_prot_table_protein_name.setText(self.protein_name.text())
        self.ui.txt_edit_prot_table_chain.setText(self.chain.text())
        self.ui.txt_edit_prot_table_seq.setText(self.sequence.text())

    def set_parameters_in_table(self):
        self.table.setItem(self.current_row, 0, QTableWidgetItem(self.ui.txt_edit_prot_table_protein_name.text()))
        self.table.setItem(self.current_row, 1, QTableWidgetItem(self.ui.txt_edit_prot_table_chain.text()))
        self.table.setItem(self.current_row, 2, QTableWidgetItem(self.ui.txt_edit_prot_table_seq.text()))
        self.close()

    def check_sequence(self):
        self.validator.validate_protein_sequence(self.ui.txt_edit_prot_table_seq,
                                                 self.ui.lbl_edit_prot_table_status_seq,
                                                 self.ui.btn_edit_prot_table_save)
