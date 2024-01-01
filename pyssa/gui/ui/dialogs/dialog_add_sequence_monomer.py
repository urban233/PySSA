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
from pymol import Qt
from PyQt5 import QtCore
from pyssa.util import input_validator, gui_utils
from pyssa.gui.ui.forms.auto_generated.auto_dialog_add_sequence_monomer import Ui_Dialog


class DialogAddSequenceMonomer(Qt.QtWidgets.QDialog):
    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.show_stage_protein_name()

        # connections
        self.ui.txt_prot_name.textChanged.connect(self.validate_protein_name)
        self.ui.txt_seq_name.textChanged.connect(self.validate_protein_sequence)
        self.ui.btn_next.clicked.connect(self.show_stage_protein_sequence)
        self.ui.btn_back.clicked.connect(self.show_stage_protein_name)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.show()

    def validate_protein_name(self):
        input_validator.InputValidator.validate_protein_name(self.ui.txt_prot_name, self.ui.lbl_prot_name_status)

    def validate_protein_sequence(self):
        input_validator.InputValidator.validate_protein_sequence(
            self.ui.txt_seq_name, self.ui.lbl_seq_name_status, self.ui.btn_add_protein
        )

    def show_stage_protein_name(self):
        gui_elements_to_hide = [
            self.ui.btn_add_protein,
            self.ui.btn_back,
            self.ui.txt_seq_name,
            self.ui.lbl_seq_name,
            self.ui.lbl_seq_name_status,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_elements_to_show = [
            self.ui.btn_next,
            self.ui.btn_cancel,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)

    def show_stage_protein_sequence(self):
        gui_elements_to_hide = [
            self.ui.btn_next,
            self.ui.btn_cancel,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.txt_prot_name.setEnabled(False)
        gui_elements_to_show = [
            self.ui.btn_add_protein,
            self.ui.btn_back,
            self.ui.txt_seq_name,
            self.ui.lbl_seq_name,
            self.ui.lbl_seq_name_status,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)


# class DialogAddSequenceMonomer(Qt.QtWidgets.QDialog):
#
#     def __init__(self, parent=None):
#         """Constructor
#
#         Args:
#         """
#         Qt.QtWidgets.QDialog.__init__(self, parent)
#         # widgets for protein name
#         self.lbl_prot_name = QtWidgets.QLabel("Protein name")
#         self.txt_prot_name = QtWidgets.QLineEdit()
#         self.lbl_prot_name_status = QtWidgets.QLabel()
#         # widgets for next
#         self.btn_next = QtWidgets.QPushButton("Next")
#         self.btn_cancel = QtWidgets.QPushButton("Cancel")
#         # widgets for protein sequence
#         self.lbl_seq_name = QtWidgets.QLabel("Protein sequence")
#         self.txt_seq_name = QtWidgets.QTextEdit()
#         self.lbl_seq_name_status = QtWidgets.QLabel()
#         # widgets for add protein or cancel
#         self.btn_add_protein = QtWidgets.QPushButton("Add")
#         self.btn_back = QtWidgets.QPushButton("Back")
#         # layouts
#         self.outer_layout = QtWidgets.QVBoxLayout()
#         self.prot_name_layout = QtWidgets.QVBoxLayout()
#         self.next_stage_layout = QtWidgets.QHBoxLayout()
#         self.seq_name_layout = QtWidgets.QVBoxLayout()
#         self.add_protein_layout = QtWidgets.QHBoxLayout()
#         self.prot_name_layout.addWidget(self.lbl_prot_name)
#         self.prot_name_layout.addWidget(self.txt_prot_name)
#         self.prot_name_layout.addWidget(self.lbl_prot_name_status)
#         self.next_stage_layout.addStretch()
#         self.next_stage_layout.addWidget(self.btn_cancel)
#         self.next_stage_layout.addWidget(self.btn_next)
#         self.seq_name_layout.addWidget(self.lbl_seq_name)
#         self.seq_name_layout.addWidget(self.txt_seq_name)
#         self.seq_name_layout.addWidget(self.lbl_seq_name_status)
#         self.add_protein_layout.addStretch()
#         self.add_protein_layout.addWidget(self.btn_back)
#         self.add_protein_layout.addWidget(self.btn_add_protein)
#         self.outer_layout.addLayout(self.prot_name_layout)
#         self.outer_layout.addLayout(self.next_stage_layout)
#         self.outer_layout.addLayout(self.seq_name_layout)
#         self.outer_layout.addLayout(self.add_protein_layout)
#         self.setLayout(self.outer_layout)
#         self.show_stage_protein_name()
#
