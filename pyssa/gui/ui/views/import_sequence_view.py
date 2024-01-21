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
"""Module for the add proteins dialog."""
import os
import pymol
from PyQt5.QtCore import pyqtSignal
from pymol import Qt
from pymol import cmd
from PyQt5 import QtCore
from PyQt5 import QtGui
from pyssa.gui.ui.forms.auto_generated import auto_import_sequence_view
from pyssa.gui.ui.styles import styles
from pyssa.util import constants, tools, gui_utils

global_var_add_model = ("", False)


class ImportSequenceView(Qt.QtWidgets.QDialog):
    """Class for a dialog to add proteins to a project."""

    """
    A pyqtsignal that is used to hand-over the protein structure information.
    """
    return_value = pyqtSignal(tuple)

    def __init__(self, parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            parent: The parent.
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = auto_import_sequence_view.Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.btn_import_sequence.setEnabled(False)
        self.ui.txt_import_sequence.setEnabled(False)
        self.ui.btn_choose_fasta_file.clicked.connect(self.load_model)
        self.ui.btn_cancel.clicked.connect(self.close_dialog)
        self.ui.btn_import_sequence.clicked.connect(self.import_sequence)
        #self.ui.txt_import_sequence.textChanged.connect(self.validate_reference_in_project)
        self.ui.lbl_status.setText("")
        self.ui.btn_choose_fasta_file.setToolTip("Click to add a .pdb file")
        self.setWindowTitle("Add an existing protein to the current project")
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self)
        styles.color_button_not_ready(self.ui.btn_import_sequence)
        # fixme: this flag needs to be set if the WhatsThat icon in the window bar should be hidden
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        # # check internet connectivity
        # if not tools.check_internet_connectivity():
        #     gui_utils.no_internet_dialog_with_custom_msg(
        #         "You do not have a working internet connection which is "
        #         "necessary for connecting to the PDB!\n"
        #         "However you can add a protein structure from "
        #         "your local filesystem.",
        #     )
        #     self.ui.txt_import_sequence.setEnabled(False)
        #     self.ui.lbl_status.setText("You cannot enter a PDB ID (no internet).")

    # @SLOT
    def validate_reference_in_project(self) -> None:
        """Checks if the entered reference protein is valid or not."""
        if len(self.ui.txt_import_sequence.text()) == 0:
            self.ui.txt_import_sequence.setStyleSheet("color: #FC5457")
            self.ui.lbl_status.setText("")
            self.ui.btn_import_sequence.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_import_sequence)
        elif len(self.ui.txt_import_sequence.text()) < 4:
            self.ui.txt_import_sequence.setStyleSheet("color: #FC5457")
            styles.color_button_not_ready(self.ui.btn_import_sequence)
            self.ui.btn_import_sequence.setEnabled(False)
            self.ui.lbl_status.setText("")
        # checks if a pdb id was entered
        elif len(self.ui.txt_import_sequence.text()) == 4:
            pdb_id = self.ui.txt_import_sequence.text().upper()
            try:
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                cmd.fetch(pdb_id, type="pdb", path=constants.SCRATCH_DIR)
                os.remove(f"{constants.SCRATCH_DIR}/{pdb_id}.pdb")
                cmd.reinitialize()
                self.ui.txt_import_sequence.setStyleSheet("color: #000000")
                self.ui.btn_import_sequence.setEnabled(True)
                styles.color_button_ready(self.ui.btn_import_sequence)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self.ui.txt_import_sequence.setStyleSheet("color: #FC5457")
                styles.color_button_not_ready(self.ui.btn_import_sequence)
                return
            except FileNotFoundError:
                self.ui.txt_import_sequence.setStyleSheet("color: #FC5457")
                self.ui.lbl_status.setText("Invalid PDB ID.")
                self.ui.btn_import_sequence.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_import_sequence)
                return
        else:
            if self.ui.txt_import_sequence.text().find("/") == -1:
                self.ui.txt_import_sequence.setStyleSheet("color: #FC5457")
                self.ui.btn_import_sequence.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_import_sequence)

            elif self.ui.txt_import_sequence.text().find("\\") == -1:
                self.ui.txt_import_sequence.setStyleSheet("color: #FC5457")
                self.ui.btn_import_sequence.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_import_sequence)
            else:
                # TODO: displays a correct filepath in red, needs to be fixed!
                self.ui.txt_import_sequence.setStyleSheet("color: #000000")
                self.ui.btn_import_sequence.setEnabled(True)
                styles.color_button_ready(self.ui.btn_import_sequence)

    def close_dialog(self) -> None:
        """Closes the dialog."""
        self.close()

    def load_model(self) -> None:
        """Loads a protein from the filesystem into the textbox."""
        try:
            # open file dialog
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(
                self,
                "Import existing sequence",
                Qt.QtCore.QDir.homePath(),
                "FASTA Files (*.fasta)",
            )
            if file_name == ("", ""):
                self.ui.lbl_status.setText("No file has been selected.")
            else:
                # display path in text box
                self.ui.txt_import_sequence.setText(str(file_name[0]))
                self.ui.btn_import_sequence.setEnabled(True)
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the protein sequence failed!")
            self.ui.lbl_status.setText("Loading the protein sequence failed!")

    def import_sequence(self) -> None:
        """Adds a protein to the global variable and closes the dialog."""
        self.return_value.emit((self.ui.txt_import_sequence.text(), len(self.ui.txt_import_sequence.text())))
        self.close()