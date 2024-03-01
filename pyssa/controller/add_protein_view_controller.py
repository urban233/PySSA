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
"""Module for the Open Dialog."""
import os
import pymol
from pymol import cmd
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.controller import interface_manager
from pyssa.internal.portal import pymol_io
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import validate_async
from pyssa.io_pyssa import bio_data
from pyssa.util import constants


class AddProteinViewController(QtCore.QObject):
    """Class for the RenameProteinViewController class"""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_add_protein_view()
        self._connect_all_ui_elements_to_slot_functions()

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_choose_protein.clicked.connect(self.__slot_load_protein_from_filesystem)
        self._view.ui.btn_add_protein.clicked.connect(self.__slot_add_protein)
        self._view.ui.txt_add_protein.textChanged.connect(self.__slot_validate_input)

    def restore_ui(self):
        self._view.ui.txt_add_protein.clear()
        self._view.ui.btn_add_protein.setEnabled(False)
        self._view.setMinimumWidth(500)

    # @SLOT
    def __slot_validate_input(self, the_entered_text) -> None:
        """Checks if the entered reference protein is valid or not."""
        self._view.ui.lbl_status.setStyleSheet(
            """QLabel {color: #ba1a1a;}"""
        )
        if len(the_entered_text) == 0:
            # empty line edit
            self._view.ui.txt_add_protein.setStyleSheet(
                """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}"""
            )
            self._view.ui.lbl_status.setText("Please enter a PDB id or choose an existing .pdb file from your filesystem!")
            self._view.ui.btn_add_protein.setEnabled(False)
        elif len(the_entered_text) < 4:
            # length of text is too small
            self._view.ui.txt_add_protein.setStyleSheet(
                """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}"""
            )
            self._view.ui.btn_add_protein.setEnabled(False)
            self._view.ui.lbl_status.setText("Please enter more characters.")
        # checks if a pdb id was entered
        else:
            self._active_task = tasks.Task(
                target=validate_async.validate_add_protein_view_input,
                args=(
                    the_entered_text, 0
                ),
                post_func=self.__await__slot_validate_input,
            )
            self._active_task.start()
            self._view.wait_spinner.start()
            self._view.ui.txt_add_protein.setStyleSheet(
                """QLineEdit {color: #000000; border-color: #DCDBE3;}"""
            )
            self._view.ui.lbl_status.setStyleSheet(
                """QLabel {color: #367AF6;}"""
            )
            self._view.ui.lbl_status.setText("Checking input ...")

        #     # pdb id is used
        #     pdb_id = self._view.ui.txt_add_protein.text().upper()
        #     tmp_filepath: str = f"{constants.SCRATCH_DIR}/{pdb_id}.pdb"
        #     bio_data.download_pdb_file(pdb_id, tmp_filepath)
        #     if os.path.exists(tmp_filepath):
        #         os.remove(tmp_filepath)
        #         self._view.ui.txt_add_protein.setStyleSheet(
        #             """QLineEdit {color: #000000; border-color: #DCDBE3;}"""
        #         )
        #         self._view.ui.btn_add_protein.setEnabled(True)
        #     else:
        #         self._view.ui.txt_add_protein.setStyleSheet(
        #             """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}"""
        #         )
        #         self._view.ui.lbl_status.setText("Invalid PDB id!")
        # else:
        #     if os.path.exists(the_entered_text):
        #         self._view.ui.txt_add_protein.setStyleSheet(
        #             """QLineEdit {color: #000000; border-color: #DCDBE3;}"""
        #         )
        #         self._view.ui.lbl_status.setText("")
        #         self._view.ui.btn_add_protein.setEnabled(True)
        #     else:
        #         self._view.ui.txt_add_protein.setStyleSheet(
        #             """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}"""
        #         )
        #         self._view.ui.lbl_status.setText("Invalid filepath!")
        #         self._view.ui.btn_add_protein.setEnabled(False)

    def __await__slot_validate_input(self, return_value: tuple):
        tmp_type: int = return_value[0]
        tmp_is_valid: bool = return_value[1]
        self._view.ui.lbl_status.setStyleSheet(
            """QLabel {color: #ba1a1a;}"""
        )
        if tmp_is_valid:
            self._view.ui.txt_add_protein.setStyleSheet(
                """QLineEdit {color: #000000; border-color: #DCDBE3;}"""
            )
            self._view.ui.lbl_status.setText("")
            self._view.ui.btn_add_protein.setEnabled(True)
        elif not tmp_is_valid and tmp_type == 1:
            self._view.ui.txt_add_protein.setStyleSheet(
                """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}"""
            )
            self._view.ui.lbl_status.setText("Invalid PDB id!")
            self._view.ui.btn_add_protein.setEnabled(False)
        elif not tmp_is_valid and tmp_type == 2:
            self._view.ui.txt_add_protein.setStyleSheet(
                """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}"""
            )
            self._view.ui.lbl_status.setText("Invalid filepath!")
            self._view.ui.btn_add_protein.setEnabled(False)
        else:
            constants.PYSSA_LOGGER.error("There is an unknown case, while validating the add protein view user input!")
        self._view.wait_spinner.stop()
        self._view.ui.txt_add_protein.setFocus()

    def __slot_load_protein_from_filesystem(self) -> None:
        """Loads a protein from the filesystem into the textbox."""
        try:
            # open file dialog
            file_name = QtWidgets.QFileDialog.getOpenFileName(
                self._view,
                "Open existing protein",
                QtCore.QDir.homePath(),
                "PDB Files (*.pdb)",
            )
            if file_name == ("", ""):
                self._view.ui.lbl_status.setText("No file has been selected.")
            else:
                # display path in text box
                self._view.ui.txt_add_protein.setText(str(file_name[0]))
                self._view.ui.btn_add_protein.setEnabled(True)
        except FileNotFoundError:
            self._view.ui.lbl_status.setText("Loading the protein structure failed!")

    def __slot_add_protein(self) -> None:
        """Adds a protein to the global variable and closes the dialog."""
        self._view.close()
        self.user_input.emit((self._view.ui.txt_add_protein.text(), len(self._view.ui.txt_add_protein.text())))
    
    def _validate_scene_name(self, text):
        new_text = ''.join(char for char in text)
        self._view.line_edit_scene_name.setText(new_text)
        if new_text in self._all_current_scenes:
            self._view.btn_add_scene.setEnabled(False)
            self._view.line_edit_scene_name.setToolTip("This scene name already exists. Please enter another name.")
            self._view.line_edit_scene_name.setStyleSheet(
                """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}"""
            )
        else:
            self._view.btn_add_scene.setEnabled(True)
            self._view.line_edit_scene_name.setToolTip("")
            self._view.line_edit_scene_name.setStyleSheet(
                """QLineEdit {color: #000000; border-color: #DCDBE3;}"""
            )
