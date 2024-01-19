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
"""Module for the Create Dialog."""
import glob
import os
import pymol
from PyQt5 import QtCore, QtWidgets
from pymol import cmd

from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.util import input_validator, gui_utils, constants, tools


class CreateProjectViewController(QtCore.QObject):
    """Class for the Create Project View Controller."""
    string_model = QtCore.QStringListModel()

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_create_view()
        self._fill_projects_list_view()
        self._connect_all_ui_elements_to_slot_functions()
        self.hide_add_protein_options()

    def _fill_projects_list_view(self) -> None:
        """Lists all projects."""
        # self._view.ui.list_create_projects_view.setModel(self._interface_manager.get_workspace_model())
        xml_pattern = os.path.join(constants.DEFAULT_WORKSPACE_PATH, '*.xml')
        self.string_model.setStringList(
            # Filters the workspace for all project files based on the xml extension
            [os.path.basename(file).replace(".xml", "") for file in glob.glob(xml_pattern)]
        )
        self._view.ui.list_create_projects_view.setModel(self.string_model)

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.txt_new_project_name.textChanged.connect(self.validate_project_name)
        self._view.ui.cb_new_add_reference.clicked.connect(self.show_add_protein_options)
        # self._view.ui.txt_new_choose_reference.textChanged.connect(self.validate_reference_in_project)
        # self._view.ui.btn_cancel(self._close)

    def hide_add_protein_options(self) -> None:
        self._view.ui.lbl_new_status_choose_reference.hide()
        self._view.ui.lbl_new_choose_reference.hide()
        self._view.ui.txt_new_choose_reference.hide()
        self._view.ui.btn_new_choose_reference.hide()

    def validate_project_name(self) -> None:
        """Validates the input of the project name in real-time."""
        input_validator.InputValidator.validate_project_name(
            self._view.ui.list_create_projects_view.model(),
            self._view.ui.txt_new_project_name,
            self._view.ui.lbl_new_status_project_name,
            self._view.ui.txt_new_choose_reference
        )

    def validate_reference_in_project(self) -> None:
        """Checks if the entered reference protein is valid or not."""
        if len(self._view.ui.txt_new_choose_reference.text()) == 0:
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
            self._view.ui.lbl_new_status_choose_reference.setText("")
            self._view.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_new_create_project)
        elif len(self._view.ui.txt_new_choose_reference.text()) < 4:
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
            styles.color_button_not_ready(self._view.ui.btn_new_create_project)
            self._view.ui.btn_new_create_project.setEnabled(False)
            self._view.ui.lbl_new_status_choose_reference.setText("")
        # checks if a pdb id was entered
        elif len(self._view.ui.txt_new_choose_reference.text()) == 4:
            pdb_id = self._view.ui.txt_new_choose_reference.text().upper()
            try:
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                cmd.fetch(pdb_id, type="pdb", path=constants.SCRATCH_DIR)
                os.remove(f"{constants.SCRATCH_DIR}/{pdb_id}.pdb")
                cmd.reinitialize()
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #000000")
                self._view.ui.btn_new_create_project.setEnabled(True)
                styles.color_button_ready(self._view.ui.btn_new_create_project)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                styles.color_button_not_ready(self._view.ui.btn_new_create_project)
                return
            except FileNotFoundError:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.lbl_new_status_choose_reference.setText("Invalid PDB ID.")
                self._view.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_new_create_project)
                return
        else:
            if self._view.ui.txt_new_choose_reference.text().find("/") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_new_create_project)

            elif self._view.ui.txt_new_choose_reference.text().find("\\") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_new_create_project)

    def show_add_protein_options(self) -> None:
        if self._view.ui.cb_new_add_reference.isChecked():
            self._view.ui.lbl_new_choose_reference.show()
            self._view.ui.txt_new_choose_reference.show()
            self._view.ui.btn_new_choose_reference.show()

    def show_add_reference(self) -> None:
        """Shows the reference input section."""
        # checkbox is checked
        self._view.ui.cb_new_add_reference.checkState()
        if self._view.ui.cb_new_add_reference.checkState() == 2:
            self._view.ui.txt_new_choose_reference.clear()
            self._view.ui.txt_new_choose_reference.setStyleSheet("background-color: white")
            self._view.ui.lbl_new_choose_reference.show()
            self._view.ui.txt_new_choose_reference.show()
            self._view.ui.btn_new_choose_reference.show()
            self._view.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_new_create_project)
            # check internet connectivity
            if not tools.check_internet_connectivity():
                gui_utils.no_internet_dialog()
                self._view.ui.txt_new_choose_reference.setEnabled(False)
                self._view.ui.lbl_new_status_choose_reference.setText("You cannot enter a PDB ID (no internet).")
                return
            self._view.ui.txt_new_choose_reference.setEnabled(True)
            self._view.ui.lbl_new_status_choose_reference.setText("")
        else:
            self._view.ui.lbl_new_choose_reference.hide()
            self._view.ui.txt_new_choose_reference.hide()
            self._view.ui.btn_new_choose_reference.hide()
            self._view.ui.lbl_new_status_choose_reference.setText("")
            self._view.ui.btn_new_create_project.setEnabled(True)
            styles.color_button_ready(self._view.ui.btn_new_create_project)

    def load_reference_in_project(self) -> None:
        """Loads a reference in a new project."""
        try:
            # open file dialog
            file_name = QtWidgets.QFileDialog.getOpenFileName(
                self._view,
                "Open Reference",
                QtCore.QDir.homePath(),
                "PDB Files (*.pdb)",
            )
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self._view.ui.txt_new_choose_reference.setText(str(file_name[0]))
            self._view.ui.txt_new_choose_reference.setEnabled(False)
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #000000")
            self._view.ui.btn_new_create_project.setEnabled(True)
            styles.color_button_ready(self._view.ui.btn_new_create_project)
        except ValueError:
            print("No file has been selected.")
