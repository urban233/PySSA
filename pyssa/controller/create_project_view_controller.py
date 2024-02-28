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
import subprocess

import pymol
from PyQt5 import QtCore, QtWidgets
from pymol import cmd
from PyQt5.QtCore import Qt
from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import input_validator, gui_utils, constants, tools


class CreateProjectViewController(QtCore.QObject):
    """Class for the Create Project View Controller."""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_create_view()
        self._restore_ui()
        self._project_names: set = self._convert_model_into_set()
        self._connect_all_ui_elements_to_slot_functions()
        self._hide_add_protein_options()

    # <editor-fold desc="Util methods">
    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.update_status_bar("Opening help center ...")
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, 0),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self):
        subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
        self._interface_manager.update_status_bar("Opening help center finished.")

    def _open_help_for_dialog(self):
        self.open_help("help/project/new_project/")

    def _restore_ui(self):
        """Restores the ui."""
        self._view.ui.list_create_projects_view.clearSelection()
        self._fill_projects_list_view()
        self._view.ui.txt_new_project_name.clear()
        self._view.ui.txt_new_project_name.setStyleSheet(
            """QLineEdit {color: #000000; border-color: #DCDBE3;}"""
        )
        self._view.ui.lbl_new_status_choose_reference.setText("")
        self._view.ui.lbl_new_status_choose_reference.setStyleSheet("color: #ba1a1a; font-size: 10px;")
        self._view.ui.lbl_new_status_project_name.setText("")
        self._view.ui.lbl_new_status_project_name.setStyleSheet("color: #ba1a1a; font-size: 10px;")
        self._view.ui.btn_new_create_project.setEnabled(False)
        self._view.ui.cb_new_add_reference.setCheckState(False)
        self._view.ui.txt_new_choose_reference.clear()
        self._hide_add_protein_options()

    def _fill_projects_list_view(self) -> None:
        """Lists all projects."""
        self._view.ui.list_create_projects_view.setModel(self._interface_manager.get_workspace_projects())

    def _convert_model_into_set(self) -> set:
        tmp_project_names = []
        for tmp_row in range(self._view.ui.list_create_projects_view.model().rowCount()):
            tmp_project_names.append(
                self._view.ui.list_create_projects_view.model().index(tmp_row, 0).data(Qt.DisplayRole)
            )
        return set(tmp_project_names)

    # </editor-fold>

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.txt_new_project_name.textChanged.connect(self._validate_project_name)
        self._view.ui.cb_new_add_reference.clicked.connect(self._show_add_protein_options)
        self._view.ui.btn_new_choose_reference.clicked.connect(self._load_reference_in_project)
        self._view.ui.btn_new_create_project.clicked.connect(self._create_new_project)
        self._view.ui.txt_new_choose_reference.textChanged.connect(self._validate_reference_in_project)

    def _create_new_project(self) -> None:
        self._view.close()
        self.user_input.emit((self._view.ui.txt_new_project_name.text(), self._view.ui.txt_new_choose_reference.text()))

    def _hide_add_protein_options(self) -> None:
        self._view.ui.lbl_new_status_choose_reference.hide()
        self._view.ui.lbl_new_choose_reference.hide()
        self._view.ui.txt_new_choose_reference.hide()
        self._view.ui.btn_new_choose_reference.hide()

    def _validate_project_name(self, the_entered_text: str) -> None:
        """Validates the input of the project name in real-time."""
        tmp_input_validator = input_validator.InputValidator(self._view.ui.txt_new_project_name)
        tmp_validate_flag, tmp_message = tmp_input_validator.validate_input_for_project_name(the_entered_text,
                                                                                             self._project_names)
        if tmp_validate_flag:
            self._view.ui.lbl_new_status_project_name.setText("")
            self._view.ui.btn_new_create_project.setEnabled(True)
        else:
            self._view.ui.lbl_new_status_project_name.setText(tmp_message)
            self._view.ui.btn_new_create_project.setEnabled(False)

    def _validate_reference_in_project(self) -> None:
        """Checks if the entered reference protein is valid or not."""
        if len(self._view.ui.txt_new_choose_reference.text()) == 0:
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
            self._view.ui.lbl_new_status_choose_reference.setText("")
            self._view.ui.btn_new_create_project.setEnabled(False)
        elif len(self._view.ui.txt_new_choose_reference.text()) < 4:
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
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

            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                return
            except FileNotFoundError:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.lbl_new_status_choose_reference.setText("Invalid PDB ID.")
                self._view.ui.btn_new_create_project.setEnabled(False)
                return
        else:
            if self._view.ui.txt_new_choose_reference.text().find("/") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)

            elif self._view.ui.txt_new_choose_reference.text().find("\\") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)

    def _show_add_protein_options(self) -> None:
        if self._view.ui.cb_new_add_reference.isChecked():
            self._view.ui.lbl_new_choose_reference.show()
            self._view.ui.lbl_new_status_choose_reference.show()
            self._view.ui.txt_new_choose_reference.show()
            self._view.ui.btn_new_choose_reference.show()
        else:
            self._view.ui.lbl_new_choose_reference.hide()
            self._view.ui.lbl_new_status_choose_reference.hide()
            self._view.ui.txt_new_choose_reference.hide()
            self._view.ui.btn_new_choose_reference.hide()

    def _load_reference_in_project(self) -> None:
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
        except ValueError:
            print("No file has been selected.")
