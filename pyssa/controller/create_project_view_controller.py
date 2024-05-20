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
"""Module for the create project view controller."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt
from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import input_validator, constants, exception
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class CreateProjectViewController(QtCore.QObject):
    """Class for the CreateProjectViewController."""
    
    user_input = QtCore.pyqtSignal(tuple)
    """Singal used to transfer data back to the previous window."""
    
    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
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
        self._view = the_interface_manager.get_create_view()
        self._restore_ui()
        self._project_names: set = self._convert_model_into_set()
        self._connect_all_ui_elements_to_slot_functions()
        self._hide_add_protein_options()

    # <editor-fold desc="Util methods">
    def open_help(self, a_page_name: str) -> None:
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        
        Raises:
            exception.IllegalArgumentError: If `a_page_name` is None.
        """
        # <editor-fold desc="Checks">
        if a_page_name is None:
            logger.error("a_page_name is None.")
            raise exception.IllegalArgumentError("a_page_name is None.")
        
        # </editor-fold>
        
        try:
            self._interface_manager.status_bar_manager.show_temporary_message("Opening help center ...")
            if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) != 1:
                self._interface_manager.documentation_window = None
            self._active_task = tasks.LegacyTask(
                target=util_async.open_documentation_on_certain_page,
                args=(a_page_name, self._interface_manager.documentation_window),
                post_func=self.__await_open_help,
            )
        except Exception as e:
            logger.error(f"Error while opening help center: {e}")
        else:
            self._active_task.start()

    def __await_open_help(self, return_value: tuple) -> None:
        """Opens the help center and performs necessary actions based on the return value.

        Args:
            return_value (tuple): The return value from opening the help center.
        """
        # <editor-fold desc="Checks">
        if return_value[0] == "":
            self._interface_manager.status_bar_manager.show_error_message("Opening help center failed!")
            return
        
        # </editor-fold>

        try:
            self._interface_manager.documentation_window = return_value[2]
            if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "The script for bringing the documentation window in front could not be found!", "Documentation",
                    custom_message_box.CustomMessageBoxIcons.ERROR.value,
                )
                tmp_dialog.exec_()
            else:
                self._interface_manager.documentation_window.restore()
                subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
                self._interface_manager.status_bar_manager.show_temporary_message("Opening help center finished.")
        except Exception as e:
            logger.error(f"Error while opening help center: {e}")
            self._interface_manager.status_bar_manager.show_error_message("Opening help center failed!")

    def _open_help_for_dialog(self) -> None:
        """Opens the help dialog."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked.")
        self.open_help("help/project/new_project/")

    def _restore_ui(self) -> None:
        """Restores the UI."""
        self._view.ui.list_create_projects_view.clearSelection()
        self._fill_projects_list_view()
        self._view.ui.txt_new_project_name.clear()
        self._view.ui.txt_new_project_name.setStyleSheet(
            """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
        )
        self._view.ui.lbl_new_status_choose_reference.setText("")
        self._view.ui.lbl_new_status_choose_reference.setStyleSheet("color: #ba1a1a; font-size: 11px;")
        self._view.ui.lbl_new_status_project_name.setText("")
        self._view.ui.lbl_new_status_project_name.setStyleSheet("color: #ba1a1a; font-size: 11px;")
        self._view.ui.btn_new_create_project.setEnabled(False)
        self._view.ui.cb_new_add_reference.setCheckState(False)
        self._view.ui.txt_new_choose_reference.clear()
        self._hide_add_protein_options()

    def _fill_projects_list_view(self) -> None:
        """Lists all projects."""
        self._view.ui.list_create_projects_view.setModel(self._interface_manager.get_workspace_projects())

    def _convert_model_into_set(self) -> set:
        """Converts the model into a set of project names.

        Returns:
            A set containing the project names.
        """
        tmp_project_names = []
        for tmp_row in range(self._view.ui.list_create_projects_view.model().rowCount()):
            tmp_project_names.append(
                self._view.ui.list_create_projects_view.model().index(tmp_row, 0).data(Qt.DisplayRole),
            )
        return set(tmp_project_names)

    # </editor-fold>

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        """Connects all UI elements to their corresponding slot functions in the class."""
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.txt_new_project_name.textChanged.connect(self._validate_project_name)
        self._view.ui.cb_new_add_reference.clicked.connect(self._show_add_protein_options)
        self._view.ui.btn_new_choose_reference.clicked.connect(self._load_reference_in_project)
        self._view.ui.btn_new_create_project.clicked.connect(self._create_new_project)
        self._view.ui.txt_new_choose_reference.textChanged.connect(self._validate_reference_in_project)

    def _create_new_project(self) -> None:
        """Creates a new project based on the user's input."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Create' button was clicked.")
        self._view.close()
        self.user_input.emit((self._view.ui.txt_new_project_name.text(), self._view.ui.txt_new_choose_reference.text()))

    def _hide_add_protein_options(self) -> None:
        """Hides the add protein options in the UI."""
        self._view.ui.lbl_new_status_choose_reference.hide()
        self._view.ui.lbl_new_choose_reference.hide()
        self._view.ui.txt_new_choose_reference.hide()
        self._view.ui.btn_new_choose_reference.hide()

    def _validate_project_name(self, the_entered_text: str) -> None:
        """Validates the input of the project name in real-time.

        Args:
            the_entered_text (str): The text entered by the user.
        """
        # <editor-fold desc="Checks">
        if the_entered_text is None:
            logger.error("the_entered_text is None.")
            raise exception.IllegalArgumentError("the_entered_text is None.")
        
        # </editor-fold>
        
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
        tmp_validate_flag, tmp_stylesheet_string, tmp_message = input_validator.validate_input_for_project_name(
            the_entered_text, self._project_names,
        )
        self._view.ui.txt_new_project_name.setStyleSheet(tmp_stylesheet_string)

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
        else:
            if self._view.ui.txt_new_choose_reference.text().find("/") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)

            elif self._view.ui.txt_new_choose_reference.text().find("\\") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)

    def _show_add_protein_options(self) -> None:
        """Show or hide add protein options based on the value of the 'cb_new_add_reference' checkbox."""
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
            logger.error("No file has been selected.")
