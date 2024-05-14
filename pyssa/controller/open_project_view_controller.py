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
"""Module for the open project view controller."""

import glob
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.controller import interface_manager, pymol_session_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.gui.ui.styles import styles
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import input_validator, constants, ui_util
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class OpenProjectViewController(QtCore.QObject):
    """Class for the Open Project View Controller."""
    return_value = QtCore.pyqtSignal(str)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_open_view()
        self._fill_projects_list_view()
        self._project_names = self._convert_model_into_set()
        self._connect_all_ui_elements_to_slot_functions()
        self.restore_default_view()

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
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

    def __await_open_help(self, return_value):
        if return_value[0] == "":
            self._interface_manager.status_bar_manager.show_error_message("Opening help center failed!")
            return

        try:
            self._interface_manager.documentation_window = return_value[2]
            if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "The script for bringing the documentation window in front could not be found!", "Documentation",
                    custom_message_box.CustomMessageBoxIcons.ERROR.value
                )
                tmp_dialog.exec_()
            else:
                self._interface_manager.documentation_window.restore()
                subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
                self._interface_manager.status_bar_manager.show_temporary_message("Opening help center finished.")
        except Exception as e:
            logger.error(f"Error while opening help center: {e}")
            self._interface_manager.status_bar_manager.show_error_message("Opening help center failed!")

    def _open_help_for_dialog(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked.")
        self.open_help("help/project/open_project/")

    def restore_default_view(self) -> None:
        self._view.ui.label_28.hide()
        self._view.ui.txt_open_search.setPlaceholderText("Search")
        self._view.ui.txt_open_search.clear()
        self._view.ui.txt_open_selected_project.clear()
        self._view.ui.btn_open_project.setEnabled(False)

    def _fill_projects_list_view(self) -> None:
        """Lists all projects."""
        self._view.ui.projects_list_view.setModel(self._interface_manager.get_workspace_model())

    def _convert_model_into_set(self) -> set:
        tmp_project_names = []
        for tmp_row in range(self._view.ui.projects_list_view.model().rowCount()):
            tmp_project_names.append(
                self._view.ui.projects_list_view.model().index(tmp_row, 0).data(Qt.DisplayRole)
            )
        return set(tmp_project_names)

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.txt_open_search.textChanged.connect(self._validate_open_search)
        self._view.ui.projects_list_view.clicked.connect(self._select_project_from_open_list)
        self._view.ui.txt_open_selected_project.textChanged.connect(self._activate_open_button)
        self._view.ui.btn_open_project.clicked.connect(self._open_selected_project)
        self._view.ui.projects_list_view.doubleClicked.connect(self._open_selected_project)
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

    def _validate_open_search(self, the_entered_text: str) -> None:
        """Validates the input of the project name in real-time."""
        ui_util.select_matching_string_in_q_list_view(
            self._view.ui.txt_open_search.text(),
            self._view.ui.projects_list_view,
            self._view.ui.txt_open_selected_project
        )

        # tmp_validate_flag, tmp_stylesheet_string, tmp_message = input_validator.check_project_names_for_given_project_name(
        #     the_entered_text, self._project_names
        # )
        #
        # self._view.ui.txt_open_search.setStyleSheet(tmp_stylesheet_string)
        # if tmp_validate_flag:
        #     self._view.ui.lbl_open_status_search.setText("")
        # else:
        #     self._view.ui.lbl_open_status_search.setText(tmp_message)

        # projects_list_view = self._view.ui.projects_list_view
        #
        # # Deselect any current item in the list view
        # if projects_list_view.currentIndex().isValid():
        #     projects_list_view.model().itemFromIndex(projects_list_view.currentIndex()).setSelected(False)
        #
        # # Assuming validate_search_input is a static method
        # input_validator.InputValidator.validate_search_input(
        #     projects_list_view.model(),
        #     self._view.ui.txt_open_search,
        #     self._view.ui.lbl_open_status_search,
        #     self._view.ui.txt_open_selected_project,
        # )

    def _select_project_from_open_list(self) -> None:
        """Sets the selected project name in the text box."""
        tmp_project_name = self._view.ui.projects_list_view.model().data(
            self._view.ui.projects_list_view.currentIndex(), Qt.DisplayRole
        )
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, f"The project '{tmp_project_name}' from the list was clicked.")
        try:
            self._view.ui.txt_open_selected_project.setText(tmp_project_name)
        except AttributeError:
            self._view.ui.txt_open_selected_project.setText("")

    def _activate_open_button(self) -> None:
        """Activates the open button."""
        if self._view.ui.txt_open_selected_project.text() == "":
            self._view.ui.btn_open_project.setEnabled(False)
            # styles.color_button_not_ready(self._view.ui.btn_open_project)
        else:
            self._view.ui.btn_open_project.setEnabled(True)
            # styles.color_button_ready(self._view.ui.btn_open_project)

    def _open_selected_project(self):
        tmp_project_name = self._view.ui.projects_list_view.model().data(
            self._view.ui.projects_list_view.currentIndex(), Qt.DisplayRole
        )
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, f"The project '{tmp_project_name}' from the list was double-clicked or the 'Open' button was clicked.")
        self._view.close()
        self.return_value.emit(self._view.ui.txt_open_selected_project.text())
