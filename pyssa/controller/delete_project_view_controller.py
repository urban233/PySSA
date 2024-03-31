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

"""Module for the Delete Dialog."""

import os
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.gui.ui.styles import styles
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import input_validator, gui_utils, constants, enums


class DeleteProjectViewController(QtCore.QObject):
    """Class for the Delete Project View Controller."""

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_delete_view()
        self._fill_projects_list_view()
        self._connect_all_ui_elements_to_slot_functions()
        self.restore_default_view()

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.status_bar_manager.show_temporary_message("Opening help center ...")
        if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) != 1:
            self._interface_manager.documentation_window = None
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, self._interface_manager.documentation_window),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self, return_value):
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

    def _open_help_for_dialog(self):
        self.open_help("help/project/delete_project/")

    def restore_default_view(self):
        self._view.ui.txt_delete_selected_projects.clear()
        self._view.ui.btn_delete_delete_project.setEnabled(False)

    def _fill_projects_list_view(self) -> None:
        """Lists all projects."""
        self._view.ui.list_delete_projects_view.setModel(self._interface_manager.get_workspace_model())

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.txt_delete_search.textChanged.connect(self.validate_delete_search)
        self._view.ui.list_delete_projects_view.clicked.connect(self.select_project_from_delete_list)
        self._view.ui.txt_delete_selected_projects.textChanged.connect(self.activate_delete_button)
        self._view.ui.btn_delete_delete_project.clicked.connect(self.delete_project)

    def validate_delete_search(self) -> None:
        """Validates the input of the project name in real-time."""
        input_validator.InputValidator.validate_project_name(
            self._view.ui.list_delete_projects_view.model(),
            self._view.ui.txt_delete_search,
            self._view.ui.lbl_delete_status_search,
            self._view.ui.txt_delete_selected_projects,
        )

    def select_project_from_delete_list(self) -> None:
        """Selects a project from the project list on the delete page."""
        try:
            self._view.ui.txt_delete_selected_projects.setText(self._view.ui.list_delete_projects_view.model().data
                                                               (self._view.ui.list_delete_projects_view.currentIndex(),
                                                                Qt.DisplayRole))

        except AttributeError:
            self._view.ui.txt_delete_selected_projects.setText("")

    def activate_delete_button(self) -> None:
        """Activates the delete button."""
        if self._view.ui.txt_delete_selected_projects.text() == "":
            self._view.ui.btn_delete_delete_project.setEnabled(False)
        else:
            self._view.ui.btn_delete_delete_project.setEnabled(True)

    def delete_project(self) -> None:
        """Deletes an existing project."""
        # popup message which warns the user that the selected project gets deleted
        tmp_dialog = custom_message_box.CustomMessageBoxDelete(
            "Are you sure you want to delete this project?",
            "Delete Project",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        response: bool = tmp_dialog.response
        tmp_index = self._view.ui.list_delete_projects_view.currentIndex()
        if response is True:
            os.remove(self._view.ui.list_delete_projects_view.model().data(tmp_index, enums.ModelEnum.FILEPATH_ROLE))  # TODO: throws permission error
            self._view.ui.list_delete_projects_view.model().removeRow(tmp_index.row())  # removes item from model
            self.restore_default_view()
        else:
            constants.PYSSA_LOGGER.info("No project has been deleted. No changes were made.")

    def _close(self):
        self._view.close()
