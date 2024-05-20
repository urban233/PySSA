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
"""Module for the delete project view controller."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import constants, enums, ui_util, exception
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class DeleteProjectViewController(QtCore.QObject):
    """Class for the DeleteProjectViewController."""

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
        self._view = the_interface_manager.get_delete_view()
        self._fill_projects_list_view()
        self._connect_all_ui_elements_to_slot_functions()
        self.restore_default_view()

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
        """Opens help dialog."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked.")
        self.open_help("help/project/delete_project/")

    def restore_default_view(self) -> None:
        """Restores the default UI."""
        self._view.ui.label_31.hide()
        self._view.ui.txt_delete_search.setPlaceholderText("Search")
        self._view.ui.txt_delete_search.clear()
        self._view.ui.txt_delete_selected_projects.clear()
        self._view.ui.btn_delete_delete_project.setEnabled(False)

    def _fill_projects_list_view(self) -> None:
        """Lists all projects."""
        self._view.ui.list_delete_projects_view.setModel(self._interface_manager.get_workspace_model())

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        """Connects all UI elements to their corresponding slot functions in the class."""
        self._view.ui.txt_delete_search.textChanged.connect(self.validate_delete_search)
        self._view.ui.list_delete_projects_view.clicked.connect(self.select_project_from_delete_list)
        self._view.ui.txt_delete_selected_projects.textChanged.connect(self.activate_delete_button)
        self._view.ui.btn_delete_delete_project.clicked.connect(self.delete_project)
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

    def validate_delete_search(self) -> None:
        """Validates the input of the project name in real-time."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
        projects_list_view = self._view.ui.list_delete_projects_view
        ui_util.select_matching_string_in_q_list_view(
            self._view.ui.txt_delete_search.text(),
            projects_list_view,
            self._view.ui.txt_delete_selected_projects,
        )

    def select_project_from_delete_list(self) -> None:
        """Selects a project from the project list on the delete page."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A project from the list of existing projects was clicked.")
        try:
            if len(self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()) == 1:
                self._view.ui.txt_delete_selected_projects.setText(self._view.ui.list_delete_projects_view.model().data
                                                                   (self._view.ui.list_delete_projects_view.currentIndex(),
                                                                    Qt.DisplayRole))
            elif len(self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()) > 1:
                self._view.ui.txt_delete_selected_projects.setText("Multiple projects selected.")
            else:
                self._view.ui.txt_delete_selected_projects.setText("")
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
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete' button was clicked.")

        if len(self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()) > 1:
            tmp_dialog = custom_message_box.CustomMessageBoxDelete(
                "Are you sure you want to delete these projects?",
                "Delete Projects",
                custom_message_box.CustomMessageBoxIcons.WARNING.value,
            )
            tmp_dialog.exec_()
            tmp_indexes = self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()
            tmp_filepaths_with_row_numbers = []
            for tmp_index in tmp_indexes:
                tmp_filepaths_with_row_numbers.append(
                    (self._view.ui.list_delete_projects_view.model().data(tmp_index, enums.ModelEnum.FILEPATH_ROLE),
                     tmp_index.row()),
                )
            tmp_filepaths_with_row_numbers.sort(reverse=True)
        elif len(self._view.ui.list_delete_projects_view.selectionModel().selectedIndexes()) == 1:
            tmp_dialog = custom_message_box.CustomMessageBoxDelete(
                "Are you sure you want to delete this project?",
                "Delete Project",
                custom_message_box.CustomMessageBoxIcons.WARNING.value,
            )
            tmp_dialog.exec_()
            tmp_filepaths_with_row_numbers = [
                (
                    self._view.ui.list_delete_projects_view.model().data(
                    self._view.ui.list_delete_projects_view.currentIndex(), enums.ModelEnum.FILEPATH_ROLE),
                    self._view.ui.list_delete_projects_view.currentIndex().row(),
                ),
            ]
        else:
            return

        response: bool = tmp_dialog.response
        if response is True:
            for tmp_filepath, tmp_row in tmp_filepaths_with_row_numbers:
                try:
                    os.remove(tmp_filepath)  # TODO: throws permission error, throws FileNotFound if more than one project gets deleted and afterwards selected
                except PermissionError:
                    tmp_dialog = custom_message_box.CustomMessageBoxOk(
                        "The project cannot be deleted, due to a permission error. Restart the application and try again.",
                        "Delete Project",
                        custom_message_box.CustomMessageBoxIcons.ERROR.value,
                    )
                    tmp_dialog.exec_()
                self._view.ui.list_delete_projects_view.model().removeRow(tmp_row)  # removes item from model
            self.restore_default_view()
        else:
            constants.PYSSA_LOGGER.info("No project has been deleted. No changes were made.")
