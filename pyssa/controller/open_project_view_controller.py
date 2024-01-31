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

import glob
import os
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.util import input_validator, constants


class OpenProjectViewController(QtCore.QObject):
    """Class for the Open Project View Controller."""
    return_value = QtCore.pyqtSignal(str)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_open_view()
        self._fill_projects_list_view()
        self._connect_all_ui_elements_to_slot_functions()
        self.restore_default_view()

    def restore_default_view(self) -> None:
        self._view.ui.txt_open_selected_project.clear()

    def _fill_projects_list_view(self) -> None:
        """Lists all projects."""
        self._view.ui.projects_list_view.setModel(self._interface_manager.get_workspace_projects())

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.txt_open_search.textChanged.connect(self._validate_open_search)
        self._view.ui.projects_list_view.clicked.connect(self._select_project_from_open_list)
        self._view.ui.txt_open_selected_project.textChanged.connect(self._activate_open_button)
        self._view.ui.btn_open_project.clicked.connect(self._open_selected_project)
        self._view.ui.projects_list_view.doubleClicked.connect(self._open_selected_project)

    def _validate_open_search(self) -> None:
        """Validates the input of the project name in real-time."""
        projects_list_view = self._view.ui.projects_list_view

        # Deselect any current item in the list view
        if projects_list_view.currentIndex().isValid():
            projects_list_view.model().itemFromIndex(projects_list_view.currentIndex()).setSelected(False)

        # Assuming validate_search_input is a static method
        input_validator.InputValidator.validate_search_input(
            projects_list_view.model(),
            self._view.ui.txt_open_search,
            self._view.ui.lbl_open_status_search,
            self._view.ui.txt_open_selected_project,
        )

    def _select_project_from_open_list(self) -> None:
        """Sets the selected project name in the text box."""
        try:
            self._view.ui.txt_open_selected_project.setText(
                self._view.ui.projects_list_view.model().data(self._view.ui.projects_list_view.currentIndex(),
                                                              Qt.DisplayRole))
        except AttributeError:
            self._view.ui.txt_open_selected_project.setText("")

    def _activate_open_button(self) -> None:
        """Activates the open button."""
        if self._view.ui.txt_open_selected_project.text() == "":
            self._view.ui.btn_open_project.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_open_project)
        else:
            self._view.ui.btn_open_project.setEnabled(True)
            styles.color_button_ready(self._view.ui.btn_open_project)

    def _open_selected_project(self):
        self._view.close()
        self.return_value.emit(self._view.ui.txt_open_selected_project.text())
