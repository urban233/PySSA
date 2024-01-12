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
"""Module for the About Dialog."""
import os
import glob


from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pyssa.gui.ui.forms.auto_generated import auto_open_project_view
from pyssa.gui.ui.styles import styles
from pyssa.util import constants, input_validator


class OpenProjectView(QtWidgets.QDialog):
    """Class representing an About dialog."""

    string_model = QtCore.QStringListModel()
    return_value = QtCore.pyqtSignal(str)

    def __init__(self, parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            parent: The parent.
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = auto_open_project_view.Ui_Dialog()
        self.ui.setupUi(self)
        self._initialize_ui()
        self._connect_all_ui_elements_to_slot_functions()

    def _initialize_ui(self) -> None:
        """Initialize the UI elements."""
        self._list_all_projects()
        self.ui.lbl_open_status_search.setText("")
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Open project")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

    def _list_all_projects(self) -> None:
        """Lists all projects."""
        xml_pattern = os.path.join(constants.DEFAULT_WORKSPACE_PATH, '*.xml')
        self.string_model.setStringList(
            [os.path.basename(file).replace(".xml", "") for file in glob.glob(xml_pattern)]  # Filters the workspace for all project files based on the xml extension
        )
        self.ui.projects_list_view.setModel(self.string_model)

    def _connect_all_ui_elements_to_slot_functions(self):
        self.ui.btn_open_project.clicked.connect(self._open_selected_project)
        self.ui.btn_cancel.clicked.connect(self._close)
        self.ui.txt_open_search.textChanged.connect(self._validate_open_search)
        self.ui.txt_open_selected_project.textChanged.connect(self._activate_open_button)
        self.ui.projects_list_view.clicked.connect(self._select_project_from_open_list)

    def _close(self):
        self.close()

    def _open_selected_project(self):
        self.return_value.emit(self.ui.txt_open_selected_project.text())
        self.close()

    def _validate_open_search(self) -> None:
        """Validates the input of the project name in real-time."""
        projects_list_view = self.ui.projects_list_view

        # Deselect any current item in the list view
        if projects_list_view.currentIndex().isValid():
            projects_list_view.model().itemFromIndex(projects_list_view.currentIndex()).setSelected(False)

        # Assuming validate_search_input is a static method
        input_validator.InputValidator.validate_search_input(
            projects_list_view.model(),
            self.ui.txt_open_search,
            self.ui.lbl_open_status_search,
            self.ui.txt_open_selected_project,
        )

    def _select_project_from_open_list(self) -> None:
        """Sets the selected project name in the text box."""
        try:
            self.ui.txt_open_selected_project.setText(self.ui.projects_list_view.model().data(self.ui.projects_list_view.currentIndex(), Qt.DisplayRole))
        except AttributeError:
            self.ui.txt_open_selected_project.setText("")

    def _activate_open_button(self) -> None:
        """Activates the open button."""
        if self.ui.txt_open_selected_project.text() == "":
            self.ui.btn_open_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_open_project)
        else:
            self.ui.btn_open_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_open_project)