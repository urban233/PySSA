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
"""Module for the use project view."""
import os
import glob

from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.gui.ui.forms.auto_generated import auto_use_project_view
from pyssa.gui.ui.styles import styles
from pyssa.util import constants, input_validator


class UseProjectView(QtWidgets.QDialog):
    """Class representing an Open dialog."""

    def __init__(self) -> None:
        """Constructor."""
        QtWidgets.QDialog.__init__(self)
        # build ui object
        self.ui = auto_use_project_view.Ui_Dialog()
        self.ui.setupUi(self)
        self._initialize_ui()
        self.resize(450, 600)
        self.setModal(True)

    def _initialize_ui(self) -> None:
        """Initialize the UI elements."""
        self.ui.lbl_use_status_search.setText("")
        styles.color_bottom_frame_button(self.ui.btn_use_create_new_project)
        self.ui.list_use_existing_projects.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")
        styles.color_bottom_frame_button(self.ui.btn_use_create_new_project)
        styles.set_stylesheet(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Use Project")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
