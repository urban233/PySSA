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
import glob

from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pyssa.gui.ui.forms.auto_generated import auto_open_project_view
from pyssa.gui.ui.styles import styles
from pyssa.util import constants, input_validator


class OpenProjectView(QtWidgets.QDialog):
    """Class representing an Open dialog."""

    def __init__(self) -> None:
        """Constructor."""
        QtWidgets.QDialog.__init__(self)
        # build ui object
        self.ui = auto_open_project_view.Ui_Dialog()
        self.ui.setupUi(self)
        self._initialize_ui()

    def _initialize_ui(self) -> None:
        """Initialize the UI elements."""
        self.ui.lbl_open_status_search.setText("")
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Open project")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)