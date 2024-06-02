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
from PyQt5.QtGui import QIcon
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants


class ProcessProgress(QtWidgets.QDialog):
    def __init__(self):
        super().__init__()

        self.progress_bar = QtWidgets.QProgressBar(self)
        self.label = QtWidgets.QLabel("Status: ")

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.progress_bar)
        layout.addWidget(self.label)

        self.setLayout(layout)
        self.setGeometry(100, 100, 600, 50)
        styles.set_stylesheet(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Progress")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

    def update_progress_bar(self, value):
        self.progress_bar.setValue(value)
