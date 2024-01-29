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
"""Module for the rename protein view."""
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5 import QtWidgets

from pyssa.gui.ui.custom_widgets import custom_line_edit
from pyssa.gui.ui.styles import styles
from pyssa.util import constants


class RenameProteinView(QtWidgets.QDialog):
    # Define a custom signal
    dialogClosed = QtCore.pyqtSignal(tuple)

    def __init__(self, parent=None) -> None:  # noqa: ANN001
        """Constructor."""
        QtWidgets.QDialog.__init__(self, parent)

        self.lbl_description = QtWidgets.QLabel("New protein name")
        self.line_edit_protein_name = custom_line_edit.CustomLineEdit()
        self.btn_rename_protein = QtWidgets.QPushButton("Rename")

        self.layout_user_input = QtWidgets.QVBoxLayout()
        self.layout_user_input.addWidget(self.lbl_description)
        self.layout_user_input.addWidget(self.line_edit_protein_name)

        self.layout_confirmation = QtWidgets.QHBoxLayout()  # Use QHBoxLayout for the button
        self.layout_confirmation.addStretch(1)  # Add stretchable space before the button
        self.layout_confirmation.addWidget(self.btn_rename_protein)

        self.layout_complete = QtWidgets.QVBoxLayout()
        self.layout_complete.addLayout(self.layout_user_input)
        self.layout_complete.addLayout(self.layout_confirmation)

        self.setLayout(self.layout_complete)

        self.setMaximumSize(600, 80)
        self.setMinimumWidth(250)
        self.resize(450, 80)

        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.setWindowTitle("Rename selected protein structure")

    def closeEvent(self, event):
        # Emit the custom signal when the window is closed
        self.dialogClosed.emit(("", False))
        event.accept()