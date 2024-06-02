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
"""Module for the add scene view."""
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5 import QtWidgets

from src.pyssa.gui.ui.custom_widgets import custom_line_edit
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants


class AddSceneView(QtWidgets.QDialog):
  """Dialog for adding a scene."""

  dialogClosed = QtCore.pyqtSignal(tuple)
  """A signal indicating that the dialog is closed."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor."""
    QtWidgets.QDialog.__init__(self, parent)

    self.lbl_description = QtWidgets.QLabel("New PyMOL scene name")
    self.line_edit_scene_name = custom_line_edit.CustomLineEdit()
    self.lbl_status = QtWidgets.QLabel("")
    self.btn_add_scene = QtWidgets.QPushButton("Add")
    self.btn_cancel = QtWidgets.QPushButton("Cancel")

    self.layout_user_input = QtWidgets.QVBoxLayout()
    self.layout_user_input.addWidget(self.lbl_description)
    self.layout_user_input.addWidget(self.line_edit_scene_name)

    self.layout_confirmation = (
        QtWidgets.QHBoxLayout()
    )  # Use QHBoxLayout for the button
    self.layout_confirmation.addWidget(self.lbl_status)
    self.layout_confirmation.addStretch(
        1
    )  # Add stretchable space before the button
    self.layout_confirmation.addWidget(self.btn_add_scene)
    self.layout_confirmation.addWidget(self.btn_cancel)

    self.layout_complete = QtWidgets.QVBoxLayout()
    self.layout_complete.addLayout(self.layout_user_input)
    self.layout_complete.addLayout(self.layout_confirmation)

    self.setLayout(self.layout_complete)

    self.setMaximumSize(600, 80)
    self.setMinimumWidth(250)
    self.resize(450, 80)

    self.btn_cancel.clicked.connect(self.close)
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    styles.set_stylesheet(self)
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )
    self.setWindowTitle("Add New PyMOL Scene")
    self.setModal(True)
