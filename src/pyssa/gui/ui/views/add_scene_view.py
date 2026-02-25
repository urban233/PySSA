#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for the add scene view."""
from src.pyssa.gui.qt import Qt
from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets

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

    self.top_frame = QtWidgets.QFrame()
    self.lbl_description = QtWidgets.QLabel("Enter a new scene name")
    self.line_edit_scene_name = custom_line_edit.CustomLineEditWithMaxLength(30)
    self.lbl_status = QtWidgets.QLabel("")
    self.bottom_frame = QtWidgets.QFrame()
    self.btn_add_scene = QtWidgets.QPushButton("Add")
    self.btn_cancel = QtWidgets.QPushButton("Cancel")

    self.layout_user_input = QtWidgets.QVBoxLayout(self.top_frame)
    self.layout_user_input.addWidget(self.lbl_description)
    self.layout_user_input.addWidget(self.line_edit_scene_name)
    self.layout_user_input.addWidget(self.lbl_status)
    self.layout_user_input.setContentsMargins(8, 8, 8, 2)

    self.layout_confirmation = QtWidgets.QHBoxLayout(self.bottom_frame)
    self.layout_confirmation.addStretch()
    self.layout_confirmation.addWidget(self.btn_add_scene)
    self.layout_confirmation.addWidget(self.btn_cancel)

    self.layout_complete = QtWidgets.QVBoxLayout()
    self.layout_complete.setContentsMargins(0, 0, 0, 0)
    self.layout_complete.addWidget(self.top_frame)
    self.layout_complete.addWidget(self.bottom_frame)

    self.setLayout(self.layout_complete)

    self.setMaximumSize(600, 70)
    self.setMinimumWidth(450)
    self.resize(450, 70)

    self.btn_cancel.clicked.connect(self.close)
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    # styles.set_stylesheet(self)
    self.bottom_frame.setStyleSheet(
      """
      QFrame {
        background-color: #f7f8fa;
        border-style: solid;
        border-width: 1px;
        border-radius: 6px;
        border-color: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #f9f9f9, stop:1 #f0f0f0);;
        border-top-color: #ebecf0;
        border-top-left-radius: 0px;
        border-top-right-radius: 0px;
    }
      """
    )
    self.lbl_status.setStyleSheet(
      """
      QLabel {
        background-color: #f3f3f3;
        border: none;
        color: #ba1a1a; 
        font-size: 11px;
      }
      """
    )
    styles.color_bottom_frame_button(self.btn_add_scene)
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )
    self.setWindowTitle("Add New PyMOL Scene")
    self.setModal(True)
