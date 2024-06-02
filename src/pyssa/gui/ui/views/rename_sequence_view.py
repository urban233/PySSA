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
"""Module for the rename sequence view."""
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from src.pyssa.gui.ui.forms.auto_generated import auto_rename_object_view
from src.pyssa.gui.ui.custom_widgets import custom_line_edit
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants


class RenameSequenceView(QtWidgets.QDialog):
  """Class for the rename sequence dialog."""

  dialogClosed = QtCore.pyqtSignal(tuple)
  """A signal indicating that the dialog is closed."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor."""
    QtWidgets.QDialog.__init__(self, parent)
    # build ui object
    self.ui = auto_rename_object_view.Ui_Dialog()
    self.ui.setupUi(self)

    self.ui.lbl_description.setText("Enter a new sequence name")
    self.ui.lbl_status.setStyleSheet("color: #ba1a1a; font-size: 11px;")
    styles.color_bottom_frame_button(self.ui.btn_rename)

    self.setMaximumSize(600, 80)
    self.setMinimumWidth(250)
    self.resize(450, 80)
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    styles.set_stylesheet(self)
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )
    self.setWindowTitle("Rename Selected Sequence")
    self.setModal(True)
