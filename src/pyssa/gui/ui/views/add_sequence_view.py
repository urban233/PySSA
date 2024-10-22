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
"""Module for the add sequence view."""
from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from src.pyssa.gui.ui.forms.auto_generated import auto_add_sequence_view
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants

global_var_add_model = ("", False)


class AddSequenceView(QtWidgets.QDialog):
  """Dialog for adding a sequence."""

  return_value = pyqtSignal(tuple)
  """A pyqtsignal that is used to hand-over the sequence information."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        parent: The parent.
    """
    QtWidgets.QDialog.__init__(self, parent)
    # build ui object
    self.ui = auto_add_sequence_view.Ui_Dialog()
    self.ui.setupUi(self)

    self._initalize_ui()

    self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
    self.ui.btn_help.setIconSize(
        self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30))
    )
    self.ui.btn_help.setText("")
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    styles.set_stylesheet(self)
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )
    self.setWindowTitle("Add Protein Sequence")
    self.setWindowModality(Qt.WindowModal)

  def _initalize_ui(self):
    self.ui.btn_next.setEnabled(False)
    self.ui.lbl_status.setText("")
    self.ui.lbl_protein_seq.hide()
    self.ui.le_protein_seq.hide()
    self.ui.btn_back.hide()
    self.ui.btn_add.setEnabled(False)
    styles.color_bottom_frame_button(self.ui.btn_add)
