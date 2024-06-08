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
"""Module for the results summary view."""
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from src.pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from src.pyssa.gui.ui.forms.auto_generated import auto_results_view
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants


class ResultsView(QtWidgets.QDialog):
  """Class representing an About dialog."""

  # string_model = QtCore.QStringListModel()
  # return_value = QtCore.pyqtSignal(str)

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        parent: The parent.
    """
    QtWidgets.QDialog.__init__(self, parent)
    # build ui object
    self.ui = auto_results_view.Ui_Dialog()
    self.ui.setupUi(self)
    self._initialize_ui()
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowModality(Qt.WindowModal)

  def _initialize_ui(self) -> None:
    """Initialize the UI elements."""
    styles.color_bottom_frame_button(self.ui.btn_view_plots)
    styles.color_bottom_frame_button(self.ui.btn_export_data)
    self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
    self.ui.btn_help.setIconSize(
        self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30))
    )
    self.ui.btn_help.setText("")
    styles.set_stylesheet(self)
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    self.setWindowTitle("Results Summary")
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )
