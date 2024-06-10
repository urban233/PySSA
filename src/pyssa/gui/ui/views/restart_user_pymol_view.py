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
"""Module for the add sequence view."""
from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from src.pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui.ui.forms.auto_generated import auto_restart_user_pymol_view
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants

global_var_add_model = ("", False)


class RestartUserPyMOLView(QtWidgets.QDialog):
  """Class for a dialog for the restart user pymol."""

  return_value = pyqtSignal(tuple)
  """A signal that transfers the message that all windows should be closed."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        parent: The parent.
    """
    QtWidgets.QDialog.__init__(self, parent)
    # build ui object
    self.ui = auto_restart_user_pymol_view.Ui_Dialog()
    self.ui.setupUi(self)
    self.ui.btn_close_pyssa.clicked.connect(self._close_pyssa)
    pixmap = QtGui.QPixmap(
        r"C:\ProgramData\IBCI\PySSA\bin\PySSA\assets\images\splash_screen_logo_002.png"
    )
    scaled_pixmap = pixmap.scaled(
        600,
        600,
        aspectRatioMode=Qt.KeepAspectRatio,
        transformMode=Qt.SmoothTransformation,
    )
    # Set the scaled pixmap to the QLabel
    self.ui.lbl_logo.setPixmap(scaled_pixmap)
    self.ui.lbl_logo.setAlignment(Qt.AlignCenter)
    self.ui.lbl_message.setText("PyMOL needs to be restarted.")
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    styles.set_stylesheet(self)
    self.ui.lbl_message.setStyleSheet("font-size: 16px;")
    self.setWindowFlags(
        self.windowFlags()
        ^ Qt.WindowContextHelpButtonHint
        ^ Qt.WindowCloseButtonHint
    )
    self.setWindowTitle("Crash")
    self.setModal(True)

  def _close_pyssa(self) -> None:
    """Sends signal to close all windows related to PySSA."""
    tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "Do you want to close all windows related to PySSA?",
        "Close PySSA",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
    )
    tmp_dialog.exec_()
    if tmp_dialog.response:
      self.return_value.emit((True, 0))
