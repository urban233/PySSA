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
"""Module for custom message boxes which are designed through subclassing QDialog."""
import enum
import logging

from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5 import QtWidgets

from src.pyssa.gui.ui.styles import styles
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class CustomMessageBoxIcons(enum.Enum):
  """Enumeration of custom message box icons."""

  INFORMATION = ":/icons/info_w200.png"
  WARNING = ":/icons/warning_w200.png"
  ERROR = ":/icons/error_w200.png"
  DANGEROUS = ":/icons/dangerous_w200.png"


class CustomMessageBox(QtWidgets.QDialog):
  """A custom message box that is based on a QDialog."""

  dialogClosed = QtCore.pyqtSignal(bool)
  """A signal indicating that the dialog is closed."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor."""
    QtWidgets.QDialog.__init__(self, parent)
    self.lbl_icon = QtWidgets.QLabel("Generic Icon Pixmap")
    self.lbl_description = QtWidgets.QLabel("Generic Message")
    self.btn_left = QtWidgets.QPushButton("Generic Left")
    self.btn_right = QtWidgets.QPushButton("Generic Right")

    # <editor-fold desc="Layouts">
    self.layout_description_msg = QtWidgets.QHBoxLayout()
    self.layout_description_msg.addWidget(self.lbl_icon)
    self.layout_description_msg.addWidget(self.lbl_description)
    self.layout_description_msg.addStretch(1)

    self.layout_buttons = (
        QtWidgets.QHBoxLayout()
    )  # Use QHBoxLayout for the button
    self.layout_buttons.addStretch(1)  # Add stretchable space before the button
    self.layout_buttons.addWidget(self.btn_left)
    self.layout_buttons.addWidget(self.btn_right)

    self.layout_complete = QtWidgets.QVBoxLayout()
    self.layout_complete.addLayout(self.layout_description_msg)
    self.layout_complete.addLayout(self.layout_buttons)

    self.setLayout(self.layout_complete)

    # </editor-fold>

    self.setMaximumSize(600, 80)
    self.setMinimumWidth(250)
    self.resize(450, 80)

    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    styles.set_stylesheet(self)
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )
    self.setWindowTitle("Generic Window Title")
    self.setModal(True)

  def closeEvent(self, event):
    """Overrides the closeEvent method of the QDialog class."""
    # Emit the custom signal when the window is closed
    self.dialogClosed.emit(False)
    event.accept()


class CustomMessageBoxDelete(CustomMessageBox):
  """A message box that contains a delete and cancel button."""

  response: bool
  """A boolean value indicating if the cancel button was clicked (set to False)."""

  def __init__(
      self, a_message: str, a_window_title: str, an_icon_path: str
  ) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        a_message (str): The message to be displayed in the dialog.
        a_window_title (str): The title of the window.
        an_icon_path (str): The path to the icon image file.

    Raises:
        exception.IllegalArgumentError: If any of the arguments is None.
    """
    # <editor-fold desc="Checks">
    if a_message is None:
      logger.error("a_message is None.")
      raise exception.IllegalArgumentError("a_message is None.")
    if a_window_title is None:
      logger.error("a_window_title is None.")
      raise exception.IllegalArgumentError("a_window_title is None.")
    if an_icon_path is None:
      logger.error("an_icon_path is None.")
      raise exception.IllegalArgumentError("an_icon_path is None.")

    # </editor-fold>

    super().__init__()

    self.response = False

    self.lbl_icon.setText("")
    self.lbl_icon.setPixmap(QtGui.QIcon(an_icon_path).pixmap(40, 40))
    self.lbl_description.setText(a_message)
    self.btn_left.setText("Delete")
    self.btn_left.setStyleSheet(
        """
            QPushButton {
                background-color: #ba1a1a; 
                color: white; 
                border: none;
            }
            QPushButton::pressed {
                background-color: #410002; 
                color: white; 
                border: none;
            }
        """
    )
    self.btn_right.setText("Cancel")
    self.setWindowTitle(a_window_title)

    self.btn_left.clicked.connect(self.__slot_left_button)
    self.btn_right.clicked.connect(self.__slot_right_button)

  def __slot_left_button(self) -> None:
    """Clicks and handles the left button.

    This method sets the response attribute to True and closes the current widget.
    """
    self.response = True
    self.close()

  def __slot_right_button(self) -> None:
    """Clicks and handles the right button.

    This method sets the response attribute to False and closes the current widget.
    """
    self.response = False
    self.close()


class CustomMessageBoxOk(CustomMessageBox):
  """A message box that contains an ok button."""

  response: bool
  """A boolean value indicating if the ok button was clicked (set to True)."""

  def __init__(
      self, a_message: str, a_window_title: str, an_icon_path: str
  ) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        a_message (str): The message to be displayed in the dialog.
        a_window_title (str): The title of the window.
        an_icon_path (str): The path to the icon image file.

    Raises:
        exception.IllegalArgumentError: If any of the arguments is None.
    """
    # <editor-fold desc="Checks">
    if a_message is None:
      logger.error("a_message is None.")
      raise exception.IllegalArgumentError("a_message is None.")
    if a_window_title is None:
      logger.error("a_window_title is None.")
      raise exception.IllegalArgumentError("a_window_title is None.")
    if an_icon_path is None:
      logger.error("an_icon_path is None.")
      raise exception.IllegalArgumentError("an_icon_path is None.")

    # </editor-fold>

    super().__init__()

    self.response = False

    self.lbl_icon.setText("")
    self.lbl_icon.setPixmap(QtGui.QIcon(an_icon_path).pixmap(40, 40))
    self.lbl_description.setText(a_message)
    self.btn_left.setText("OK")
    styles.color_bottom_frame_button(self.btn_left)
    self.btn_right.hide()
    self.setWindowTitle(a_window_title)

    self.btn_left.clicked.connect(self.__slot_left_button)

  def __slot_left_button(self) -> None:
    """Clicks and handles the left button.

    This method sets the response attribute to True and closes the current widget.
    """
    self.response = True
    self.close()


class CustomMessageBoxYesNo(CustomMessageBox):
  """A message box that contains a yes and no button."""

  response: bool
  """A boolean value indicating if the no button was clicked (set to False)."""

  def __init__(
      self, a_message: str, a_window_title: str, an_icon_path: str
  ) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        a_message (str): The message to be displayed in the dialog.
        a_window_title (str): The title of the window.
        an_icon_path (str): The path to the icon image file.

    Raises:
        exception.IllegalArgumentError: If any of the arguments is None.
    """
    # <editor-fold desc="Checks">
    if a_message is None:
      logger.error("a_message is None.")
      raise exception.IllegalArgumentError("a_message is None.")
    if a_window_title is None:
      logger.error("a_window_title is None.")
      raise exception.IllegalArgumentError("a_window_title is None.")
    if an_icon_path is None:
      logger.error("an_icon_path is None.")
      raise exception.IllegalArgumentError("an_icon_path is None.")

    # </editor-fold>

    super().__init__()

    self.response = False

    self.lbl_icon.setText("")
    self.lbl_icon.setPixmap(QtGui.QIcon(an_icon_path).pixmap(40, 40))
    self.lbl_description.setText(a_message)
    self.btn_left.setText("Yes")
    styles.color_bottom_frame_button(self.btn_left)
    self.btn_right.setText("No")
    self.setWindowTitle(a_window_title)

    self.btn_left.clicked.connect(self.__slot_left_button)
    self.btn_right.clicked.connect(self.__slot_right_button)

  def __slot_left_button(self) -> None:
    """Clicks and handles the left button.

    This method sets the response attribute to True and closes the current widget.
    """
    self.response = True
    self.close()

  def __slot_right_button(self) -> None:
    """Clicks and handles the right button.

    This method sets the response attribute to False and closes the current widget.
    """
    self.response = False
    self.close()
