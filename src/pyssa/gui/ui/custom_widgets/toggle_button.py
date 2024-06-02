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
"""Module for the custom toggle button widget."""
import logging
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt, pyqtSignal

from src.pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ToggleButton(QtWidgets.QCheckBox):
  """A custom toggle button widget based on a QCheckBox."""

  def __init__(
      self,
      width: int = 40,
      bg_color: str = "#d5d5d5",
      circle_color: str = "#5a5a5a",
      active_circle_color: str = "#fff",
      active_color: str = "#367AF6",
  ) -> None:
    """Initializes the checkbox with the given parameters.

    Args:
        width (int): The width of the checkbox.
        bg_color (str): The background color of the checkbox.
        circle_color (str): The color of the empty circle.
        active_circle_color (str): The color of the filled circle when checkbox is active.
        active_color (str): The color of the checkbox when active.
    """
    QtWidgets.QCheckBox.__init__(self)
    # default parameters
    self.setFixedSize(width, 19)
    self.setCursor(Qt.PointingHandCursor)

    # colors
    self._bg_color = bg_color
    self._circle_color = circle_color
    self._active_circle_color = active_circle_color
    self._active_color = active_color

  def hitButton(self, pos) -> bool:
    """Overrides the QCheckBox.hitButton method.

    Args:
        pos: QPoint object representing the position of the button press

    Returns:
        True if the button was hit, False otherwise.
    """
    return self.contentsRect().contains(pos)

  def paintEvent(self, e) -> None:
    """Overrides the QCheckBox.paintEvent method.

    Args:
        e: QPaintEvent object representing the event that triggered the painting

    This method is called automatically whenever the widget needs to be repainted. It is used to paint the custom appearance of the widget.
    """
    # set painter
    painter = QtGui.QPainter(self)
    painter.setRenderHint(QtGui.QPainter.Antialiasing)
    painter.setPen(Qt.NoPen)
    rectangle = QtCore.QRect(0, 0, self.width(), self.height())

    if not self.isChecked():
      painter.setBrush(QtGui.QColor(self._bg_color))
      painter.drawRoundedRect(
          0,
          0,
          rectangle.width(),
          self.height(),
          self.height() / 2,
          self.height() / 2,
      )
      painter.setBrush(QtGui.QColor(self._circle_color))
      painter.drawEllipse(3, 2, 15, 15)
    else:
      painter.setBrush(QtGui.QColor(self._active_color))
      painter.drawRoundedRect(
          0,
          0,
          rectangle.width(),
          self.height(),
          self.height() / 2,
          self.height() / 2,
      )
      painter.setBrush(QtGui.QColor(self._active_circle_color))
      painter.drawEllipse(self.width() - 18, 2, 15, 15)
    painter.end()


class ToggleWidget(QtWidgets.QWidget):
  """A widget that combines the toggle button with an on/off label."""

  toggleChanged = pyqtSignal(bool)
  """A custom signal that should be used instead of the stateChanged signal."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.toggle_button = ToggleButton()
    self.toggle_label = QtWidgets.QLabel("Off  ")

    self._layout = QtWidgets.QHBoxLayout()
    self._layout.addWidget(self.toggle_label)
    self._layout.addWidget(self.toggle_button)
    self.setLayout(self._layout)

    self.toggle_button.stateChanged.connect(self.switch_toggle_label_text)

  def switch_toggle_label_text(self) -> None:
    """Toggles the text of the label based on the state of the toggle button.

    If the toggle button is checked, the label text will be set to "On",
    otherwise it will be set to "Off". It also emits the toggleChanged signal
    with the appropriate state.
    """
    if self.toggle_button.isChecked():
      self.toggle_label.setText("On  ")
      self.toggleChanged.emit(True)
    else:
      self.toggle_label.setText("Off  ")
      self.toggleChanged.emit(False)
