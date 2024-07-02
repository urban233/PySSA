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
"""Module for the custom label widget."""
from PyQt5 import QtWidgets
from PyQt5 import QtCore


class PermanentMessageLabel(QtWidgets.QLabel):
  """A custom label widget that can send a signal if the text changes."""

  textChanged = QtCore.pyqtSignal(str)
  """A signal indicating that the text changed."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor."""
    super().__init__(parent)

  def setText(self, text) -> None:
    """Overrides the setText method of the QLabel class.

    Args:
        text (str): The text to set.
    """
    super().setText(text)
    self.textChanged.emit(text)
