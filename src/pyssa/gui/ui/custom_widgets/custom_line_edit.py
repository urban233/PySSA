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
"""Module for the custom line edit widget."""
from PyQt5 import QtWidgets
from PyQt5 import QtCore


class CustomLineEdit(QtWidgets.QLineEdit):
  """A custom line edit widget that allows only a subset of chars."""

  def keyPressEvent(self, event) -> None:
    """Overrides keyPressEvent of QLineEdit class."""
    # Get the key code
    key_text = event.text()

    # Check if the key is allowed, Backspace is allowed, or if it's an empty string (allowing empty input)
    allowed_chars = set(
        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-"
    )
    if (
        key_text not in allowed_chars
        and key_text != ""
        and event.key() != QtCore.Qt.Key_Backspace
    ):
      # Ignore the key event
      return

    # Call the base class implementation to handle other keys
    super(CustomLineEdit, self).keyPressEvent(event)


# class CustomLineEditForEnteringNumbers(QtWidgets.QLineEdit):
#     def keyPressEvent(self, event):
#         print("Hello")
#         # Get the key code
#         key_text = event.text()
#
#         # Check if the key is allowed, Backspace is allowed, or if it's an empty string (allowing empty input)
#         allowed_chars = set("0123456789")
#         if key_text not in allowed_chars and key_text != "" and event.key() != QtCore.Qt.Key_Backspace:
#             # Ignore the key event
#             return
#
#         # Call the base class implementation to handle other keys
#         super(CustomLineEditForEnteringNumbers, self).keyPressEvent(event)
