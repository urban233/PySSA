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
"""Module for sequence table delegates."""
import logging

from PyQt5 import QtWidgets

from src.pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class InputCheckDelegate(QtWidgets.QStyledItemDelegate):
  """Custom input delegate to validate a certain input style."""

  def __init__(self, table_widget: QtWidgets.QTableWidget) -> None:
    """Constructor.

    Args:
        table_widget (QtWidgets.QTableWidget): The table widget to add the delegate to.
    """
    super().__init__()
    self.table_widget: QtWidgets.QTableWidget = table_widget
    self._first_text = ""

  def createEditor(self, parent, option, index) -> QtWidgets.QLineEdit:
    """Overrides the createEditor method of the QStyledItemDelegate class.

    Args:
        parent (QtWidgets.QWidget): The parent widget to which the editor will be added.
        option (bool): The option for the editor.
        index (QtCore.QModelIndex): The index for the editor.

    Returns:
        QtWidgets.QLineEdit: The created editor.
    """
    self.editor = QtWidgets.QLineEdit(parent)
    tmp_column = index.column()
    self._first_text = self.table_widget.currentItem().text()
    self.editor.textChanged.connect(
        lambda text, a_column=tmp_column: self.validateInput(text, a_column)
    )
    return self.editor

  def validateInput(self, text, a_column) -> None:
    """Validates the input based on the given text and column.

    Args:
        text: The input text to be validated.
        a_column: The column index used to determine the validation rules.
    """
    if a_column == 0:
      # validate name
      allowed_chars = (
          "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-"
      )
      tmp_corrected_text = "".join(
          char for char in text if char in allowed_chars
      )
    elif a_column == 1:
      # validate chain letter
      full_name = self.table_widget.item(
          self.table_widget.currentRow(), 0
      ).text()
      base_name = full_name.rsplit("_", 1)[0]
      not_allowed_chars: set = self.pool_chains_by_base_name(self.table_widget)[
          base_name
      ]
      not_allowed_chars.remove(self._first_text)
      tmp_corrected_text = "".join(
          char for char in text if char not in not_allowed_chars
      )
      tmp_corrected_text = tmp_corrected_text[:1]
    elif a_column == 2:
      # validate protein sequence
      allowed_chars = {
          "C",
          "D",
          "S",
          "Q",
          "K",
          "I",
          "P",
          "T",
          "F",
          "N",
          "G",
          "H",
          "L",
          "R",
          "W",
          "A",
          "V",
          "E",
          "Y",
          "M",
      }
      tmp_corrected_text = "".join(
          char for char in text if char in allowed_chars
      )
    else:
      tmp_corrected_text = ""
    self.editor.setText(tmp_corrected_text)

  def pool_chains_by_base_name(
      self, a_table_widget: QtWidgets.QTableWidget
  ) -> dict:
    """Connects the chains with their base names.

    Args:
        a_table_widget (QtWidgets.QTableWidget): The table widget containing the data.

    Returns:
        dict: A dictionary with the base names as keys and sets of chains as values.

    Description:
        This method takes a table widget as input and iterates through each row of the widget.
        For each row, it extracts the full name, base name, and chain from the respective columns.
        It then uses the base name as a key and adds the chain to a set of chains in a dictionary.
        If the base name already exists in the dictionary, the chain is added to the existing set,
        otherwise a new set is created for the base name.
    """
    pool = {}
    rows = a_table_widget.rowCount()
    for i in range(rows):
      full_name = a_table_widget.item(i, 0).text()
      base_name = full_name.rsplit("_", 1)[0]
      chain = a_table_widget.item(i, 1).text()
      if base_name not in pool:
        pool[base_name] = set()
      pool[base_name].add(chain)
    return pool
