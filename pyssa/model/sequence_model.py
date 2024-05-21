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
"""Module contains the sequence model."""
import logging

from PyQt5 import QtGui

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class SequenceModel(QtGui.QStandardItemModel):
  """Contains the sequences of the project in form of a QStandardItemModel."""

  def __init__(self) -> None:
    """Constructor."""
    super(SequenceModel, self).__init__()

  def add_sequence(self, a_sequence: str) -> None:
    """Adds a sequence to the model.

    Args:
        a_sequence (str): The amino acid sequence to be added.

    Raises:
        exception.IllegalArgumentError: If a_sequence is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_sequence is None or a_sequence == "":
      logger.error("a_sequence is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_sequence is either None or an empty string."
      )

    # </editor-fold>

    i = self.rowCount()
    j = 0
    for tmp_amino_acid in a_sequence:
      self.setItem(i, j, QtGui.QStandardItem(tmp_amino_acid))
      logger.debug(self.item(i, j).text())
