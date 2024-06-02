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
"""Module for all task workers which run in separate threads."""
import logging

from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5 import QtWidgets
from src.pyssa.internal.data_structures.data_classes import prediction_protein_info
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import prediction_util
from src.pyssa.internal.prediction_engines import esmfold

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

"""
This module is currently not in use!
It is an experimental module for using and integrating ESMFold with PySSA. 
"""


class EsmFoldWorker(QObject):
  """Class for a worker which runs the ESM fold prediction."""

  finished = pyqtSignal()
  progress = pyqtSignal(int)
  return_value = pyqtSignal(list)
  table: QtWidgets.QTableWidget

  def __init__(self, table_prot_to_predict: QtWidgets.QTableWidget) -> None:
    """Constructor.

    Args:
        table_prot_to_predict: a QTableWidget containing the proteins to predict.
    """
    super().__init__()
    self.table = table_prot_to_predict

  def run(self) -> None:
    """This function is a reimplementation of the QObject run method."""
    predictions: list[prediction_protein_info.PredictionProteinInfo] = (
        prediction_util.get_prediction_name_and_seq_from_table(self.table)
    )
    output = esmfold.EsmFold(predictions).run_prediction()
    self.return_value.emit(output)
    self.finished.emit()
