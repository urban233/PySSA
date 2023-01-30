#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""This module contains helper function for the prediction process."""
import pathlib
import logging
from pyssa.util import constants
from pyssa.io_pyssa import filesystem_helpers
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def get_prediction_name_and_seq_from_table(table) -> list[tuple[str, str]]:
    """This function gets the names and sequences of the table which stores the predictions to run

    Args:
        table:
            pyqt table which contains the proteins to predict
    Returns:
        list of tuples with name and sequence
    """
    # list which consists of tuples of the protein name and protein sequence
    predictions: list[tuple[str, str]] = []
    for i in range(table.rowCount()):
        tmp_name = table.verticalHeaderItem(i).text()
        tmp_seq = table.item(i, 1).text()
        predictions.append((tmp_name, tmp_seq))
    return predictions


def get_relaxed_rank_1_pdb_file(protein_name_seq_tuples) -> list[tuple]:
    """This function gets the prediction models which were relaxed and ranked number one

    Args:
        protein_name_seq_tuples:
            list of tuples which consists of the protein name and sequence
    Returns:
        a list of tuples with the name of the modelled protein and the actual filename
    """
    # <editor-fold desc="Checks">
    if protein_name_seq_tuples[0][0] == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if protein_name_seq_tuples[0][1] == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    prediction_results = filesystem_helpers.create_generic_dictionary_from_directory(pathlib.Path(constants.PREDICTION_PDB_DIR))
    filenames = []
    for tmp_prediction in protein_name_seq_tuples:
        filename = [key for key, value in prediction_results.items() if value == f"{tmp_prediction[0]}_relaxed_rank_1"]
        if len(filename) == 1:
            filenames.append((tmp_prediction, filename[0]))
    return filenames
