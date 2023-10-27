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
import os
import pathlib
import logging
from pyssa.util import constants
from pyssa.internal.data_structures.data_classes import prediction_protein_info
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def get_prediction_name_and_seq_from_table(table) -> list[prediction_protein_info.PredictionProteinInfo]:
    """This function gets the names and sequences of the table which stores the predictions to run.

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
    # create an empty dictionary
    groups = {}

    # loop over the tuples in the list
    for t in predictions:
        key = t[0]
        value = t[1]
        # if the key is not already in the dictionary, add it with an empty list
        if key not in groups:
            groups[key] = []
        # add the tuple's second element to the list associated with the key
        groups[key].append(value)
    prediction_runs = []
    for tmp_prot_name, tmp_seqs in groups.items():
        prediction_runs.append(prediction_protein_info.PredictionProteinInfo(tmp_prot_name, tmp_seqs))
    return prediction_runs


def get_relaxed_rank_1_pdb_file(proteins_to_predict: list[prediction_protein_info.PredictionProteinInfo]) -> list[tuple]:
    """This function gets the prediction models which were relaxed and ranked number one.

    Args:
        proteins_to_predict:
            list of tuples which consists of the protein name and sequence
    Returns:
        a list of tuples with the name of the modelled protein and the actual filename
    """
    # <editor-fold desc="Checks">
    if proteins_to_predict[0].name == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not proteins_to_predict[0].sequences:
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    filenames = []
    logger.debug(proteins_to_predict)
    for tmp_protein_to_predict in proteins_to_predict:
        logger.debug(tmp_protein_to_predict)
        logger.info(os.listdir(pathlib.Path(constants.PREDICTION_PDB_DIR)))
        for tmp_filename in os.listdir(pathlib.Path(constants.PREDICTION_PDB_DIR)):
            if tmp_filename.find(f"{tmp_protein_to_predict.name}_relaxed_rank_001") != -1:
                filenames.append((tmp_protein_to_predict, tmp_filename))
                logger.debug(tmp_filename)
    logger.debug(filenames)
    return filenames
