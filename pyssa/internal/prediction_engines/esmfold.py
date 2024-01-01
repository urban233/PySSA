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
"""This module contains the esmfold class."""
import os
import logging
import pathlib
import shutil

import requests

from pyssa.internal.data_structures.data_classes import prediction_protein_info
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class EsmFold:
    """This class contains information about the ESMFold prediction engine."""

    """
    a list of protein sequences which should get predicted
    """
    protein_prediction_infos: list[prediction_protein_info.PredictionProteinInfo]

    def __init__(self, protein_prediction_infos):
        """Constructor.

        Args:
            protein_prediction_infos: list of valid protein sequences with their protein names

        Raises:
            ValueError: if protein sequences list is None
            ValueError: if protein sequences list is empty
        """
        self.protein_prediction_infos = protein_prediction_infos

    def predict_single_protein(self, sequence, pdb_filepath):
        """This function predicts a protein structure with the ESMFold API

        Args:
            sequence: valid protein sequence which gets predicted
            pdb_filepath: filepath of the pdb file which should be created through the prediction

        Raises:
            ValueError: if the filepath is only a filename
            NotADirectoryError: if the directory where the file should be created is not found
            ConnectionError: if the ESMFold server throws the error INTERNAL SERVER ERROR
        """

        # <editor-fold desc="Checks">
        if os.path.dirname(pdb_filepath) == "":
            logger.error("Invalid filepath, add dirname to complete the filepath!")
            raise ValueError("Invalid filepath, add dirname to complete the filepath!")
        if not os.path.exists(os.path.dirname(pdb_filepath)):
            logger.error(f"There is no directory called: {os.path.dirname(pdb_filepath)}")
            raise NotADirectoryError(f"There is no directory called: {os.path.dirname(pdb_filepath)}")

        # </editor-fold>

        url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
        response = requests.post(url, data=sequence)
        print(response.text)  # only for debug purpose
        error_count = 0
        if response.text == "INTERNAL SERVER ERROR":
            while response.text != "INTERNAL SERVER ERROR" or error_count > 10:
                logger.debug(f"No. of errors: {error_count}")
                response = requests.post(url, data=sequence)
                error_count += 1
        if response.text != "INTERNAL SERVER ERROR":
            with open(pdb_filepath, "w") as f:
                f.write(response.text)
        else:
            raise ConnectionError("Too many internal server errors.")

    def run_prediction(self) -> list:
        """This function starts the esm fold prediction process.

        Returns:
            a list with prediction protein infos which could not be predicted.
        """

        # <editor-fold desc="Checks">
        if self.protein_prediction_infos is None:
            logger.error("Protein sequences must not be None")
            raise ValueError("Protein sequences must not be None")
        if len(self.protein_prediction_infos) == 0:
            logger.error("Protein sequences list is empty!")
            raise ValueError("Protein sequences list is empty!")

        if not os.path.exists(constants.ESMFOLD_DIR):
            os.mkdir(constants.ESMFOLD_DIR)
        else:
            shutil.rmtree(constants.ESMFOLD_DIR)
            os.mkdir(constants.ESMFOLD_DIR)
        if not os.path.exists(constants.ESMFOLD_PDB_DIR):
            os.mkdir(constants.ESMFOLD_PDB_DIR)
        else:
            shutil.rmtree(constants.ESMFOLD_PDB_DIR)
            os.mkdir(constants.ESMFOLD_PDB_DIR)

        # </editor-fold>

        i = 1
        failed_predictions = []
        for tmp_protein_prediction_info in self.protein_prediction_infos:
            try:
                self.predict_single_protein(
                    tmp_protein_prediction_info.sequences[0],
                    str(pathlib.Path(f"{constants.ESMFOLD_PDB_DIR}/{tmp_protein_prediction_info.name}.pdb")),
                )
            except ConnectionError:
                failed_predictions.append(tmp_protein_prediction_info.sequences[0])
                logger.error(
                    f"Prediction with protein name: {tmp_protein_prediction_info.name} failed. "
                    f"Prediction will be retried one more time after all others finished."
                )
            finally:
                logger.info(
                    f"Prediction with protein name: {tmp_protein_prediction_info.name} finished. "
                    f"Run {i} of {len(self.protein_prediction_infos)} finished."
                )
                i += 1
        failed_multiple_attempts = []
        if len(failed_predictions) > 0:
            for tmp_protein_prediction_info in failed_predictions:
                logger.info("Retry previously failed predictions ...")
                try:
                    self.predict_single_protein(
                        tmp_protein_prediction_info.sequences[0],
                        str(pathlib.Path(f"{constants.ESMFOLD_PDB_DIR}/{tmp_protein_prediction_info.name}.pdb")),
                    )
                except ConnectionError:
                    failed_multiple_attempts.append(tmp_protein_prediction_info)
                    logger.error(
                        f"Prediction with protein name: {tmp_protein_prediction_info.name} failed. "
                        f"Prediction cannot be completed at the moment."
                    )
            logger.debug(f"In the first run these predictions failed: {failed_predictions}")
            logger.debug(f"In the second run these predictions failed: {failed_multiple_attempts}")

        return failed_multiple_attempts
