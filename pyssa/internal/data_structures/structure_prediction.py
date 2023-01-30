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
"""Module for structure prediction class"""
import os
import pathlib
import shutil
import pymol
import logging
from pyssa.internal.data_structures import protein
from pyssa.internal.portal import pymol_io
from pyssa.internal.data_processing import data_transformer
from pyssa.internal.prediction_engines import colabbatch
from pyssa.util import constants
from pyssa.util import types
from pyssa.util import prediction_util
from pyssa.logging_pyssa import log_handlers
from pyssa.logging_pyssa import loggers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class StructurePrediction:
    """This class is used to organize the structure prediction process."""

    # <editor-fold desc="Class attributes">
    """
    a list of tuples of protein name and sequence 
    """
    predictions: list[tuple[str, str]]
    """
    the configuration settings for the prediction
    """
    prediction_configuration: types.PREDICTION_CONFIG
    """
    the current project in use 
    """
    project: types.PROJECT

    # </editor-fold>

    def __init__(self,
                 predictions: list[tuple[str, str]],
                 prediction_config: types.PREDICTION_CONFIG,
                 current_project: types.PROJECT) -> None:
        self.predictions = predictions
        self.prediction_configuration = prediction_config
        self.project = current_project
        loggers.log_multiple_variable_values(logger, "Constructor", [("predictions", self.predictions),
                                                                     ("prediction_configuration", self.prediction_configuration),
                                                                     ("project", self.project)])

    @staticmethod
    def create_tmp_directories() -> None:
        """This function creates tmp directories in a scratch folder to organize prediction inputs and outputs

        """
        # TODO: is there a more elegant way to do it?
        if not os.path.exists(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions")):
            os.mkdir(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
        if not os.path.exists(constants.PREDICTION_FASTA_DIR):
            os.mkdir(constants.PREDICTION_FASTA_DIR)
        if not os.path.exists(constants.PREDICTION_PDB_DIR):
            os.mkdir(constants.PREDICTION_PDB_DIR)

    def create_fasta_files_for_prediction(self) -> None:
        """This function creates fasta file based on the predictions in the object

        """
        protein_sequences: list[types.PROTEIN_SEQUENCE] = data_transformer.transform_protein_name_seq_tuple_to_sequence_obj(self.predictions)
        logger.debug(f"Variable: protein_sequences; Value: {protein_sequences} in function create_fasta_files_for_prediction")
        for tmp_seq_to_predict in protein_sequences:
            tmp_seq_to_predict.write_fasta_file(constants.PREDICTION_FASTA_DIR)
        if len(os.listdir(constants.PREDICTION_FASTA_DIR)) == 0:
            logger.error("No fasta files were created!")
            raise FileNotFoundError

    def run_prediction(self) -> None:
        """This function runs a structure prediction.

        """
        colabbatch.Colabbatch(self.prediction_configuration).run_prediction()

    def move_best_prediction_models(self):
        """This function moves the best prediction model(s) to the project directory

        """
        best_prediction_models: list[tuple] = prediction_util.get_relaxed_rank_1_pdb_file(self.predictions)
        for tmp_prediction in best_prediction_models:
            src = pathlib.Path(f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[1]}")
            dest = pathlib.Path(f"{self.project.get_pdb_path()}/{tmp_prediction[1]}")
            shutil.copy(src, dest)
            os.rename(f"{self.project.get_pdb_path()}/{tmp_prediction[1]}",
                      f"{self.project.get_pdb_path()}/{tmp_prediction[0]}.pdb")
            self.project.add_existing_protein(protein.Protein(tmp_prediction[0], pathlib.Path(self.project.get_pdb_path())))
        shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
        try:
            pymol_io.load_protein(self.project.proteins[0].filepath, self.project.proteins[0].filename,
                                  self.project.proteins[0].molecule_object)
        except pymol.CmdException:
            print("Loading the model failed.")
            return
        except FileNotFoundError:
            print("Prediction was unsuccessful")
            return
