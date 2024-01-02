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
"""Module for structure prediction class."""
import os
import pathlib
import shutil
import logging

from pyssa.internal.data_structures.data_classes import prediction_protein_info
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.internal.data_processing import data_transformer
from pyssa.internal.prediction_engines import colabbatch
from pyssa.util import constants
from pyssa.util import prediction_util
from pyssa.util import exception
from pyssa.io_pyssa import path_util
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class StructurePrediction:
    """This class is used to organize the structure prediction process."""

    # <editor-fold desc="Class attributes">
    """
    a list of tuples of protein name and sequence
    """
    predictions: list[prediction_protein_info.PredictionProteinInfo]
    """
    the configuration settings for the prediction
    """
    prediction_config: prediction_configuration.PredictionConfiguration
    """
    the current project in use
    """
    current_project: project.Project

    # </editor-fold>

    def __init__(
        self,
        predictions: list[prediction_protein_info.PredictionProteinInfo],
        prediction_config: prediction_configuration.PredictionConfiguration,
        current_project: "project.Project",
    ) -> None:
        """Constructor.

        Args:
            predictions: a list of prediction protein infos.
            prediction_config: the configuration settings for the prediction.
            current_project: the current project in use.
        """
        self.predictions: list[prediction_protein_info.PredictionProteinInfo] = predictions
        self.prediction_configuration = prediction_config
        self.project = current_project

    @staticmethod
    def create_tmp_directories() -> None:
        """This function creates tmp directories in a scratch folder to organize prediction inputs and outputs."""
        # TODO: is there a more elegant way to do it?
        if not os.path.exists(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions")):
            os.mkdir(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
        if not os.path.exists(constants.PREDICTION_FASTA_DIR):
            os.mkdir(constants.PREDICTION_FASTA_DIR)
        if not os.path.exists(constants.PREDICTION_PDB_DIR):
            os.mkdir(constants.PREDICTION_PDB_DIR)

    def create_fasta_files_for_prediction(self) -> None:
        """This function creates fasta file based on the predictions in the object.

        Raises:
            FastaFilesCouldNotBeCreatedError: If a fasta file could not be created.
            FastaFilesNotFoundError: If a fasta file could not be found.
        """
        try:
            protein_objects: list[protein.Protein] = data_transformer.transform_protein_name_seq_tuple_to_sequence_obj(
                self.predictions,
            )
        except exception.IllegalArgumentError:
            logger.error("Invalid argument")
            raise exception.FastaFilesNotCreatedError("")
        except Exception as e:
            logger.error("Unexpected error:", e)
            raise exception.FastaFilesNotCreatedError("")

        logger.debug(
            f"Variable: protein_sequences; Value: {protein_objects} in function create_fasta_files_for_prediction",
        )

        for tmp_protein_to_predict in protein_objects:
            try:
                tmp_protein_to_predict.write_fasta_file(constants.PREDICTION_FASTA_DIR)
            except exception.IllegalArgumentError:
                logger.error("Invalid argument")
                raise exception.FastaFilesNotCreatedError("")
            except exception.DirectoryNotFoundError:
                logger.error("Directory does not exists of given filepath.")
                raise exception.FastaFilesNotCreatedError("")
            except Exception as e:
                logger.error("Unexpected error:", e)
                raise exception.FastaFilesNotCreatedError("")

        if len(os.listdir(constants.PREDICTION_FASTA_DIR)) == 0:
            logger.error("No fasta files were created!")
            raise exception.FastaFilesNotFoundError("")

    def run_prediction(self) -> None:
        """This function runs a structure prediction.

        Raises:
            PredictionEndedWithError: If the prediction ended with an error.
        """
        try:
            colabbatch.Colabbatch(self.prediction_configuration).run_prediction()
        except exception.IllegalArgumentError:
            logger.error(
                f"The Prediction configuration: ({self.prediction_configuration.amber_force_field}, "
                f"{self.prediction_configuration.amber_force_field}) is invalid!",
            )
            raise exception.PredictionEndedWithError("")
        except exception.PredictionEndedWithError:
            logger.error("Prediction ended with errors!")
            raise exception.PredictionEndedWithError("")

    def move_best_prediction_models(self) -> None:
        """This function moves the best prediction model(s) to the project directory.

        Raises:
            UnableToFindColabfoldModelError: If the prediction rank 1 model could not be found.
            FileNotFoundError: If the filepath of the rank 1 model could not be found.
        """
        logger.debug(self.predictions)
        best_prediction_models: list[
            tuple[prediction_protein_info.PredictionProteinInfo, str]
        ] = prediction_util.get_relaxed_rank_1_pdb_file(self.predictions)
        logger.debug(f"These are the models created by ColabFold: {best_prediction_models}")
        if len(best_prediction_models) == 0:
            logger.error("The prediction process finished but no rank 1 model could be found.")
            shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
            raise exception.UnableToFindColabfoldModelError("")

        for tmp_prediction in best_prediction_models:
            try:
                src = path_util.FilePath(f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[1]}")
            except FileNotFoundError:
                logger.error(
                    "This path does not exists: %s",
                    path_util.FilePath(
                        f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[1]}",
                    ).get_filepath(),
                )
                raise FileNotFoundError()
            dest = pathlib.Path(f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[0].name}.pdb")
            os.rename(src.get_filepath(), dest)
            logger.debug(tmp_prediction[0].name)
            self.project.add_existing_protein(
                protein.Protein(tmp_prediction[0].name, pdb_filepath=path_util.FilePath(dest)),
            )
            logger.debug(self.project.proteins)
        try:
            shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
        except Exception as e:
            logger.warning("Scratch dir could not be deleted: %s", e)
        else:
            logger.info("Scratch dir for prediction process deleted successfully.")
