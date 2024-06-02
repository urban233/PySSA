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
"""Module for structure prediction class."""
import os
import pathlib
import shutil
import logging
from typing import Optional

import zmq
from PyQt5 import QtCore

from src.auxiliary_pymol import auxiliary_pymol_client
from src.pyssa.controller import database_manager
from src.pyssa.internal.data_structures.data_classes import prediction_protein_info
from src.pyssa.internal.data_structures import protein, job
from src.pyssa.internal.data_structures import project
from src.pyssa.internal.data_structures.data_classes import prediction_configuration
from src.pyssa.internal.data_processing import data_transformer
from src.pyssa.internal.prediction_engines import colabbatch
from src.pyssa.util import constants
from src.pyssa.util import prediction_util
from src.pyssa.util import exception
from src.pyssa.io_pyssa import path_util, bio_data
from src.pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class StructurePrediction:
  """This class is used to organize the structure prediction process."""

  # <editor-fold desc="Class attributes">
  predictions: list["prediction_protein_info.PredictionProteinInfo"]
  """A list of tuples of protein name and sequence."""

  prediction_config: "prediction_configuration.PredictionConfiguration"
  """The configuration settings for the prediction."""

  current_project: "project.Project"
  """The current project in use."""

  # </editor-fold>

  def __init__(
      self,
      predictions: list["prediction_protein_info.PredictionProteinInfo"],
      prediction_config: "prediction_configuration.PredictionConfiguration",
  ) -> None:
    """Constructor.

    Args:
        predictions (prediction_protein_info.PredictionProteinInfo): A list of prediction protein infos.
        prediction_config (prediction_configuration.PredictionConfiguration): The configuration settings for the prediction.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if predictions is None:
      logger.error("predictions is None.")
      raise exception.IllegalArgumentError("predictions is None.")
    if prediction_config is None:
      logger.error("prediction_config is None.")
      raise exception.IllegalArgumentError("prediction_config is None.")

    # </editor-fold>

    self.predictions: list[prediction_protein_info.PredictionProteinInfo] = (
        predictions
    )
    self.prediction_configuration = prediction_config

  @staticmethod
  def create_tmp_directories() -> None:
    """Creates tmp directories in a scratch folder to organize prediction inputs and outputs."""
    if os.path.exists(
        pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions")
    ):
      shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
      os.mkdir(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
    else:
      os.mkdir(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
    if os.path.exists(constants.PREDICTION_FASTA_DIR):
      shutil.rmtree(constants.PREDICTION_FASTA_DIR)
      os.mkdir(constants.PREDICTION_FASTA_DIR)
    else:
      os.mkdir(constants.PREDICTION_FASTA_DIR)
    if os.path.exists(constants.PREDICTION_PDB_DIR):
      shutil.rmtree(constants.PREDICTION_PDB_DIR)
      os.mkdir(constants.PREDICTION_PDB_DIR)
    else:
      os.mkdir(constants.PREDICTION_PDB_DIR)

  def create_fasta_files_for_prediction(self) -> None:
    """Creates fasta file based on the predictions in the object.

    Raises:
        FastaFilesCouldNotBeCreatedError: If a fasta file could not be created.
        FastaFilesNotFoundError: If a fasta file could not be found.
    """
    try:
      protein_objects: list[protein.Protein] = (
          data_transformer.transform_protein_name_seq_tuple_to_sequence_obj(
              self.predictions,
          )
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
    """Runs a structure prediction.

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

  def move_best_prediction_models(
      self,
  ) -> Optional[
      list[tuple[prediction_protein_info.PredictionProteinInfo, str]]
  ]:
    """Moves the best prediction model(s) to the project directory.

    Raises:
        UnableToFindColabfoldModelError: If the prediction rank 1 model could not be found.
        FileNotFoundError: If the filepath of the rank 1 model could not be found.
    """
    logger.debug(self.predictions)
    best_prediction_models: list[
        tuple[prediction_protein_info.PredictionProteinInfo, str]
    ] = prediction_util.get_relaxed_rank_1_pdb_file(self.predictions)
    logger.debug(
        f"These are the models created by ColabFold: {best_prediction_models}"
    )
    if len(best_prediction_models) == 0:
      logger.error(
          "The prediction process finished but no rank 1 model could be found."
      )
      shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
      raise exception.UnableToFindColabfoldModelError("")

    for tmp_prediction in best_prediction_models:
      try:
        src = path_util.FilePath(
            pathlib.Path(
                f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[1]}"
            )
        )
      except FileNotFoundError:
        logger.error(
            "This path does not exists: %s",
            path_util.FilePath(
                pathlib.Path(
                    f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[1]}"
                ),
            ).get_filepath(),
        )
        raise FileNotFoundError()
      dest = pathlib.Path(
          f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[0].name}.pdb"
      )
      os.rename(src.get_filepath(), dest)
      logger.debug(tmp_prediction[0].name)
      return best_prediction_models
    return None

  def add_proteins_to_project(
      self,
      the_main_socket: zmq.Socket,
      a_socket: zmq.Socket,
      the_general_purpose_socket: zmq.Socket,
      best_prediction_models: list,
      a_project: "project.Project",
      the_project_lock: QtCore.QMutex,
  ) -> None:
    """Add proteins to a project.

    Args:
        the_main_socket (zmq.Socket): The main socket.
        a_socket (zmq.Socket): A socket.
        the_general_purpose_socket (zmq.Socket): The general purpose socket.
        best_prediction_models (list): The best prediction models.
        a_project (project.Project): A project.
        the_project_lock (QtCore.QMutex): The project lock.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_main_socket is None:
      logger.error("the_main_socket is None.")
      raise exception.IllegalArgumentError("the_main_socket is None.")
    if a_socket is None:
      logger.error("a_socket is None.")
      raise exception.IllegalArgumentError("a_socket is None.")
    if the_general_purpose_socket is None:
      logger.error("the_general_purpose_socket is None.")
      raise exception.IllegalArgumentError(
          "the_general_purpose_socket is None."
      )
    if best_prediction_models is None:
      logger.error("best_prediction_models is None.")
      raise exception.IllegalArgumentError("best_prediction_models is None.")
    if a_project is None:
      logger.error("a_project is None.")
      raise exception.IllegalArgumentError("a_project is None.")
    if the_project_lock is None:
      logger.error("the_project_lock is None.")
      raise exception.IllegalArgumentError("the_project_lock is None.")

    # </editor-fold>

    the_project_lock.lock()
    for tmp_prediction in best_prediction_models:
      tmp_protein = protein.Protein(tmp_prediction[0].name)
      tmp_protein.add_protein_structure_data_from_local_pdb_file(
          pathlib.Path(
              f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction[0].name}.pdb"
          ),
          the_main_socket,
          the_general_purpose_socket,
      )
      pdb_filepath = pathlib.Path(
          f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein.get_molecule_object()}.pdb"
      )
      try:
        bio_data.build_pdb_file(tmp_protein.get_pdb_data(), pdb_filepath)
      except exception.IllegalArgumentError:
        logger.error(
            f"The argument pdb data is not usable: {tmp_protein.get_pdb_data}."
        )
        raise exception.UnableToCreatePdbFileError("")
      except exception.DirectoryNotFoundError:
        logger.error(f"The argument pdb_filepath is illegal: {pdb_filepath}!")
        raise exception.UnableToCreatePdbFileError("")
      except PermissionError:
        logger.error(f"The argument pdb_filepath is illegal: {pdb_filepath}!")
        raise exception.UnableToCreatePdbFileError("")
      except exception.UnableToOpenFileError:
        logger.error("pdb file could not be opened for writing.")
        raise exception.UnableToOpenFileError("")

      tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
          the_main_socket,
          a_socket,
          job.PredictionJobDescription(pdb_filepath),
      )
      tmp_protein.pymol_session = tmp_reply["data"][0]

      with database_manager.DatabaseManager(
          str(a_project.get_database_filepath())
      ) as db_manager:
        logger.info(
            f"Inserting {tmp_protein.get_molecule_object()} into current project, from prediction thread."
        )
        tmp_protein.db_project_id = a_project.get_id()
        tmp_protein.add_id_to_all_chains(
            db_manager.get_next_id_of_chain_table()
        )
        tmp_protein.set_id(db_manager.insert_new_protein(tmp_protein))
      a_project.add_existing_protein(
          tmp_protein,
      )
    the_project_lock.unlock()
    try:
      shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
    except Exception as e:
      logger.warning("Scratch dir could not be deleted: %s", e)
    else:
      logger.info("Scratch dir for prediction process deleted successfully.")
