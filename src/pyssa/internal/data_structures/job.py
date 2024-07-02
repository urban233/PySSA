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
"""Module for the job classes."""
import copy
import logging
import subprocess
from typing import Optional

import zmq
from PyQt5 import QtCore

from src.auxiliary_pymol import auxiliary_pymol_client
from src.pyssa.controller import database_manager
from src.pyssa.gui.ui.custom_widgets import job_entry
from src.pyssa.internal.data_structures import project, structure_prediction, structure_analysis, protein_pair
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import enums, exception, analysis_util

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Job(QtCore.QObject):
  """Base class for all jobs."""

  # <editor-fold desc="Class attributes">
  frozen_project: "project.Project"
  """A copied version of the project before the job."""

  type: "enums.JobType"
  """The type of job."""

  # </editor-fold>


class PredictionJob(Job):
  """Job for any type of structure prediction."""

  # <editor-fold desc="Class attributes">
  update_job_entry_signal = QtCore.pyqtSignal(tuple)
  """A signal to update the job entry with the form (job entry widget, progress description, progress value)."""

  cancel_job_signal = QtCore.pyqtSignal(tuple)
  """A signal to cancel the job."""

  # </editor-fold>

  def __init__(
      self,
      the_main_socket: zmq.Socket,
      a_socket: zmq.Socket,
      the_general_purpose_socket: zmq.Socket,
      a_project: "project.Project",
      the_prediction_protein_infos,
      the_prediction_configuration,
      the_project_lock: QtCore.QMutex,
  ) -> None:
    """Constructor.

    Args:
        the_main_socket (zmq.Socket): The main socket for communication.
        a_socket (zmq.Socket): A socket for communication.
        the_general_purpose_socket (zmq.Socket): The general purpose socket for communication.
        a_project (project.Project): The project object.
        the_prediction_protein_infos (fixme: Needs to be checked): The prediction protein infos.
        the_prediction_configuration (fixme: Needs to be checked): The prediction configuration.
        the_project_lock (QtCore.QMutex): The project lock object.

    Notes:
        The project needs to be a pointer!

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
    if a_project is None:
      logger.error("a_project is None.")
      raise exception.IllegalArgumentError("a_project is None.")
    if the_prediction_protein_infos is None:
      logger.error("the_prediction_protein_infos is None.")
      raise exception.IllegalArgumentError(
          "the_prediction_protein_infos is None."
      )
    if the_prediction_configuration is None:
      logger.error("the_prediction_configuration is None.")
      raise exception.IllegalArgumentError(
          "the_prediction_configuration is None."
      )
    if the_project_lock is None:
      logger.error("the_project_lock is None.")
      raise exception.IllegalArgumentError("the_project_lock is None.")

    # </editor-fold>

    super().__init__()
    self.frozen_project = copy.deepcopy(a_project)
    self.type = enums.JobType.PREDICTION
    self._main_socket = the_main_socket
    self._socket = a_socket
    self._general_purpose_socket = the_general_purpose_socket
    self.prediction_protein_infos = the_prediction_protein_infos
    self.prediction_configuration = the_prediction_configuration
    self.project_lock = the_project_lock
    self.job_entry_widget: Optional["job_entry.JobEntryWidget"] = None

  def run_job(self) -> None:
    """Defines how the prediction will be run."""
    structure_prediction_obj = structure_prediction.StructurePrediction(
        self.prediction_protein_infos,
        self.prediction_configuration,
    )
    self.update_job_entry_signal.emit(
        (self.job_entry_widget, "Creating temp directories ...", 5)
    )
    structure_prediction_obj.create_tmp_directories()
    logger.info("Tmp directories were created.")

    # <editor-fold desc="Creates fasta files for prediction">
    self.update_job_entry_signal.emit(
        (self.job_entry_widget, "Creating FASTA files ...", 10)
    )
    try:
      structure_prediction_obj.create_fasta_files_for_prediction()
    except exception.FastaFilesNotCreatedError:
      tmp_msg = "Fasta files were not created."
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    except exception.FastaFilesNotFoundError:
      tmp_msg = "Fasta files were not found."
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    except Exception as e:
      tmp_msg = f"Unexpected error: {e}"
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    else:
      logger.info("Fasta files were successfully created.")

    # </editor-fold>

    # <editor-fold desc="Runs structure prediction">
    self.update_job_entry_signal.emit(
        (self.job_entry_widget, "Running ColabFold prediction ...", 25)
    )  # 25 must only be used for this exact step
    try:
      structure_prediction_obj.run_prediction()
    except exception.PredictionEndedWithError:
      tmp_msg = "Prediction ended with error."
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    else:
      logger.info("Prediction process finished.")

    # </editor-fold>

    # <editor-fold desc="Saves predicted protein to project">
    self.update_job_entry_signal.emit(
        (self.job_entry_widget, "Saving best prediction results ...", 85)
    )
    try:
      tmp_best_prediction_models = (
          structure_prediction_obj.move_best_prediction_models()
      )
    except exception.UnableToFindColabfoldModelError:
      tmp_msg = "Could not move rank 1 model, because it does not exists."
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    except FileNotFoundError:
      tmp_msg = "Could not move rank 1 model, because it does not exists."
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    except Exception as e:
      logger.error(f"Unexpected error: {e}")
      tmp_msg = "Could not move rank 1 model, because it does not exists."
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    else:
      structure_prediction_obj.add_proteins_to_project(
          self._main_socket,
          self._socket,
          self._general_purpose_socket,
          tmp_best_prediction_models,
          self.frozen_project,
          self.project_lock,
      )
      subprocess.run(["wsl", "--shutdown"])
      logger.info("WSL gets shutdown.")
      if (
          self.job_entry_widget.job_base_information.job_type
          == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
      ):
        self.job_entry_widget.job_base_information.job_progress = (
            enums.JobProgress.RUNNING
        )
      else:
        self.job_entry_widget.job_base_information.job_progress = (
            enums.JobProgress.FINISHED
        )
      self.update_job_entry_signal.emit(
          (self.job_entry_widget, "A structure prediction job finished.", 100)
      )
    # </editor-fold>

    # self._interface_manager.status_bar_manager.hide_progress_bar()
    # tmp_exit_code = result[0]
    # tmp_exit_code_description = result[1]
    # if tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
    #
    #     # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #     #     "Prediction failed because there was an error writing the fasta file(s)!",
    #     #     "Structure Prediction",
    #     #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #     # )
    #     # tmp_dialog.exec_()
    #     constants.PYSSA_LOGGER.error(
    #         f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #     )
    #     self._interface_manager.status_bar_manager.show_error_message(
    #         "Prediction failed because there was an error writing the fasta file(s)!")
    # elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
    #     # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #     #     "Prediction failed because the fasta file(s) could not be found!",
    #     #     "Structure Prediction",
    #     #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #     # )
    #     # tmp_dialog.exec_()
    #     constants.PYSSA_LOGGER.error(
    #         f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #     )
    #     self._interface_manager.status_bar_manager.show_error_message(
    #         "Prediction failed because the fasta file(s) could not be found!")
    # elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
    #     # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #     #     "Prediction failed because a subprocess failed!",
    #     #     "Structure Prediction",
    #     #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #     # )
    #     # tmp_dialog.exec_()
    #     constants.PYSSA_LOGGER.error(
    #         f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #     )
    #     self._interface_manager.status_bar_manager.show_error_message(
    #         "Prediction failed because a subprocess failed!")
    # elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
    #     # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #     #     "Prediction failed because of an unknown error!",
    #     #     "Structure Prediction",
    #     #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #     # )
    #     # tmp_dialog.exec_()
    #     constants.PYSSA_LOGGER.error(
    #         f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #     )
    #     self._interface_manager.status_bar_manager.show_error_message(
    #         "Prediction failed because of an unknown error!")
    # elif tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
    #     # Prediction was successful
    #     tmp_proteins_to_add = self._main_view_state.get_not_matching_proteins(
    #         self._interface_manager.get_current_project().proteins
    #     )
    #     for tmp_protein in tmp_proteins_to_add:
    #         self._interface_manager.add_protein_to_proteins_model(tmp_protein)
    #
    #     self._interface_manager.main_tasks_manager.prediction_task = None
    #     self._interface_manager.refresh_main_view()
    #
    #     if self.active_custom_message_box is not None:
    #         self.active_custom_message_box.close()
    #
    #     self._main_view_state.restore_main_view_state()
    #
    #     constants.PYSSA_LOGGER.info("All structure predictions are done.")
    #     self._interface_manager.status_bar_manager.show_temporary_message("All structure predictions are done.")
    #     self._interface_manager.stop_wait_cursor()
    # else:
    #     self._interface_manager.status_bar_manager.show_error_message(
    #         "Prediction failed because of an unknown case!")

  def cancel_job(self) -> None:
    """Emits the `cancel_job_signal` signal to cancel the job."""
    self.cancel_job_signal.emit(
        (self.type, self.job_entry_widget, self, self.prediction_protein_infos)
    )


class DistanceAnalysisJob(Job):
  """Job for a distance analysis."""

  # <editor-fold desc="Class attributes">
  update_job_entry_signal = QtCore.pyqtSignal(tuple)
  """A signal to update the job entry with the form (job entry widget, progress description, progress value)."""

  # </editor-fold>

  def __init__(
      self,
      the_main_socket: zmq.Socket,
      a_socket: zmq.Socket,
      a_project: "project.Project",
      the_project_lock: QtCore.QMutex,
      a_list_with_analysis_names: list,
      a_cutoff: float,
      cycles: int,
  ) -> None:
    """Constructor.

    Args:
        the_main_socket (zmq.Socket): The main socket used for communication.
        a_socket (zmq.Socket): The specific socket used for communication.
        a_project (project.Project): The project to be analyzed.
        the_project_lock (QtCore.QMutex): Lock object for thread synchronization.
        a_list_with_analysis_names (list): List of names of the analysis to be performed.
        a_cutoff (float): The cutoff value for the analysis.
        cycles (int): The number of cycles to be performed in the analysis.

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
    if a_project is None:
      logger.error("a_project is None.")
      raise exception.IllegalArgumentError("a_project is None.")
    if the_project_lock is None:
      logger.error("the_project_lock is None.")
      raise exception.IllegalArgumentError("the_project_lock is None.")
    if a_list_with_analysis_names is None:
      logger.error("a_list_with_analysis_names is None.")
      raise exception.IllegalArgumentError(
          "a_list_with_analysis_names is None."
      )
    if a_cutoff is None:
      logger.error("a_cutoff is None.")
      raise exception.IllegalArgumentError("a_cutoff is None.")
    if cycles is None:
      logger.error("cycles is None.")
      raise exception.IllegalArgumentError("cycles is None.")

    # </editor-fold>

    super().__init__()
    self.type = enums.JobType.DISTANCE_ANALYSIS
    self._main_socket = the_main_socket
    self._socket = a_socket
    self.frozen_project = copy.deepcopy(a_project)
    self.list_with_analysis_names = a_list_with_analysis_names
    self.cutoff = a_cutoff
    self.cycles = cycles
    self.project_lock = the_project_lock
    self.job_entry_widget: Optional["job_entry.JobEntryWidget"] = None

  def run_job(self) -> None:
    """Defines how the distance analysis will be run."""
    try:
      self.update_job_entry_signal.emit(
          (self.job_entry_widget, "Transforming input ...", 15)
      )
      analysis_runs = structure_analysis.Analysis(self.frozen_project)
      analysis_runs.analysis_list = (
          analysis_util.transform_gui_input_to_practical_data(
              self.list_with_analysis_names,
              self.frozen_project,
              self.cutoff,
              self.cycles,
          )
      )
      logger.debug(
          f"Analysis runs before actual analysis: {analysis_runs.analysis_list}"
      )
      self.update_job_entry_signal.emit(
          (self.job_entry_widget, "Running the distance analysis ...", 50)
      )
      analysis_runs.run_analysis(
          "distance", False, self._main_socket, self._socket
      )
      logger.debug(
          f"Analysis runs after actual analysis: {analysis_runs.analysis_list}"
      )
      self.update_job_entry_signal.emit(
          (self.job_entry_widget, "Saving results ...", 80)
      )
      for tmp_protein_pair in analysis_runs.analysis_list:
        tmp_protein_pair.db_project_id = self.frozen_project.get_id()
        copy_tmp_protein_pair: "protein_pair.ProteinPair" = copy.deepcopy(
            tmp_protein_pair
        )
        with database_manager.DatabaseManager(
            str(self.frozen_project.get_database_filepath())
        ) as db_manager:
          copy_tmp_protein_pair.set_id(
              db_manager.insert_new_protein_pair(copy_tmp_protein_pair)
          )
        # Protein pair gets added to "self.frozen_project" of this class
        self.frozen_project.add_protein_pair(copy_tmp_protein_pair)
    except exception.UnableToSetupAnalysisError:
      tmp_msg = "Setting up the analysis runs failed therefore the distance analysis failed."
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    except exception.UnableToDoAnalysisError as e:
      tmp_msg = f"Error: {e}"
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
        enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    except Exception as e:
      tmp_msg = f"Unknown error: {e}"
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    else:
      if (
          self.job_entry_widget.job_base_information.job_type
          == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
      ):
        self.job_entry_widget.job_base_information.job_progress = (
            enums.JobProgress.RUNNING
        )
      else:
        self.job_entry_widget.job_base_information.job_progress = (
            enums.JobProgress.FINISHED
        )
      self.update_job_entry_signal.emit(
          (self.job_entry_widget, "A distance analysis job finished.", 100)
      )


class PredictionAndDistanceAnalysisJob(Job):
  """Job for a prediction and distance analysis."""

  update_job_entry_signal = QtCore.pyqtSignal(tuple)
  """A signal to update the job entry with the form (job entry widget, progress description, progress value)."""

  def __init__(
      self,
      a_prediction_job: "PredictionJob",
      a_distance_analysis_job: "DistanceAnalysisJob",
  ) -> None:
    """Constructor.

    Args:
        a_prediction_job (PredictionJob): An instance of the PredictionJob class.
        a_distance_analysis_job (DistanceAnalysisJob): An instance of the DistanceAnalysisJob class.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_prediction_job is None:
      logger.error("a_prediction_job is None.")
      raise exception.IllegalArgumentError("a_prediction_job is None.")
    if a_distance_analysis_job is None:
      logger.error("a_distance_analysis_job is None.")
      raise exception.IllegalArgumentError("a_distance_analysis_job is None.")

    # </editor-fold>

    super().__init__()
    self.prediction_job: "PredictionJob" = a_prediction_job
    self.distance_analysis_job: "DistanceAnalysisJob" = a_distance_analysis_job
    self.type = enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
    self.job_entry_widget: Optional["job_entry.JobEntryWidget"] = None

  def run_job(self) -> None:
    """Defines how the prediction and distance analysis will be run."""
    self.job_entry_widget.job_base_information.job_progress = (
        enums.JobProgress.RUNNING
    )
    self.update_job_entry_signal.emit(
        (self.job_entry_widget, "Running ColabFold prediction ...", 33)
    )
    self.prediction_job.run_job()
    self.distance_analysis_job.frozen_project = (
        self.prediction_job.frozen_project
    )
    self.update_job_entry_signal.emit(
        (self.job_entry_widget, "Running a distance analysis ...", 66)
    )
    self.distance_analysis_job.run_job()
    self.job_entry_widget.job_base_information.job_progress = (
        enums.JobProgress.FINISHED
    )
    self.update_job_entry_signal.emit(
        (
            self.job_entry_widget,
            "A ColabFold prediction and distance analysis job finished.",
            100,
        )
    )


class RayTracingJob(Job):
  """Job for any type of ray-tracing."""

  update_job_entry_signal = QtCore.pyqtSignal(tuple)
  """A signal to update the job entry with the form (job entry widget, progress description, progress value)."""

  def __init__(
      self,
      the_main_socket: zmq.Socket,
      a_socket: zmq.Socket,
      the_destination_image_filepath: str,
      the_cached_session_filepath: str,
      image_ray_trace_mode: int,
      image_ray_texture: int,
      image_renderer: str,
  ) -> None:
    """Constructor.

    Args:
        the_main_socket (zmq.Socket): The main socket used for communication.
        a_socket (zmq.Socket): A secondary socket used for communication.
        the_destination_image_filepath (str): The filepath to save the rendered image.
        the_cached_session_filepath (str): The filepath to the cached session file.
        image_ray_trace_mode (int): The mode for the ray tracing.
        image_ray_texture (int): The texture for the ray tracing.
        image_renderer (str): The renderer used for rendering.

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
    if the_destination_image_filepath is None:
      logger.error("the_destination_image_filepath is None.")
      raise exception.IllegalArgumentError(
          "the_destination_image_filepath is None."
      )
    if the_cached_session_filepath is None:
      logger.error("the_cached_session_filepath is None.")
      raise exception.IllegalArgumentError(
          "the_cached_session_filepath is None."
      )
    if image_ray_trace_mode is None:
      logger.error("image_ray_trace_mode is None.")
      raise exception.IllegalArgumentError("image_ray_trace_mode is None.")
    if image_ray_texture is None:
      logger.error("image_ray_texture is None.")
      raise exception.IllegalArgumentError("image_ray_texture is None.")
    if image_renderer is None:
      logger.error("image_renderer is None.")
      raise exception.IllegalArgumentError("image_renderer is None.")

    # </editor-fold>

    super().__init__()
    self.type = enums.JobType.RAY_TRACING
    self._main_socket = the_main_socket
    self._socket = a_socket
    self.dest_image_filepath = the_destination_image_filepath
    self.cached_session_filepath = the_cached_session_filepath
    self.image_ray_trace_mode = image_ray_trace_mode
    self.image_ray_texture = image_ray_texture
    self.image_renderer = image_renderer
    self.job_entry_widget: "job_entry.JobEntryWidget" = None

  def run_job(self) -> None:
    """Defines how the ray-tracing will be run."""
    self.update_job_entry_signal.emit(
        (self.job_entry_widget, "Starting rendering process ...", 33)
    )
    try:
      auxiliary_pymol_client.send_request_to_auxiliary_pymol(
          self._main_socket,
          self._socket,
          RayTracingJobDescription(
              self.dest_image_filepath,
              self.cached_session_filepath,
              self.image_ray_trace_mode,
              self.image_ray_texture,
              self.image_renderer,
          ),
      )
    except Exception as e:
      tmp_msg = f"Unknown error: {e}"
      logger.error(tmp_msg)
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FAILED
      )
      self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
    else:
      self.job_entry_widget.job_base_information.job_progress = (
          enums.JobProgress.FINISHED
      )
      self.update_job_entry_signal.emit(
          (self.job_entry_widget, "Create ray-tracing image job finished.", 100)
      )


class PredictionJobDescription:
  """Job description for a prediction job."""

  def __init__(self, a_pdb_filepath: str) -> None:
    """Constructor.

    Args:
        a_pdb_filepath (str): The filepath of the PDB file.

    Raises:
        exception.IllegalArgumentError: If a_pdb_filepath is None.
    """
    # <editor-fold desc="Checks">
    if a_pdb_filepath is None:
      logger.error("a_pdb_filepath is None.")
      raise exception.IllegalArgumentError("a_pdb_filepath is None.")

    # </editor-fold>

    self.type: "enums.JobType" = enums.JobType.PREDICTION
    self.pdb_filepath: str = str(a_pdb_filepath)

  def get_dict(self) -> dict:
    """Gets a dictionary representation of the object.

    Returns:
        A dictionary containing the job type and PDB filepath.
    """
    return {
        enums.JobDescriptionKeys.JOB_TYPE.value: self.type.value,
        enums.JobDescriptionKeys.PDB_FILEPATH.value: self.pdb_filepath,
    }


class DistanceAnalysisJobDescription:
  """Job description for a distance analysis job."""

  def __init__(
      self,
      the_protein_pair_name: str,
      a_protein_1_pdb_cache_filepath: str,
      a_protein_2_pdb_cache_filepath: str,
      a_protein_1_pymol_selection_string: str,
      a_protein_2_pymol_selection_string: str,
      a_cutoff: float,
      the_cycles: int,
  ) -> None:
    """Constructor.

    Args:
        the_protein_pair_name (str): The name of the protein pair.
        a_protein_1_pdb_cache_filepath (str): The file path of the PDB cache for protein 1.
        a_protein_2_pdb_cache_filepath (str): The file path of the PDB cache for protein 2.
        a_protein_1_pymol_selection_string (str): The PyMOL selection string for protein 1.
        a_protein_2_pymol_selection_string (str): The PyMOL selection string for protein 2.
        a_cutoff (float): The cutoff value.
        the_cycles (int): The number of cycles.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_protein_pair_name is None:
      logger.error("the_protein_pair_name is None.")
      raise exception.IllegalArgumentError("the_protein_pair_name is None.")
    if a_protein_1_pdb_cache_filepath is None:
      logger.error("a_protein_1_pdb_cache_filepath is None.")
      raise exception.IllegalArgumentError(
          "a_protein_1_pdb_cache_filepath is None."
      )
    if a_protein_2_pdb_cache_filepath is None:
      logger.error("a_protein_2_pdb_cache_filepath is None.")
      raise exception.IllegalArgumentError(
          "a_protein_2_pdb_cache_filepath is None."
      )
    if a_protein_1_pymol_selection_string is None:
      logger.error("a_protein_1_pymol_selection_string is None.")
      raise exception.IllegalArgumentError(
          "a_protein_1_pymol_selection_string is None."
      )
    if a_protein_2_pymol_selection_string is None:
      logger.error("a_protein_2_pymol_selection_string is None.")
      raise exception.IllegalArgumentError(
          "a_protein_2_pymol_selection_string is None."
      )
    if a_cutoff is None:
      logger.error("a_cutoff is None.")
      raise exception.IllegalArgumentError("a_cutoff is None.")
    if the_cycles is None:
      logger.error("the_cycles is None.")
      raise exception.IllegalArgumentError("the_cycles is None.")

    # </editor-fold>

    self.type: "enums.JobType" = enums.JobType.DISTANCE_ANALYSIS
    self.the_protein_pair_name = the_protein_pair_name
    self.a_protein_1_pdb_cache_filepath = str(a_protein_1_pdb_cache_filepath)
    self.a_protein_2_pdb_cache_filepath = str(a_protein_2_pdb_cache_filepath)
    self.a_protein_1_pymol_selection_string = a_protein_1_pymol_selection_string
    self.a_protein_2_pymol_selection_string = a_protein_2_pymol_selection_string
    self.cutoff = a_cutoff
    self.cycles = the_cycles

  def get_dict(self) -> dict:
    """Gets a dictionary representation of the object.

    Returns:
        A dictionary containing the instance attributes.
    """
    return {
        enums.JobDescriptionKeys.JOB_TYPE.value: self.type.value,
        enums.JobDescriptionKeys.PROTEIN_PAIR_NAME.value: self.the_protein_pair_name,
        enums.JobDescriptionKeys.PROTEIN_1_PDB_CACHE_FILEPATH.value: self.a_protein_1_pdb_cache_filepath,
        enums.JobDescriptionKeys.PROTEIN_2_PDB_CACHE_FILEPATH.value: self.a_protein_2_pdb_cache_filepath,
        enums.JobDescriptionKeys.PROTEIN_1_PYMOL_SELECTION_STRING.value: self.a_protein_1_pymol_selection_string,
        enums.JobDescriptionKeys.PROTEIN_2_PYMOL_SELECTION_STRING.value: self.a_protein_2_pymol_selection_string,
        enums.JobDescriptionKeys.CUTOFF.value: self.cutoff,
        enums.JobDescriptionKeys.CYCLES.value: self.cycles,
    }


class RayTracingJobDescription:
  """Job description for a ray-tracing job."""

  def __init__(
      self,
      a_destination_image_filepath: str,
      a_cached_session_filepath: str,
      an_image_ray_trace_mode: int,
      an_image_ray_texture: int,
      an_image_renderer: str,
  ) -> None:
    """Constructor.

    Args:
        a_destination_image_filepath (str): The filepath of the destination image.
        a_cached_session_filepath (str): The filepath of the cached session.
        an_image_ray_trace_mode (int): The ray trace mode for the image.
        an_image_ray_texture (int): The ray texture for the image.
        an_image_renderer (str): The renderer for the image.
    """
    # <editor-fold desc="Checks">
    if a_destination_image_filepath is None:
      logger.error("a_destination_image_filepath is None.")
      raise exception.IllegalArgumentError(
          "a_destination_image_filepath is None."
      )
    if a_cached_session_filepath is None:
      logger.error("a_cached_session_filepath is None.")
      raise exception.IllegalArgumentError("a_cached_session_filepath is None.")
    if an_image_ray_trace_mode is None:
      logger.error("an_image_ray_trace_mode is None.")
      raise exception.IllegalArgumentError("an_image_ray_trace_mode is None.")
    if an_image_ray_texture is None:
      logger.error("an_image_ray_texture is None.")
      raise exception.IllegalArgumentError("an_image_ray_texture is None.")
    if an_image_renderer is None:
      logger.error("an_image_renderer is None.")
      raise exception.IllegalArgumentError("an_image_renderer is None.")

    # </editor-fold>

    self.type: "enums.JobType" = enums.JobType.RAY_TRACING
    self.dest_image_filepath: str = str(a_destination_image_filepath)
    self.cached_session_filepath: str = str(a_cached_session_filepath)
    self.image_ray_trace_mode = an_image_ray_trace_mode
    self.image_ray_texture = an_image_ray_texture
    self.image_renderer = an_image_renderer

  def get_dict(self) -> dict:
    """Gets a dictionary representation of the object.

    Returns:
        A dictionary containing the instance attributes.
    """
    return {
        enums.JobDescriptionKeys.JOB_TYPE.value: self.type.value,
        enums.JobDescriptionKeys.IMAGE_DESTINATION_FILEPATH.value: self.dest_image_filepath,
        enums.JobDescriptionKeys.CACHED_SESSION_FILEPATH.value: self.cached_session_filepath,
        enums.JobDescriptionKeys.RAY_TRACE_MODE.value: self.image_ray_trace_mode,
        enums.JobDescriptionKeys.RAY_TEXTURE.value: self.image_ray_texture,
        enums.JobDescriptionKeys.RAY_TRACING_RENDERER.value: self.image_renderer,
    }


class GeneralPurposeJobDescription:
  """Job description for a general purpose job."""

  def __init__(
      self, a_job_short_description: "enums.JobShortDescription"
  ) -> None:
    """Constructor.

    Args:
        a_job_short_description (enums.JobShortDescription): The short description of the job.

    Raises:
        exception.IllegalArgumentError: If a_job_short_description is None.
    """
    # <editor-fold desc="Checks">
    if a_job_short_description is None:
      logger.error("a_job_short_description is None.")
      raise exception.IllegalArgumentError("a_job_short_description is None.")

    # </editor-fold>

    self.type: "enums.JobType" = enums.JobType.GENERAL_PURPOSE
    self.job_short_description: "enums.JobShortDescription" = (
        a_job_short_description
    )
    self.job_information = {
        enums.JobDescriptionKeys.JOB_TYPE.value: self.type.value,
        enums.JobDescriptionKeys.JOB_SHORT_DESCRIPTION.value: self.job_short_description.value,
    }

  def setup_dict(self, an_argument_dict: dict) -> None:
    """Update the job information dictionary with the given argument dictionary.

    Args:
        an_argument_dict (dict): A dictionary containing the information to be added or updated in the job information dictionary.

    Raises:
        exception.IllegalArgumentError: If an_argument_dict is None.
    """
    # <editor-fold desc="Checks">
    if an_argument_dict is None:
      logger.error("an_argument_dict is None.")
      raise exception.IllegalArgumentError("an_argument_dict is None.")

    # </editor-fold>

    self.job_information.update(an_argument_dict)

  def get_dict(self) -> dict:
    """Gets a dictionary representation of the object.

    Returns:
        A dictionary containing the instance attributes.
    """
    return self.job_information
