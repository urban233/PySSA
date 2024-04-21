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
"""Module for the job classes."""
import copy
import logging
import subprocess

from PyQt5 import QtCore

from pyssa.controller import database_manager
from pyssa.gui.ui.custom_widgets import job_entry
from pyssa.internal.data_structures import project, structure_prediction, structure_analysis
from pyssa.internal.portal import auxiliary_pymol
from pyssa.logging_pyssa import log_handlers
from pyssa.util import enums, exception, exit_codes, analysis_util

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Job(QtCore.QObject):
    """Base class for all jobs."""
    frozen_project: "project.Project"
    type: "enums.JobType"


class PredictionJob(Job):
    """Job for any type of structure prediction."""

    update_job_entry_signal = QtCore.pyqtSignal(tuple)  # (job entry widget, progress description, progress value)
    cancel_job_signal = QtCore.pyqtSignal(tuple)

    def __init__(self,
                 a_project: "project.Project",
                 the_prediction_protein_infos,
                 the_prediction_configuration,
                 the_project_lock: QtCore.QMutex):
        """

        Notes:
            The project needs to be a pointer!
        """
        super().__init__()
        self.frozen_project = copy.deepcopy(a_project)
        self.type = enums.JobType.PREDICTION
        self.prediction_protein_infos = the_prediction_protein_infos
        self.prediction_configuration = the_prediction_configuration
        self.project_lock = the_project_lock
        self.job_entry_widget: "job_entry.JobEntryWidget" = None

    def run_job(self):
        """Defines how the prediction will be run."""
        structure_prediction_obj = structure_prediction.StructurePrediction(
            self.prediction_protein_infos,
            self.prediction_configuration,
            self.frozen_project,
        )
        self.update_job_entry_signal.emit((self.job_entry_widget, "Creating temp directories ...", 5))
        structure_prediction_obj.create_tmp_directories()
        logger.info("Tmp directories were created.")

        # <editor-fold desc="Creates fasta files for prediction">
        self.update_job_entry_signal.emit((self.job_entry_widget, "Creating FASTA files ...", 10))
        try:
            structure_prediction_obj.create_fasta_files_for_prediction()
        except exception.FastaFilesNotCreatedError:
            tmp_msg = "Fasta files were not created."
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        except exception.FastaFilesNotFoundError:
            tmp_msg = "Fasta files were not found."
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        except Exception as e:
            tmp_msg = f"Unexpected error: {e}"
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        else:
            logger.info("Fasta files were successfully created.")

        # </editor-fold>

        # <editor-fold desc="Runs structure prediction">
        self.update_job_entry_signal.emit((self.job_entry_widget, "Running ColabFold prediction ...", 25))  # 25 must only be used for this exact step
        try:
            structure_prediction_obj.run_prediction()
        except exception.PredictionEndedWithError:
            tmp_msg = "Prediction ended with error."
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        else:
            logger.info("Prediction process finished.")

        # </editor-fold>

        # <editor-fold desc="Saves predicted protein to project">
        self.update_job_entry_signal.emit((self.job_entry_widget, "Saving best prediction results ...", 85))
        try:
            tmp_best_prediction_models = structure_prediction_obj.move_best_prediction_models()
        except exception.UnableToFindColabfoldModelError:
            tmp_msg = "Could not move rank 1 model, because it does not exists."
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        except FileNotFoundError:
            tmp_msg = "Could not move rank 1 model, because it does not exists."
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            tmp_msg = "Could not move rank 1 model, because it does not exists."
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        else:
            structure_prediction_obj.add_proteins_to_project(tmp_best_prediction_models,
                                                             self.frozen_project,
                                                             self.project_lock)
            subprocess.run(["wsl", "--shutdown"])
            logger.info("WSL gets shutdown.")
            if self.job_entry_widget.job_base_information.job_type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS:
                self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.RUNNING
            else:
                self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FINISHED
            self.update_job_entry_signal.emit((self.job_entry_widget, "A structure prediction job finished.", 100))
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

    def cancel_job(self):
        self.cancel_job_signal.emit((self.type, self.job_entry_widget, self, self.prediction_protein_infos))


class DistanceAnalysisJob(Job):
    """Job for a distance analysis."""

    update_job_entry_signal = QtCore.pyqtSignal(tuple)  # (job entry widget, progress description, progress value)

    def __init__(self,
                 a_project,
                 the_project_lock: QtCore.QMutex,
                 a_list_with_analysis_names,
                 a_cutoff: float,
                 cycles: int
                 ):
        super().__init__()
        self.frozen_project = copy.deepcopy(a_project)
        self.type = enums.JobType.DISTANCE_ANALYSIS
        self.list_with_analysis_names = a_list_with_analysis_names
        self.cutoff = a_cutoff
        self.cycles = cycles
        self.project_lock = the_project_lock
        self.job_entry_widget: "job_entry.JobEntryWidget" = None

    def run_job(self):
        try:
            self.update_job_entry_signal.emit((self.job_entry_widget, "Transforming input ...", 15))
            analysis_runs = structure_analysis.Analysis(self.frozen_project)
            analysis_runs.analysis_list = analysis_util.transform_gui_input_to_practical_data(
                self.list_with_analysis_names,
                self.frozen_project,
                self.cutoff,
                self.cycles,
            )
            logger.debug(f"Analysis runs before actual analysis: {analysis_runs.analysis_list}")
            self.update_job_entry_signal.emit((self.job_entry_widget, "Running the distance analysis ...", 50))
            analysis_runs.run_analysis("distance", False)
            logger.debug(f"Analysis runs after actual analysis: {analysis_runs.analysis_list}")
            self.update_job_entry_signal.emit((self.job_entry_widget, "Saving results ...", 80))
            for tmp_protein_pair in analysis_runs.analysis_list:
                tmp_protein_pair.db_project_id = self.frozen_project.get_id()
                copy_tmp_protein_pair = copy.deepcopy(tmp_protein_pair)
                with database_manager.DatabaseManager(str(self.frozen_project.get_database_filepath())) as db_manager:
                    db_manager.open_project_database()
                    copy_tmp_protein_pair.set_id(db_manager.insert_new_protein_pair(copy_tmp_protein_pair))
                    db_manager.close_project_database()
                # Protein pair gets added to "self.frozen_project" of this class
                self.frozen_project.add_protein_pair(copy_tmp_protein_pair)
        except exception.UnableToSetupAnalysisError:
            tmp_msg = "Setting up the analysis runs failed therefore the distance analysis failed."
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        except Exception as e:
            tmp_msg = f"Unknown error: {e}"
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        else:
            if self.job_entry_widget.job_base_information.job_type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS:
                self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.RUNNING
            else:
                self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FINISHED
            self.update_job_entry_signal.emit((self.job_entry_widget, "A distance analysis job finished.", 100))


class PredictionAndDistanceAnalysisJob(Job):
    """Job for a prediction and distance analysis."""

    update_job_entry_signal = QtCore.pyqtSignal(tuple)  # (job entry widget, progress description, progress value)

    def __init__(self,
                 a_prediction_job,
                 a_distance_analysis_job):
        super().__init__()
        self.prediction_job: "PredictionJob" = a_prediction_job
        self.distance_analysis_job: "DistanceAnalysisJob" = a_distance_analysis_job
        self.type = enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
        self.job_entry_widget: "job_entry.JobEntryWidget" = None

    def run_job(self):
        self.update_job_entry_signal.emit((self.job_entry_widget, "Running ColabFold prediction ...", 33))
        self.prediction_job.run_job()
        self.distance_analysis_job.frozen_project = self.prediction_job.frozen_project
        self.update_job_entry_signal.emit((self.job_entry_widget, "Running a distance analysis ...", 66))
        self.distance_analysis_job.run_job()
        self.update_job_entry_signal.emit((self.job_entry_widget, "A ColabFold prediction and distance analysis job finished.", 100))


class RayTracingJob(Job):
    """Job for any type of ray-tracing."""

    update_job_entry_signal = QtCore.pyqtSignal(tuple)  # (job entry widget, progress description, progress value)

    def __init__(self,
                 the_destination_image_filepath,
                 the_cached_session_filepath,
                 image_ray_trace_mode,
                 image_ray_texture,
                 image_renderer):
        super().__init__()
        self.type = enums.JobType.RAY_TRACING
        self.dest_image_filepath = the_destination_image_filepath
        self.cached_session_filepath = the_cached_session_filepath
        self.image_ray_trace_mode = image_ray_trace_mode
        self.image_ray_texture = image_ray_texture
        self.image_renderer = image_renderer
        self.job_entry_widget: "job_entry.JobEntryWidget" = None

    def run_job(self):
        self.update_job_entry_signal.emit((self.job_entry_widget, "Starting rendering process ...", 33))
        try:
            auxiliary_pymol.AuxiliaryPyMOL.create_ray_traced_image(
                self.dest_image_filepath,
                self.cached_session_filepath,
                self.image_ray_trace_mode,
                self.image_ray_texture,
                self.image_renderer
            )
        except Exception as e:
            tmp_msg = f"Unknown error: {e}"
            logger.error(tmp_msg)
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FAILED
            self.update_job_entry_signal.emit((self.job_entry_widget, tmp_msg, 100))
        else:
            self.job_entry_widget.job_base_information.job_progress = enums.JobProgress.FINISHED
            self.update_job_entry_signal.emit((self.job_entry_widget, "Create ray-tracing image job finished.", 100))
