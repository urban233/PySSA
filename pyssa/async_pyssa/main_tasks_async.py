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
"""Module for all asynchronous functions used in the main tasks."""
import copy
import logging
import subprocess
import time

from pyssa.controller import interface_manager, database_manager
from pyssa.internal.data_structures import structure_prediction, structure_analysis
from pyssa.internal.thread.async_pyssa import custom_signals, locks
from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception, exit_codes, analysis_util

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def predict_protein_with_colabfold(
    the_prediction_protein_infos: list["prediction_protein_info.PredictionProteinInfo"],
    the_prediction_configuration: "prediction_configuration.PredictionConfiguration",
    a_project: "project.Project",
    the_custom_progress_signal: "custom_signals.ProgressSignal",
    the_pymol_lock: "locks.PyMOL_LOCK",
    the_disable_pymol_signal: "custom_signals.DisablePyMOLSignal"
) -> tuple:
    """Runs structure prediction for a monomeric protein.

    Args:
        the_prediction_protein_infos: a list with protein names and sequences to predict.
        the_prediction_configuration: a prediction configuration for the monomeric protein prediction.
        a_project: a project object that contains all the proteins.

    Returns:
        a tuple with the exit code and exit code description
    """
    # TODO: needs checks
    # TODO: needs tests!
    structure_prediction_obj = structure_prediction.StructurePrediction(
        the_prediction_protein_infos,
        the_prediction_configuration,
        a_project,
    )
    the_custom_progress_signal.emit_signal("Creating temp directories ...", 5)
    structure_prediction_obj.create_tmp_directories()
    logger.info("Tmp directories were created.")

    # <editor-fold desc="Creates fasta files for prediction">
    the_custom_progress_signal.emit_signal("Creating FASTA files ...", 10)
    try:
        structure_prediction_obj.create_fasta_files_for_prediction()
    except exception.FastaFilesNotCreatedError:
        logger.error("Fasta files were not created.")
        return (exit_codes.ERROR_WRITING_FASTA_FILES[0], exit_codes.ERROR_WRITING_FASTA_FILES[1])
    except exception.FastaFilesNotFoundError:
        logger.error("Fasta files were not found.")
        return (exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0], exit_codes.ERROR_FASTA_FILES_NOT_FOUND[1])
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        logger.info("Fasta files were successfully created.")

    # </editor-fold>

    # <editor-fold desc="Runs structure prediction">
    the_custom_progress_signal.emit_signal("Running ColabFold prediction ...", 25)
    try:
        structure_prediction_obj.run_prediction()
    except exception.PredictionEndedWithError:
        logger.error("Prediction ended with error.")
        return (exit_codes.ERROR_PREDICTION_FAILED[0], exit_codes.ERROR_PREDICTION_FAILED[1])
    else:
        logger.info("Prediction process finished.")

    # </editor-fold>

    # <editor-fold desc="Saves predicted protein to project">
    the_disable_pymol_signal.emit_signal("ColabFold Prediction")
    while the_pymol_lock.is_locked() is False:
        print("Waiting in separate thread for PyMOL LOCK ...")
        time.sleep(1)
    print("LOCK acquired.")

    the_custom_progress_signal.emit_signal("Saving best prediction results ...", 85)
    try:
        structure_prediction_obj.move_best_prediction_models()
        logger.info("Saved predicted pdb file into XML file.")
    except exception.UnableToFindColabfoldModelError:
        logger.error("Could not move rank 1 model, because it does not exists.")
        return (
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[0],
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[1],
        )
    except FileNotFoundError:
        logger.error("Could not move rank 1 model, because it does not exists.")
        return (
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[0],
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[1],
        )
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        logger.error("Could not move rank 1 model, because it does not exists.")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        subprocess.run(["wsl", "--shutdown"])
        logger.info("WSL gets shutdown.")
        return (exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])

    # </editor-fold>


def run_distance_analysis(
    a_list_with_analysis_names: list,
    a_project: "project.Project",
    the_settings: "settings.Settings",
    an_make_images_flag: bool,
    the_custom_progress_signal: "custom_signals.ProgressSignal",
    the_pymol_lock: "locks.PyMOL_LOCK",
    the_disable_pymol_signal: "custom_signals.DisablePyMOLSignal") -> tuple:
    """Runs the distance analysis for all protein pairs in the given job.

    Args:
        a_list_with_analysis_names: a list of the raw QListWidget analysis names.
        a_project: the current project.
        the_settings: the application settings.
        an_make_images_flag: True if images should be generated and False if not.

    Returns:
        a tuple containing exit codes.
    """
    # TODO: checks are needed
    logger.info("Running distance analysis in QThread using the Task class.")
    the_disable_pymol_signal.emit_signal("Distance Analysis")
    while the_pymol_lock.is_locked() is False:
        print("Waiting in separate thread for PyMOL LOCK ...")
        time.sleep(1)
    print("LOCK acquired.")
    try:
        analysis_runs = structure_analysis.Analysis(a_project)
        analysis_runs.analysis_list = analysis_util.transform_gui_input_to_practical_data(
            a_list_with_analysis_names,
            a_project,
            the_settings,
        )
        logger.debug(f"Analysis runs before actual analysis: {analysis_runs.analysis_list}")
        analysis_runs.run_analysis("distance", an_make_images_flag)
        logger.debug(f"Analysis runs after actual analysis: {analysis_runs.analysis_list}")
        for tmp_protein_pair in analysis_runs.analysis_list:
            tmp_protein_pair.db_project_id = a_project.get_id()
            copy_tmp_protein_pair = copy.deepcopy(tmp_protein_pair)
            with database_manager.DatabaseManager(str(a_project.get_database_filepath())) as db_manager:
                copy_tmp_protein_pair.set_id(db_manager.insert_new_protein_pair(copy_tmp_protein_pair))
            # Protein pair gets added to "a_project" argument of this function
            a_project.add_protein_pair(copy_tmp_protein_pair)

    except exception.UnableToSetupAnalysisError:
        logger.error("Setting up the analysis runs failed therefore the distance analysis failed.")
        return (
            exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[0],
            exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[1],
        )
    except Exception as e:
        logger.error(f"Unknown error: {e}")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        return (exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1], analysis_runs.analysis_list)
