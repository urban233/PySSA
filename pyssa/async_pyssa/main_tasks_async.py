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
"""Module for all asynchronous functions used in the main tasks."""
import logging
import subprocess

from pyssa.internal.data_structures import structure_prediction
from pyssa.internal.thread.async_pyssa import custom_signals
from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception, exit_codes

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def predict_protein_with_colabfold(
    the_prediction_protein_infos: list["prediction_protein_info.PredictionProteinInfo"],
    the_prediction_configuration: "prediction_configuration.PredictionConfiguration",
    a_project: "project.Project",
    the_custom_progress_signal: "custom_signals.ProgressSignal"
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
