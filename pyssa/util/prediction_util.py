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
import shutil
from pyssa.util import constants
from pyssa.util import exception
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


def get_relaxed_rank_1_pdb_file(
    proteins_to_predict: list[prediction_protein_info.PredictionProteinInfo],
) -> list[tuple]:
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


def delete_pyssa_colabfold_directory_in_wsl2() -> bool:
    """Deletes the pyssa_colabfold directory in WSL2.

    Raises:
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    tmp_pyssa_colabfold_path: str = r"\\wsl$\almaColabfold9\home\rhel_user\pyssa_colabfold"
    try:
        if os.path.exists(tmp_pyssa_colabfold_path):
            shutil.rmtree(tmp_pyssa_colabfold_path)
            return True
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToDeleteDirectoryError("Could not remove the pyssa_colabfold dir.")
    return False

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "root",
    #          "rm", "-r", "/home/rhel_user/pyssa_colabfold", "||", "true",  # fixme: don't know if this works
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to delete the pyssa_colabfold directory in WSL2!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully deleted the pyssa_colabfold directory.")


def delete_scratch_directory_in_wsl2() -> bool:
    """Deletes the scratch directory in WSL2.

    Raises:
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    tmp_scratch_path: str = r"\\wsl$\almaColabfold9\home\rhel_user\scratch"
    try:
        if os.path.exists(tmp_scratch_path):
            shutil.rmtree(tmp_scratch_path)
            return True
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToDeleteDirectoryError("Could not remove the scratch dir.")
    return False

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "root",
    #          "rm", "-r", "/home/rhel_user/scratch", "||", "true",
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to delete the scratch directory in WSL2!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully deleted the scratch directory.")


def copy_pyssa_colabfold_directory_to_wsl2() -> None:
    """Copies the pyssa_colabfold directory to WSL2.

    Raises:
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    tmp_pyssa_colabfold_wsl_path: str = r"\\wsl$\almaColabfold9\home\rhel_user\pyssa_colabfold"
    tmp_pyssa_colabfold_windows_path: str = str(pathlib.Path(f"{constants.PLUGIN_PATH}/pyssa_colabfold"))
    if not os.path.exists(tmp_pyssa_colabfold_windows_path):
        logger.error(f"The variable 'tmp_pyssa_colabfold_windows_path' is illegal: {tmp_pyssa_colabfold_windows_path}!")
        raise exception.DirectoryNotFoundError("An variable is illegal.")
    try:
        shutil.copytree(tmp_pyssa_colabfold_windows_path, tmp_pyssa_colabfold_wsl_path)
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToCopyDirectoryError("Could not copy the pyssa_colabfold dir.")
    else:
        logger.info("Successfully copied the pyssa_colabfold directory.")

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "rhel_user",
    #          "cp", "-r",
    #          f"{constants.PLUGIN_PATH_WSL_NOTATION}/pyssa_colabfold",
    #          "/home/rhel_user/",
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to copy the pyssa_colabfold directory in WSL2!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully copied the pyssa_colabfold directory.")


def delete_original_batch_py_file() -> bool:
    """Deletes the original batch.py file of colabfold.

    Raises:
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    # Remove original batch.py of Colabfold
    tmp_batch_py_filepath: str = r"\\wsl$\almaColabfold9/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py"
    try:
        if os.path.exists(tmp_batch_py_filepath):
            os.remove(tmp_batch_py_filepath)
            return True
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToDeleteDirectoryError("Could not remove the original batch.py file.")
    return False

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "root",
    #          "rm",
    #          "/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py", "||", "true",
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to delete the batch.py file!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully deleted the batch.py file.")


def copy_modified_batch_py_file() -> None:
    """Copies the modified batch.py file to the WSL2.

    Raises:
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    tmp_batch_py_filepath_wsl: str = r"\\wsl$\almaColabfold9/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py"
    tmp_batch_py_filepath_windows: str = f"{constants.PLUGIN_PATH}/pyssa_colabfold/colabfold_sub/batch.py"
    try:
        shutil.copy2(tmp_batch_py_filepath_windows, tmp_batch_py_filepath_wsl)
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToCopyFileError("")
    else:
        logger.info("Successfully copied the modified batch.py file.")

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "root",
    #          "cp",
    #          f"{constants.PLUGIN_PATH_WSL_NOTATION}/pyssa_colabfold/colabfold_sub/batch.py",
    #          "/home/rhel_user/localcolabfold/colabfold-conda/lib/python3.10/site-packages/colabfold/batch.py",
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to copy the modified batch.py file!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully copied the modified batch.py file.")


def create_fasta_directory_in_wsl2(the_fasta_path: str) -> None:
    """Creates the fasta directory inside the WSL2.

    Raises:
        IllegalArgumentError: If the argument is None.
        SubprocessExecutionError: If return code of subprocess is non-zero.
    """
    # <editor-fold desc="Checks">
    if the_fasta_path is None:
        logger.error("The argument filename is illegal.")
        raise exception.IllegalArgumentError("")

    # </editor-fold>
    try:
        os.makedirs(the_fasta_path)
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToCreateDirectoryError("")
    else:
        logger.info("Successfully created the fasta directory.")

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "rhel_user",
    #          "mkdir", "-p", the_fasta_path,
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to create the fasta directory in WSL2!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully created the fasta directory.")


def create_pdb_directory_in_wsl2(the_pdb_path: str) -> None:
    """Creates the pdb directory inside the WSL2.

    Raises:
        IllegalArgumentError: If the argument is None.
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    # <editor-fold desc="Checks">
    if the_pdb_path is None:
        logger.error("The argument filename is illegal.")
        raise exception.IllegalArgumentError("")

    # </editor-fold>

    try:
        os.makedirs(the_pdb_path)
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToCreateDirectoryError("")
    else:
        logger.info("Successfully created the fasta directory.")

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "rhel_user",
    #          "mkdir", "-p", the_pdb_path,
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to create the pdb directory in WSL2!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully created the pdb directory.")


def copy_fasta_files_from_windows_to_wsl2(the_fasta_path: str) -> None:
    """Copies fasta files from Windows host to WSL2.

    Raises:
        IllegalArgumentError: If an argument is None.
        SubprocessExecutionError: If return code of subprocess is non-zero
    """
    # <editor-fold desc="Checks">
    if the_fasta_path is None:
        logger.error("The argument filename is illegal.")
        raise exception.IllegalArgumentError("")

    # </editor-fold>

    tmp_fasta_files_path: str = f"{constants.PREDICTION_FASTA_DIR}/*.fasta"
    try:
        shutil.copytree(tmp_fasta_files_path, the_fasta_path)
    except Exception as e:
        logger.error("Unexpected error!", e)
        raise exception.UnableToCopyFileError("")
    else:
        logger.info("Successfully copied the fasta files.")

    # try:
    #     subprocess.run(
    #         ["wsl", "-d", "almaColabfold9",
    #          "-u", "rhel_user",
    #          "cp", "-r", f"{settings_dir_unix_notation}/scratch/local_predictions/fasta/*.fasta", the_fasta_path,
    #          ], check=True,
    #     )
    # except subprocess.CalledProcessError:
    #     message: str = "An error occurred while trying to copy the fasta files from Windows host to WSL2 filesystem!"
    #     logger.error(message)
    #     raise exception.SubprocessExecutionError(message)
    # except Exception as e:
    #     logger.error("Unexpected error!", e)
    #     raise exception.SubprocessExecutionError(f"Unexpected error!: {e}")
    # else:
    #     logger.info("Successfully copied the fasta files.")
