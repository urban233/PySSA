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
"""This module contains the colabbatch class."""
import os
import logging
import subprocess
import shutil
import pathlib
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Colabbatch:
    """This class contains information about the colabbatch prediction engine."""

    # <editor-fold desc="Class attributes">
    """
    name of the user who has opened the pyssa
    """
    user_name: str = os.getlogin()
    """
    path where the fasta files will be stored, in unix path format
    """
    fasta_path: str = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/fasta"
    """
    path where the pdb files will be stored, in unix path format
    """
    pdb_path: str = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/pdb"
    """
    the configuration settings for the prediction
    """
    prediction_config: prediction_configuration.PredictionConfiguration

    # </editor-fold>

    def __init__(self, prediction_configuration: prediction_configuration.PredictionConfiguration) -> None:
        """Constructor.

        Args:
            prediction_configuration:
                the configuration settings for the prediction

        Raises:
            ValueError: raised if an argument is illegal
        """
        # <editor-fold desc="Checks">
        if prediction_configuration.amber_force_field is None:
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if prediction_configuration.templates == "":
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")

        # </editor-fold>

        self.prediction_configuration = prediction_configuration

    def run_prediction(self):
        """This function starts the wsl and runs a prediction."""
        #subprocess.run([constants.POWERSHELL_EXE, str(constants.CONVERT_DOS_TO_UNIX_EXE), str(constants.COLABFOLD_PREDICT_SCRIPT_WIN)])
        #subprocess.run([constants.POWERSHELL_EXE, "-ExecutionPolicy", "Bypass", "-File", constants.CONVERT_DOS_TO_UNIX])
        #subprocess.run([str(constants.CONVERT_DOS_TO_UNIX)])
        os.startfile(str(constants.CONVERT_DOS_TO_UNIX))
        # running prediction script
        logger.debug(f"templates: {self.prediction_configuration.templates}, amber: {self.prediction_configuration.amber_force_field}")
        if self.prediction_configuration.templates == "none" and self.prediction_configuration.amber_force_field is True:
            logger.info("Run prediction with no templates and with amber force field correction.")
            try:
                subprocess.run(["wsl", "--set-default", constants.WSL_DISTRO_NAME])
                subprocess.run(["wsl", constants.COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT,
                                self.fasta_path, self.pdb_path])
                subprocess.run(["wsl", "--shutdown"])
            except OSError:
                logger.error("Something went wrong during the prediction process.")
                shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
                return
        elif self.prediction_configuration.templates == "none" and self.prediction_configuration.amber_force_field is False:
            logger.info("Run prediction with no templates and with no amber force field correction.")
            try:
                subprocess.run(["wsl", "--set-default", constants.WSL_DISTRO_NAME])
                subprocess.run(["wsl", constants.COLABFOLD_PREDICT_NO_AMBER_AND_TEMPLATES_SCRIPT,
                                self.fasta_path, self.pdb_path])
                subprocess.run(["wsl", "--shutdown"])
            except OSError:
                logger.error("Something went wrong during the prediction process.")
                shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
                return
        elif self.prediction_configuration.templates == "pdb70" and self.prediction_configuration.amber_force_field is False:
            logger.info("Run prediction with default pdb70 templates and with no amber force field correction.")
            try:
                subprocess.run(["wsl", "--set-default", constants.WSL_DISTRO_NAME])
                subprocess.run(["wsl", constants.COLABFOLD_PREDICT_NO_AMBER_SCRIPT,
                                self.fasta_path, self.pdb_path])
                subprocess.run(["wsl", "--shutdown"])
            except OSError:
                logger.error("Something went wrong during the prediction process.")
                shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
                return
        elif self.prediction_configuration.templates == "pdb70" and self.prediction_configuration.amber_force_field is True:
            logger.info("Run prediction with default pdb70 templates and with amber force field correction.")
            # new idea --------
            shell_result = subprocess.run(["wsl", "--set-default", constants.WSL_DISTRO_NAME], stdout=subprocess.PIPE, text=True)
            if shell_result.returncode == 0:
                print("Command executed successfully")
                output = shell_result.stdout
                # Process the output as needed
                lines = output.split('\n')
                for line in lines:
                    print(line)
            else:
                print(f"Command failed with return code {shell_result.returncode}")
            shell_result = subprocess.run(
                ["wsl", "/home/ubuntu_colabfold/.pyssa/localcolabfold/colabfold-conda/bin/colabfold_batch", self.fasta_path, self.pdb_path, "--amber", "--templates"], stdout=subprocess.PIPE, text=True)
            if shell_result.returncode == 0:
                print("Colabfold_batch command executed successfully")
            else:
                print(f"Command failed with return code {shell_result.returncode}")
            with open(constants.COLABFOLD_LOG_FILE_PATH) as tmp_file:
                contents = tmp_file.read()
                logger.debug(contents)
            subprocess.run(["wsl", "--shutdown"])

            # old idea --------
            # try:
            #     subprocess.run(["wsl", "--set-default", constants.WSL_DISTRO_NAME])
            #     subprocess.run(["wsl", constants.COLABFOLD_PREDICT_SCRIPT,
            #                     self.fasta_path, self.pdb_path])
            #     subprocess.run(["wsl", "--shutdown"])
            # except OSError:
            #     logger.error("Something went wrong during the prediction process.")
            #     shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
            #     return
        else:
            logger.error(f"Invalid prediction configuration. templates: {self.prediction_configuration.templates}, amber: {self.prediction_configuration.amber_force_field}")
            raise ValueError("Invalid prediction configuration.")
