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
import sys

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
    fasta_path: str = f"/home/rhel_user/scratch/local_predictions/fasta"
    """
    path where the pdb files will be stored, in unix path format
    """
    pdb_path: str = f"/home/rhel_user/scratch/local_predictions/pdb"
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

    def run_prediction(self) -> None:
        """This function starts the wsl2, podman machine, podman container and runs a prediction.

        Raises:
            RuntimeError: if powershell command fails
        """
        # string conversion
        str_conversion_1 = constants.SETTINGS_DIR.replace("\\", "/")
        str_conversion_2 = str_conversion_1.replace(":", "")
        str_conversion_3 = str_conversion_2.replace("C", "c")
        settings_dir_unix_notation = f"/mnt/{str_conversion_3}"
        # new way
        # list of paths of the fasta dir and pdb dir which is needed by colabbatch
        colabbatch_paths = [self.fasta_path, self.pdb_path]

        # <editor-fold desc="Linux pre-process">
        # lists with podman commands which need to be used for the prediction process

        podman_machine_start_command = ["podman", "machine", "start"]
        podman_machine_stop_command = ["podman", "machine", "stop"]
        podman_container_delete_command = ["podman", "container", "rm", constants.CONTAINER_NAME]
        podman_container_stop_command = ["podman", "container", "stop", constants.CONTAINER_NAME]
        podman_container_create_command = ["podman", "run", "-itd", "-v", f"{settings_dir_unix_notation}:/home/ubuntu_colabfold/{os.getlogin()}", "--name", constants.CONTAINER_NAME, constants.IMAGE_NAME]
        podman_container_start_command = ["podman", "container", "start", constants.CONTAINER_NAME]
        podman_container_exec_colabbatch_command = ["podman", "container", "exec", constants.CONTAINER_NAME,
                                                    "/home/ubuntu_colabfold/localcolabfold/colabfold-conda/bin/colabfold_batch"]

        if sys.platform.startswith("linux"):
            # podman virtual machine gets started to be able to run containers
            powershell_result = subprocess.run(podman_machine_start_command)
            if powershell_result.returncode == 0:
                logger.info("Podman machine started.")
            elif powershell_result.returncode == 125:
                logger.warning("Podman machine is already running.")
            else:
                logger.error(f"Podman machine could not be started! "
                             f"Command exited with code: {powershell_result.returncode}")
                raise RuntimeError("Podman machine could not be started!")

            subprocess.run(podman_container_stop_command)
            subprocess.run(podman_container_delete_command)
            powershell_result = subprocess.run(podman_container_create_command)
            if powershell_result.returncode != 0:
                logger.error(f"Podman container could not be created! "
                             f"Command exited with code: {powershell_result.returncode}")
                raise RuntimeError("Podman container could not be created!")
        # </editor-fold>

        # checks cases for args to use in colabbatch command
        if self.prediction_configuration.templates == "pdb70" and self.prediction_configuration.amber_force_field is True:
            logger.info("Run prediction with default pdb70 templates and with amber force field correction.")
            colabbatch_args = ["--amber", "--templates"]
        elif self.prediction_configuration.templates == "none" and self.prediction_configuration.amber_force_field is True:
            logger.info("Run prediction with no templates and with amber force field correction.")
            colabbatch_args = ["--amber"]
        elif self.prediction_configuration.templates == "pdb70" and self.prediction_configuration.amber_force_field is False:
            logger.info("Run prediction with default pdb70 templates and with no amber force field correction.")
            colabbatch_args = ["--templates"]
        elif self.prediction_configuration.templates == "none" and self.prediction_configuration.amber_force_field is False:
            logger.info("Run prediction with no templates and with no amber force field correction.")
            colabbatch_args = []
        else:
            logger.error(
                f"Invalid prediction configuration. templates: {self.prediction_configuration.templates}, "
                f"amber: {self.prediction_configuration.amber_force_field}",
            )
            raise ValueError("Invalid prediction configuration.")

        if sys.platform.startswith("linux"):
            # starts localcolabfold-container, which runs the prediction
            powershell_result = subprocess.run(podman_container_start_command)
            if powershell_result.returncode != 0:
                logger.error(f"Podman localcolabfold-container could not be started! "
                             f"Command exited with code: {powershell_result.returncode}")
                raise RuntimeError("Podman localcolabfold-container could not be started!")

            # concatenating lists to create the complete prediction command run in the container
            powershell_result = subprocess.run(podman_container_exec_colabbatch_command + colabbatch_paths + colabbatch_args)
            if powershell_result.returncode != 0:
                logger.error(f"Podman exec command which runs colabbatch failed! "
                             f"Command exited with code: {powershell_result.returncode}")
                raise RuntimeError("Podman exec command which runs colabbatch failed!")

            powershell_result = subprocess.run(podman_container_stop_command)
            if powershell_result.returncode != 0:
                logger.error(f"Podman container could not be stopped! "
                             f"Command exited with code: {powershell_result.returncode}")
                raise RuntimeError("Podman container could not be stopped!")
            powershell_result = subprocess.run(podman_container_delete_command)
            # if powershell_result.returncode != 0:
            #     logger.error(f"Podman container could not be deleted! "
            #                  f"Command exited with code: {powershell_result.returncode}")
            #     raise RuntimeError("Podman container could not be deleted")

            # stops podman virtual machine after completing prediction process
            powershell_result = subprocess.run(podman_machine_stop_command)
            if powershell_result.returncode != 0:
                logger.error(f"Podman machine could not be shutdown! "
                             f"Command exited with code: {powershell_result.returncode}")
                raise RuntimeError("Podman machine could not be shutdown!")
        elif sys.platform.startswith("win32"):
            if not os.path.exists(constants.WSL_SCRATCH_DIR):
                subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch"])
                subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch/local_predictions"])
                subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch/local_predictions/fasta"])
                subprocess.run(["wsl", "-d", "almaColabfold9", "mkdir", "/home/rhel_user/scratch/local_predictions/pdb"])

            print(settings_dir_unix_notation)
            subprocess.run(["wsl", "-d", "almaColabfold9", "cp", "-r", f"{settings_dir_unix_notation}/scratch/local_predictions/fasta", "/home/rhel_user/scratch/local_predictions/"])
            # concatenating lists to create the complete prediction command run in the container
            wsl_almacolabfold_exec_colabbatch_command = ["wsl", "-d", "almaColabfold9", "/home/rhel_user/localcolabfold/colabfold-conda/bin/colabfold_batch"]
            complete_command = wsl_almacolabfold_exec_colabbatch_command + colabbatch_paths + colabbatch_args
            print(complete_command)
            powershell_result = subprocess.run(
                wsl_almacolabfold_exec_colabbatch_command + colabbatch_paths + colabbatch_args
            )
            if powershell_result.returncode != 0:
                logger.error(f"WSL2 almaColabfold9 exec command which runs colabbatch failed! "
                             f"Command exited with code: {powershell_result.returncode}")
                raise RuntimeError("WSL2 almaColabfold9 exec command which runs colabbatch failed!")
            elif powershell_result.returncode == 0:
                logger.info("Prediction process was successful, copying results ...")
                subprocess.run(
                    [
                        "wsl", "-d", "almaColabfold9", "cp", "-r", "/home/rhel_user/scratch/local_predictions/pdb",
                        f"{settings_dir_unix_notation}/scratch/local_predictions"
                     ]
                )
        else:
            print("Unsupported operating system detected.")
        # shuts down the WSL2 environment used by the podman virtual machine
        powershell_result = subprocess.run(["wsl", "--shutdown"])
        if powershell_result.returncode != 0:
            logger.error(f"WSL2 could not be shutdown! "
                         f"Command exited with code: {powershell_result.returncode}")
            raise RuntimeError("WSL2 could not be shutdown!")
