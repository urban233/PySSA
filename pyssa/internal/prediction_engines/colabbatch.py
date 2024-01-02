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
import zmq
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.util import constants, prediction_util
from pyssa.util import exception

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
    fasta_path: str = "/home/rhel_user/scratch/local_predictions/fasta"
    """
    path where the pdb files will be stored, in unix path format
    """
    pdb_path: str = "/home/rhel_user/scratch/local_predictions/pdb"
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

        str_conversion_1 = constants.SETTINGS_DIR.replace("\\", "/")
        str_conversion_2 = str_conversion_1.replace(":", "")
        str_conversion_3 = str_conversion_2.replace("C", "c")
        self.settings_dir_unix_notation = f"/mnt/{str_conversion_3}"

    def setup_prediction_service(self) -> None:
        """Sets up the structure prediction service in WSL2 almaColabfold9 distro.

        Raises:
            Wsl2PreparationFailedError: If an error occurs during prediction service setup
        """
        try:
            if prediction_util.delete_pyssa_colabfold_directory_in_wsl2():
                logger.info("The pyssa_colabfold directory has been successfully deleted.")
            else:
                logger.warning("The pyssa_colabfold directory was not found.")
            prediction_util.copy_pyssa_colabfold_directory_to_wsl2()
            prediction_util.move_modified_batch_file_to_wsl2()
        except exception.SubprocessExecutionError:
            logger.error("An error occurred during subprocess execution!")
            raise exception.Wsl2PreparationFailedError("")
        except Exception as e:
            logger.error("Unexpected error!", e)
            raise exception.Wsl2PreparationFailedError(f"Unexpected error!: {e}")
        else:
            logger.info("Prediction service setup was successful.")

    def send_prediction_request(self) -> bool:
        """Sends structure prediction request based on custom arguments.

        Raises:
            PredictionEndedWithError: if prediction ended with any kind of error
        """
        # Start service process in WSL2
        service_process = subprocess.Popen(
            [
                "wsl",
                "-d",
                "almaColabfold9",
                "-u",
                "rhel_user",
                "/home/rhel_user/localcolabfold/colabfold-conda/bin/python3",
                "/home/rhel_user/pyssa_colabfold/service.py",
            ],
        )

        context = zmq.Context()
        socket = context.socket(zmq.REQ)
        # Connect to the WSL2 IP address and port
        socket.connect("tcp://127.0.0.1:7016")

        message = {
            "fasta_dir": str(constants.PREDICTION_FASTA_DIR),
            "pdb_dir": str(constants.PREDICTION_PDB_DIR),
            "use_amber": self.prediction_configuration.amber_force_field,
            "use_templates": self.prediction_configuration.templates,
        }
        socket.send_json(message)

        result = socket.recv_json()
        # Check results for any errors
        for tmp_key in result.keys():
            logger.debug(f"Response from service: {result}")
            if tmp_key == "error":
                service_process.terminate()
                raise exception.PredictionEndedWithError(f"Prediction run into the error: {result[tmp_key]}!")
        tmp_output: bool = True
        service_process.terminate()
        return tmp_output

    def run_prediction(self) -> None:
        """This function starts the wsl2, podman machine, podman container and runs a prediction.

        Raises:
            PredictionEndedWithError: if prediction ended with any kind of error
        """
        logger.info("Setup of prediction service.")
        try:
            self.setup_prediction_service()
        except exception.Wsl2PreparationFailedError:
            logger.error("Prediction service setup failed!")
            raise exception.PredictionEndedWithError("")

        logger.info("Send prediction request to WSL2.")
        try:
            tmp_output = self.send_prediction_request()
        except exception.PredictionEndedWithError:
            logger.error("Prediction ended with error.")
            subprocess.run(["wsl", "--shutdown"])
            raise exception.PredictionEndedWithError("")
        else:
            logger.info("Received success message of prediction service.")
            subprocess.run(["wsl", "--shutdown"])
            if tmp_output:
                logger.info("Prediction process finished, copying results ...")
                try:
                    subprocess.run(
                        [
                            "wsl",
                            "-d",
                            "almaColabfold9",
                            "-u",
                            "rhel_user",
                            "cp",
                            "-r",
                            "/home/rhel_user/scratch/local_predictions/pdb",
                            f"{self.settings_dir_unix_notation}/scratch/local_predictions",
                        ],
                        check=True,
                    )
                except subprocess.CalledProcessError:
                    logger.error("Could not copy prediction results to Windows host!")
                    raise exception.PredictionEndedWithError("Could not copy prediction results to Windows host!")
                else:
                    logger.info("Copying prediction results to Windows host was successful.")
            else:
                logger.error("Prediction finished with errors! No results were copied to the Windows host!")
                raise exception.PredictionEndedWithError("")
        return
