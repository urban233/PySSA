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
"""Module for structure prediction class"""

import os
import pathlib
import subprocess
import shutil
import pymol
import logging
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures.data_classes import prediction_list
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.internal.portal import pymol_io
from pyssa.util import constants
from pyssa.logging_pyssa import loggers
from pyssa.logging_pyssa import generic_messages


class StructurePrediction:

    def __init__(self,
                 predictions: list[tuple[str, str]],
                 prediction_config: prediction_configuration.PredictionConfiguration,
                 current_project: project.Project):
        self.predictions = predictions
        self.prediction_configuration = prediction_config
        self.project = current_project
        messages = generic_messages.get_variables_values(__file__,
                                                         "constructor", [("predictions", self.predictions),
                                                                         ("prediction_configuration", self.prediction_configuration),
                                                                         ("project", self.project)])
        loggers.log_multiple_messages(loggers.prediction_worker, logging.DEBUG, messages)

    @staticmethod
    def create_tmp_directories():
        """This function creates tmp directories in a scratch folder to organize prediction inputs and outputs

        """
        # TODO: is there a more elegant way to do it?
        if not os.path.exists(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions")):
            os.mkdir(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
        if not os.path.exists(constants.PREDICTION_FASTA_DIR):
            os.mkdir(constants.PREDICTION_FASTA_DIR)
        if not os.path.exists(constants.PREDICTION_PDB_DIR):
            os.mkdir(constants.PREDICTION_PDB_DIR)

    def create_fasta_files_for_prediction(self):
        """This function creates fasta file based on the predictions in the object

        """
        prot_entries = []
        last_header = self.predictions[0][0]
        pred_list = prediction_list.PredictionList("", [])
        if len(self.predictions) == 1:
            pred_list.protein_name = self.predictions[0][0]
            pred_list.protein_sequence.append(self.predictions[0][1])
            prot_entries.append(pred_list)
        else:
            for tmp_prediction in self.predictions:
                current_header = tmp_prediction[0]
                if last_header == current_header:
                    pred_list.protein_name = tmp_prediction[0]
                    pred_list.protein_sequence.append(tmp_prediction[1])
                else:
                    prot_entries.append(pred_list)
                    pred_list = prediction_list.PredictionList("", [])
                    last_header = current_header
        loggers.prediction_worker.debug(generic_messages.get_basic_variable_value(__file__,
                                                                                  "create_fasta_files_for_prediction",
                                                                                  "pred_list", pred_list))
        loggers.prediction_worker.debug(generic_messages.get_basic_variable_value(__file__,
                                                                                  "create_fasta_files_for_prediction",
                                                                                  "prot_entries", prot_entries))
        for tmp_prot_to_predict in prot_entries:
            tmp_prot_to_predict.write_fasta_file()
        if len(os.listdir(constants.PREDICTION_FASTA_DIR)) == 0:
            loggers.prediction_worker.critical("No fasta files were created!!!")
            raise FileNotFoundError
        else:
            for tmp_file in os.listdir(constants.PREDICTION_FASTA_DIR):
                loggers.prediction_worker.debug(generic_messages.get_basic_variable_value(__file__, "create_fasta_files_for_prediction",
                                                                                          "filename", tmp_file))

    def run_prediction(self):
        user_name = os.getlogin()
        fasta_path = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/fasta"
        pdb_path = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/pdb"
        # running prediction script
        if self.prediction_configuration.templates == "none":
            try:
                subprocess.run([constants.POWERSHELL_EXE, constants.CONVERT_DOS_TO_UNIX])
                subprocess.run(["wsl", constants.COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT,
                                fasta_path, pdb_path])
                subprocess.run(["wsl", "--shutdown"])
            except OSError:
                shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
                return
        else:
            try:
                subprocess.run([constants.POWERSHELL_EXE, constants.CONVERT_DOS_TO_UNIX])
                subprocess.run(["wsl", constants.COLABFOLD_PREDICT_SCRIPT,
                                fasta_path, pdb_path])
                subprocess.run(["wsl", "--shutdown"])
            except OSError:
                shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
                return

    def move_best_prediction_models(self):
        """This function moves the best prediction model(s) to the project directory

        """
        prediction_results: list[str] = os.listdir(pathlib.Path(constants.PREDICTION_PDB_DIR))
        for tmp_prediction in self.predictions:
            for filename in prediction_results:
                check = filename.find(f"{tmp_prediction[0]}_relaxed_rank_1")
                if check != -1:
                    src = pathlib.Path(f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{filename}")
                    dest = pathlib.Path(
                        f"{self.project.get_pdb_path()}/{filename}")
                    shutil.copy(src, dest)
                    os.rename(f"{self.project.get_pdb_path()}/{filename}",
                              f"{self.project.get_pdb_path()}/{tmp_prediction[0]}.pdb")
                    tmp_protein = protein.Protein(tmp_prediction[0], pathlib.Path(self.project.get_pdb_path()))
                    self.project.add_existing_protein(tmp_protein)
                    break
        shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
        try:
            pymol_io.load_protein(self.project.proteins[0])
        except pymol.CmdException:
            print("Loading the model failed.")
            return
        except FileNotFoundError:
            print("Prediction was unsuccessful")
            return
