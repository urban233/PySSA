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
"""This module contains worker classes for all processes done in the threadpool"""
import os
import subprocess
import logging
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pymol import cmd
from pyssa.internal.data_structures import structure_prediction
from pyssa.internal.data_structures import structure_analysis
from pyssa.internal.data_processing import data_transformer
from pyssa.internal.data_structures.data_classes import protein_analysis_info
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.io_pyssa import safeguard
from pyssa.util import prediction_util
from pyssa.util import constants
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.logging_pyssa import loggers
from pyssa.logging_pyssa import log_handlers
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project
    from pyssa.internal.data_structures import settings

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class WorkerSignals(QtCore.QObject):
    """This class is used to define signals for all worker classes"""

    # <editor-fold desc="Class attributes">
    """
    signal which returns no data, and should be emitted if process has successfully finished
    """
    finished = QtCore.pyqtSignal()
    """
    signal which returns tuple (exctype, value, traceback.format_exc()), and should be emitted if an error occurs
    """
    error = QtCore.pyqtSignal(tuple)
    """
    signal which returns object data returned from processing, and should be emitted if results are generated
    """
    result = QtCore.pyqtSignal(object)
    """
    signal which returns int progress, and should be emitted if the progress should get displayed
    """
    progress = QtCore.pyqtSignal(int)

    # </editor-fold>


class PredictionWorkerPool(QtCore.QRunnable):
    """This class is a worker class for the prediction process.

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    """

    # <editor-fold desc="Class attributes">
    """
    pyqt table which contains the proteins to predict
    """
    table: QtWidgets.QTableWidget
    """
    the configuration settings for the prediction
    """
    prediction_config: prediction_configuration.PredictionConfiguration
    """
    the current project in use
    """
    app_project: project.Project
    """
    the signals to use, for the worker
    """
    signals = WorkerSignals()

    # </editor-fold>

    def __init__(self,
                 table_prot_to_predict: QtWidgets.QTableWidget,
                 prediction_config: prediction_configuration.PredictionConfiguration,
                 app_project: 'project.Project') -> None:
        """Constructor.

        Args:
            table_prot_to_predict:
                pyqt table which contains the proteins to predict
            prediction_config:
                the prediction_config object
            app_project:
                current project

        Raises:
            ValueError: raised if an argument is illegal
        """
        super(PredictionWorkerPool, self).__init__()
        # <editor-fold desc="Checks">
        if not safeguard.Safeguard.check_if_value_is_not_none(table_prot_to_predict):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(prediction_config):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(app_project):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")

        # </editor-fold>

        self.table = table_prot_to_predict
        self.prediction_configuration = prediction_config
        self.app_project = app_project

    def run(self):
        """This function is a reimplementation of the QRunnable run method. It does the structure prediction.

        """
        predictions: list[tuple[str, str]] = prediction_util.get_prediction_name_and_seq_from_table(self.table)
        loggers.log_single_variable_value(constants.PREDICTION_WORKER_LOGGER, "run", "predictions", predictions)
        structure_prediction_obj = structure_prediction.StructurePrediction(predictions, self.prediction_configuration,
                                                                            self.app_project)
        structure_prediction_obj.create_tmp_directories()
        structure_prediction_obj.create_fasta_files_for_prediction()
        structure_prediction_obj.run_prediction()
        structure_prediction_obj.move_best_prediction_models()
        subprocess.run(["wsl", "--shutdown"])
        self.signals.finished.emit()


class AnalysisWorkerPool(QtCore.QRunnable):
    """This class is a worker class for the analysis process.

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    """

    # <editor-fold desc="Class attributes">
    """
    the list where all analysis runs are stored
    """
    list_analysis_overview: QtWidgets.QListWidget
    """
    the checkbox if images should be taken or not
    """
    cb_analysis_images: QtWidgets.QCheckBox
    """
    the status bar of the main window
    """
    status_bar: QtWidgets.QStatusBar
    """
    the current project in use
    """
    app_project: project.Project
    """
    the settings of pyssa
    """
    app_settings: settings.Settings
    """
    the signals to use, for the worker
    """
    signals = WorkerSignals()

    # </editor-fold>

    def __init__(self,
                 list_analysis_overview: QtWidgets.QListWidget,
                 cb_analysis_images: QtWidgets.QCheckBox,
                 status_bar: QtWidgets.QStatusBar,
                 app_project: 'project.Project',
                 app_settings: 'settings.Settings',
                 _init_batch_analysis_page) -> None:
        """Constructor

        Args:
            list_analysis_overview:
                list where all analysis runs are stored
            cb_analysis_images:
                checkbox which controls whether images should be created or not
            status_bar:
                status bar object of the main window
            app_project:
                project of the main window
            app_settings:
                settings of the main window
            _init_batch_analysis_page:
                function which clears the page, without parenthesis!!!

        Raises:
            ValueError: raised if an argument is illegal
        """
        super(AnalysisWorkerPool, self).__init__()
        # <editor-fold desc="Checks">
        if not safeguard.Safeguard.check_if_value_is_not_none(list_analysis_overview):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(cb_analysis_images):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(status_bar):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(app_project):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(app_settings):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")

        # </editor-fold>

        self.list_analysis_overview = list_analysis_overview
        self.cb_analysis_images = cb_analysis_images
        self.status_bar = status_bar
        self.app_project = app_project
        self.app_settings = app_settings
        self._init_batch_analysis_page = _init_batch_analysis_page

    def get_analysis_runs(self):
        """This function creates a data format which is used for the analysis

        """
        loggers.log_single_variable_value(constants.PREDICTION_WORKER_LOGGER,
                                          "get_analysis_runs",
                                          "self.app_project.proteins", self.app_project.proteins)
        loggers.log_single_variable_value(constants.PREDICTION_WORKER_LOGGER,
                                          "get_analysis_runs",
                                          "self.list_analysis_overview", self.list_analysis_overview)
        batch_analysis = []
        for row_no in range(self.list_analysis_overview.count()):
            tmp_batch_analysis = self.list_analysis_overview.item(row_no).text()
            separator_index = tmp_batch_analysis.find("_vs_")
            prot_1 = tmp_batch_analysis[:separator_index]

            if prot_1.find(";") != -1:
                prot_1_name = prot_1[:prot_1.find(";")]
                prot_1_chains = prot_1[prot_1.find(";") + 1:].split(",")
            else:
                prot_1_name = prot_1
                prot_1_chains = None
            prot_2 = tmp_batch_analysis[separator_index + 4:]
            if prot_2.find(";") != -1:
                prot_2_name = prot_2[:prot_2.find(";")]
                prot_2_chains = prot_2[prot_2.find(";") + 1:].split(",")
            else:
                prot_2_name = prot_2
                prot_2_chains = None

            loggers.log_multiple_variable_values(constants.PREDICTION_WORKER_LOGGER, "get_analysis_runs", [
                ("prot_1_name", prot_1_name),
                ("prot_2_name", prot_2_name),
                ("prot_1_chains", prot_1_chains),
                ("prot_2_chains", prot_2_chains),
            ])
            tmp_prot_1 = protein_analysis_info.ProteinAnalysisInfo(prot_1_name, prot_1_chains, tmp_batch_analysis)
            tmp_prot_2 = protein_analysis_info.ProteinAnalysisInfo(prot_2_name, prot_2_chains, tmp_batch_analysis)
            batch_analysis.append((tmp_prot_1, tmp_prot_2))

        transformer = data_transformer.DataTransformer(self.app_project, batch_analysis)
        # contains analysis-ready data format: list(tuple(prot_1, prot_2, export_dir, name), ...)
        return transformer.transform_data_for_analysis()

    def run_analysis_process(self, batch_data):
        """This function runs the analysis process.

        Args:
            batch_data:
                data in the analysis-ready format
        """
        analyzer = structure_analysis.Analysis()
        analyzer.run_analysis()

        loggers.log_single_variable_value(constants.PREDICTION_WORKER_LOGGER, "run_analysis_process",
                                          "batch_data", batch_data)
        for analysis_data in batch_data:
            if not os.path.exists(analysis_data[2]):
                os.mkdir(analysis_data[2])
            else:
                basic_boxes.ok("Structure Analysis", f"The structure analysis: {analysis_data[3]} already exists!",
                               QtWidgets.QMessageBox.Critical)
                constants.PREDICTION_WORKER_LOGGER.warning(f"The structure analysis: {analysis_data[3]} already exists!")
                self._init_batch_analysis_page()
                return

            cmd.reinitialize()
            structure_analysis_obj = structure_analysis.StructureAnalysis(
                reference_protein=[analysis_data[0]], model_proteins=[analysis_data[1]],
                ref_chains=analysis_data[0].chains, model_chains=analysis_data[1].chains,
                export_dir=analysis_data[2], cycles=self.app_settings.get_cycles(),
                cutoff=self.app_settings.get_cutoff(),
            )
            if self.cb_analysis_images.isChecked():
                structure_analysis_obj.response_create_images = True
            structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.ref_chains,
                                                                 structure_analysis_obj.reference_protein)
            structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.model_chains,
                                                                 structure_analysis_obj.model_proteins)
            protein_pairs = structure_analysis_obj.create_protein_pairs()
            structure_analysis_obj.do_analysis_in_pymol(protein_pairs, self.status_bar)
            protein_pairs[0].name = analysis_data[3]
            protein_pairs[0].cutoff = self.app_settings.cutoff
            self.app_project.add_protein_pair(protein_pairs[0])
            protein_pairs[0].serialize_protein_pair(self.app_project.get_objects_protein_pairs_path())
            self.app_project.serialize_project(self.app_project.project_path, "project")

    def run(self):
        """This function is a reimplementation of the QRunnable run method.

        """
        self.run_analysis_process(self.get_analysis_runs())
        self.signals.finished.emit()
