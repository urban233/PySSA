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
import subprocess
import logging
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from internal.data_structures.data_classes import prediction_protein_info
from pyssa.internal.data_structures import structure_prediction
from pyssa.internal.data_structures import structure_analysis
from pyssa.internal.data_processing import data_transformer
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.io_pyssa import safeguard
from pyssa.util import prediction_util
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
    app_project: 'project.Project'
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
        predictions: list[prediction_protein_info.PredictionProteinInfo] = prediction_util.get_prediction_name_and_seq_from_table(self.table)
        structure_prediction_obj = structure_prediction.StructurePrediction(predictions, self.prediction_configuration,
                                                                            self.app_project)
        structure_prediction_obj.create_tmp_directories()
        logger.info("Tmp directories were created.")
        structure_prediction_obj.create_fasta_files_for_prediction()
        logger.info("Fasta files were created.")
        structure_prediction_obj.run_prediction()
        logger.info("Prediction ran successfully.")
        structure_prediction_obj.move_best_prediction_models()
        logger.info("Saved predicted pdb file into XML file.")
        subprocess.run(["wsl", "--shutdown"])
        logger.info("WSL gets shutdown.")
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
    app_project: 'project.Project'
    """
    the settings of pyssa
    """
    app_settings: 'settings.Settings'
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

    def transform_gui_input_to_practical_data(self) -> list:
        """This function transforms the input from the gui to a practical data basis which can be used to setup
        analysis runs

        """
        distance_analysis_runs = []
        for row_no in range(self.list_analysis_overview.count()):
            distance_analysis_runs.append(
                data_transformer.DistanceAnalysisDataTransformer(
                    self.list_analysis_overview.item(row_no).text(),
                    self.app_project,
                    self.app_settings
                ).transform_gui_input_to_distance_analysis_object()
            )
        return distance_analysis_runs

    def set_up_analysis_runs(self) -> 'structure_analysis.Analysis':
        """This function creates protein pairs and distance analysis objects for the analysis runs.

        """
        analysis_runs = structure_analysis.Analysis(self.app_project)
        for tmp_distance_analysis in self.transform_gui_input_to_practical_data():
            analysis_runs.analysis_list.append(tmp_distance_analysis)
        return analysis_runs

    def run_analysis(self) -> None:
        self.set_up_analysis_runs().run_analysis()

    def run(self):
        """This function is a reimplementation of the QRunnable run method.

        """
        # # transform raw data: get data in a form to create protein pair and analysis objects
        # self.transform_gui_input_to_practical_data()
        # # setup data for analysis: create protein_pair and distance analysis objects
        # self.set_up_analysis_runs()
        logger.debug(f"Memory address of worker {self}")
        # do the analysis runs
        self.run_analysis()
        # emit finish signal
        self.signals.finished.emit()

    def __del__(self):
        logger.debug("Worker has been destroyed!")
