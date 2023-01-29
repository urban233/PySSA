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
import logging
import os
import subprocess
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pymol import cmd
from pyssa.internal.data_structures import structure_prediction
from pyssa.internal.data_structures import structure_analysis
from pyssa.internal.data_processing import data_transformer
from pyssa.internal.data_structures.data_classes import protein_analysis_info
from pyssa.util import gui_utils
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.logging_pyssa import loggers
from pyssa.logging_pyssa import generic_messages


class WorkerSignals(QtCore.QObject):
    """This class is used to define signals for all worker classes

    Available signals:
        finished
            No data
        error
            tuple (exctype, value, traceback.format_exc() )
        result
            object data returned from processing, anything
        progress
            int indicating % progress
    """
    finished = QtCore.pyqtSignal()
    error = QtCore.pyqtSignal(tuple)
    result = QtCore.pyqtSignal(object)
    progress = QtCore.pyqtSignal(int)


class PredictionWorkerPool(QtCore.QRunnable):
    """This class is a worker class for the prediction process

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    """
    def __init__(self,
                 table_prot_to_predict,
                 prediction_configuration,
                 app_project):
        super(PredictionWorkerPool, self).__init__()
        self.table = table_prot_to_predict
        self.prediction_configuration = prediction_configuration
        self.app_project = app_project
        self.signals = WorkerSignals()

    def run(self):
        predictions: list[tuple[str, str]] = gui_utils.get_prediction_name_and_seq_from_table(self.table)
        loggers.prediction_worker.debug(generic_messages.get_basic_variable_value(__file__, "run", "predictions", predictions))
        structure_prediction_obj = structure_prediction.StructurePrediction(predictions, self.prediction_configuration,
                                                                            self.app_project)
        structure_prediction_obj.create_tmp_directories()
        structure_prediction_obj.create_fasta_files_for_prediction()
        structure_prediction_obj.run_prediction()
        structure_prediction_obj.move_best_prediction_models()
        subprocess.run(["wsl", "--shutdown"])
        self.signals.finished.emit()


class AnalysisWorkerPool(QtCore.QRunnable):
    def __init__(self,
                 list_analysis_overview,
                 cb_analysis_images,
                 status_bar,
                 app_project,
                 app_settings,
                 _init_batch_analysis_page):
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
        """
        super(AnalysisWorkerPool, self).__init__()
        self.list_analysis_overview = list_analysis_overview
        self.cb_analysis_images = cb_analysis_images
        self.status_bar = status_bar
        self.app_project = app_project
        self.app_settings = app_settings
        self._init_batch_analysis_page = _init_batch_analysis_page
        self.signals = WorkerSignals()

    def get_analysis_runs(self):
        loggers.analysis_worker.debug(f"The proteins of the project loaded are: {self.app_project.proteins}")
        batch_analysis = []
        loggers.analysis_worker.debug(f"list_analysis_overview: {self.list_analysis_overview.count()}")
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

            messages = generic_messages.get_variables_values("get_analysis_runs", [
                ("prot_1_name", prot_1_name),
                ("prot_2_name", prot_2_name),
                ("prot_1_chains", prot_1_chains),
                ("prot_2_chains", prot_2_chains),
            ])
            loggers.log_multiple_messages(loggers.analysis_worker, logging.DEBUG, messages)
            tmp_prot_1 = protein_analysis_info.ProteinAnalysisInfo(prot_1_name, prot_1_chains, tmp_batch_analysis)
            tmp_prot_2 = protein_analysis_info.ProteinAnalysisInfo(prot_2_name, prot_2_chains, tmp_batch_analysis)
            batch_analysis.append((tmp_prot_1, tmp_prot_2))

        transformer = data_transformer.DataTransformer(self.app_project, batch_analysis)
        # contains analysis-ready data format: list(tuple(prot_1, prot_2, export_dir, name), ...)
        return transformer.transform_data_for_analysis()

    def run_analysis_process(self, batch_data):
        loggers.analysis_worker.debug(f"Step: run_analysis_process; batch_data: {batch_data}")
        for analysis_data in batch_data:
            if not os.path.exists(analysis_data[2]):
                os.mkdir(analysis_data[2])
            else:
                basic_boxes.ok("Structure Analysis", f"The structure analysis: {analysis_data[3]} already exists!",
                               QtWidgets.QMessageBox.Critical)
                loggers.analysis_worker.warning(f"The structure analysis: {analysis_data[3]} already exists!")
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
        self.run_analysis_process(self.get_analysis_runs())
        self.signals.finished.emit()


class PredictionWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal()
    progress = QtCore.pyqtSignal(int)

    def __init__(self,
                 table_prot_to_predict,
                 prediction_configuration,
                 app_project):
        super().__init__()
        self.table = table_prot_to_predict
        self.prediction_configuration = prediction_configuration
        self.app_project = app_project

    def update_attributes(self,
                          table_prot_to_predict,
                          prediction_configuration,
                          app_project):
        self.table = table_prot_to_predict
        self.prediction_configuration = prediction_configuration
        self.app_project = app_project
        messages = generic_messages.get_variables_values("update_attributes", [
            ("table", self.table),
            ("prediction_configuration", self.prediction_configuration),
            ("app_project", self.app_project),
        ])
        loggers.log_multiple_messages(loggers.analysis_worker, logging.DEBUG, messages)

    def run(self):
        predictions: list[tuple[str, str]] = gui_utils.get_prediction_name_and_seq_from_table(self.table)
        structure_prediction_obj = structure_prediction.StructurePrediction(predictions, self.prediction_configuration, self.app_project)
        structure_prediction_obj.create_tmp_directories()
        structure_prediction_obj.create_fasta_files_for_prediction()
        structure_prediction_obj.run_prediction()
        structure_prediction_obj.move_best_prediction_models()
        self.finished.emit()


class AnalysisWorker(QtCore.QObject):
    finished = QtCore.pyqtSignal()
    progress = QtCore.pyqtSignal(int)

    def __init__(self,
                 list_analysis_overview,
                 cb_analysis_images,
                 status_bar,
                 app_project,
                 app_settings,
                 _init_batch_analysis_page):
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
        """
        super().__init__()
        self.list_analysis_overview = list_analysis_overview
        self.cb_analysis_images = cb_analysis_images
        self.status_bar = status_bar
        self.app_project = app_project
        self.app_settings = app_settings
        self._init_batch_analysis_page = _init_batch_analysis_page

    def update_attributes(self,
                          list_analysis_overview,
                          cb_analysis_images,
                          status_bar,
                          app_project,
                          app_settings):
        self.list_analysis_overview = list_analysis_overview
        self.cb_analysis_images = cb_analysis_images
        self.status_bar = status_bar
        self.app_project = app_project
        self.app_settings = app_settings
        messages = generic_messages.get_variables_values("update_attributes", [
            ("list_analysis_overview", self.list_analysis_overview),
            ("cb_analysis_images", self.cb_analysis_images),
            ("status_bar", self.status_bar),
            ("app_project", self.app_project),
            ("app_settings", self.app_settings),
        ])
        loggers.log_multiple_messages(loggers.analysis_worker, logging.DEBUG, messages)

    def get_analysis_runs(self):
        loggers.analysis_worker.debug(f"The proteins of the project loaded are: {self.app_project.proteins}")
        batch_analysis = []
        loggers.analysis_worker.debug(f"list_analysis_overview: {self.list_analysis_overview.count()}")
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

            messages = generic_messages.get_variables_values("get_analysis_runs", [
                ("prot_1_name", prot_1_name),
                ("prot_2_name", prot_2_name),
                ("prot_1_chains", prot_1_chains),
                ("prot_2_chains", prot_2_chains),
            ])
            loggers.log_multiple_messages(loggers.analysis_worker, logging.DEBUG, messages)
            tmp_prot_1 = protein_analysis_info.ProteinAnalysisInfo(prot_1_name, prot_1_chains, tmp_batch_analysis)
            tmp_prot_2 = protein_analysis_info.ProteinAnalysisInfo(prot_2_name, prot_2_chains, tmp_batch_analysis)
            batch_analysis.append((tmp_prot_1, tmp_prot_2))

        transformer = data_transformer.DataTransformer(self.app_project, batch_analysis)
        # contains analysis-ready data format: list(tuple(prot_1, prot_2, export_dir, name), ...)
        return transformer.transform_data_for_analysis()

    def run_analysis_process(self, batch_data):
        loggers.analysis_worker.debug(f"Step: run_analysis_process; batch_data: {batch_data}")
        for analysis_data in batch_data:
            if not os.path.exists(analysis_data[2]):
                os.mkdir(analysis_data[2])
            else:
                basic_boxes.ok("Structure Analysis", f"The structure analysis: {analysis_data[3]} already exists!",
                               QtWidgets.QMessageBox.Critical)
                loggers.analysis_worker.warning(f"The structure analysis: {analysis_data[3]} already exists!")
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
        self.run_analysis_process(self.get_analysis_runs())
        self.finished.emit()
