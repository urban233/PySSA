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
import shutil
import subprocess
import logging
import pathlib


from pyssa.external_modules.mega import mega
from PyQt5.QtWidgets import QMessageBox
from pymol import cmd
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.internal.data_structures.data_classes import prediction_protein_info
from pyssa.io_pyssa import path_util
from pyssa.internal.data_structures import structure_prediction
from pyssa.internal.data_structures import structure_analysis
from pyssa.internal.data_processing import data_transformer
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.io_pyssa import safeguard
from pyssa.util import prediction_util
from pyssa.logging_pyssa import log_handlers
from typing import TYPE_CHECKING
from pyssa.util import constants

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
        logger.debug(distance_analysis_runs)
        return distance_analysis_runs

    def set_up_analysis_runs(self) -> 'structure_analysis.Analysis':
        """This function creates protein pairs and distance analysis objects for the analysis runs.

        """
        analysis_runs = structure_analysis.Analysis(self.app_project)
        for tmp_distance_analysis in self.transform_gui_input_to_practical_data():
            analysis_runs.analysis_list.append(tmp_distance_analysis)
        return analysis_runs

    def run_analysis(self) -> None:
        self.set_up_analysis_runs().run_analysis(self.cb_analysis_images)

    def run(self):
        """This function is a reimplementation of the QRunnable run method.

        """
        logger.debug(f"Memory address of worker {self}")
        # do the analysis runs
        self.run_analysis()
        # emit finish signal
        self.signals.finished.emit()

    def __del__(self):
        logger.debug("Worker has been destroyed!")


class BatchImageWorkerPool(QtCore.QRunnable):
    """This class is a worker class for the analysis process.

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    """

    # <editor-fold desc="Class attributes">
    """
    the list where all analysis runs are stored
    """
    list_analysis_images: QtWidgets.QListWidget
    """
    the list where all analysis runs are stored for which images should be created
    """
    list_analysis_for_image_creation_overview: QtWidgets.QListWidget
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
                 list_analysis_images: QtWidgets.QListWidget,
                 list_analysis_for_image_creation_overview: QtWidgets.QListWidget,
                 status_bar: QtWidgets.QStatusBar,
                 app_project: 'project.Project') -> None:
        """Constructor

        Args:
            list_analysis_images:
                list where all analysis runs are stored
            status_bar:
                status bar object of the main window
            app_project:
                project of the main window
            app_settings:
                settings of the main window

        Raises:
            ValueError: raised if an argument is illegal
        """
        super(BatchImageWorkerPool, self).__init__()
        # <editor-fold desc="Checks">
        if not safeguard.Safeguard.check_if_value_is_not_none(status_bar):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(app_project):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")

        # </editor-fold>

        self.list_analysis_images = list_analysis_images
        self.list_analysis_for_image_creation_overview = list_analysis_for_image_creation_overview
        self.status_bar = status_bar
        self.app_project = app_project

    def run(self):
        for i in range(self.list_analysis_for_image_creation_overview.count()):
            tmp_protein_pair = self.app_project.search_protein_pair(self.list_analysis_for_image_creation_overview.item(i).text())
            cmd.reinitialize()
            tmp_protein_pair.load_pymol_session()
            if not os.path.exists(constants.SCRATCH_DIR_IMAGES):
                os.mkdir(constants.SCRATCH_DIR_IMAGES)
            if not os.path.exists(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR):
                os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR)
            if not os.path.exists(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)
            tmp_protein_pair.distance_analysis.take_image_of_protein_pair(
                filename=f"structure_aln_{tmp_protein_pair.name}",
                representation="cartoon")
            tmp_protein_pair.distance_analysis.analysis_results.set_structure_aln_image(
                path_util.FilePath(
                    pathlib.Path(
                        f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{tmp_protein_pair.name}.png")
                )
            )
            logger.debug(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image[0])
            tmp_protein_pair.distance_analysis.take_image_of_interesting_regions(
                tmp_protein_pair.distance_analysis.cutoff,
                f"interesting_reg_{tmp_protein_pair.name}")
            interesting_region_filepaths = []
            for tmp_filename in os.listdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                interesting_region_filepaths.append(
                    path_util.FilePath(
                        pathlib.Path(
                            f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{tmp_filename}")
                    )
                )
            tmp_protein_pair.distance_analysis.analysis_results.set_interesting_region_images(
                interesting_region_filepaths)
            shutil.rmtree(constants.SCRATCH_DIR_IMAGES)
            # emit finish signal
            self.signals.finished.emit()


class ColabfoldInstallerWorkerPool(QtCore.QRunnable):
    """This class is a worker class for the analysis process.

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    """

    # <editor-fold desc="Class attributes">
    """
    the list where all analysis runs are stored
    """
    list_analysis_images: QtWidgets.QListWidget
    """
    the list where all analysis runs are stored for which images should be created
    """
    list_analysis_for_image_creation_overview: QtWidgets.QListWidget
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

    def __init__(self, install) -> None:
        """Constructor

        Args:
            install:
                set true if Colabfold should be installed

        Raises:
            ValueError: raised if an argument is illegal
        """
        super(ColabfoldInstallerWorkerPool, self).__init__()
        self.install = install
        self.local_install: tuple[bool, str] = False, ""

    def run(self):
        if self.install is False:
            # if localcolab should be uninstalled
            try:
                subprocess.run([constants.POWERSHELL_EXE, constants.REMOVE_WSL_POWERSHELL])
            except:
                basic_boxes.ok("Local Colabfold removal",
                               "The uninstallation failed. Please re-run the process or consult the documentation.",
                               QMessageBox.Critical)
                return
        else:
            # logical message: the user wants to install local colabfold
            if not os.path.exists(pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/wsl/")):
                os.mkdir(pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/wsl/"))
            if not os.path.exists(constants.WSL_STORAGE_PATH):
                os.mkdir(constants.WSL_STORAGE_PATH)

            def download_file_with_curl(destination):
                mega_handler = mega.Mega()
                # Download the file
                user = mega_handler.login()
                # Download the file
                user.download_url(url=constants.DISTRO_DOWNLOAD_URL, dest_path=destination)

            if self.local_install[0] is False:
                download_file_with_curl(str(constants.WSL_DISTRO_IMPORT_PATH))
            else:
                shutil.move(self.local_install[1], str(constants.WSL_DISTRO_IMPORT_PATH))
            subprocess.run([constants.POWERSHELL_EXE, constants.CONVERT_DOS_TO_UNIX])
            subprocess.run(["wsl", "--import", constants.WSL_DISTRO_NAME, str(constants.WSL_STORAGE_PATH),
                            str(constants.WSL_DISTRO_IMPORT_PATH)])
            subprocess.run(["wsl", "--set-default", constants.WSL_DISTRO_NAME])
            if not os.path.exists(constants.WSL_DISK_PATH):
                basic_boxes.ok("Colabfold installation", "Installation failed, please try again.", QMessageBox.Critical)

        # emit finish signal
        self.signals.finished.emit()

# class ResultsWorkerPool(QtCore.QRunnable):
#
#     """
#     the signals to use, for the worker
#     """
#     signals = WorkerSignals()
#
#     def __init__(self, app_project, ui, show_analysis_results_options, show_results_interactions):
#         """Constructor
#
#                 Args:
#
#
#                 Raises:
#                     ValueError: raised if an argument is illegal
#                 """
#         super(ResultsWorkerPool, self).__init__()
#         self.app_project = app_project
#         self.ui = ui
#         self.show_analysis_results_options = show_analysis_results_options
#         self.show_results_interactions = show_results_interactions
#
#     def load_results(self):
#         shutil.rmtree(constants.CACHE_DIR)
#         os.mkdir(constants.CACHE_DIR)
#         os.mkdir(constants.CACHE_IMAGES)
#         self.results_name = self.ui.cb_results_analysis_options.currentText()
#
#         if self.results_name == "":
#             self.show_analysis_results_options()
#             return
#         gui_elements_to_hide = []
#         # TODO: implement image check for xml format
#         # if not os.path.exists(pathlib.Path(f"{current_results_path}/images")):
#         #     # no images where made
#         #     gui_elements_to_hide.append(self.ui.lbl_results_structure_alignment)
#         #     gui_elements_to_hide.append(self.ui.btn_view_struct_alignment)
#         #     gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
#         #     gui_elements_to_hide.append(self.ui.list_results_interest_regions)
#         #     gui_elements_to_hide.append(self.ui.btn_view_interesting_region)
#         # elif os.path.exists(pathlib.Path(f"{current_results_path}/images")):
#         #     if not os.path.exists(pathlib.Path(f"{current_results_path}/images/structure_alignment.png")):
#         #         gui_elements_to_hide.append(self.ui.lbl_results_structure_alignment)
#         #         gui_elements_to_hide.append(self.ui.btn_view_struct_alignment)
#         #     elif not os.path.exists(pathlib.Path(f"{current_results_path}/images/interesting_")):
#         #         gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
#         #         gui_elements_to_hide.append(self.ui.list_results_interest_regions)
#         #         gui_elements_to_hide.append(self.ui.btn_view_interesting_region)
#
#         tmp_protein_pair = self.app_project.search_protein_pair(self.results_name)
#         distance_data: dict[str, np.ndarray] = tmp_protein_pair.distance_analysis.analysis_results.distance_data
#         distance_data_array = np.array([distance_data["index"], distance_data["ref_chain"], distance_data["ref_pos"],
#                                         distance_data["ref_resi"], distance_data["model_chain"],
#                                         distance_data["model_pos"],
#                                         distance_data["model_resi"], distance_data["distance"]])
#
#         filesystem_io.XmlDeserializer(self.app_project.get_project_xml_path()).deserialize_analysis_images(
#             tmp_protein_pair.name, tmp_protein_pair.distance_analysis.analysis_results)
#         tmp_protein_pair.distance_analysis.analysis_results.create_image_png_files_from_base64()
#
#         distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
#
#         # check if histogram can be created
#         distance_list.sort()
#         x, y = np.histogram(distance_list, bins=np.arange(0, distance_list[len(distance_list) - 1], 0.25))
#         if x.size != y.size:
#             x = np.resize(x, (1, y.size))
#         # this conversion is needed for the pyqtgraph library!
#         x = x.tolist()
#         try:
#             x = x[0]
#         except IndexError:
#             # histogram could not be created
#             gui_elements_to_hide.append(self.ui.lbl_results_distance_histogram)
#             gui_elements_to_hide.append(self.ui.btn_view_distance_histogram)
#
#         if gui_elements_to_hide:
#             self.show_results_interactions(gui_elements_to_hide=gui_elements_to_hide)
#         else:
#             self.show_results_interactions()
#
#         self.ui.list_results_interest_regions.clear()
#         for tmp_filename in os.listdir(constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
#             self.ui.list_results_interest_regions.addItem(tmp_filename)
#         self.ui.list_results_interest_regions.sortItems()
#         self.ui.txt_results_rmsd.setText(str(tmp_protein_pair.distance_analysis.analysis_results.rmsd))
#         self.ui.txt_results_aligned_residues.setText(
#             str(tmp_protein_pair.distance_analysis.analysis_results.aligned_aa))
#         cmd.reinitialize()
#         tmp_protein_pair.load_pymol_session()
#
#     def run(self):
#         self.load_results()
#         self.signals.finished.emit()