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
import copy
import logging
import os
import pathlib
import shutil
import subprocess

from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5 import QtWidgets
from pymol import cmd
from pyssa.internal.data_processing import data_transformer
from pyssa.internal.data_structures import protein, protein_pair, structure_prediction, structure_analysis
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.io_pyssa import filesystem_io, safeguard, path_util
from pyssa.io_pyssa.xml_pyssa import element_names, attribute_names
from pyssa.logging_pyssa import log_handlers
from pyssa.util import tools, constants, prediction_util, exception, exit_codes
from pyssa.internal.prediction_engines import esmfold
from pyssa.internal.data_structures import project
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import settings

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def setup_worker_for_work(tmp_thread, tmp_worker, return_value_func):
    tmp_worker.moveToThread(tmp_thread)
    tmp_thread.started.connect(tmp_worker.run)
    tmp_worker.finished.connect(tmp_thread.quit)
    tmp_worker.finished.connect(tmp_worker.deleteLater)
    tmp_thread.finished.connect(tmp_thread.deleteLater)
    tmp_worker.return_value.connect(return_value_func)
    return tmp_thread


class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(tuple)
    workspace_path: str

    def __init__(self, workspace):
        super().__init__()
        self.workspace_path = workspace

    def run(self):
        protein_dict, protein_names = tools.scan_workspace_for_non_duplicate_proteins(pathlib.Path(self.workspace_path))
        self.return_value.emit((protein_dict, protein_names))
        self.finished.emit()


class CreateUseProjectWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(list)
    workspace_path: str
    proteins_to_copy: list

    def __init__(self, workspace, proteins_to_copy):
        super().__init__()
        self.workspace_path = workspace
        self.proteins_to_copy = proteins_to_copy

    def run(self):
        proteins_for_new_project = []
        protein_infos = (tools.scan_workspace_for_non_duplicate_proteins(pathlib.Path(self.workspace_path)))[0]
        for tmp_protein in self.proteins_to_copy:
            for tmp_protein_info in protein_infos:
                if tmp_protein_info.name == tmp_protein:
                    """Var: project_proteins is a list which contains all proteins from a single project"""
                    xml_deserializer = filesystem_io.XmlDeserializer(
                        pathlib.Path(f"{self.workspace_path}/{tmp_protein_info.project_name}.xml"))
                    for xml_protein in xml_deserializer.xml_root.iter(element_names.PROTEIN):
                        if xml_protein.attrib[attribute_names.ID] == tmp_protein_info.id:
                            basic_information = xml_protein.attrib
                            pdb_lines = []
                            session_data_base64 = ""
                            for tmp_data in xml_protein:
                                if tmp_data.tag == "pdb_data":
                                    for tmp_atom in tmp_data.findall("atom"):
                                        pdb_lines.append(tmp_atom.text)
                                elif tmp_data.tag == "session_data":
                                    session_data_base64 = tmp_data.attrib[attribute_names.PROTEIN_SESSION]
                                else:
                                    raise ValueError
                            tmp_protein_obj = protein.Protein(
                                molecule_object=basic_information[attribute_names.PROTEIN_MOLECULE_OBJECT],
                                pdb_xml_string=xml_protein)
                            tmp_protein_obj.set_all_attributes(basic_information, pdb_lines, session_data_base64)
            proteins_for_new_project.append(tmp_protein_obj)

        self.return_value.emit(proteins_for_new_project)
        self.finished.emit()


class LoadResultsWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(str)
    protein_pair_of_results: protein_pair.ProteinPair
    app_project_xml_filepath: str
    image_type: str

    def __init__(self, protein_pair_of_results, app_project_xml_filepath):
        super().__init__()
        self.protein_pair_of_results = protein_pair_of_results
        self.image_type = constants.IMAGES_NONE
        self.app_project_xml_filepath = app_project_xml_filepath

    def run(self):
        filesystem_io.XmlDeserializer(self.app_project_xml_filepath).deserialize_analysis_images(
            self.protein_pair_of_results.name, self.protein_pair_of_results.distance_analysis.analysis_results)
        if len(self.protein_pair_of_results.distance_analysis.analysis_results.structure_aln_image) != 0 and len(
                self.protein_pair_of_results.distance_analysis.analysis_results.interesting_regions_images) != 0:
            # if both image types were made during analysis
            self.protein_pair_of_results.distance_analysis.analysis_results.create_image_png_files_from_base64()
            self.image_type = constants.IMAGES_ALL
        elif len(self.protein_pair_of_results.distance_analysis.analysis_results.structure_aln_image) != 0 and len(
                self.protein_pair_of_results.distance_analysis.analysis_results.interesting_regions_images) == 0:
            # only struct align image were made
            self.protein_pair_of_results.distance_analysis.analysis_results.create_image_png_files_from_base64()
            self.image_type = constants.IMAGES_STRUCT_ALN_ONLY
        else:
            # no images were made
            self.image_type = constants.IMAGES_NONE

        self.return_value.emit(self.image_type)
        self.finished.emit()


class BatchImageWorker(QObject):
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
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(str)

    # </editor-fold>

    def __init__(self,
                 list_analysis_images: QtWidgets.QListWidget,
                 list_analysis_for_image_creation_overview: QtWidgets.QListWidget,
                 status_bar: QtWidgets.QStatusBar,
                 app_project: 'project.Project') -> None:
        """Constructor.

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
        super().__init__()
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
                        f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{tmp_protein_pair.name}.png"),
                ),
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
                            f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{tmp_filename}"),
                    ),
                )
            tmp_protein_pair.distance_analysis.analysis_results.set_interesting_region_images(
                interesting_region_filepaths)
            shutil.rmtree(constants.SCRATCH_DIR_IMAGES)
            # emit finish signal
            self.finished.emit()


class EsmFoldWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(list)
    table: QtWidgets.QTableWidget

    def __init__(self, table_prot_to_predict: QtWidgets.QTableWidget):
        super().__init__()
        self.table = table_prot_to_predict

    def run(self):
        predictions: list[prediction_protein_info.PredictionProteinInfo] = prediction_util.get_prediction_name_and_seq_from_table(
            self.table)
        output = esmfold.EsmFold(predictions).run_prediction()
        self.return_value.emit(output)
        self.finished.emit()


class ColabfoldWorker(QObject):
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
    finished = pyqtSignal(int, str)
    progress = pyqtSignal(int)
    return_value = pyqtSignal(list)
    """
    a signal which gets used if a process failed due to an error
    """
    non_zero_exit_code = pyqtSignal()

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
        super().__init__()
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

    def run(self) -> None:
        """This function is a reimplementation of the QObject run method. It does the structure prediction."""
        predictions: list[
            prediction_protein_info.PredictionProteinInfo
        ] = prediction_util.get_prediction_name_and_seq_from_table(self.table)
        structure_prediction_obj = structure_prediction.StructurePrediction(predictions,
                                                                            self.prediction_configuration,
                                                                            self.app_project)
        structure_prediction_obj.create_tmp_directories()
        logger.info("Tmp directories were created.")

        # Create fasta files for prediction
        try:
            structure_prediction_obj.create_fasta_files_for_prediction()
        except exception.FastaFilesNotCreatedError:
            logger.error("Fasta files were not created.")
            self.finished.emit(exit_codes.ERROR_WRITING_FASTA_FILES[0], exit_codes.ERROR_WRITING_FASTA_FILES[1])
            return
        except exception.FastaFilesNotFoundError:
            logger.error("Fasta files were not found.")
            self.finished.emit(exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0], exit_codes.ERROR_FASTA_FILES_NOT_FOUND[1])
            return
        except Exception as e:
            logger.error("Unexpected error:", e)
            self.finished.emit(exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
            return
        else:
            logger.info("Fasta files were successfully created.")

        # Run structure prediction
        try:
            structure_prediction_obj.run_prediction()
        except exception.PredictionEndedWithError:
            logger.error("Prediction ended with error.")
            self.finished.emit(exit_codes.ERROR_PREDICTION_FAILED[0], exit_codes.ERROR_PREDICTION_FAILED[1])
            return
        else:
            logger.info("Prediction process finished.")

        try:
            structure_prediction_obj.move_best_prediction_models()
            logger.info("Saved predicted pdb file into XML file.")
        except exception.UnableToFindColabfoldModelError:
            logger.error("Could not move rank 1 model, because it does not exists.")
            self.finished.emit(exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[0], exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[1])
            return
        except FileNotFoundError:
            logger.error("Could not move rank 1 model, because it does not exists.")
            self.finished.emit(exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[0], exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[1])
            return
        except Exception as e:
            logger.error("Unexpected error:", e)
            logger.error("Could not move rank 1 model, because it does not exists.")
            self.finished.emit(exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
            return
        subprocess.run(["wsl", "--shutdown"])
        logger.info("WSL gets shutdown.")
        self.finished.emit(exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])


class DistanceAnalysisWorker(QObject):
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
    finished = pyqtSignal(int, str)
    progress = pyqtSignal(int)
    return_value = pyqtSignal(list)

    # </editor-fold>

    def __init__(self,
                 list_analysis_overview: QtWidgets.QListWidget,
                 cb_analysis_images: QtWidgets.QCheckBox,
                 status_bar: QtWidgets.QStatusBar,
                 app_project: 'project.Project',
                 app_settings: 'settings.Settings',
                 _init_batch_analysis_page) -> None:
        """Constructor.

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
        super().__init__()
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
        """Transforms the input from the gui to a practical data basis that can be used to set up analysis runs.

        Raises:
            UnableToTransformDataForAnalysis: If the transformation process failed.
        """
        distance_analysis_runs = []
        logger.debug(f"Count of distance analysis runs: {self.list_analysis_overview.count()}")

        try:
            for row_no in range(self.list_analysis_overview.count()):
                input_transformer = data_transformer.DistanceAnalysisDataTransformer(
                    self.list_analysis_overview.item(row_no).text(),
                    self.app_project,
                    self.app_settings,
                )
                protein_pair_for_analysis = input_transformer.transform_gui_input_to_distance_analysis_object()
                new_protein_pair = copy.deepcopy(protein_pair_for_analysis)
                distance_analysis_runs.append(new_protein_pair)
            logger.debug(f"These are the distance analysis runs, after the data transformation: {distance_analysis_runs}.")
        except exception.IllegalArgumentError:
            logger.error("Transformation of data failed because an argument was illegal.")
            raise exception.UnableToTransformDataForAnalysisError("")
        except exception.UnableToTransformDataForAnalysisError:
            logger.error("Transformation of data failed because the transformation process failed.")
            raise exception.UnableToTransformDataForAnalysisError("")
        except Exception as e:
            logger.error(f"Unknown error: {e}")
            raise exception.UnableToTransformDataForAnalysisError("")
        return distance_analysis_runs

    def set_up_analysis_runs(self) -> 'structure_analysis.Analysis':
        """This function creates protein pairs and distance analysis objects for the analysis runs.

        Raises:
            UnableToSetupAnalysisError: If the analysis runs setup failed.
        """
        try:
            analysis_runs = structure_analysis.Analysis(self.app_project)
            analysis_runs.analysis_list = self.transform_gui_input_to_practical_data()
            logger.debug(analysis_runs.analysis_list[
                             0].distance_analysis.get_protein_pair().protein_1.pymol_selection.selection_string)
        except exception.UnableToTransformDataForAnalysisError:
            logger.error("Setting up the analysis runs failed.")
            raise exception.UnableToSetupAnalysisError("")
        except Exception as e:
            logger.error(f"Unknown error: {e}")
            raise exception.UnableToSetupAnalysisError("")
        else:
            return analysis_runs

    def run(self) -> None:
        """This function is a reimplementation of the QRunnable run method."""
        logger.debug(f"Memory address of worker {self}")
        try:
            self.set_up_analysis_runs().run_analysis(self.cb_analysis_images)
        except exception.UnableToSetupAnalysisError:
            logger.error("Setting up the analysis runs failed therefore the distance analysis failed.")
            self.finished.emit(exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[0], exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[1])
        except Exception as e:
            logger.error(f"Unknown error: {e}")
            self.finished.emit(exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
        else:
            self.finished.emit(exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])


class PreviewRayImageWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(tuple)
    renderer: str

    def __init__(self, renderer):
        super().__init__()
        self.renderer = renderer

    def run(self):
        cmd.ray(2400, 2400, renderer=int(self.renderer))
        self.finished.emit()


class SaveRayImageWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(tuple)
    renderer: str
    full_file_name: str

    def __init__(self, renderer, full_file_name):
        super().__init__()
        self.renderer = renderer
        self.full_file_name = full_file_name

    def run(self):
        cmd.ray(2400, 2400, renderer=int(self.renderer))
        cmd.png(self.full_file_name, dpi=300)
        self.finished.emit()


class OpenProjectWorker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(project.Project)
    workspace_path: pathlib.Path
    project_name: str
    app_settings: 'settings.Settings'

    def __init__(self, workspace_path, project_name, app_settings) -> None:
        """Constructor."""
        super().__init__()
        self.workspace_path = workspace_path
        self.project_name = project_name
        self.app_settings = app_settings

    def run(self) -> None:
        """Logic of the open project process."""
        # show project management options in side menu
        tmp_project_path = pathlib.Path(f"{self.workspace_path}/{self.project_name}")
        tmp_app_project = project.Project.deserialize_project(tmp_project_path, self.app_settings)
        self.return_value.emit(tmp_app_project)
        self.finished.emit()
