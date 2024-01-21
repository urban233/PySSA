import logging
import os
import pathlib
import shutil
import subprocess
import platform

import pymol

from pymol import cmd
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from Bio import SeqRecord
from Bio import SeqIO
from xml import sax

from pyssa.controller import results_view_controller
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import predict_monomer_view, delete_project_view
from pyssa.gui.ui.dialogs import dialog_startup, dialog_settings_global, dialog_tutorial_videos, dialog_about
from pyssa.internal.data_structures import project, settings, protein
from pyssa.internal.data_structures.data_classes import prediction_protein_info
from pyssa.internal.portal import graphic_operations, pymol_io
from pyssa.internal.thread import tasks, task_workers
from pyssa.io_pyssa import safeguard, filesystem_io, path_util
from pyssa.logging_pyssa import log_handlers
from pyssa.presenter import main_presenter_async
from pyssa.util import constants, enums, exception, main_window_util, exit_codes, prediction_util, gui_utils, tools
from pyssa.gui.ui.views import main_view
from pyssa.gui.ui.views import open_project_view
from pyssa.model import application_model
from pyssa.controller import interface_manager, distance_analysis_view_controller, predict_monomer_view_controller, \
    delete_project_view_controller, create_project_view_controller, open_project_view_controller, database_manager
from pyssa.util import globals

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class MainViewController:
    """Class for main presenter of the pyssa plugin."""

    """
    The main view of the pyssa plugin.
    """
    _view: "main_view.MainView"

    """
    
    """
    _app_model: "application_model.ApplicationModel"

    """
    The path to the active workspace.
    """
    _workspace_path: pathlib.Path

    """
    The application settings object.
    """
    _application_settings: "settings.Settings"

    """
    A watcher for the project state.
    """
    _interface_manager: "interface_manager.InterfaceManager"

    """
    The active task of the application.
    """
    _active_task: tasks.Task

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        """Constructor.

        Args:
            a_view: the main view of the pyssa plugin.

        Raises:
            IllegalArgumentException: If an argument is illegal.
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(the_interface_manager, logger)

        # </editor-fold>

        self._view: "main_view.MainView" = the_interface_manager.get_main_view()
        self._interface_manager: "interface_manager.InterfaceManager" = the_interface_manager
        self._database_manager = database_manager.DatabaseManager("")
        self._app_model: "application_model.ApplicationModel" = application_model.ApplicationModel(project.Project())
        self._external_view = None
        self._protein_model = QtGui.QStandardItemModel()
        self._workspace_path = constants.DEFAULT_WORKSPACE_PATH
        self._workspace_status = f"Current workspace: {str(self._workspace_path)}"
        self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")

        self._setup_statusbar()
        self._connect_all_ui_elements_with_slot_functions()

    def _setup_statusbar(self) -> None:
        """Sets up the status bar and fills it with the current workspace."""
        self._interface_manager.get_main_view().setStatusBar(self._interface_manager.get_main_view().status_bar)

    def _connect_all_ui_elements_with_slot_functions(self):
        self._view.ui.action_new_project.triggered.connect(self._create_project)
        self._view.ui.project_tab_widget.currentChanged.connect(self._update_tab)
        self._view.ui.action_open_project.triggered.connect(self._open_project)
        self._view.ui.action_delete_project.triggered.connect(self._delete_project)
        self._view.ui.actionImport.triggered.connect(self.import_project)
        self._view.ui.actionExport.triggered.connect(self.export_current_project)
        self._view.ui.action_close_project.triggered.connect(self._close_project)

        self._view.ui.action_results_summary.triggered.connect(self._results_summary)
        self._view.ui.actionPreview.triggered.connect(self.preview_image)

        self._view.ui.action_edit_settings.triggered.connect(self.open_settings_global)
        self._view.ui.action_restore_settings.triggered.connect(self.restore_settings)
        self._view.ui.actionShow_log_in_explorer.triggered.connect(self.open_logs)
        self._view.ui.actionClear_Logs.triggered.connect(self.clear_all_log_files)
        self._view.ui.actionDocumentation.triggered.connect(self.open_documentation)
        self._view.ui.actionTutorials.triggered.connect(self.open_tutorial)
        self._view.ui.action_about.triggered.connect(self.open_about)
        self._view.ui.action_predict_monomer.triggered.connect(self._predict_monomer)
        self._view.ui.action_distance_analysis.triggered.connect(self._distance_analysis)

        # proteins tab
        self._view.ui.proteins_tree_view.clicked.connect(self._show_protein_information)
        self._view.ui.btn_create_protein_scene.clicked.connect(self.save_scene)
        self._view.ui.btn_update_protein_scene.clicked.connect(self.update_scene)
        self._view.cb_chain_color.currentIndexChanged.connect(self._change_chain_color_proteins)
        self._view.cb_chain_representation.currentIndexChanged.connect(self._change_chain_representation_proteins)
        self._view.ui.btn_import_protein.clicked.connect(self._import_protein_structure)
        # seqs tab
        self._view.ui.seqs_list_view.clicked.connect(self._show_sequence_information)
        self._view.ui.pushButton.clicked.connect(self._add_sequence)
        self._view.ui.btn_import_seq.clicked.connect(self._import_sequence)

        # protein pairs tab
        self._view.ui.protein_pairs_tree_view.clicked.connect(self._show_protein_information_of_protein_pair)
        self._view.ui.btn_create_protein_pair_scene.clicked.connect(self.save_scene)
        self._view.ui.btn_update_protein_pair_scene.clicked.connect(self.update_scene)
        self._view.ui.protein_pairs_tree_view.clicked.connect(self._check_for_results)

    def _update_tab(self) -> None:
        # if self._view.ui.project_tab_widget.currentIndex() == 1:
        #     print("Changed proteins")
        #     self._interface_manager.refresh_protein_model()
        # elif self._view.ui.project_tab_widget.currentIndex() == 2:
        #     print("Changed protein pairs")
        #     self._interface_manager.refresh_protein_pair_model()
        pass

    def _close_project(self):
        """Closes the current project"""
        self._interface_manager.set_new_project(project.Project())
        self._interface_manager.refresh_main_view()

    def _create_project(self) -> None:
        self._external_controller = create_project_view_controller.CreateProjectViewController(self._interface_manager)
        self._external_controller.user_input.connect(self._post_create_project)
        self._interface_manager.get_create_view().show()

    def _post_create_project(self, user_input: tuple) -> None:
        tmp_project_name, tmp_protein_name = user_input
        self._database_manager.set_database_filepath(
            str(pathlib.Path(f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}")))
        self._database_manager.build_new_database()
        self._database_manager.open_project_database()
        tmp_project = project.Project(tmp_project_name,
                                      self._interface_manager.get_application_settings().workspace_path)
        tmp_project.set_id(self._database_manager.insert_new_project(tmp_project.get_project_name(),
                                                                     platform.system()))
        if len(tmp_protein_name) == 4:
            tmp_ref_protein = protein.Protein(tmp_protein_name.upper())
            tmp_ref_protein.db_project_id = tmp_project.get_id()
            tmp_ref_protein.add_protein_structure_data_from_pdb_db(tmp_protein_name.upper())
            tmp_ref_protein.create_new_pymol_session()
            tmp_ref_protein.save_pymol_session_as_base64_string()
            tmp_project.add_existing_protein(tmp_ref_protein)
            tmp_ref_protein.db_project_id = self._database_manager.insert_new_protein(tmp_ref_protein)
            constants.PYSSA_LOGGER.info("Create project finished with protein from the PDB.")
        elif len(tmp_protein_name) > 0:
            # local pdb file as input
            pdb_filepath = pathlib.Path(tmp_protein_name)
            graphic_operations.setup_default_session_graphic_settings()
            tmp_ref_protein = protein.Protein(
                pdb_filepath.name.replace(".pdb","")
            )
            tmp_ref_protein.db_project_id = tmp_project.get_id()
            tmp_ref_protein.add_protein_structure_data_from_local_pdb_file(pathlib.Path(tmp_protein_name))
            tmp_ref_protein.create_new_pymol_session()
            tmp_ref_protein.save_pymol_session_as_base64_string()
            tmp_project.add_existing_protein(tmp_ref_protein)
            tmp_ref_protein.db_project_id = self._database_manager.insert_new_protein(tmp_ref_protein)
            constants.PYSSA_LOGGER.info("Create project finished with protein from local filesystem.")
        else:
            constants.PYSSA_LOGGER.info("Create empty project finished.")
        self._interface_manager.set_new_project(tmp_project)
        self._interface_manager.refresh_main_view()

    def _open_project(self) -> None:
        self._external_controller = open_project_view_controller.OpenProjectViewController(self._interface_manager)
        self._external_controller.return_value.connect(self._post_open_project)
        self._interface_manager.get_open_view().show()

    def _post_open_project(self, return_value: str):
        self._interface_manager.start_wait_spinner()
        self._interface_manager.update_status_bar("Opening existing project ...")
        tmp_project_name = return_value
        tmp_project_database_filepath = str(
            pathlib.Path(
                f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}"
            )
        )
        self._database_manager.set_database_filepath(tmp_project_database_filepath)
        self._database_manager.open_project_database()
        print(self._database_manager.get_number_of_sequences())
        print(self._database_manager.get_number_of_proteins())
        print(self._database_manager.get_number_of_protein_pairs())
        tmp_project = self._database_manager.get_project_as_object(
            tmp_project_name, self._interface_manager.get_application_settings().workspace_path
        )
        self._interface_manager.set_new_project(tmp_project)
        self._interface_manager.refresh_main_view()
        self._interface_manager.stop_wait_spinner()
        self._interface_manager.update_status_bar("Opening existing project finished.")
        # tmp_filepath = pathlib.Path(f"{return_value}.xml")
        #
        # self._active_task = tasks.Task(
        #     target=main_presenter_async.open_project,
        #     args=(
        #         self._workspace_path,
        #         tmp_filepath,
        #         self._interface_manager.get_application_settings(),
        #     ),
        #     post_func=self.__await_open_project,
        # )
        # self._active_task.start()

    #
    # def __await_open_project(self, a_result: tuple) -> None:
    #     self._interface_manager.set_new_project(a_result[1])
    #     self._interface_manager.update_status_bar(self._workspace_status)
    #     self._interface_manager.refresh_main_view()
    #     self._interface_manager.stop_wait_spinner()

    def _delete_project(self) -> None:
        self._external_controller = delete_project_view_controller.DeleteProjectViewController(self._interface_manager)
        self._interface_manager.get_delete_view().show()

    def _show_protein_information(self) -> None:
        tmp_type = self._view.ui.proteins_tree_view.model().data(
            self._view.ui.proteins_tree_view.currentIndex(), enums.ModelEnum.TYPE_ROLE
        )
        if tmp_type == "protein":
            tmp_first_chain = self._view.ui.proteins_tree_view.currentIndex().child(0, 0)
            self._interface_manager.show_chain_pymol_parameters(tmp_first_chain)
        elif tmp_type == "chain":
            self._interface_manager.show_chain_pymol_parameters(self._view.ui.proteins_tree_view.currentIndex())
        else:
            raise ValueError("Unknown type!")
        self._view.status_bar.showMessage(f"Active protein structure: {tmp_type}")

    def _show_protein_information_of_protein_pair(self) -> None:
        tmp_type = self._view.ui.protein_pairs_tree_view.model().data(
            self._view.ui.protein_pairs_tree_view.currentIndex(), enums.ModelEnum.TYPE_ROLE
        )
        if tmp_type == "protein":
            tmp_first_chain = self._view.ui.protein_pairs_tree_view.currentIndex().child(0, 0)
            self._interface_manager.show_chain_pymol_parameter_for_protein_pairs(tmp_first_chain)
        elif tmp_type == "chain":
            self._interface_manager.show_chain_pymol_parameter_for_protein_pairs(self._view.ui.protein_pairs_tree_view.currentIndex())
        elif tmp_type == "protein_pair":
            pass
        else:
            print(tmp_type)
            raise ValueError("Unknown type!")
        self._view.status_bar.showMessage(f"Active protein structure: {tmp_type}")

    def _change_chain_color_proteins(self) -> None:
        if self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE) == "chain":
            tmp_protein: "protein.Protein" = self._view.ui.proteins_tree_view.currentIndex().parent().data(enums.ModelEnum.OBJECT_ROLE)
            tmp_chain = self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
            tmp_protein_chain = tmp_protein.get_chain_by_letter(tmp_chain.chain_letter)
            tmp_protein_chain.pymol_parameters["chain_color"] = self._view.cb_chain_color.currentText()
        elif self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE) == "protein":
            tmp_protein = self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
            tmp_protein.chains[0].pymol_parameters["chain_color"] = self._view.cb_chain_color.currentText()
        self._interface_manager.get_current_project().serialize_project(self._interface_manager.get_current_project().get_project_xml_path())

    def _change_chain_representation_proteins(self) -> None:
        if self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE) == "chain":
            tmp_protein: "protein.Protein" = self._view.ui.proteins_tree_view.currentIndex().parent().data(enums.ModelEnum.OBJECT_ROLE)
            tmp_chain = self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
            tmp_protein_chain = tmp_protein.get_chain_by_letter(tmp_chain.chain_letter)
            tmp_protein_chain.pymol_parameters["chain_representation"] = self._view.cb_chain_representation.currentText()
        elif self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE) == "protein":
            tmp_protein = self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
            tmp_protein.chains[0].pymol_parameters["chain_representation"] = self._view.cb_chain_representation.currentText()

    def _add_sequence(self):
        pass

    def _post_add_sequence(self, return_value: tuple):
        tmp_seq_name = return_value[0]
        tmp_sequence = return_value[1]
        tmp_seq_record = SeqRecord.SeqRecord(tmp_sequence, name=tmp_seq_name)
        self._interface_manager.get_current_project().sequences.append(tmp_seq_record)
        self._interface_manager.get_current_project().serialize_project(self._interface_manager.get_current_project().get_project_xml_path())

    def _show_sequence_information(self):
        self._interface_manager.show_sequence_parameters(
            self._view.ui.seqs_list_view.currentIndex()
        )

    def update_status(self, message: str) -> None:
        """Updates the status bar of the main view with a custom message."""
        self._view.status_bar.showMessage(message)

    def _predict_monomer(self):
        self._external_controller = predict_monomer_view_controller.PredictMonomerViewController(
            self._interface_manager
        )
        self._external_controller.job_input.connect(self._post_predict_monomer)
        self._interface_manager.get_predict_monomer_view().show()

    # <editor-fold desc="Sequences tab methods">
    def _import_sequence(self) -> None:
        self._interface_manager.get_import_sequence_view().return_value.connect(self._post_import_sequence)
        self._interface_manager.get_import_sequence_view().show()

    def _post_import_sequence(self, return_value: tuple):
        tmp_fasta_filepath, _ = return_value
        tmp_record = SeqIO.read(tmp_fasta_filepath, "fasta")
        self._interface_manager.get_current_project().sequences.append(tmp_record)
        self._interface_manager.refresh_sequence_model()
        self._interface_manager.refresh_main_view()

    # </editor-fold>

    # <editor-fold desc="Proteins tab methods">
    def _import_protein_structure(self):
        self._interface_manager.get_add_protein_view().return_value.connect(self._post_import_protein_structure)
        self._interface_manager.get_add_protein_view().show()

    def _post_import_protein_structure(self, return_value: tuple):
        tmp_protein_name, tmp_name_len = return_value
        if tmp_name_len == 4:
            tmp_ref_protein = protein.Protein(tmp_protein_name.upper())
            tmp_ref_protein.db_project_id = self._interface_manager.get_current_project().get_id()
            tmp_ref_protein.add_protein_structure_data_from_pdb_db(tmp_protein_name.upper())
            tmp_ref_protein.create_new_pymol_session()
            tmp_ref_protein.save_pymol_session_as_base64_string()
            self._interface_manager.get_current_project().add_existing_protein(tmp_ref_protein)
            tmp_ref_protein.db_project_id = self._database_manager.insert_new_protein(tmp_ref_protein)
            constants.PYSSA_LOGGER.info("Create project finished with protein from the PDB.")
        elif tmp_name_len > 0:
            # local pdb file as input
            pdb_filepath = pathlib.Path(tmp_protein_name)
            graphic_operations.setup_default_session_graphic_settings()
            tmp_ref_protein = protein.Protein(
                pdb_filepath.name.replace(".pdb","")
            )
            tmp_ref_protein.db_project_id = self._interface_manager.get_current_project().get_id()
            tmp_ref_protein.add_protein_structure_data_from_local_pdb_file(pathlib.Path(tmp_protein_name))
            tmp_ref_protein.create_new_pymol_session()
            tmp_ref_protein.save_pymol_session_as_base64_string()
            self._interface_manager.get_current_project().add_existing_protein(tmp_ref_protein)
            tmp_ref_protein.db_project_id = self._database_manager.insert_new_protein(tmp_ref_protein)
            constants.PYSSA_LOGGER.info("Create project finished with protein from local filesystem.")
        self._interface_manager.refresh_protein_model()
        self._interface_manager.refresh_main_view()

    @staticmethod
    def update_scene() -> None:
        """Updates the current selected PyMOL scene."""
        cmd.scene(key="auto", action="update")

    def save_scene(self) -> None:
        """Saves the current view as a new PyMOL scene."""
        # returns tuple with (name, bool)
        scene_name = QtWidgets.QInputDialog.getText(self._view, "Save Scene", "Enter scene name:")
        if scene_name[1]:
            cmd.scene(key=scene_name[0], action="append")

    # </editor-fold>

    # <editor-fold desc="Logic">
    def _post_predict_monomer(self, result: tuple) -> None:
        """Sets up the worker for the prediction of the proteins."""
        self._view.wait_spinner.start()

        # <editor-fold desc="Check if WSL2 and ColabFold are installed">
        if globals.g_os == "win32":
            constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
            if not dialog_settings_global.is_wsl2_installed():
                constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
                self._interface_manager.get_application_settings().wsl_install = 0
                basic_boxes.ok(
                    "Prediction",
                    "Prediction failed because the WSL2 environment is not installed!",
                    QtWidgets.QMessageBox.Critical,
                )
                return
            constants.PYSSA_LOGGER.info("Checking if Local Colabfold is installed ...")
            if not dialog_settings_global.is_local_colabfold_installed():
                constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")
                self._interface_manager.get_application_settings().local_colabfold = 0
                basic_boxes.ok(
                    "Prediction",
                    "Prediction failed because the ColabFold is not installed!",
                    QtWidgets.QMessageBox.Critical,
                )
                return

        # </editor-fold>

        self.prediction_type = constants.PREDICTION_TYPE_PRED_MONO_ANALYSIS
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        self._interface_manager.update_status_bar("Begin prediction process ...")
        if result[3] is True:
            constants.PYSSA_LOGGER.info("Running prediction with subsequent analysis.")
            # Analysis should be run after the prediction
            self._active_task = tasks.Task(
                target=main_presenter_async.predict_protein_with_colabfold,
                args=(
                    result[1],
                    result[2],
                    self._interface_manager.get_current_project(),
                ),
                post_func=self.__await_monomer_prediction_for_subsequent_analysis,
            )
            self._active_task.start()
        else:
            constants.PYSSA_LOGGER.info("Running only a prediction.")
            # No analysis after prediction
            self._active_task = tasks.Task(
                target=main_presenter_async.predict_protein_with_colabfold,
                args=(
                    result[1],
                    result[2],
                    self._interface_manager.get_current_project(),
                ),
                post_func=self.__await_predict_protein_with_colabfold,
            )
            self._active_task.start()

        self._view.status_bar.showMessage("A prediction is currently running ...")
        self.block_box_prediction = QtWidgets.QMessageBox()
        self.block_box_prediction.setIcon(QtWidgets.QMessageBox.Information)
        self.block_box_prediction.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self.block_box_prediction)
        self.block_box_prediction.setWindowTitle("Structure Prediction")
        self.block_box_prediction.setText("A prediction is currently running.")
        btn_abort = self.block_box_prediction.addButton("Abort", QtWidgets.QMessageBox.ActionRole)
        self.block_box_prediction.exec_()
        if self.block_box_prediction.clickedButton() == btn_abort:
            self.abort_prediction()
            self.block_box_prediction.close()
            self._view.wait_spinner.stop()
        else:
            self.block_box_prediction.close()
            self._view.wait_spinner.stop()
        # self._interface_manager.update_status_bar("")

    def abort_prediction(self) -> None:
        """Aborts the running prediction."""
        constants.PYSSA_LOGGER.info("Structure prediction process was aborted manually.")
        subprocess.run(["wsl", "--shutdown"])
        constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
        filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
        constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
        basic_boxes.ok("Abort prediction", "The structure prediction was aborted.", QtWidgets.QMessageBox.Information)
        self._interface_manager.refresh_main_view()

    def __await_monomer_prediction_for_subsequent_analysis(self, result: tuple) -> None:
        tmp_exit_code = result[0]
        tmp_exit_code_description = [1]
        if tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            self.block_box_prediction.destroy(True)
            constants.PYSSA_LOGGER.info("All structure predictions are done.")
            self._interface_manager.update_status_bar("All structure predictions are done.")
            constants.PYSSA_LOGGER.info("Begin analysis process.")
            self._interface_manager.update_status_bar("Running analysis process ...")

            tmp_raw_analysis_run_names: list = []
            for row_no in range(self._interface_manager.get_predict_monomer_view().ui.list_pred_analysis_mono_overview.count()):
                tmp_raw_analysis_run_names.append(
                    self._interface_manager.get_predict_monomer_view().ui.list_pred_analysis_mono_overview.item(row_no).text())

            self._database_manager.close_project_database()
            self._active_task = tasks.Task(
                target=main_presenter_async.run_distance_analysis,
                args=(
                    tmp_raw_analysis_run_names,
                    self._interface_manager.get_current_project(),
                    self._interface_manager.get_application_settings(),
                    self._interface_manager.get_predict_monomer_view().ui.cb_pred_analysis_mono_images.isChecked(),
                ),
                post_func=self.post_analysis_process,
            )
            self._active_task.start()
            if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
                os.mkdir(constants.SCRATCH_DIR_ANALYSIS)

        elif tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because there was an error writing the fasta file(s)!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because the fasta file(s) could not be found!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because a subprocess failed!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because of an unknown error!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        self._view.wait_spinner.stop()
        self._interface_manager.refresh_main_view()

    def post_analysis_process(self, an_exit_code: tuple[int, str, list]) -> None:
        """Post process after the analysis thread finished."""
        constants.PYSSA_LOGGER.debug("post_analysis_process() started ...")
        if an_exit_code[0] == exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[0]:
            basic_boxes.ok(
                "Distance analysis",
                "Distance analysis failed because there was an error during the analysis!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
            )
        elif an_exit_code[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            basic_boxes.ok(
                "Distance analysis",
                "Distance analysis failed because of an unknown error!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
            )
        elif an_exit_code[0] == exit_codes.EXIT_CODE_ZERO[0]:
            # self._interface_manager.get_current_project().serialize_project(
            #     self._interface_manager.get_current_project().get_project_xml_path()
            # )
            constants.PYSSA_LOGGER.info("Project has been saved to project database.")
            basic_boxes.ok(
                "Structure analysis",
                "All structure analysis' are done. Go to results to check the new results.",
                QtWidgets.QMessageBox.Information,
            )
            constants.PYSSA_LOGGER.info("All structure analysis' are done.")
        self._database_manager.open_project_database()
        self._interface_manager.refresh_main_view()
        self._interface_manager.stop_wait_spinner()

    def __await_predict_protein_with_colabfold(self, result: tuple) -> None:
        """Process which runs after each prediction job."""
        print(result)
        tmp_exit_code = result[0]
        tmp_exit_code_description = result[1]
        if tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because there was an error writing the fasta file(s)!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because the fasta file(s) could not be found!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because a subprocess failed!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because of an unknown error!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            self._interface_manager.get_current_project().serialize_project(self._interface_manager.get_current_project().get_project_xml_path())
            constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Structure prediction",
                "All structure predictions are done. Go to View to check the new proteins.",
                QtWidgets.QMessageBox.Information,
            )
            constants.PYSSA_LOGGER.info("All structure predictions are done.")
        else:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because of an unknown case!",
                QtWidgets.QMessageBox.Critical,
            )
        self._view.wait_spinner.stop()

    # </editor-fold>

    def _distance_analysis(self):
        self._external_controller = distance_analysis_view_controller.DistanceAnalysisViewController(
            self._interface_manager
        )
        self._external_controller.job_input.connect(self.start_process_batch)
        self._interface_manager.get_distance_analysis_view().show()

    def start_process_batch(self, job_input: tuple) -> None:
        """Sets up the worker for the analysis task."""
        constants.PYSSA_LOGGER.info("Begin analysis process.")

        _, tmp_raw_analysis_run_names, tmp_checkbox_state = job_input

        self._active_task = tasks.Task(
            target=main_presenter_async.run_distance_analysis,
            args=(
                tmp_raw_analysis_run_names,
                self._interface_manager.get_current_project(),
                self._interface_manager.get_application_settings(),
                tmp_checkbox_state,
            ),
            post_func=self.post_analysis_process,
        )
        self._active_task.start()

        if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
            os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
        #self.block_box_analysis.exec_()

    def _check_for_results(self) -> None:
        if self._view.ui.protein_pairs_tree_view.model().data(self._view.ui.protein_pairs_tree_view.currentIndex(), Qt.DisplayRole).find("_vs_") != -1:
            self._view.ui.action_results_summary.setEnabled(True)
        else:
            self._view.ui.action_results_summary.setEnabled(False)

    def _results_summary(self) -> None:
        tmp_protein_pair = self._view.ui.protein_pairs_tree_view.model().data(self._view.ui.protein_pairs_tree_view.currentIndex(),
                                                                              enums.ModelEnum.OBJECT_ROLE)
        self._external_controller = results_view_controller.ResultsViewController(self._interface_manager,
                                                                                  tmp_protein_pair)
        self._interface_manager.get_results_view().show()

    # <editor-fold desc="Settings menu methods">
    def open_settings_global(self) -> None:
        """Opens the dialog for the global settings."""
        dialog = dialog_settings_global.DialogSettingsGlobal()
        dialog.exec_()
        self._interface_manager.update_settings()
        self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")

    def restore_settings(self) -> None:
        """Restores the settings.xml file to the default values."""
        out = gui_utils.warning_dialog_restore_settings("Are you sure you want to restore all settings?")
        if out:
            tools.restore_default_settings(self._application_settings)
            self._view.status_bar.showMessage("Settings were successfully restored.")
            logging.info("Settings were successfully restored.")
        else:
            self._view.status_bar.showMessage("Settings were not modified.")
            logging.info("Settings were not modified.")

    # </editor-fold>

    # <editor-fold desc="About menu methods">
    def open_logs(self) -> None:
        """Opens a file explorer with all log files and can open a log file in the default application."""
        file_dialog = QtWidgets.QFileDialog()
        log_path = str(constants.LOG_PATH)
        file_dialog.setDirectory(log_path)
        file_path, _ = file_dialog.getOpenFileName(self._view, "Select a log file to open", "", "LOG File (*.log)")
        if file_path:
            os.startfile(file_path)

    @staticmethod
    def clear_all_log_files() -> None:
        """Clears all log files generated under .pyssa/logs."""
        response = basic_boxes.yes_or_no(
            "Clear log files",
            "Are you sure you want to delete all log files?",
            QtWidgets.QMessageBox.Information,
        )
        if response:
            try:
                shutil.rmtree(str(constants.LOG_PATH))
            except PermissionError:
                print("The active log file was not deleted.")
            if len(os.listdir(str(constants.LOG_PATH))) == 1:
                basic_boxes.ok("Clear log files", "All log files could be deleted.", QtWidgets.QMessageBox.Information)
                constants.PYSSA_LOGGER.info("All log files were deleted.")
            else:
                basic_boxes.ok("Clear log files", "Not all log files could be deleted.", QtWidgets.QMessageBox.Warning)
                constants.PYSSA_LOGGER.warning("Not all log files were deleted!")

    @staticmethod
    def open_tutorial() -> None:
        """Opens the official tutorial pdf file."""
        tmp_dialog = dialog_tutorial_videos.TutorialVideosDialog()
        tmp_dialog.exec_()

    @staticmethod
    def open_documentation() -> None:
        """Opens the official plugin documentation as PDF."""
        os.startfile(constants.DOCS_PATH)

    @staticmethod
    def open_about() -> None:
        """Opens the About dialog."""
        dialog = dialog_about.DialogAbout()
        dialog.exec_()

    # </editor-fold>

    # <editor-fold desc="Image menu methods">
    # TODO: images need to be reimplemented
    def post_preview_image(self) -> None:
        """Hides the block box of the preview process."""
        #self.block_box_uni.hide()
        #self.block_box_uni.destroy(True)
        self._view.status_bar.showMessage("Finished preview of ray-traced image.")
        QtWidgets.QApplication.restoreOverrideCursor()

    def preview_image(self) -> None:
        """Previews the image."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        if self._view.ui.cb_ray_tracing.isChecked():
            self._view.status_bar.showMessage("Preview ray-traced image ...")
            # <editor-fold desc="Worker setup">
            # --Begin: worker setup
            self.tmp_thread = QtCore.QThread()
            self.tmp_worker = task_workers.PreviewRayImageWorker(self.renderer)
            self.tmp_thread = task_workers.setup_worker_for_work(
                self.tmp_thread,
                self.tmp_worker,
                self.display_view_page,
            )
            self.tmp_worker.finished.connect(self.post_preview_image)
            self.tmp_thread.start()
            # --End: worker setup

            # </editor-fold>
            gui_utils.setup_standard_block_box(
                self.block_box_uni,
                "Preview ray-trace image",
                "Creating preview for the ray-traced image ...",
            )
            #self.block_box_uni.exec_()
        else:
            self._view.status_bar.showMessage("Preview draw image ...")
            cmd.draw(2400, 2400)
            self._view.status_bar.showMessage("Finished preview of drawn image.")
            QtWidgets.QApplication.restoreOverrideCursor()

    def post_save_image(self) -> None:
        """Displays a message box which informs that the process has finished."""
        self.block_box_uni.hide()
        self.block_box_uni.destroy(True)
        self._view.status_bar.showMessage("Finished image creation.")
        QtWidgets.QApplication.restoreOverrideCursor()
        basic_boxes.ok("Finished image creation", "The image has been created.", QtWidgets.QMessageBox.Information)

    def save_image(self) -> None:
        """Saves the image as a png file."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        if self._view.ui.cb_ray_tracing.isChecked():
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display(
                        "info",
                        "No file has been selected.",
                        self._view.status_bar,
                        "No file has been selected.",
                    )
                    return
                self._view.status_bar.showMessage("Creating ray-traced image ...")

                # <editor-fold desc="Worker setup">
                # --Begin: worker setup
                self.tmp_thread = QtCore.QThread()
                self.tmp_worker = task_workers.SaveRayImageWorker(self.renderer, full_file_name[0])
                self.tmp_thread = task_workers.setup_worker_for_work(
                    self.tmp_thread,
                    self.tmp_worker,
                    self.display_view_page,
                )
                self.tmp_worker.finished.connect(self.post_save_image)
                self.tmp_thread.start()
                # --End: worker setup

                # </editor-fold>
                gui_utils.setup_standard_block_box(
                    self.block_box_uni,
                    "Save ray-trace image",
                    "Creating the ray-traced image ...",
                )
                self.block_box_uni.exec_()

                # cmd.ray(2400, 2400, renderer=int(self.renderer))
                # cmd.png(full_file_name[0], dpi=300)

            except FileExistsError:
                tools.quick_log_and_display(
                    "error",
                    "File exists already.",
                    self._view.status_bar,
                    "File exists already.",
                )
            except pymol.CmdException:
                tools.quick_log_and_display(
                    "error",
                    "Unexpected Error from PyMOL while saving the " "an image",
                    self._view.status_bar,
                    "Unexpected Error from PyMOL",
                )
        else:
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display(
                        "info",
                        "No file has been selected.",
                        self._view.status_bar,
                        "No file has been selected.",
                    )
                    return
                self._view.status_bar.showMessage("Creating draw image ...")
                cmd.draw(2400, 2400)
                cmd.png(full_file_name[0], dpi=300)
                self._view.status_bar.showMessage("Finished image creation.")
                basic_boxes.ok(
                    "Finished image creation",
                    "The image has been created.",
                    QtWidgets.QMessageBox.Information,
                )
            except FileExistsError:
                tools.quick_log_and_display(
                    "error",
                    "File exists already.",
                    self._view.status_bar,
                    "File exists already.",
                )
            except pymol.CmdException:
                tools.quick_log_and_display(
                    "error",
                    "Unexpected Error from PyMOL while saving the " "an image",
                    self._view.status_bar,
                    "Unexpected Error from PyMOL",
                )
            finally:
                QtWidgets.QApplication.restoreOverrideCursor()

    # </editor-fold>

    # <editor-fold desc="Import, Export functions">
    def import_project(self) -> None:
        """Imports a project.xml into the current workspace."""
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getOpenFileName(
            self._view,
            "Select a project file to import",
            "",
            "XML Files (*.xml)",
        )
        if file_path:
            file = QtCore.QFile(file_path)
            if not file.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text):
                print("Error: Cannot open file for reading")
                return

            tmp_project = project.Project()
            handler = filesystem_io.ProjectParserHandler(tmp_project, self._interface_manager.get_application_settings())
            parser = sax.make_parser()
            parser.setContentHandler(handler)
            parser.parse(file_path)
            file.close()
            tmp_project = handler.get_project()

            tmp_project.set_workspace_path(self._workspace_path)
            if len(tmp_project.proteins) <= 1:
                if self._interface_manager.get_application_settings().wsl_install == 0:
                    basic_boxes.ok(
                        "Create new project",
                        "Please install local colabfold to import this project!",
                        QtWidgets.QMessageBox.Warning,
                    )
                    return
                elif self._interface_manager.get_application_settings().local_colabfold == 0:  # noqa: RET505
                    basic_boxes.ok(
                        "Create new project",
                        "Please install local colabfold to import this project!",
                        QtWidgets.QMessageBox.Warning,
                    )
                    return
            new_filepath = pathlib.Path(f"{self._workspace_path}/{tmp_project.get_project_name()}.xml")
            tmp_project.serialize_project(new_filepath)
            self._interface_manager.set_new_project(
                self._interface_manager.get_current_project().deserialize_project(
                    new_filepath, self._interface_manager.get_application_settings()
                )
            )
            constants.PYSSA_LOGGER.info(
                f"Opening the project {self._interface_manager.get_current_project().get_project_name()}.")
            self._view.ui.lbl_project_name.setText(self._interface_manager.get_current_project().get_project_name())
            self._interface_manager.refresh_main_view()
            basic_boxes.ok(
                "Import Project",
                "The project was successfully imported.",
                QtWidgets.QMessageBox.Information,
            )

    def export_current_project(self) -> None:
        """Exports the current project to an importable format."""
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(self._view, "Save current project", "", "XML Files (*.xml)")
        if file_path:
            self._interface_manager.get_current_project().serialize_project(pathlib.Path(file_path))
            basic_boxes.ok(
                "Export Project",
                "The project was successfully exported.",
                QtWidgets.QMessageBox.Information,
            )

    # </editor-fold>
