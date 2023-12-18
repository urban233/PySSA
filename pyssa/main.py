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
"""Module for the main window of the pyssa plugin."""
from typing import TYPE_CHECKING
from urllib.request import urlopen
from urllib.error import URLError
import logging
import os
import shutil
import subprocess
import sys
import threading
import pathlib
import csv
import requests

import numpy as np
import pymol
from pymol import cmd
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt  # pylint: disable=no-name-in-module
from PyQt5 import QtMultimedia

from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures import project_watcher
from pyssa.internal.data_structures import settings
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.internal.data_structures.data_classes import stage
from pyssa.internal.data_structures.data_classes import current_session
from pyssa.internal.portal import graphic_operations, pymol_io
from pyssa.internal.thread import workers, task_workers
from pyssa.gui.ui.forms.auto_generated.auto_main_window import Ui_MainWindow
from pyssa.gui.ui.dialogs import dialog_settings_global, dialog_add_model, dialog_display_docs, dialog_help
from pyssa.gui.ui.dialogs import dialog_startup
from pyssa.gui.ui.dialogs import dialog_distance_plot
from pyssa.gui.ui.dialogs import dialog_distance_histogram
from pyssa.gui.ui.dialogs import dialog_about
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.dialogs.dialog_tutorial_videos import TutorialVideosDialog
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.io_pyssa import safeguard, bio_data
from pyssa.io_pyssa import filesystem_io
from pyssa.io_pyssa import path_util
from pyssa.util import pyssa_keys, exit_codes
from pyssa.util import globals
from pyssa.util import protein_pair_util, session_util, exception
from pyssa.util import constants, input_validator, gui_page_management, tools, gui_utils


if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair


class MainWindow(QtWidgets.QMainWindow):
    """This class contains all information about the MainWindow in the application.

    Args:
        app_project:
            The active project.Project object which contains all information about the currently loaded project.
        _project_watcher:
            The project watcher object to determine the status of the project.
        scratch_path:
            The path where all files are stored temporarily before being copied in the project directory.
        workspace_path:
            The path where all projects are stored.
        workspace:
            A QLabel which holds information from workspace_path.
        status_bar:
            The status_bar object for the main window.
        no_of_selected_chains:
            Holds the number of selected chains from the single analysis page.
        plot_dialog:
            Is a Qt dialog window which is used in combination with pyqtgraph.
        view_box:
            The view box object from pyqtgraph, to manipulate the graphs view.
        local_pred_monomer_management:
            The object which contains information in form of stages of the page.
        local_pred_multimer_management:
            The object which contains information in form of stages of the page.
        single_analysis_management:
            The object which contains information in form of stages of the page.
    """

    def __init__(self) -> None:
        """Constructor."""
        # <editor-fold desc="Initialize the ui building process">
        super().__init__()
        # build ui object
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.lbl_page_title.setText("Home")
        self.setMinimumWidth(580)
        self.setMinimumHeight(200)

        constants.PYSSA_LOGGER.info("Successful initialization of basic UI.")
        # </editor-fold>

        # <editor-fold desc="OS Check">
        if sys.platform.startswith("win32"):
            # Windows path
            globals.g_os = "win32"
            globals.g_plugin_path = pathlib.Path(f"C:\\ProgramData\\pyssa\\mambaforge_pyssa\\pyssa-mamba-env\\Lib\\site-packages\\pymol\\pymol_path\\data\\startup\\{constants.PLUGIN_NAME}")
        elif sys.platform.startswith("linux"):
            # Linux path
            globals.g_os = "linux"
            globals.g_plugin_path = f"/home/{os.getlogin()}/.local/pyssa/pyssa-mamba-env/lib/python3.10/site-packages/pmg_tk/startup/{constants.PLUGIN_NAME}"

        constants.PYSSA_LOGGER.info(f"Started on platform: {globals.g_os}")
        # </editor-fold>

        globals.g_settings = settings.Settings("", "")
        pixmapi = QtWidgets.QStyle.SP_MessageBoxQuestion
        icon = self.style().standardIcon(pixmapi)
        self.ui.btn_info.setIcon(icon)
        self.ui.btn_info.setText("")
        self.ui.btn_info.setFixedWidth(50)

        # <editor-fold desc="Program directory check">
        if not os.path.exists(str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa"))):
            os.mkdir(str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa")))
        if not os.path.exists(str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/logs"))):
            os.mkdir(str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/logs")))

        constants.PYSSA_LOGGER.info(f"Checked program and logs directory.")
        # </editor-fold>

        # <editor-fold desc="Setup App Settings">
        self.app_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        if not os.path.exists(constants.SETTINGS_FULL_FILEPATH):
            constants.PYSSA_LOGGER.info(f"Settings file not found, open configuration dialog.")
            # Configuration dialog to setup setting file
            dialog = dialog_startup.DialogStartup()
            dialog.exec_()

            # checks if the cancel button was pressed
            if dialog_startup.global_var_terminate_app == 1:
                os.remove(constants.SETTINGS_FULL_FILEPATH)
                constants.PYSSA_LOGGER.info(f"Configuration dialog closed, and removed new settings file.")
                sys.exit()

            self.app_settings.app_launch = 1
            self.app_settings.workspace_path = pathlib.Path(dialog_startup.global_var_startup_workspace)

            constants.PYSSA_LOGGER.info(f"Demo projects are getting downloaded and extracted ...")
            import zipfile
            with zipfile.ZipFile(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"), 'r') as zip_ref:
                zip_ref.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
            constants.PYSSA_LOGGER.info(f"Demo projects are downloaded and extracted.\n Import of demo projects started ...")

            path_of_demo_projects = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
            tmp_project = project.Project("", dialog_startup.global_var_startup_workspace)
            for tmp_filename in os.listdir(path_of_demo_projects):
                try:
                    tmp_project = tmp_project.deserialize_project(
                        pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}"),
                        self.app_settings,
                    )
                except exception.IllegalArgumentError:
                    constants.PYSSA_LOGGER.warning("The workspace path does not exist on this system, "
                                                   "but this is due to the demo projects.")
                tmp_project.set_workspace_path(dialog_startup.global_var_startup_workspace)
                new_filepath = pathlib.Path(f"{dialog_startup.global_var_startup_workspace}/{tmp_filename}")
                tmp_project.serialize_project(new_filepath)
            constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
            try:
                os.remove(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"))
            except FileNotFoundError:
                constants.PYSSA_LOGGER.warning("Zip archive of demo projects could not be found!")
            constants.PYSSA_LOGGER.info("Serialize settings ...")
            self.app_settings.serialize_settings()
            constants.PYSSA_LOGGER.info("Serialize settings finished.")

            QtWidgets.QApplication.restoreOverrideCursor()

        try:
            self.app_settings = self.app_settings.deserialize_settings()
        except ValueError:
            constants.PYSSA_LOGGER.warning("The settings file is damaged.")
            gui_utils.error_dialog_settings(
                "The settings file is damaged. You have to restore the settings to use PySSA!", "", self.app_settings,
            )
        if globals.g_os == "win32":
            constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
            if dialog_settings_global.is_wsl2_installed():
                self.app_settings.wsl_install = 1
                constants.PYSSA_LOGGER.info("WSL2 is installed.")
            else:
                self.app_settings.wsl_install = 0
                constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
        else:
            self.app_settings.wsl_install = 1

        constants.PYSSA_LOGGER.info("Checking if Local Colabfold is installed ...")
        if dialog_settings_global.is_local_colabfold_installed():
            self.app_settings.local_colabfold = 1
            constants.PYSSA_LOGGER.info("Local Colabfold is installed.")
        else:
            self.app_settings.local_colabfold = 0
            constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")

        globals.g_settings = self.app_settings
        # </editor-fold>

        constants.PYSSA_LOGGER.info("Checking version of localhost with remote ...")
        url = "https://w-hs.sciebo.de/s/0kR11n8PkLDo1gB/download"
        try:
            response = requests.get(url)
            print(f"Current version: {constants.VERSION_NUMBER[1:]}")
            print(f"Status code: {response.status_code}")
            if response.status_code == 503:
                constants.PYSSA_LOGGER.warning("The connection to the update servers failed. If could not determine if a new version is avaliable.")
                basic_boxes.ok(
                    "No connection",
                    f"The connection to the update servers failed. If could not determine if a new version is avaliable. You can start the PySSA now.",
                    QtWidgets.QMessageBox.Information,
                )
            else:
                print(f"Latest version: {response.text}")
                if response.text != constants.VERSION_NUMBER[1:]:
                    constants.PYSSA_LOGGER.info("here is a new version of PySSA available.")
                    basic_boxes.ok(
                        "New version",
                        f"There is a new version for your PySSA!\nTo install the latest version {response.text}, open the PySSA Installer and click on update.",
                        QtWidgets.QMessageBox.Information,
                    )
        except requests.exceptions.RequestException as e:
            constants.PYSSA_LOGGER.error("Downloading the version file of the remote failed!")
            print("Failed to download the file:", e)
        constants.PYSSA_LOGGER.info("Checking version of localhost with remote finished.")

        # <editor-fold desc="Class attributes">
        self.app_project = project.Project("", pathlib.Path(""))
        self._project_watcher = project_watcher.ProjectWatcher(self.app_project, no_of_pdb_files=None)
        self.scratch_path = constants.SCRATCH_DIR
        self.workspace_path = self.app_settings.workspace_path
        self.workspace = QtWidgets.QLabel(f"Current Workspace: {self.workspace_path}")
        self.status_bar = QtWidgets.QStatusBar()
        self.prediction_configuration = prediction_configuration.PredictionConfiguration(True, "pdb70")
        self.results_name = ""
        self.no_of_selected_chains = 0
        self.plot_dialog = QtWidgets.QDialog(self)
        self.view_box = None
        self.block_box_expert_install = basic_boxes.no_buttons(
            "Local Colabfold installation",
            "An installation process is currently running.",
            QtWidgets.QMessageBox.Information,
        )
        self.prediction_type = 0
        self.current_session = current_session.CurrentSession("", "", "")
        self.is_distance_plot_open = False
        self.distance_plot_dialog = None

        constants.PYSSA_LOGGER.info("Setup class attributes finished.")
        # </editor-fold>

        # sets up the status bar
        self._setup_statusbar()
        tools.create_directory(constants.SETTINGS_DIR, "scratch")
        self._setup_default_configuration()

        constants.PYSSA_LOGGER.info("Setup rest of GUI related elements ...")

        # <editor-fold desc="GUI page management">
        # -- Gui page management vars
        self.local_pred_monomer_management: gui_page_management.GuiPageManagement
        self.local_pred_multimer_management: gui_page_management.GuiPageManagement
        self.single_analysis_management: gui_page_management.GuiPageManagement
        self.batch_analysis_management: gui_page_management.GuiPageManagement
        self.results_management: gui_page_management.GuiPageManagement
        self.monomer_prediction_analysis_management: gui_page_management.GuiPageManagement
        self.multi_prediction_analysis_management: gui_page_management.GuiPageManagement

        # management functions
        self._create_local_pred_monomer_management()
        self._create_local_pred_multimer_management()
        self._create_single_analysis_management()
        self._create_batch_analysis_management()
        self._create_results_management()
        self._create_monomer_prediction_analysis_management()
        self._create_multimer_prediction_analysis_management()

        # </editor-fold>

        # <editor-fold desc="Worker definitions">
        self.worker_prediction = workers.PredictionWorkerPool(
            self.ui.table_pred_analysis_mono_prot_to_predict, self.prediction_configuration, self.app_project,
        )
        self.worker_prediction.signals.finished.connect(self.post_prediction_process)

        self.worker_analysis = workers.AnalysisWorkerPool(
            self.ui.list_analysis_batch_overview,
            self.ui.cb_analysis_images,
            self.status_bar,
            self.app_project,
            self.app_settings,
            self._init_batch_analysis_page,
        )
        self.worker_analysis.signals.finished.connect(self.post_analysis_process)
        self.worker_analysis.setAutoDelete(True)

        self.worker_image_creation = workers.BatchImageWorkerPool(
            self.ui.list_analysis_images_struct_analysis,
            self.ui.list_analysis_images_creation_struct_analysis,
            self.status_bar,
            self.app_project,
        )
        self.worker_image_creation.signals.finished.connect(self.post_image_creation_process)
        self.worker_image_creation.setAutoDelete(True)

        # </editor-fold>

        # <editor-fold desc="Block box definitions">
        self.block_box_analysis = basic_boxes.no_buttons(
            "Analysis", "An analysis is currently running, please wait.", QtWidgets.QMessageBox.Information,
        )
        self.block_box_prediction: QtWidgets.QMessageBox = QtWidgets.QMessageBox()
        self.block_box_images = basic_boxes.no_buttons(
            "Analysis Images", "Images getting created, please wait.", QtWidgets.QMessageBox.Information,
        )
        self.block_box_uni = basic_boxes.no_buttons("Generic", "Generic", QtWidgets.QMessageBox.Information)

        # </editor-fold>

        # configure gui element properties
        self.ui.txt_results_aligned_residues.setAlignment(QtCore.Qt.AlignRight)
        self.ui.table_pred_mono_prot_to_predict.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.ui.table_pred_mono_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self.ui.table_pred_multi_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self.ui.list_pred_analysis_multi_ref_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.list_pred_analysis_multi_model_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)

        # helper attributes
        self.pymol_session_specs = {
            pyssa_keys.SESSION_SPEC_PROTEIN: [0, ""],
            pyssa_keys.SESSION_SPEC_COLOR: [0, ""],
            pyssa_keys.SESSION_SPEC_REPRESENTATION: [0, ""],
            pyssa_keys.SESSION_SPEC_BG_COLOR: [0, ""],
        }

        # setup defaults for pages
        self._init_fill_combo_boxes()
        self._init_new_page()
        self._init_use_page()
        self._init_local_pred_mono_page()
        self._init_local_pred_multi_page()
        # self._init_sequence_vs_pdb_page()
        self._init_single_analysis_page()
        self._init_batch_analysis_page()
        self.ui.action_toggle_notebook_visibility.setVisible(False)
        self.ui.action_settings_model_w_off_colab_notebook.setVisible(False)
        self.last_sidebar_button = QtWidgets.QPushButton()
        self.ui.table_pred_mono_prot_to_predict.setEditTriggers(self.ui.table_pred_mono_prot_to_predict.NoEditTriggers)
        self.ui.table_pred_multi_prot_to_predict.setEditTriggers(
            self.ui.table_pred_multi_prot_to_predict.NoEditTriggers,
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setEditTriggers(
            self.ui.table_pred_analysis_mono_prot_to_predict.NoEditTriggers,
        )
        self.ui.table_pred_analysis_multi_prot_to_predict.setEditTriggers(
            self.ui.table_pred_analysis_multi_prot_to_predict.NoEditTriggers,
        )

        # connections
        self._connect_all_gui_elements()
        # create tooltips
        self._create_all_tooltips()
        self._project_watcher.show_valid_options(self.ui)
        self._project_watcher.check_workspace_for_projects(self.workspace_path, self.ui)
        self.project_scanner = filesystem_io.ProjectScanner(self.app_project)
        # sets threadpool
        self.threadpool = QtCore.QThreadPool()
        # create scratch and cache dir
        if not os.path.exists(constants.SCRATCH_DIR):
            os.mkdir(constants.SCRATCH_DIR)
        if not os.path.exists(constants.CACHE_DIR):
            os.mkdir(constants.CACHE_DIR)

        if self.app_settings.wsl_install == 1 and self.app_settings.local_colabfold == 0:
            self.ui.action_install_from_file.setVisible(True)
        else:
            self.ui.action_install_from_file.setVisible(False)
        
        # temp gui changes (need to be changed in the designer)
        self.ui.lbl_hotspots_resi_show.setText("Residue(s) as sticks")
        self.ui.lbl_hotspots_resi_hide.setText("Residue(s) as sticks")
        
        # fixme: should the pdf documentation be accessible through the pyssa gui?
        self.ui.action_help_docs_pdf.setText("Documentation")
        self.ui.action_help_docs_pdf.setVisible(True)
        self.ui.action_help_docs.setText("Tutorials")
        self.ui.action_help_docs.setVisible(True)
        # sets additional parameters
        self.ui.lbl_logo.setPixmap(QtGui.QPixmap(str(pathlib.Path(f"{constants.PLUGIN_ROOT_PATH}/assets/images/pyssa_logo.png"))))
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("PySSA")
        constants.PYSSA_LOGGER.info(f"PySSA started with version {constants.VERSION_NUMBER}.")

    # <editor-fold desc="GUI page management functions">
    def _create_local_pred_monomer_management(self) -> None:
        """Creates a list of gui stages."""
        # gui element management
        tmp_stages = [
            # add protein stage
            stage.Stage(
                {
                    "label_proteins_to_predict": self.ui.lbl_pred_mono_prot_to_predict,
                    "table_proteins_to_predict": self.ui.table_pred_mono_prot_to_predict,
                },
                {
                    "remove_button": self.ui.btn_pred_mono_seq_to_predict_remove,
                    "next_button": self.ui.btn_pred_mono_seq_to_predict,
                },
            ),
            # protein name stage
            stage.Stage(
                {
                    "label_protein_name": self.ui.lbl_pred_mono_prot_name,
                    "text_field_protein_name": self.ui.txt_pred_mono_prot_name,
                    "label_protein_name_status": self.ui.lbl_pred_mono_prot_name_status,
                },
                {
                    "back_button": self.ui.btn_pred_mono_back,
                    "next_button": self.ui.btn_pred_mono_next,
                },
            ),
            # protein sequence stage
            stage.Stage(
                {
                    "label_protein_sequence": self.ui.lbl_pred_mono_seq_name,
                    "text_field_protein_sequence": self.ui.txt_pred_mono_seq_name,
                    "label_protein_sequence_status": self.ui.lbl_pred_mono_seq_name_status,
                },
                {
                    "back_button": self.ui.btn_pred_mono_back_2,
                    "next_button": self.ui.btn_pred_mono_add_protein,
                },
            ),
            # prediction stage (with advanced configurations)
            stage.Stage(
                {
                    "label_advanced_config": self.ui.lbl_pred_mono_advanced_config,
                    "button_advanced_config": self.ui.btn_pred_mono_advanced_config,
                },
                {
                    "predict_button": self.ui.btn_pred_mono_predict,
                },
            ),
        ]
        self.local_pred_monomer_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_local_pred_multimer_management(self) -> None:
        """Creates a list of gui stages."""
        # gui element management
        tmp_stages = [
            # add protein stage
            stage.Stage(
                {
                    "label_proteins_to_predict": self.ui.lbl_pred_multi_prot_to_predict,
                    "table_proteins_to_predict": self.ui.table_pred_multi_prot_to_predict,
                },
                {
                    "remove_button": self.ui.btn_pred_multi_prot_to_predict_remove,
                    "next_button": self.ui.btn_pred_multi_prot_to_predict_add,
                },
            ),
            # protein name stage
            stage.Stage(
                {
                    "label_protein_name": self.ui.lbl_pred_multi_prot_name,
                    "text_field_protein_name": self.ui.txt_pred_multi_prot_name,
                    "label_protein_name_status": self.ui.lbl_pred_multi_prot_name_status,
                },
                {
                    "back_button": self.ui.btn_pred_multi_back,
                    "next_button": self.ui.btn_pred_multi_next,
                },
            ),
            # protein sequence stage
            stage.Stage(
                {
                    "label_protein_sequence": self.ui.lbl_pred_multi_prot_seq,
                    "text_field_protein_sequence": self.ui.txt_pred_multi_prot_seq,
                    "label_protein_sequence_status": self.ui.lbl_pred_multi_prot_seq_status,
                    "label_protein_sequence_add": self.ui.lbl_pred_multi_prot_seq_add,
                    "button_protein_sequence_add": self.ui.btn_pred_multi_prot_seq_add,
                    "label_protein_sequence_overview": self.ui.lbl_pred_multi_prot_seq_overview,
                    "list_protein_sequence_overview": self.ui.list_pred_multi_prot_seq_overview,
                    "button_protein_sequence_overview_remove": self.ui.btn_pred_multi_prot_seq_overview_remove,
                    "label_protein_to_predict": self.ui.lbl_pred_multi_prot_to_predict_2,
                },
                {
                    "back_button": self.ui.btn_pred_multi_back_2,
                    "next_button": self.ui.btn_pred_multi_prot_to_predict_add_2,
                },
            ),
            # prediction stage (with advanced configurations)
            stage.Stage(
                {
                    "label_advanced_config": self.ui.lbl_pred_multi_advanced_config,
                    "button_advanced_config": self.ui.btn_pred_multi_advanced_config,
                },
                {
                    "predict_button": self.ui.btn_pred_multi_predict,
                },
            ),
        ]
        self.local_pred_multimer_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_single_analysis_management(self) -> None:
        """Creates a list of gui stages."""
        # gui element management
        tmp_stages = [
            # choose protein structures: stage 0
            stage.Stage(
                {
                    "label_protein_structure_1": self.ui.lbl_analysis_prot_struct_1,
                    "box_protein_structure_1": self.ui.box_analysis_prot_struct_1,
                    "label_vs": self.ui.lbl_analysis_vs,
                    "label_protein_structure_2": self.ui.lbl_analysis_prot_struct_2,
                    "box_protein_structure_2": self.ui.box_analysis_prot_struct_2,
                },
                {
                    "next_button": self.ui.btn_analysis_next,
                },
            ),
            # choose chains from prot structure 1: stage 1
            stage.Stage(
                {
                    "label_protein_structure_1_chains": self.ui.lbl_analysis_ref_chains,
                    "list_protein_structure_1_chains": self.ui.list_analysis_ref_chains,
                },
                {
                    "back_button": self.ui.btn_analysis_back,
                    "next_button": self.ui.btn_analysis_next_2,
                },
            ),
            # choose chains from prot structure 2: stage 2
            stage.Stage(
                {
                    "label_protein_structure_2_chains": self.ui.lbl_analysis_model_chains,
                    "list_protein_structure_2_chains": self.ui.list_analysis_model_chains,
                    "label_images": self.ui.lbl_analysis_images,
                    "checkbox_images": self.ui.cb_analysis_images,
                },
                {
                    "back_button": self.ui.btn_analysis_back_2,
                    "start_button": self.ui.btn_analysis_start,
                },
            ),
        ]
        self.single_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_batch_analysis_management(self) -> None:
        """Creates a list of gui stages."""
        # gui element management
        tmp_stages = [
            # add a prot analysis: stage 0
            stage.Stage(
                {
                    "label_batch_analysis_overview": self.ui.lbl_analysis_batch_overview,
                    "box_protein_structure_1": self.ui.list_analysis_batch_overview,
                },
                {
                    "add_button": self.ui.btn_analysis_batch_add,
                    "remove_button": self.ui.btn_analysis_batch_remove,
                },
            ),
            # choose protein structures: stage 1
            stage.Stage(
                {
                    "label_protein_structure_1": self.ui.lbl_analysis_batch_prot_struct_1,
                    "box_protein_structure_1": self.ui.box_analysis_batch_prot_struct_1,
                    "label_vs": self.ui.lbl_analysis_batch_vs,
                    "label_protein_structure_2": self.ui.lbl_analysis_batch_prot_struct_2,
                    "box_protein_structure_2": self.ui.box_analysis_batch_prot_struct_2,
                },
                {
                    "next_button": self.ui.btn_analysis_batch_next,
                    "back_button": self.ui.btn_analysis_batch_back,
                },
            ),
            # choose chains from prot structure 1: stage 2
            stage.Stage(
                {
                    "label_protein_structure_1_chains": self.ui.lbl_analysis_batch_ref_chains,
                    "list_protein_structure_1_chains": self.ui.list_analysis_batch_ref_chains,
                },
                {
                    "back_button": self.ui.btn_analysis_batch_back_2,
                    "next_button": self.ui.btn_analysis_batch_next_2,
                },
            ),
            # choose chains from prot structure 2: stage 3
            stage.Stage(
                {
                    "label_protein_structure_2_chains": self.ui.lbl_analysis_batch_model_chains,
                    "list_protein_structure_2_chains": self.ui.list_analysis_batch_model_chains,
                },
                {
                    "back_button": self.ui.btn_analysis_batch_back_3,
                    "next_button": self.ui.btn_analysis_batch_next_3,
                },
            ),
            # start batch run: stage 4
            stage.Stage(
                {
                    "label_images": self.ui.lbl_analysis_batch_images,
                    "checkbox_images": self.ui.cb_analysis_batch_images,
                },
                {
                    "start_button": self.ui.btn_analysis_batch_start,
                },
            ),
        ]
        self.batch_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_results_management(self) -> None:
        """Creates a list of gui stages."""
        # gui element management
        tmp_stages = [
            # choose protein structures: stage 0
            stage.Stage(
                {
                    "label_analysis_options": self.ui.lbl_results_analysis_options,
                    "box_results_analysis_options": self.ui.cb_results_analysis_options,
                },
                {
                    "": None,
                },
            ),
            # choose chains from prot structure 1: stage 1
            stage.Stage(
                {
                    "label_results_rmsd": self.ui.lbl_results_rmsd,
                    "text_results_rmsd": self.ui.txt_results_rmsd,
                    "label_color_rmsd": self.ui.lbl_color_rmsd,
                    "button_color_rmsd": self.ui.btn_color_rmsd,
                    "label_results_aligned_residues": self.ui.lbl_results_aligned_residues,
                    "text_results_aligned_residues": self.ui.txt_results_aligned_residues,
                    "label_results_distance_plot": self.ui.lbl_results_distance_plot,
                    "button_view_distance_plot": self.ui.btn_view_distance_plot,
                    "label_results_distance_histogram": self.ui.lbl_results_distance_histogram,
                    "button_view_distance_histogram": self.ui.btn_view_distance_histogram,
                    "label_results_distance_table": self.ui.lbl_results_distance_table,
                    "button_view_distance_table": self.ui.btn_view_distance_table,
                    "label_results_structure_alignment": self.ui.lbl_results_structure_alignment,
                    "button_view_struct_alignment": self.ui.btn_view_struct_alignment,
                    "label_results_interest_regions": self.ui.lbl_results_interest_regions,
                    "list_results_interest_regions": self.ui.list_results_interest_regions,
                    "button_results_interest_regions": self.ui.btn_view_interesting_region,
                },
                {
                    "": None,
                },
            ),
        ]
        self.results_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_monomer_prediction_analysis_management(self) -> None:
        """Creates a list of gui stages."""
        tmp_stages = [
            # add protein stage 0
            stage.Stage(
                {
                    "label_proteins_to_predict": self.ui.lbl_pred_analysis_mono_prot_to_predict,
                    "table_proteins_to_predict": self.ui.table_pred_analysis_mono_prot_to_predict,
                },
                {
                    "remove_button": self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                    "next_button": self.ui.btn_pred_analysis_mono_seq_to_predict,
                },
            ),
            # protein name stage 1
            stage.Stage(
                {
                    "label_protein_name": self.ui.lbl_pred_analysis_mono_prot_name,
                    "text_field_protein_name": self.ui.txt_pred_analysis_mono_prot_name,
                    "label_protein_name_status": self.ui.lbl_pred_analysis_mono_prot_name_status,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_mono_back,
                    "next_button": self.ui.btn_pred_analysis_mono_next,
                },
            ),
            # protein sequence stage 2
            stage.Stage(
                {
                    "label_protein_sequence": self.ui.lbl_pred_analysis_mono_seq_name,
                    "text_field_protein_sequence": self.ui.txt_pred_analysis_mono_seq_name,
                    "label_protein_sequence_status": self.ui.lbl_pred_analysis_mono_seq_name_status,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_mono_back_2,
                    "next_button": self.ui.btn_pred_analysis_mono_add_protein,
                },
            ),
            # prediction stage (with advanced configurations) 3
            stage.Stage(
                {
                    "label_advanced_config": self.ui.lbl_pred_mono_advanced_config_2,
                    "button_advanced_config": self.ui.btn_pred_mono_advanced_config_2,
                },
                {
                    "label_go_to_analysis": self.ui.lbl_pred_analysis_mono_to_analysis_setup,
                    "button_go_to_analysis": self.ui.btn_pred_analysis_mono_go_analysis_setup,
                },
            ),
            # add a prot analysis: stage 4
            stage.Stage(
                {
                    "label_batch_analysis_overview": self.ui.lbl_pred_analysis_mono_overview,
                    "box_protein_structure_1": self.ui.list_pred_analysis_mono_overview,
                },
                {
                    "add_button": self.ui.btn_pred_analysis_mono_add,
                    "remove_button": self.ui.btn_pred_analysis_mono_remove,
                },
            ),
            # choose protein structures: stage 5
            stage.Stage(
                {
                    "label_protein_structure_1": self.ui.lbl_pred_analysis_mono_prot_struct_1,
                    "box_protein_structure_1": self.ui.box_pred_analysis_mono_prot_struct_1,
                    "label_vs": self.ui.lbl_analysis_batch_vs_2,
                    "label_protein_structure_2": self.ui.lbl_pred_analysis_mono_prot_struct_2,
                    "box_protein_structure_2": self.ui.box_pred_analysis_mono_prot_struct_2,
                },
                {
                    "next_button": self.ui.btn_pred_analysis_mono_next_2,
                    "back_button": self.ui.btn_pred_analysis_mono_back_3,
                },
            ),
            # choose chains from prot structure 1: stage 6
            stage.Stage(
                {
                    "label_protein_structure_1_chains": self.ui.lbl_pred_analysis_mono_ref_chains,
                    "list_protein_structure_1_chains": self.ui.list_pred_analysis_mono_ref_chains,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_mono_back_4,
                    "next_button": self.ui.btn_pred_analysis_mono_next_3,
                },
            ),
            # choose chains from prot structure 2: stage 7
            stage.Stage(
                {
                    "label_protein_structure_2_chains": self.ui.lbl_pred_analysis_mono_model_chains,
                    "list_protein_structure_2_chains": self.ui.list_pred_analysis_mono_model_chains,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_mono_back_5,
                    "next_button": self.ui.btn_pred_analysis_mono_next_4,
                },
            ),
            # start batch run: stage 8
            stage.Stage(
                {
                    "label_images": self.ui.lbl_pred_analysis_mono_images,
                    "checkbox_images": self.ui.cb_pred_analysis_mono_images,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_mono_back_pred_setup,
                    "start_button": self.ui.btn_pred_analysis_mono_start,
                },
            ),
        ]
        self.monomer_prediction_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_multimer_prediction_analysis_management(self) -> None:
        """Creates a list of gui stages."""
        tmp_stages = [
            # add protein stage
            stage.Stage(
                {
                    "label_proteins_to_predict": self.ui.lbl_pred_analysis_multi_prot_to_predict,
                    "table_proteins_to_predict": self.ui.table_pred_analysis_multi_prot_to_predict,
                },
                {
                    "remove_button": self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                    "next_button": self.ui.btn_pred_analysis_multi_prot_to_predict_add,
                },
            ),
            # protein name stage
            stage.Stage(
                {
                    "label_protein_name": self.ui.lbl_pred_analysis_multi_prot_name,
                    "text_field_protein_name": self.ui.txt_pred_analysis_multi_prot_name,
                    "label_protein_name_status": self.ui.lbl_pred_analysis_multi_prot_name_status,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_multi_back,
                    "next_button": self.ui.btn_pred_analysis_multi_next,
                },
            ),
            # protein sequence stage
            stage.Stage(
                {
                    "label_protein_sequence": self.ui.lbl_pred_analysis_multi_prot_seq,
                    "text_field_protein_sequence": self.ui.txt_pred_analysis_multi_prot_seq,
                    "label_protein_sequence_status": self.ui.lbl_pred_analysis_multi_prot_seq_status,
                    "label_protein_sequence_add": self.ui.lbl_pred_multi_prot_seq_add_2,
                    "button_protein_sequence_add": self.ui.btn_pred_analysis_multi_prot_seq_add,
                    "label_protein_sequence_overview": self.ui.lbl_pred_analysis_multi_prot_seq_overview,
                    "list_protein_sequence_overview": self.ui.list_pred_analysis_multi_prot_seq_overview,
                    "button_protein_sequence_overview_remove": self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                    "label_protein_to_predict": self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_multi_back_2,
                    "next_button": self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
                },
            ),
            # prediction stage (with advanced configurations)
            stage.Stage(
                {
                    "label_advanced_config": self.ui.lbl_pred_analysis_multi_advanced_config,
                    "button_advanced_config": self.ui.btn_pred_analysis_multi_advanced_config,
                },
                {
                    "label_go_to_analysis": self.ui.lbl_pred_analysis_multi_to_analysis_setup,
                    "button_go_to_analysis": self.ui.btn_pred_analysis_multi_go_analysis_setup,
                },
            ),
            # add a prot analysis: stage 0
            stage.Stage(
                {
                    "label_batch_analysis_overview": self.ui.lbl_pred_analysis_multi_overview,
                    "box_protein_structure_1": self.ui.list_pred_analysis_multi_overview,
                },
                {
                    "add_button": self.ui.btn_pred_analysis_multi_add,
                    "remove_button": self.ui.btn_pred_analysis_multi_remove,
                },
            ),
            # choose protein structures: stage 1
            stage.Stage(
                {
                    "label_protein_structure_1": self.ui.lbl_pred_analysis_multi_prot_struct_1,
                    "box_protein_structure_1": self.ui.box_pred_analysis_multi_prot_struct_1,
                    "label_vs": self.ui.lbl_analysis_batch_vs_3,
                    "label_protein_structure_2": self.ui.lbl_pred_analysis_multi_prot_struct_2,
                    "box_protein_structure_2": self.ui.box_pred_analysis_multi_prot_struct_2,
                },
                {
                    "next_button": self.ui.btn_pred_analysis_multi_next_2,
                    "back_button": self.ui.btn_pred_analysis_multi_back_3,
                },
            ),
            # choose chains from prot structure 1: stage 2
            stage.Stage(
                {
                    "label_protein_structure_1_chains": self.ui.lbl_pred_analysis_multi_ref_chains,
                    "list_protein_structure_1_chains": self.ui.list_pred_analysis_multi_ref_chains,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_multi_back_4,
                    "next_button": self.ui.btn_pred_analysis_multi_next_3,
                },
            ),
            # choose chains from prot structure 2: stage 3
            stage.Stage(
                {
                    "label_protein_structure_2_chains": self.ui.lbl_pred_analysis_multi_model_chains,
                    "list_protein_structure_2_chains": self.ui.list_pred_analysis_multi_model_chains,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_multi_back_5,
                    "next_button": self.ui.btn_pred_analysis_multi_next_4,
                },
            ),
            # start batch run: stage 4
            stage.Stage(
                {
                    "label_images": self.ui.lbl_pred_analysis_multi_images,
                    "checkbox_images": self.ui.cb_pred_analysis_multi_images,
                },
                {
                    "back_button": self.ui.btn_pred_analysis_multi_back_pred_setup,
                    "start_button": self.ui.btn_pred_analysis_multi_start,
                },
            ),
        ]
        self.multimer_prediction_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    # </editor-fold>

    def _setup_statusbar(self) -> None:
        """Sets up the status bar and fills it with the current workspace."""
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage(str(self.workspace.text()))

    def _setup_default_configuration(self) -> None:
        """Sets up the default values for specific gui elements."""
        self.ui.lbl_current_project_name.setText("")
        # menu
        # side menu

        # new project page
        self.ui.btn_new_create_project.setEnabled(False)
        self.ui.cb_new_add_reference.setCheckable(False)
        self.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
        # open project page
        self.ui.lbl_open_status_search.setText("")
        self.ui.btn_open_open_project.setEnabled(False)
        # delete project page
        self.ui.lbl_delete_status_search.setText("")
        self.ui.btn_delete_delete_project.setEnabled(False)
        # edit project page

        # view project page

        # use project page

        # new sequence page
        # sequence vs .pdb page
        self.ui.btn_s_v_p_start.setEnabled(False)
        self.ui.list_s_v_p_ref_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        # single analysis page
        self.ui.lbl_analysis_model_chains.hide()
        self.ui.list_analysis_model_chains.hide()
        self.ui.btn_analysis_back.hide()
        self.ui.btn_analysis_start.hide()
        # batch analysis page

        # results page

        # image page

    def _connect_all_gui_elements(self) -> None:
        """Connects all gui elements with their corresponding slots."""
        # <editor-fold desc="Menu">
        self.ui.action_file_quit.triggered.connect(self.quit_app)
        self.ui.action_file_restore_settings.triggered.connect(self.restore_settings)
        self.ui.action_settings_edit_all.triggered.connect(self.open_settings_global)
        self.ui.action_help_docs.triggered.connect(self.open_tutorial)
        self.ui.action_help_docs_pdf.triggered.connect(self.open_documentation)
        self.ui.action_help_about.triggered.connect(self.open_about)
        self.ui.action_settings_open_logs.triggered.connect(self.open_logs)
        self.ui.action_settings_clear_logs.triggered.connect(self.clear_all_log_files)
        # </editor-fold>

        self.ui.btn_info.clicked.connect(self.open_page_information)

        # <editor-fold desc="Side Menu">
        self.ui.btn_new_page.clicked.connect(self.display_new_page)
        self.ui.btn_open_page.clicked.connect(self.display_open_page)
        self.ui.btn_delete_page.clicked.connect(self.display_delete_page)
        self.ui.btn_save_project.clicked.connect(self.save_project)
        self.ui.btn_edit_page.clicked.connect(self.display_edit_page)
        self.ui.btn_view_page.clicked.connect(self.display_view_page)
        self.ui.btn_use_page.clicked.connect(self.start_display_use_page)
        self.ui.btn_import_project.clicked.connect(self.import_project)
        self.ui.btn_export_project.clicked.connect(self.export_current_project)
        self.ui.btn_close_project.clicked.connect(self.close_project)
        self.ui.btn_pred_cloud_monomer_page.clicked.connect(self.display_esm_pred_mono)
        self.ui.btn_pred_local_monomer_page.clicked.connect(self.display_local_pred_mono)
        self.ui.btn_pred_local_multimer_page.clicked.connect(self.display_local_pred_multi)
        self.ui.btn_prediction_abort.clicked.connect(self.abort_prediction)
        self.ui.btn_pred_analysis_monomer_page.clicked.connect(self.display_monomer_pred_analysis)
        self.ui.btn_pred_analysis_multimer_page.clicked.connect(self.display_multimer_pred_analysis)
        self.ui.btn_single_analysis_page.clicked.connect(self.display_single_analysis_page)
        self.ui.btn_batch_analysis_page.clicked.connect(self.display_job_analysis_page)
        self.ui.btn_image_analysis_page.clicked.connect(self.display_image_analysis_page)
        self.ui.btn_results_page.clicked.connect(self.display_results_page)
        # self.ui.btn_analysis_abort.clicked.connect(self.abort_analysis)
        self.ui.btn_manage_session.clicked.connect(self.display_manage_pymol_session)
        self.ui.btn_image_page.clicked.connect(self.display_image_page)
        self.ui.btn_hotspots_page.clicked.connect(self.display_hotspots_page)

        # </editor-fold>

        # <editor-fold desc="New project page">
        self.ui.btn_new_choose_reference.clicked.connect(self.load_reference_in_project)
        self.ui.txt_new_project_name.textChanged.connect(self.validate_project_name)
        self.ui.txt_new_choose_reference.textChanged.connect(self.validate_reference_in_project)
        self.ui.cb_new_add_reference.stateChanged.connect(self.show_add_reference)
        self.ui.btn_new_create_project.clicked.connect(self.create_new_project)

        # </editor-fold>

        # <editor-fold desc="Open project page">
        self.ui.btn_open_open_project.clicked.connect(self.pre_open_project)
        self.ui.list_open_projects.doubleClicked.connect(self.pre_open_project)
        self.ui.txt_open_search.textChanged.connect(self.validate_open_search)
        self.ui.txt_open_selected_project.textChanged.connect(self.activate_open_button)
        self.ui.list_open_projects.currentItemChanged.connect(self.select_project_from_open_list)

        # </editor-fold>

        # <editor-fold desc="Delete project page">
        self.ui.btn_delete_delete_project.clicked.connect(self.delete_project)
        self.ui.txt_delete_search.textChanged.connect(self.validate_delete_search)
        self.ui.txt_delete_selected_projects.textChanged.connect(self.activate_delete_button)
        self.ui.list_delete_projects.currentItemChanged.connect(self.select_project_from_delete_list)

        # </editor-fold>

        # <editor-fold desc="Edit project page">
        self.ui.btn_edit_page.clicked.connect(self.display_edit_page)
        self.ui.list_edit_project_proteins.currentItemChanged.connect(self.check_for_cleaning)
        self.ui.btn_edit_clean_new_prot.clicked.connect(self.clean_protein_new)
        self.ui.btn_edit_clean_update_prot.clicked.connect(self.clean_protein_update)
        self.ui.btn_edit_project_delete.clicked.connect(self.delete_protein)
        self.ui.btn_edit_existing_protein_struct.clicked.connect(self.add_existing_protein)
        self.ui.btn_edit_project_save.clicked.connect(self.save_selected_protein_structure_as_pdb_file)
        # </editor-fold>

        # <editor-fold desc="View project page">
        self.ui.btn_view_project_show.clicked.connect(self.view_sequence)
        self.ui.btn_view_project_show_structure.clicked.connect(self.view_structure)
        self.ui.list_view_project_proteins.doubleClicked.connect(self.view_sequence)
        self.ui.list_view_project_proteins.itemClicked.connect(self.view_show_options)
        # </editor-fold>

        # <editor-fold desc="Use project page">
        self.ui.txt_use_project_name.textChanged.connect(self.validate_use_project_name)
        self.ui.btn_use_next.clicked.connect(self.show_protein_selection_for_use)
        self.ui.txt_use_search.textChanged.connect(self.validate_use_search)
        self.ui.btn_use_add_available_protein_structures.clicked.connect(self.add_protein_structure_to_new_project)
        self.ui.list_use_available_protein_structures.doubleClicked.connect(self.add_protein_structure_to_new_project)
        self.ui.list_use_available_protein_structures.itemClicked.connect(self.use_enable_add)
        self.ui.btn_use_remove_selected_protein_structures.clicked.connect(self.remove_protein_structure_to_new_project)
        self.ui.list_use_selected_protein_structures.doubleClicked.connect(self.remove_protein_structure_to_new_project)
        self.ui.list_use_selected_protein_structures.itemClicked.connect(self.use_enable_remove)
        self.ui.btn_use_back.clicked.connect(self.hide_protein_selection_for_use)
        self.ui.btn_use_create_new_project.clicked.connect(self.pre_create_use_project)

        # </editor-fold>

        # <editor-fold desc="ESMFold Monomer Prediction page">
        self.ui.btn_esm_seq_to_predict.clicked.connect(self.cloud_esm_add_seq_to_predict)
        self.ui.btn_esm_seq_to_predict_remove.clicked.connect(self.cloud_esm_remove)
        self.ui.btn_esm_next.clicked.connect(self.cloud_esm_next)
        self.ui.btn_esm_back.clicked.connect(self.cloud_esm_back)
        self.ui.btn_esm_next_2.clicked.connect(self.cloud_esm_add_protein)
        self.ui.btn_esm_back_2.clicked.connect(self.cloud_esm_back_2)
        self.ui.txt_esm_prot_name.textChanged.connect(self.cloud_esm_validate_protein_name)
        self.ui.txt_esm_prot_seq.textChanged.connect(self.cloud_esm_validate_protein_sequence)
        self.ui.btn_esm_predict.clicked.connect(self.predict_esm_monomer)

        self.ui.table_esm_prot_to_predict.itemSelectionChanged.connect(self.cloud_esm_item_changed)
        # </editor-fold>

        # <editor-fold desc="Monomer local prediction page">
        self.ui.btn_pred_mono_seq_to_predict.clicked.connect(self.local_pred_mono_add_seq_to_predict)
        self.ui.btn_pred_mono_seq_to_predict_remove.clicked.connect(self.local_pred_mono_remove)
        self.ui.btn_pred_mono_next.clicked.connect(self.local_pred_mono_next)
        self.ui.btn_pred_mono_back.clicked.connect(self.local_pred_mono_back)
        self.ui.btn_pred_mono_add_protein.clicked.connect(self.local_pred_mono_add_protein)
        self.ui.btn_pred_mono_back_2.clicked.connect(self.local_pred_mono_back_2)
        self.ui.txt_pred_mono_prot_name.textChanged.connect(self.local_pred_mono_validate_protein_name)
        self.ui.txt_pred_mono_seq_name.textChanged.connect(self.local_pred_mono_validate_protein_sequence)
        self.ui.btn_pred_mono_advanced_config.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_mono_predict.clicked.connect(self.predict_local_monomer)

        self.ui.table_pred_mono_prot_to_predict.itemSelectionChanged.connect(self.local_pred_mono_item_changed)
        # </editor-fold>

        # <editor-fold desc="Multimer prediction page">
        self.ui.btn_pred_multi_prot_to_predict_add.clicked.connect(self.local_pred_multi_add)
        self.ui.btn_pred_multi_prot_to_predict_remove.clicked.connect(self.local_pred_multi_remove)
        self.ui.btn_pred_multi_next.clicked.connect(self.local_pred_multi_next)
        self.ui.btn_pred_multi_back.clicked.connect(self.local_pred_multi_back)
        self.ui.btn_pred_multi_prot_to_predict_add_2.clicked.connect(self.local_pred_multi_prot_to_predict_add_2)
        self.ui.btn_pred_multi_back_2.clicked.connect(self.local_pred_multi_back_2)
        self.ui.txt_pred_multi_prot_name.textChanged.connect(self.local_pred_multi_validate_protein_name)
        self.ui.txt_pred_multi_prot_seq.textChanged.connect(self.local_pred_multi_validate_protein_sequence)
        self.ui.btn_pred_multi_prot_seq_add.clicked.connect(self.local_pred_multi_add_sequence_to_list)
        self.ui.btn_pred_multi_prot_seq_overview_remove.clicked.connect(self.local_pred_multi_remove_sequence_to_list)
        self.ui.btn_pred_multi_advanced_config.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_multi_predict.clicked.connect(self.predict_local_multimer)

        self.ui.list_pred_multi_prot_seq_overview.itemClicked.connect(
            self.local_pred_multi_prot_seq_overview_item_changed,
        )
        self.ui.table_pred_multi_prot_to_predict.itemSelectionChanged.connect(
            self.local_pred_multi_prot_to_predict_item_changed,
        )
        # self.ui.btn_local_pred_multi_single.clicked.connect(self.show_local_pred_multi_stage_protein_name)
        # self.ui.btn_local_pred_multi_back_prediction_mode.clicked.connect(self.show_local_pred_multi_stage_prediction_mode)
        # # single connections
        # self.ui.btn_local_pred_multi_next.clicked.connect(self.show_local_pred_multi_stage_protein_sequence_single)
        # self.ui.btn_local_pred_multi_back.clicked.connect(self.show_local_pred_multi_stage_protein_name)
        # self.ui.btn_local_pred_multi_next_2.clicked.connect(self.show_local_pred_multi_stage_prediction_single)
        # self.ui.btn_local_pred_multi_back_2.clicked.connect(self.show_local_pred_multi_stage_protein_sequence_single)
        # # batch connections
        # self.ui.btn_local_pred_multi_batch.clicked.connect(self.show_local_pred_multi_stage_protein_sequence_batch)
        # #self.ui.btn_local_pred_multi_back_3.clicked.connect(self.hide_protein_sequence_stage_batch)
        # # text fields
        # self.ui.txt_local_pred_multi_protein_name.textChanged.connect(self.validate_local_pred_multi)
        # self.ui.txt_local_pred_multi_prot_seq.textChanged.connect(self.validate_local_pred_multi)
        # single analysis page
        # self.ui.btn_analysis_next.clicked.connect(self.show_single_analysis_stage_1)
        # self.ui.btn_analysis_next_2.clicked.connect(self.show_single_analysis_stage_2)
        # self.ui.btn_analysis_back.clicked.connect(self.show_single_analysis_stage_0)
        # self.ui.btn_analysis_back_2.clicked.connect(self.show_single_analysis_stage_1)
        # self.ui.btn_analysis_start.clicked.connect(self.start_process)
        # self.ui.box_analysis_prot_struct_1.currentIndexChanged.connect(self.check_if_prot_structs_are_filled)
        # self.ui.box_analysis_prot_struct_2.currentIndexChanged.connect(self.check_if_prot_structs_are_filled)
        # self.ui.list_analysis_ref_chains.itemSelectionChanged.connect(self.count_selected_chains_for_prot_struct_1)
        # self.ui.list_analysis_model_chains.itemSelectionChanged.connect(self.check_if_same_no_of_chains_selected)
        # </editor-fold>

        # <editor-fold desc="Monomer Prediction + Analysis page">
        # <editor-fold desc="Prediction section">
        self.ui.btn_pred_analysis_mono_seq_to_predict.clicked.connect(self.mono_pred_analysis_add_seq_to_predict)
        self.ui.btn_pred_analysis_mono_seq_to_predict_remove.clicked.connect(
            self.mono_pred_analysis_remove_protein_to_predict,
        )
        self.ui.btn_pred_analysis_mono_next.clicked.connect(self.mono_pred_analysis_next)
        self.ui.btn_pred_analysis_mono_back.clicked.connect(self.mono_pred_analysis_back)
        self.ui.btn_pred_analysis_mono_add_protein.clicked.connect(self.mono_pred_analysis_add_protein)
        self.ui.btn_pred_analysis_mono_back_2.clicked.connect(self.mono_pred_analysis_back_2)
        self.ui.txt_pred_analysis_mono_prot_name.textChanged.connect(self.mono_pred_analysis_validate_protein_name)
        self.ui.txt_pred_analysis_mono_seq_name.textChanged.connect(self.mono_pred_analysis_validate_protein_sequence)
        self.ui.btn_pred_mono_advanced_config_2.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_analysis_mono_go_analysis_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self.ui.table_pred_analysis_mono_prot_to_predict.itemSelectionChanged.connect(
            self.mono_pred_analysis_prediction_overview_item_clicked,
        )

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self.ui.btn_pred_analysis_mono_add.clicked.connect(self.mono_pred_analysis_structure_analysis_add)
        self.ui.btn_pred_analysis_mono_remove.clicked.connect(self.remove_mono_pred_analysis_analysis_run)
        self.ui.btn_pred_analysis_mono_back_3.clicked.connect(self.mono_pred_analysis_structure_analysis_back_3)
        self.ui.btn_pred_analysis_mono_next_2.clicked.connect(self.mono_pred_analysis_structure_analysis_next_2)
        self.ui.btn_pred_analysis_mono_back_4.clicked.connect(self.mono_pred_analysis_structure_analysis_back_4)
        self.ui.btn_pred_analysis_mono_next_3.clicked.connect(self.mono_pred_analysis_structure_analysis_next_3)
        self.ui.btn_pred_analysis_mono_back_5.clicked.connect(self.mono_pred_analysis_structure_analysis_back_5)
        self.ui.btn_pred_analysis_mono_next_4.clicked.connect(self.mono_pred_analysis_structure_analysis_next_4)
        self.ui.box_pred_analysis_mono_prot_struct_1.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled,
        )
        self.ui.box_pred_analysis_mono_prot_struct_2.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled,
        )
        self.ui.list_pred_analysis_mono_ref_chains.itemSelectionChanged.connect(
            self.count_mono_pred_analysis_selected_chains_for_prot_struct_1,
        )
        self.ui.list_pred_analysis_mono_model_chains.itemSelectionChanged.connect(
            self.check_mono_pred_analysis_if_same_no_of_chains_selected,
        )
        self.ui.btn_pred_analysis_mono_back_pred_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self.ui.btn_pred_analysis_mono_start.clicked.connect(self.start_monomer_prediction_analysis)
        self.ui.list_pred_analysis_mono_overview.clicked.connect(
            self.mono_pred_analysis_structure_analysis_overview_clicked,
        )

        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Multimer Prediction + Analysis page">
        # <editor-fold desc="Prediction section">
        self.ui.btn_pred_analysis_multi_prot_to_predict_add.clicked.connect(self.multi_pred_analysis_add)
        self.ui.btn_pred_analysis_multi_prot_to_predict_remove.clicked.connect(
            self.multi_pred_analysis_remove_protein_to_predict,
        )
        self.ui.btn_pred_analysis_multi_next.clicked.connect(self.multi_pred_analysis_next)
        self.ui.btn_pred_analysis_multi_back.clicked.connect(self.multi_pred_analysis_back)
        self.ui.btn_pred_analysis_multi_prot_seq_add.clicked.connect(self.multi_pred_analysis_add_sequence_to_list)
        self.ui.btn_pred_analysis_multi_prot_seq_overview_remove.clicked.connect(
            self.multi_pred_analysis_remove_sequence_to_list,
        )
        self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.clicked.connect(
            self.multi_pred_analysis_add_protein_to_predict,
        )
        self.ui.btn_pred_analysis_multi_back_2.clicked.connect(self.multi_pred_analysis_back_2)

        self.ui.txt_pred_analysis_multi_prot_name.textChanged.connect(self.multi_pred_analysis_validate_protein_name)
        self.ui.txt_pred_analysis_multi_prot_seq.textChanged.connect(self.multi_pred_analysis_validate_protein_sequence)
        self.ui.btn_pred_analysis_multi_advanced_config.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_analysis_multi_go_analysis_setup.clicked.connect(self.switch_multimer_pred_analysis_tab)

        self.ui.list_pred_analysis_multi_prot_seq_overview.clicked.connect(
            self.multi_pred_analysis_prot_seq_overview_item_changed,
        )
        self.ui.table_pred_analysis_multi_prot_to_predict.itemSelectionChanged.connect(
            self.multi_pred_analysis_prot_to_predict_item_changed,
        )

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self.ui.btn_pred_analysis_multi_add.clicked.connect(self.multi_pred_analysis_structure_analysis_add)
        self.ui.btn_pred_analysis_multi_remove.clicked.connect(self.remove_multi_pred_analysis_analysis_run)
        self.ui.btn_pred_analysis_multi_back_3.clicked.connect(self.multi_pred_analysis_structure_analysis_back_3)
        self.ui.btn_pred_analysis_multi_next_2.clicked.connect(self.multi_pred_analysis_structure_analysis_next_2)
        self.ui.btn_pred_analysis_multi_back_4.clicked.connect(self.multi_pred_analysis_structure_analysis_back_4)
        self.ui.btn_pred_analysis_multi_next_3.clicked.connect(self.multi_pred_analysis_structure_analysis_next_3)
        self.ui.btn_pred_analysis_multi_back_5.clicked.connect(self.multi_pred_analysis_structure_analysis_back_5)
        self.ui.btn_pred_analysis_multi_next_4.clicked.connect(self.multi_pred_analysis_structure_analysis_next_4)
        self.ui.box_pred_analysis_multi_prot_struct_1.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled,
        )
        self.ui.box_pred_analysis_multi_prot_struct_2.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled,
        )
        self.ui.list_pred_analysis_multi_ref_chains.itemSelectionChanged.connect(
            self.count_multi_pred_analysis_selected_chains_for_prot_struct_1,
        )
        self.ui.list_pred_analysis_multi_model_chains.itemSelectionChanged.connect(
            self.check_multi_pred_analysis_if_same_no_of_chains_selected,
        )
        self.ui.btn_pred_analysis_multi_back_pred_setup.clicked.connect(self.switch_multimer_pred_analysis_tab)
        self.ui.btn_pred_analysis_multi_start.clicked.connect(self.start_multimer_prediction_analysis)
        self.ui.list_pred_analysis_multi_overview.clicked.connect(
            self.multi_pred_analysis_structure_analysis_overview_clicked,
        )
        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Batch analysis page">
        self.ui.btn_analysis_batch_add.clicked.connect(self.structure_analysis_add)
        self.ui.btn_analysis_batch_remove.clicked.connect(self.remove_analysis_run)
        self.ui.btn_analysis_batch_back.clicked.connect(self.structure_analysis_back)
        self.ui.btn_analysis_batch_next.clicked.connect(self.structure_analysis_next)
        self.ui.btn_analysis_batch_back_2.clicked.connect(self.structure_analysis_back_2)
        self.ui.btn_analysis_batch_next_2.clicked.connect(self.structure_analysis_next_2)
        self.ui.btn_analysis_batch_back_3.clicked.connect(self.structure_analysis_back_3)
        self.ui.btn_analysis_batch_next_3.clicked.connect(self.structure_analysis_next_3)
        self.ui.box_analysis_batch_prot_struct_1.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch,
        )
        self.ui.box_analysis_batch_prot_struct_2.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch,
        )
        self.ui.list_analysis_batch_ref_chains.itemSelectionChanged.connect(
            self.count_batch_selected_chains_for_prot_struct_1,
        )
        self.ui.list_analysis_batch_model_chains.itemSelectionChanged.connect(
            self.check_if_same_no_of_chains_selected_batch,
        )
        self.ui.btn_analysis_batch_start.clicked.connect(self.start_process_batch)
        self.ui.list_analysis_batch_overview.clicked.connect(self.structure_analysis_overview_clicked)

        # </editor-fold>

        # <editor-fold desc="Analysis images page">
        self.ui.btn_add_analysis_images_struct_analysis.clicked.connect(self.add_protein_pair_to_image_creation_queue)
        self.ui.list_analysis_images_struct_analysis.doubleClicked.connect(
            self.add_protein_pair_to_image_creation_queue,
        )
        self.ui.list_analysis_images_struct_analysis.clicked.connect(self.analysis_images_enable_add)
        self.ui.btn_remove_analysis_images_creation_struct_analysis.clicked.connect(
            self.remove_protein_pair_from_image_creation_queue,
        )
        self.ui.list_analysis_images_creation_struct_analysis.doubleClicked.connect(
            self.remove_protein_pair_from_image_creation_queue,
        )
        self.ui.list_analysis_images_creation_struct_analysis.clicked.connect(self.analysis_images_enable_remove)
        self.ui.btn_start_automatic_image_creation.clicked.connect(self.start_automatic_image_creation)

        # </editor-fold>

        # <editor-fold desc="Results page">
        self.ui.cb_results_analysis_options.currentIndexChanged.connect(self.pre_load_results)
        self.ui.btn_color_rmsd.clicked.connect(self.color_protein_pair_by_rmsd)
        self.ui.btn_view_struct_alignment.clicked.connect(self.display_structure_alignment)
        self.ui.btn_view_distance_plot.clicked.connect(self.display_distance_plot)
        self.ui.btn_view_distance_histogram.clicked.connect(self.display_distance_histogram)
        self.ui.btn_view_interesting_region.clicked.connect(self.display_interesting_region)
        self.ui.btn_view_distance_table.clicked.connect(self.display_distance_table)

        # </editor-fold>

        # <editor-fold desc="Manage">
        self.ui.box_manage_choose_protein.activated.connect(self.choose_manage_open_protein)
        self.ui.box_manage_choose_color.activated.connect(self.choose_manage_color_selected_protein)
        self.ui.box_manage_choose_representation.activated.connect(self.choose_manage_representation)
        self.ui.box_manage_choose_bg_color.activated.connect(self.choose_manage_bg_color)
        self.ui.btn_disulfid_bond_show.clicked.connect(self.show_disulfid_bonds_as_sticks)
        self.ui.btn_disulfid_bond_hide.clicked.connect(self.hide_disulfid_bonds_as_sticks)
        
        # </editor-fold>

        # <editor-fold desc="Image page">
        self.ui.btn_update_scene.clicked.connect(self.update_scene)
        self.ui.btn_save_scene.clicked.connect(self.save_scene)
        self.ui.btn_save_image.clicked.connect(self.save_image)
        self.ui.btn_preview_image.clicked.connect(self.preview_image)
        self.ui.box_representation.activated.connect(self.show_representation)
        self.ui.box_bg_color.activated.connect(self.choose_bg_color)
        self.ui.box_renderer.activated.connect(self.choose_renderer)
        self.ui.box_ray_trace_mode.activated.connect(self.choose_ray_trace_mode)
        self.ui.box_ray_texture.activated.connect(self.choose_ray_texture)
        self.ui.cb_transparent_bg.stateChanged.connect(self.decide_transparent_bg)
        self.ui.cb_ray_tracing.stateChanged.connect(self.decide_ray_tracing)

        # </editor-fold>

        # <editor-fold desc="Hotspots page">
        self.ui.list_hotspots_choose_protein.currentItemChanged.connect(self.open_protein)
        self.ui.btn_hotspots_resi_show.clicked.connect(self.show_resi_sticks)
        self.ui.btn_hotspots_resi_hide.clicked.connect(self.hide_resi_sticks)
        self.ui.btn_hotspots_resi_zoom.clicked.connect(self.zoom_resi_position)

        # </editor-fold>

    def _create_all_tooltips(self) -> None:
        """Creates all tooltips for the gui elements."""
        self.status_bar.setToolTip("Status information: Current process")
        # new project page
        self.ui.btn_new_choose_reference.setToolTip("Click to add a .pdb file")
        # sidebar
        self.ui.lbl_current_project_name.setToolTip("Name of the current project")
        # edit page
        self.ui.btn_edit_project_save.setToolTip("Save as a .pdb file")
        # view page
        # fixme: is this important? self.ui.list_view_project_proteins.setToolTip("Proteins of the current project")
        self.ui.txtedit_view_sequence.setToolTip("Protein sequence of the selected protein")
        # use page
        self.ui.txt_use_search.setToolTip("Enter a protein name to search in your current workspace")
        # prediction Monomer
        self.ui.table_pred_mono_prot_to_predict.setToolTip("Protein monomers which get predicted")
        self.ui.btn_pred_mono_seq_to_predict.setToolTip("Set up a protein which can be used for a prediction")
        self.ui.table_pred_analysis_mono_prot_to_predict.setToolTip("Protein monomers which get predicted")
        self.ui.list_pred_analysis_mono_overview.setToolTip("Protein pairs which get analyzed")
        self.ui.btn_pred_analysis_mono_seq_to_predict.setToolTip("Set up a protein which can be used for a prediction")
        # prediction Multimer
        self.ui.table_pred_multi_prot_to_predict.setToolTip("Protein multimers which get predicted")
        self.ui.btn_pred_multi_prot_to_predict_add.setToolTip("Set up a protein which can be used for a prediction")
        self.ui.table_pred_analysis_multi_prot_to_predict.setToolTip("Protein multimers which get predicted")
        self.ui.list_pred_analysis_multi_overview.setToolTip("Protein pairs which get analyzed")
        self.ui.btn_pred_analysis_multi_prot_to_predict_add.setToolTip(
            "Set up a protein which can be used for a prediction",
        )
        # image page
        self.ui.btn_save_scene.setToolTip("Create new PyMOL scene")
        self.ui.btn_update_scene.setToolTip("Overwrite current scene")
        self.ui.btn_preview_image.setToolTip("Preview current viewpoint")
        self.ui.btn_save_image.setToolTip("Save current viewpoint as png file")
        self.ui.cb_ray_tracing.setToolTip("Enable ray-tracing")
        self.ui.cb_transparent_bg.setToolTip("Enable transparent background")
        self.ui.box_representation.setToolTip("Choose a representation")
        self.ui.box_bg_color.setToolTip("Choose a background color")
        self.ui.box_renderer.setToolTip("Choose a ray-tracing renderer")
        self.ui.box_ray_trace_mode.setToolTip("Choose a ray-trace mode")

    # <editor-fold desc="Page init functions">
    def _init_fill_combo_boxes(self) -> None:
        """Fills all combo boxes of the plugin."""
        item_list_representation = [
            "",
            "cartoon",
            "ribbon",
        ]
        gui_utils.fill_combo_box(self.ui.box_representation, item_list_representation)
        gui_utils.fill_combo_box(self.ui.box_manage_choose_representation, item_list_representation)
        # combo box BgColor
        item_list_bg_color = [
            "",
            "black",
            "white",
        ]
        gui_utils.fill_combo_box(self.ui.box_bg_color, item_list_bg_color)
        gui_utils.fill_combo_box(self.ui.box_manage_choose_bg_color, item_list_bg_color)
        # combo box Renderer
        item_list_renderer = [
            "",
            "default renderer",
            "PyMOL internal renderer",
        ]
        gui_utils.fill_combo_box(self.ui.box_renderer, item_list_renderer)
        # combo box RayTraceMode
        item_list_ray_trace_mode = [
            "",
            "normal color",
            "normal color + black outline",
            "black outline only",
            "quantized color + black outline",
        ]
        gui_utils.fill_combo_box(self.ui.box_ray_trace_mode, item_list_ray_trace_mode)
        # combo box Ray Texture
        item_list_ray_texture = [
            "",
            "None",
            "Matte 1",
            "Matte 2",
            "Swirl 1",
            "Fiber",
        ]
        gui_utils.fill_combo_box(self.ui.box_ray_texture, item_list_ray_texture)
        # combo box protein Colors
        gui_utils.fill_combo_box(self.ui.box_manage_choose_color, constants.PYMOL_COLORS)

    def _init_new_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self.ui.txt_new_project_name.clear()
        self.ui.txt_new_choose_reference.clear()
        self.ui.lbl_new_status_project_name.setText("")
        self.ui.lbl_new_status_choose_reference.setText("")
        self.ui.cb_new_add_reference.setCheckState(0)
        self.ui.btn_new_create_project.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_new_create_project)

    def _init_use_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        gui_elements = [
            self.ui.lbl_use_search,
            self.ui.lbl_use_status_search,
            self.ui.txt_use_search,
            self.ui.btn_use_add_available_protein_structures,
            self.ui.lbl_use_available_protein_structures,
            self.ui.list_use_available_protein_structures,
            self.ui.btn_use_remove_selected_protein_structures,
            self.ui.lbl_use_selected_protein_structures,
            self.ui.list_use_selected_protein_structures,
            self.ui.btn_use_back,
            self.ui.btn_use_create_new_project,
            self.ui.lbl_use_new_project,
        ]
        gui_utils.hide_gui_elements(gui_elements)
        self.ui.txt_use_project_name.clear()
        self.ui.lbl_use_status_project_name.setText("")
        self.ui.txt_use_search.clear()
        self.ui.lbl_use_status_search.setText("")
        self.ui.list_use_available_protein_structures.clear()
        self.ui.list_use_selected_protein_structures.clear()
        self.ui.list_use_existing_projects.clear()
        self.ui.btn_use_next.setEnabled(False)
        self.hide_protein_selection_for_use()

    def _init_edit_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self.ui.list_edit_project_proteins.clear()
        gui_elements_to_hide = [
            self.ui.lbl_edit_clean_new_prot,
            self.ui.btn_edit_clean_new_prot,
            self.ui.lbl_edit_clean_update_prot,
            self.ui.btn_edit_clean_update_prot,
            self.ui.label_12,
            self.ui.btn_edit_project_delete,
            self.ui.btn_edit_project_save,
            self.ui.label_13,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_edit_project_proteins)
        # self.project_scanner.scan_project_for_valid_proteins(self.ui.list_edit_project_proteins)

    def _init_sequence_vs_pdb_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self.ui.list_s_v_p_ref_chains.clear()
        # # sets up defaults: Prediction + Analysis
        #
        # #self.ui.label_18.hide()
        # #self.ui.txt_prediction_project_name.hide()
        # #self.ui.lbl_prediction_status_project_name.hide()
        # # stage 1
        # self.ui.lbl_prediction_status_project_name.setText("")
        #
        # #self.ui.list_widget_projects.hide()
        # #self.ui.btn_prediction_next_1.hide()
        # # stage 2
        # self.ui.lbl_prediction_load_reference.hide()
        # self.ui.txt_prediction_load_reference.clear()
        # self.ui.txt_prediction_load_reference.hide()
        # self.ui.lbl_prediction_status_load_reference.setText("")
        # self.ui.btn_prediction_back_2.setEnabled(False)
        # self.ui.btn_prediction_back_2.hide()
        # self.ui.btn_prediction_next_2.setEnabled(False)
        # self.ui.btn_prediction_next_2.hide()
        #
        # # stage 3
        # self.ui.lbl_prediction_ref_chains.hide()
        # self.ui.list_widget_ref_chains.hide()
        # # TO-DO: stylesheet needs to be integrated into the styles.css
        # self.ui.list_widget_ref_chains.setStyleSheet("border-style: solid;"
        #                                              "border-width: 2px;"
        #                                              "border-radius: 8px;"
        #                                              "border-color: #DCDBE3;")
        # self.ui.list_widget_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        # self.ui.btn_prediction_back_3.setEnabled(False)
        # self.ui.btn_prediction_back_3.hide()
        # self.ui.btn_prediction_start.setEnabled(False)
        # self.ui.btn_prediction_start.hide()
        # # TO-DO: needs to be removed if model chains are implemented at the right spot
        # self.ui.lbl_prediction_model_chains.hide()
        # self.ui.txt_prediction_chain_model.hide()

    def _init_single_analysis_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self.single_analysis_management.show_stage_x(0)
        self.single_analysis_management.disable_all_next_buttons()
        self.single_analysis_management.set_empty_string_in_label()
        self.ui.lbl_analysis_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_analysis_prot_struct_2.setText("Protein structure 2")
        self.ui.box_analysis_prot_struct_1.setCurrentIndex(0)
        self.ui.box_analysis_prot_struct_2.setCurrentIndex(0)

    def _init_results_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        # stage 1
        self.ui.list_results_interest_regions.clear()
        self.ui.txt_results_rmsd.clear()
        self.ui.txt_results_aligned_residues.clear()

    def _init_analysis_image_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self.ui.list_analysis_images_struct_analysis.setEnabled(True)
        self.ui.list_analysis_images_creation_struct_analysis.setEnabled(True)
        self.display_image_analysis_page()

    def _init_image_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self.ui.box_representation.setCurrentIndex(0)
        self.ui.box_bg_color.setCurrentIndex(0)
        self.ui.box_renderer.setCurrentIndex(0)
        self.ui.box_ray_trace_mode.setCurrentIndex(0)
        self.ui.box_ray_texture.setCurrentIndex(0)
        self.ui.cb_ray_tracing.setChecked(False)
        self.ui.cb_transparent_bg.setChecked(False)
        self.ui.label_10.hide()
        self.ui.box_ray_trace_mode.hide()
        self.ui.label_14.hide()
        self.ui.box_ray_texture.hide()

    def _init_batch_analysis_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        # sets up defaults: Batch
        self.batch_analysis_management.show_stage_x(0)
        self.ui.list_analysis_batch_overview.clear()
        self.ui.btn_analysis_batch_remove.hide()

    def _init_all_pages(self) -> None:
        """Collection of all init methods to reset all page defaults."""
        self._init_local_pred_mono_page()
        self._init_local_pred_multi_page()
        self._init_mono_pred_analysis_page()
        self._init_multi_pred_analysis_page()
        self._init_sequence_vs_pdb_page()
        self._init_single_analysis_page()
        self._init_batch_analysis_page()
        self._init_analysis_image_page()
        self._init_results_page()
        self._init_image_page()

    # </editor-fold>

    # <editor-fold desc="Display page functions">
    def display_home_page(self) -> None:
        """Displays the homepage of the plugin."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 0, "Home")

    def display_single_analysis_page(self) -> None:
        """Displays the single analysis work area."""
        self.fill_protein_structure_boxes()
        self.ui.list_analysis_ref_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.list_analysis_model_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        # regular work area opening
        self._init_single_analysis_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 3, "Single Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_single_analysis_page,
        )

    def display_job_analysis_page(self) -> None:
        """Displays the job analysis work area."""
        self.ui.list_analysis_batch_ref_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.list_analysis_batch_model_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        # regular work area opening
        self._init_batch_analysis_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 4, "Structure Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_batch_analysis_page,
        )

    def display_results_page(self) -> None:
        """Displays the results work area."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        results = []
        results.insert(0, "")
        i = 0
        for tmp_protein_pair in self.app_project.protein_pairs:
            results.append(tmp_protein_pair.name)
            if tmp_protein_pair.name == self.results_name:
                pass
            i += 1
        self.ui.cb_results_analysis_options.clear()
        gui_utils.fill_combo_box(self.ui.cb_results_analysis_options, results)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 5, "Results")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_results_page)
        # self.ui.cb_results_analysis_options.setCurrentIndex(current_results_index + 1)

    def display_image_page(self) -> None:
        """Displays the image work area."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self._init_image_page()
        if self.ui.box_renderer.currentText() == "":
            self.ui.cb_ray_tracing.hide()
            self.ui.label_26.hide()
        else:
            self.ui.cb_ray_tracing.show()
            self.ui.label_26.show()
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_image_page)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 6, "Image")

    def display_new_page(self) -> None:
        """Displays the new project work area."""
        self._init_new_page()
        self.ui.list_new_projects.clear()
        # pre-process
        gui_elements_to_hide = [
            self.ui.lbl_new_choose_reference,
            self.ui.txt_new_choose_reference,
            self.ui.btn_new_choose_reference,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.status_bar.showMessage(self.workspace.text())
        tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_new_projects)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 7, "Create new project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_new_page)

    def display_open_page(self) -> None:
        """Displays the open project work area."""
        self.ui.txt_open_search.clear()
        self.ui.txt_open_selected_project.clear()
        if safeguard.Safeguard.check_filepath(self.workspace_path):
            self.ui.list_open_projects.clear()
            # pre-process
            self.status_bar.showMessage(self.workspace.text())
            try:
                tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_open_projects)
            except PermissionError:
                gui_utils.error_dialog_settings(
                    "The settings file is corrupted. Please restore the settings!", "", self.app_settings,
                )
                self.display_home_page()
                return
            tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 8, "Open existing project")
            self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_open_page)
        else:
            gui_utils.error_dialog_settings(
                "The settings file is corrupted. Please restore the settings!", "", self.app_settings,
            )
            self.display_home_page()

    def display_delete_page(self) -> None:
        """Displays the "delete" project work area."""
        self.ui.txt_delete_search.clear()
        self.ui.txt_delete_selected_projects.clear()
        self.ui.list_delete_projects.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_delete_projects)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 9, "Delete existing project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_delete_page)

    def display_edit_page(self) -> None:
        """Displays the edit project page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        self._init_edit_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 13, "Edit proteins of current project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_edit_page)
        if len(self.ui.list_edit_project_proteins) > 2:
            self.ui.btn_edit_existing_protein_struct.hide()
            self.ui.lbl_edit_existing_protein_struct.hide()
        else:
            self.ui.btn_edit_existing_protein_struct.show()
            self.ui.lbl_edit_existing_protein_struct.show()

    def start_display_use_page(self) -> None:
        """Displays the use project page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
    
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)

        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.Worker(self.workspace_path)
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_use_page)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        self._init_use_page()
        self.ui.list_use_available_protein_structures.clear()
        self.ui.list_use_selected_protein_structures.clear()
        tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_use_existing_projects)
        gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_use_selected_protein_structures)

    def display_use_page(self, return_value) -> None:
        """Displays the use project page, after cpu intense task (post thread method).

        Args:
            return_value: the value which gets returned from the thread process
        """
        # this for-loop is necessary for eliminating all proteins which are in the current project from the ones which
        # are available
        for i in range(self.ui.list_use_selected_protein_structures.count()):
            self.ui.list_use_selected_protein_structures.setCurrentRow(i)
            tmp_prot_name = self.ui.list_use_selected_protein_structures.currentItem().text()
            if tmp_prot_name in return_value[1]:
                return_value[1].remove(tmp_prot_name)

        self.ui.list_use_available_protein_structures.addItems(return_value[1])

        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 14, "Use existing project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_use_page)
        QtWidgets.QApplication.restoreOverrideCursor()
        # tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 14, "Use existing project")
        # self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
        #                                                         self.ui.btn_use_page)
        # QtWidgets.QApplication.restoreOverrideCursor()

        # QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        # self._init_use_page()
        # self.ui.list_use_available_protein_structures.clear()
        # self.ui.list_use_selected_protein_structures.clear()
        # valid_projects = tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_use_existing_projects)
        # # filesystem operations
        # gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_use_selected_protein_structures)
        # # self.project_scanner.scan_project_for_valid_proteins(self.ui.list_use_selected_protein_structures)
        # protein_dict, protein_names = tools.scan_workspace_for_non_duplicate_proteins(self.workspace_path)
        # global_variables.global_var_workspace_proteins = protein_dict
        # # this for-loop is necessary for eliminating all proteins which are in the current project from the ones which
        # # are available
        # for i in range(self.ui.list_use_selected_protein_structures.count()):
        #     self.ui.list_use_selected_protein_structures.setCurrentRow(i)
        #     tmp_prot_name = self.ui.list_use_selected_protein_structures.currentItem().text()
        #     if tmp_prot_name in protein_names:
        #         protein_names.remove(tmp_prot_name)
        #
        # self.ui.list_use_available_protein_structures.addItems(protein_names)
        # tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 14, "Use existing project")
        # self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
        #                                                         self.ui.btn_use_page)
        # QtWidgets.QApplication.restoreOverrideCursor()

    def display_hotspots_page(self) -> None:
        """Displays the hotspots page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        try:
            print(self.ui.list_hotspots_choose_protein.currentItem().text())
        except AttributeError:
            self.ui.list_hotspots_choose_protein.clear()
            gui_elements_to_hide = [
                self.ui.lbl_hotspots_resi_no,
                self.ui.sp_hotspots_resi_no,
                self.ui.lbl_hotspots_resi_show,
                self.ui.btn_hotspots_resi_show,
                self.ui.lbl_hotspots_resi_hide,
                self.ui.btn_hotspots_resi_hide,
                self.ui.lbl_hotspots_resi_zoom,
                self.ui.btn_hotspots_resi_zoom,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_hotspots_choose_protein)
            gui_utils.fill_list_view_with_protein_pair_names(self.app_project, self.ui.list_hotspots_choose_protein)
        # self.project_scanner.scan_project_for_valid_proteins(self.ui.list_hotspots_choose_protein)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 18, "Hotspots")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_hotspots_page)

    def display_manage_pymol_session(self) -> None:
        """Displays the manage pymol session page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self.ui.box_manage_choose_protein.clear()
        pymol_objs = cmd.get_object_list()
        pymol_objs.insert(0, "")
        for tmp_object in pymol_objs:
            self.ui.box_manage_choose_protein.addItem(tmp_object)
        self.ui.box_manage_choose_protein.setCurrentIndex(self.pymol_session_specs[pyssa_keys.SESSION_SPEC_PROTEIN][0])
        self.ui.box_manage_choose_color.setCurrentIndex(self.pymol_session_specs[pyssa_keys.SESSION_SPEC_COLOR][0])
        self.ui.box_manage_choose_representation.setCurrentIndex(
            self.pymol_session_specs[pyssa_keys.SESSION_SPEC_REPRESENTATION][0],
        )
        self.ui.box_manage_choose_bg_color.setCurrentIndex(
            self.pymol_session_specs[pyssa_keys.SESSION_SPEC_BG_COLOR][0],
        )
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 24, "Manage PyMOL session")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_manage_session)
        gui_elements_to_hide = [
            self.ui.lbl_manage_choose_color,
            self.ui.box_manage_choose_color,
            self.ui.lbl_manage_choose_representation,
            self.ui.box_manage_choose_representation,
            self.ui.lbl_manage_choose_bg_color,
            self.ui.box_manage_choose_bg_color,
            self.ui.lbl_disulfid_bond_1,
            self.ui.lbl_disulfid_bond_2,
            self.ui.btn_disulfid_bond_show,
            self.ui.btn_disulfid_bond_hide,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    # </editor-fold>

    def restore_settings(self) -> None:
        """Restores the settings.xml file to the default values."""
        out = gui_utils.warning_dialog_restore_settings("Are you sure you want to restore all settings?")
        if out:
            tools.restore_default_settings(self.app_settings)
            self.status_bar.showMessage("Settings were successfully restored.")
            logging.info("Settings were successfully restored.")
        else:
            self.status_bar.showMessage("Settings were not modified.")
            logging.info("Settings were not modified.")

    def quit_app(self) -> None:
        """Closes the entire plugin."""
        self.close()

    def open_settings_global(self) -> None:
        """Opens the dialog for the global settings."""
        dialog = dialog_settings_global.DialogSettingsGlobal()
        dialog.exec_()
        self.app_settings = self.app_settings.deserialize_settings()
        globals.g_settings = self.app_settings
        self.workspace_path = globals.g_settings.workspace_path
        self.workspace = QtWidgets.QLabel(f"Current Workspace: {self.workspace_path}")
        self._setup_statusbar()

    def open_logs(self):
        file_dialog = QtWidgets.QFileDialog()
        log_path = str(constants.LOG_PATH)
        file_dialog.setDirectory(log_path)
        file_path, _ = file_dialog.getOpenFileName(self, "Select a log file to open", "", "LOG File (*.log)")
        if file_path:
            os.startfile(file_path)

    def clear_all_log_files(self):
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
                basic_boxes.ok("Clear log files",
                               "All log files could be deleted.",
                               QtWidgets.QMessageBox.Information)
                constants.PYSSA_LOGGER.info("All log files were deleted.")
            else:
                basic_boxes.ok("Clear log files",
                               "Not all log files could be deleted.",
                               QtWidgets.QMessageBox.Warning)
                constants.PYSSA_LOGGER.warning("Not all log files were deleted!")

    @staticmethod
    def open_tutorial() -> None:
        """Opens the official tutorial pdf file."""
        # os.startfile(constants.TUTORIAL_PATH)
        tmp_dialog = TutorialVideosDialog()
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

    def open_page_information(self) -> None:
        """Opens the message box, to display extra information based on the page."""
        if self.ui.lbl_page_title.text() == "Home":
            with open(f'{constants.HELP_HOME_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Create new project":
            with open(f'{constants.HELP_CREATE_NEW_PROJECT_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Open existing project":
            with open(f'{constants.HELP_OPEN_EXISTING_PROJECT_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Delete existing project":
            with open(f'{constants.HELP_DELETE_PROJECT_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Edit proteins of current project":
            with open(f'{constants.HELP_EDIT_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "View proteins of current project":
            with open(f'{constants.HELP_VIEW_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Use existing project":
            with open(f'{constants.HELP_USE_EXISTING_PROJECT_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Structure Analysis":
            with open(f'{constants.HELP_STRUCTURE_ANALYSIS_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Results":
            with open(f'{constants.HELP_RESULTS_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Analysis Images":
            with open(f'{constants.HELP_ANALYSIS_IMAGES_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Image":
            with open(f'{constants.HELP_IMAGE_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Hotspots":
            with open(f'{constants.HELP_HOTSPOTS_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Manage PyMOL session":
            with open(f'{constants.HELP_MANAGE_PYMOL_SESSION_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Local Monomer Prediction":
            with open(f'{constants.HELP_LOCAL_MONOMER_PREDICTION_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Local Multimer Prediction":
            with open(f'{constants.HELP_LOCAL_MULTIMER_PREDICTION_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Monomer Prediction + Analysis":
            with open(f'{constants.HELP_MONOMER_PREDICTION_ANALYSIS_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        elif self.ui.lbl_page_title.text() == "Multimer Prediction + Analysis":
            with open(f'{constants.HELP_MULTIMER_PREDICTION_ANALYSIS_HTML_PATH}', 'r', encoding='utf-8') as file:
                html_content = file.read()
                file.close()
            tmp_dialog = dialog_help.DialogHelp(html_content)
            tmp_dialog.exec_()
        return

    # <editor-fold desc="New project page functions">
    def show_add_reference(self) -> None:
        """Shows the reference input section."""
        # checkbox is checked
        self.ui.cb_new_add_reference.checkState()
        if self.ui.cb_new_add_reference.checkState() == 2:
            self.ui.txt_new_choose_reference.clear()
            self.ui.txt_new_choose_reference.setStyleSheet("background-color: white")
            self.ui.lbl_new_choose_reference.show()
            self.ui.txt_new_choose_reference.show()
            self.ui.btn_new_choose_reference.show()
            self.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_new_create_project)
            # check internet connectivity
            if not tools.check_internet_connectivity():
                gui_utils.no_internet_dialog()
                self.ui.txt_new_choose_reference.setEnabled(False)
                self.ui.lbl_new_status_choose_reference.setText("You cannot enter a PDB ID (no internet).")
                return
            self.ui.txt_new_choose_reference.setEnabled(True)
            self.ui.lbl_new_status_choose_reference.setText("")
        else:
            self.ui.lbl_new_choose_reference.hide()
            self.ui.txt_new_choose_reference.hide()
            self.ui.btn_new_choose_reference.hide()
            self.ui.lbl_new_status_choose_reference.setText("")
            self.ui.btn_new_create_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_new_create_project)

    def load_reference_in_project(self) -> None:
        """Loads a reference in a new project."""
        try:
            # open file dialog
            file_name = QtWidgets.QFileDialog.getOpenFileName(
                self, "Open Reference", QtCore.QDir.homePath(), "PDB Files (*.pdb)",
            )
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self.ui.txt_new_choose_reference.setText(str(file_name[0]))
            self.ui.txt_new_choose_reference.setEnabled(False)
            self.ui.txt_new_choose_reference.setStyleSheet("color: #000000")
            self.ui.btn_new_create_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_new_create_project)
        except ValueError:
            print("No file has been selected.")

    def validate_reference_in_project(self) -> None:
        """Checks if the entered reference protein is valid or not."""
        if len(self.ui.txt_new_choose_reference.text()) == 0:
            self.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
            self.ui.lbl_new_status_choose_reference.setText("")
            self.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_new_create_project)
        elif len(self.ui.txt_new_choose_reference.text()) < 4:
            self.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
            styles.color_button_not_ready(self.ui.btn_new_create_project)
            self.ui.btn_new_create_project.setEnabled(False)
            self.ui.lbl_new_status_choose_reference.setText("")
        # checks if a pdb id was entered
        elif len(self.ui.txt_new_choose_reference.text()) == 4:
            pdb_id = self.ui.txt_new_choose_reference.text().upper()
            try:
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                cmd.fetch(pdb_id, type="pdb", path=self.scratch_path)
                os.remove(f"{self.scratch_path}/{pdb_id}.pdb")
                cmd.reinitialize()
                self.ui.txt_new_choose_reference.setStyleSheet("color: #000000")
                self.ui.btn_new_create_project.setEnabled(True)
                styles.color_button_ready(self.ui.btn_new_create_project)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                styles.color_button_not_ready(self.ui.btn_new_create_project)
                return
            except FileNotFoundError:
                self.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self.ui.lbl_new_status_choose_reference.setText("Invalid PDB ID.")
                self.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_new_create_project)
                return
        else:
            if self.ui.txt_new_choose_reference.text().find("/") == -1:
                self.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_new_create_project)

            elif self.ui.txt_new_choose_reference.text().find("\\") == -1:
                self.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_new_create_project)

    def validate_project_name(self) -> None:
        """Validates the input of the project name in real-time."""
        input_validator.InputValidator.validate_project_name(
            self.ui.list_new_projects,
            self.ui.txt_new_project_name,
            self.ui.lbl_new_status_project_name,
            self.ui.btn_new_create_project,
            self.ui.cb_new_add_reference,
        )

    def create_new_project(self) -> None:
        """Creates a new project with the content of the new page."""
        # <editor-fold desc="Checks">
        if self.app_settings.wsl_install == 0:
            basic_boxes.ok(
                "Create new project",
                "Please install local colabfold to create a project!",
                QtWidgets.QMessageBox.Warning,
            )
            return
        if self.app_settings.local_colabfold == 0:
            basic_boxes.ok(
                "Create new project",
                "Please install local colabfold to create a project!",
                QtWidgets.QMessageBox.Warning,
            )
            return

        # </editor-fold>

        new_project_name = self.ui.txt_new_project_name.text()
        self.app_project = project.Project(new_project_name, self.workspace_path)
        if self.ui.cb_new_add_reference.checkState() == 2:
            if len(self.ui.txt_new_choose_reference.text()) == 4:
                tmp_ref_protein = pymol_io.get_protein_from_pdb(self.ui.txt_new_choose_reference.text().upper())
            else:
                # local pdb file as input
                pdb_filepath = path_util.FilePath(pathlib.Path(self.ui.txt_new_choose_reference.text()))
                graphic_operations.setup_default_session_graphic_settings()
                tmp_ref_protein = protein.Protein(
                    molecule_object=pdb_filepath.get_filename(), pdb_filepath=pdb_filepath,
                )
            # fixme: should the disulfid-bonds be displayed globally
            # cmd.select(name="disulfides", selection="byres (resn CYS and name SG) within 2 of (resn CYS and name SG)")
            # cmd.color(color="atomic", selection="disulfides and not elem C")
            # cmd.set("valence", 0)  # this needs to be better implemented
            # cmd.show("sticks", "disulfides")
            self.app_project.add_existing_protein(tmp_ref_protein)
        self.ui.cb_new_add_reference.setCheckState(0)

        self.app_project.serialize_project(
            pathlib.Path(f"{self.workspace_path}/{self.app_project.get_project_name()}.xml"),
        )
        constants.PYSSA_LOGGER.info(f"Created the project {self.app_project.get_project_name()}.")
        protein_names = []
        for tmp_protein in self.app_project.proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        constants.PYSSA_LOGGER.debug(f"These are the proteins {protein_names}.")
        self._project_watcher.current_project = self.app_project
        self.project_scanner.project = self.app_project
        constants.PYSSA_LOGGER.info(
            f"{self._project_watcher.current_project.get_project_name()} is the current project.",
        )
        self._project_watcher.on_home_page = False
        # update gui
        self._project_watcher.show_valid_options(self.ui)
        self.ui.lbl_current_project_name.setText(new_project_name)
        self.status_bar.showMessage(f"Current project path: {self.workspace_path}/{new_project_name}")
        self.display_view_page()

    # </editor-fold>

    # <editor-fold desc="Open project page functions">
    def validate_open_search(self) -> None:
        """Validates the input of the project name in real-time."""
        if self.ui.list_open_projects.currentItem() is not None:
            self.ui.list_open_projects.currentItem().setSelected(False)
        # set color for lineEdit
        input_validator.InputValidator.validate_search_input(
            self.ui.list_open_projects,
            self.ui.txt_open_search,
            self.ui.lbl_open_status_search,
            self.ui.txt_open_selected_project,
        )

    def select_project_from_open_list(self) -> None:
        """Sets the selected project name in the text box."""
        try:
            self.ui.txt_open_selected_project.setText(self.ui.list_open_projects.currentItem().text())
        except AttributeError:
            self.ui.txt_open_selected_project.setText("")

    def activate_open_button(self) -> None:
        """Activates the open button."""
        if self.ui.txt_open_selected_project.text() == "":
            self.ui.btn_open_open_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_open_open_project)
        else:
            self.ui.btn_open_open_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_open_open_project)

    def thread_func_open_project(self) -> None:
        """Logic of the open project process."""
        # show project management options in side menu
        tmp_project_path = pathlib.Path(f"{self.workspace_path}/{self.ui.list_open_projects.currentItem().text()}")
        self.app_project = project.Project.deserialize_project(tmp_project_path, self.app_settings)
        constants.PYSSA_LOGGER.info(f"Opening the project {self.app_project.get_project_name()}.")
        self._project_watcher.current_project = self.app_project
        self.project_scanner.project = self.app_project
        constants.PYSSA_LOGGER.info(
            f"{self._project_watcher.current_project.get_project_name()} is the current project.",
        )
        self.ui.lbl_current_project_name.setText(self.app_project.get_project_name())
        self._project_watcher.on_home_page = False
        self._project_watcher.show_valid_options(self.ui)

    def pre_open_project(self) -> None:
        """Opens an existing project."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)

        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.OpenProjectWorker(self.workspace_path, 
                                                         self.ui.list_open_projects.currentItem().text(),
                                                         self.app_settings)
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.open_project)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>
        
        constants.PYSSA_LOGGER.info(f"Opening the project {self.app_project.get_project_name()}.")

    def open_project(self, project_obj: project.Project):
        self.app_project = project_obj
        self._project_watcher.current_project = self.app_project
        self.project_scanner.project = self.app_project
        constants.PYSSA_LOGGER.info(
            f"{self._project_watcher.current_project.get_project_name()} is the current project.",
        )
        self.ui.lbl_current_project_name.setText(self.app_project.get_project_name())
        self._project_watcher.on_home_page = False
        self._project_watcher.show_valid_options(self.ui)

        # tmp_thread = threading.Thread(target=self.thread_func_open_project, daemon=True)
        # tmp_thread.start()
        # tmp_thread.join()
        cmd.reinitialize()
        #self.ui.list_hotspots_choose_protein.clear()
        # try:
        #     self.ui.list_hotspots_choose_protein.currentItem().setText(None)
        # except AttributeError:
        #     print("")
        self.ui.btn_manage_session.hide()
        self.display_view_page()
        QtWidgets.QApplication.restoreOverrideCursor()
    
    # </editor-fold>

    # <editor-fold desc="Delete project page functions">
    def select_project_from_delete_list(self) -> None:
        """Selects a project from the project list on the delete page."""
        try:
            self.ui.txt_delete_selected_projects.setText(self.ui.list_delete_projects.currentItem().text())
        except AttributeError:
            self.ui.txt_delete_selected_projects.setText("")

    def activate_delete_button(self) -> None:
        """Activates the delete button."""
        if self.ui.txt_delete_selected_projects.text() == "":
            self.ui.btn_delete_delete_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_delete_delete_project)
        else:
            self.ui.btn_delete_delete_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_delete_delete_project)

    def validate_delete_search(self) -> None:
        """Validates the input of the project name in real-time."""
        if self.ui.list_delete_projects.currentItem() is not None:
            self.ui.list_delete_projects.currentItem().setSelected(False)
        # set color for lineEdit
        input_validator.InputValidator.validate_search_input(
            self.ui.list_delete_projects,
            self.ui.txt_delete_search,
            self.ui.lbl_delete_status_search,
            self.ui.txt_delete_selected_projects,
        )

    def delete_project(self) -> None:
        """Deletes an existing project."""
        # popup message which warns the user that the selected project gets deleted
        response: bool = gui_utils.warning_message_project_gets_deleted()
        tmp_project_name = self.ui.txt_delete_selected_projects.text()
        if response is True:
            os.remove(pathlib.Path(f"{self.workspace_path}/{self.ui.txt_delete_selected_projects.text()}"))
            if self.ui.txt_delete_selected_projects.text() == self.ui.lbl_current_project_name.text():
                self.ui.lbl_current_project_name.clear()
            self.ui.txt_delete_selected_projects.clear()
            # update list
            self.ui.list_delete_projects.clear()
            # pre-process
            self.status_bar.showMessage(self.workspace.text())
            self.ui.list_delete_projects.clear()
            self.status_bar.showMessage(self.workspace.text())
            tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_delete_projects)
            constants.PYSSA_LOGGER.info(f"The project {tmp_project_name} was successfully deleted.")
            if self.ui.list_delete_projects.count() == 0:
                self.display_home_page()
                self._project_watcher.check_workspace_for_projects(self.workspace_path, self.ui)
        else:
            constants.PYSSA_LOGGER.info("No project has been deleted. No changes were made.")
            return

    # </editor-fold>

    # <editor-fold desc="Save project functions">
    def save_project(self) -> None:
        """Saves the project.xml."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_save_project)
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        QtWidgets.QApplication.restoreOverrideCursor()
        basic_boxes.ok("Save Project", "The project was successfully saved.", QtWidgets.QMessageBox.Information)

    # </editor-fold>

    # <editor-fold desc="Edit project page functions">
    def check_for_cleaning(self) -> None:
        """Checks if the selected protein can be cleaned."""
        try:
            protein_name = self.ui.list_edit_project_proteins.currentItem().text()
        except AttributeError:
            return
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        cmd.reinitialize()
        self.ui.btn_manage_session.show()
        tmp_protein = self.app_project.search_protein(protein_name.replace(".pdb", ""))
        tmp_protein.load_protein_in_pymol()
        # check if selected protein contains any organic or solvent molecules which can be removed
        if cmd.select("organic") > 0 or cmd.select("solvent") > 0:
            gui_elements_to_show = [
                self.ui.lbl_edit_clean_new_prot,
                self.ui.btn_edit_clean_new_prot,
                self.ui.lbl_edit_clean_update_prot,
                self.ui.btn_edit_clean_update_prot,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
        else:
            gui_elements_to_hide = [
                self.ui.lbl_edit_clean_new_prot,
                self.ui.btn_edit_clean_new_prot,
                self.ui.lbl_edit_clean_update_prot,
                self.ui.btn_edit_clean_update_prot,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        # check if selected protein is in any existing protein pair
        if self.app_project.check_if_protein_is_in_any_protein_pair(protein_name) is True:
            gui_elements_to_hide = [
                self.ui.label_12,
                self.ui.btn_edit_project_delete,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            gui_elements_to_show = [
                self.ui.label_12,
                self.ui.btn_edit_project_delete,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
        self.ui.btn_edit_project_save.show()
        self.ui.label_13.show()

    def clean_protein_new(self) -> None:
        """Cleans the selected protein structure and creates a new cleaned structure."""
        tmp_protein = self.app_project.search_protein(
            self.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""),
        )
        clean_tmp_protein = tmp_protein.clean_protein(new_protein=True)
        constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", clean_tmp_protein.get_molecule_object())
        self.app_project.add_existing_protein(clean_tmp_protein)
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        self._init_edit_page()

    def clean_protein_update(self) -> None:
        """Cleans the selected protein structure."""
        tmp_protein = self.app_project.search_protein(
            self.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""),
        )
        if basic_boxes.yes_or_no(
            "Clean protein",
            "Are you sure you want to clean this protein?\n" "This will remove all organic and solvent components!",
            QtWidgets.QMessageBox.Information,
        ):
            tmp_protein.clean_protein()
            constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", tmp_protein.get_molecule_object())
            self.app_project.serialize_project(self.app_project.get_project_xml_path())
            self._init_edit_page()
        else:
            constants.PYSSA_LOGGER.info("No protein has been cleaned.")
            return

    def delete_protein(self) -> None:
        """Deletes the selected protein structure."""
        protein_name = self.ui.list_edit_project_proteins.currentItem().text()
        response = gui_utils.warning_message_protein_gets_deleted()
        if response:
            try:
                self.app_project.delete_specific_protein(protein_name)
            except ValueError:
                constants.PYSSA_LOGGER.error(
                    "The protein %s could not be deleted, because it is not in the project.", protein_name,
                )
            self.ui.list_edit_project_proteins.clear()
            gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_edit_project_proteins)
            self._project_watcher.show_valid_options(self.ui)
        else:
            constants.PYSSA_LOGGER.info("No protein was deleted.")

    def add_existing_protein(self) -> None:
        tmp_dialog = dialog_add_model.DialogAddModel()
        tmp_dialog.exec_()
        if len(dialog_add_model.global_var_add_model[0]) == 4 and dialog_add_model.global_var_add_model[1] is True:
            tmp_ref_protein = pymol_io.get_protein_from_pdb(dialog_add_model.global_var_add_model[0].upper())
            self.app_project.add_existing_protein(tmp_ref_protein)
        elif len(dialog_add_model.global_var_add_model[0]) == 0 and dialog_add_model.global_var_add_model[1] is False:
            print("Dialog closed.")
        elif len(dialog_add_model.global_var_add_model[0]) == 4 and dialog_add_model.global_var_add_model[1] is False:
            print("Unexpected Error.")
        else:
            # local pdb file as input
            pdb_filepath = path_util.FilePath(pathlib.Path(dialog_add_model.global_var_add_model[0]))
            graphic_operations.setup_default_session_graphic_settings()
            tmp_ref_protein = protein.Protein(
                molecule_object=pdb_filepath.get_filename(), pdb_filepath=pdb_filepath,
            )
            self.app_project.add_existing_protein(tmp_ref_protein)
        self._project_watcher.show_valid_options(self.ui)
        self.display_edit_page()

    def save_selected_protein_structure_as_pdb_file(self) -> None:
        """Saves selected protein as pdb file."""
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(self, "Save protein structure", "", "Protein Data Bank File (*.pdb)")
        if file_path:
            tmp_protein = self.app_project.search_protein(self.ui.list_edit_project_proteins.currentItem().text())
            try:
                bio_data.convert_xml_string_to_pdb_file(bio_data.convert_pdb_data_list_to_xml_string(tmp_protein.get_pdb_data()), pathlib.Path(file_path))
            except:
                basic_boxes.ok(
                    "Save protein structure",
                    "Saving the protein as .pdb file failed!", QtWidgets.QMessageBox.Error,
                )
            else:
                basic_boxes.ok(
                    "Save protein structure",
                    "The protein was successfully saved as .pdb file.", QtWidgets.QMessageBox.Information,
                )
    # </editor-fold>

    # <editor-fold desc="View project page functions">
    def display_view_page(self) -> None:
        """Displays the edit project page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self.ui.list_view_project_proteins.clear()
        self.ui.txtedit_view_sequence.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        # list all proteins from pdb directory
        gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_view_project_proteins)
        # self.project_scanner.scan_project_for_valid_proteins(list_view_project_proteins=self.ui.list_view_project_proteins)

        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 11, "View proteins of current project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_view_page)
        gui_elements_to_hide = [
            self.ui.btn_view_project_show,
            self.ui.btn_view_project_show_structure,
            self.ui.txtedit_view_sequence,
            self.ui.label_9,
            self.ui.label_11,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def view_show_options(self) -> None:
        """Controls which gui elements get shown."""
        gui_elements_to_show = [
            self.ui.btn_view_project_show,
            self.ui.btn_view_project_show_structure,
            self.ui.txtedit_view_sequence,
            self.ui.label_9,
            self.ui.label_11,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)

    def view_sequence(self) -> None:
        """Displays the sequence of the selected protein in a text box."""
        tmp_protein_basename = self.ui.list_view_project_proteins.currentItem().text()
        tmp_protein_sequences = self.app_project.search_protein(tmp_protein_basename).get_protein_sequences()
        self.ui.txtedit_view_sequence.clear()
        for tmp_sequence in tmp_protein_sequences:
            self.ui.txtedit_view_sequence.append("".join(tmp_sequence.sequence))
        # fixme: experimental sequence viewer gui
        # dialog = dialog_sequence_viewer.SequenceViewer(tmp_protein_sequences, tmp_protein_filename)
        # dialog.exec_()

    def view_structure(self) -> None:
        """Displays the structure of the selected protein in pymol."""
        protein_name = self.ui.list_view_project_proteins.currentItem().text()
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        cmd.reinitialize()
        self.ui.btn_manage_session.show()
        try:
            self.app_project.search_protein(protein_name).load_protein_pymol_session()
            constants.PYSSA_LOGGER.info("Loaded PyMOL session of protein %s", protein_name)
        except pymol.CmdException:
            constants.PYSSA_LOGGER.error("Error while loading protein in PyMOL!")
            return
        self.current_session = current_session.CurrentSession(
            "protein", protein_name, self.app_project.search_protein(protein_name).pymol_session,
        )
        print(self.current_session)

    # </editor-fold>

    # <editor-fold desc="Use project page functions">
    def validate_use_project_name(self) -> None:
        """Validates the input of the project name in real-time."""
        input_validator.InputValidator.validate_project_name(
            self.ui.list_use_existing_projects,
            self.ui.txt_use_project_name,
            self.ui.lbl_use_status_project_name,
            self.ui.btn_use_next,
        )

    def validate_use_search(self) -> None:
        """Validates the input of the protein name in real-time."""
        message = "Protein structure does not exists."
        input_validator.InputValidator.validate_search_input(
            self.ui.list_use_available_protein_structures,
            self.ui.txt_use_search,
            self.ui.lbl_use_status_search,
            status_message=message,
        )

    def add_protein_structure_to_new_project(self) -> None:
        """Adds the selected protein to the list which is used to create the new project."""
        prot_to_add = self.ui.list_use_available_protein_structures.currentItem().text()
        self.ui.list_use_selected_protein_structures.addItem(prot_to_add)
        self.ui.list_use_available_protein_structures.takeItem(
            self.ui.list_use_available_protein_structures.currentRow(),
        )
        self.ui.btn_use_add_available_protein_structures.setEnabled(False)
        if self.ui.list_use_available_protein_structures.count() > 0:
            try:
                self.ui.list_use_available_protein_structures.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in use available proteins list on Use page.")

    def remove_protein_structure_to_new_project(self) -> None:
        """Removes the selected protein from the list which is used to create the new project."""
        prot_to_remove = self.ui.list_use_selected_protein_structures.currentItem()
        self.ui.list_use_selected_protein_structures.takeItem(self.ui.list_use_selected_protein_structures.currentRow())
        self.ui.list_use_available_protein_structures.addItem(prot_to_remove)
        self.ui.btn_use_remove_selected_protein_structures.setEnabled(False)
        if self.ui.list_use_selected_protein_structures.count() > 0:
            try:
                self.ui.list_use_selected_protein_structures.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in use selected proteins list on Use page.")

    def show_protein_selection_for_use(self) -> None:
        """Shows the two lists for the protein selection."""
        gui_elements_to_show = [
            self.ui.lbl_use_search,
            self.ui.lbl_use_status_search,
            self.ui.txt_use_search,
            self.ui.btn_use_add_available_protein_structures,
            self.ui.lbl_use_available_protein_structures,
            self.ui.list_use_available_protein_structures,
            self.ui.btn_use_remove_selected_protein_structures,
            self.ui.lbl_use_selected_protein_structures,
            self.ui.list_use_selected_protein_structures,
            self.ui.btn_use_back,
            self.ui.btn_use_create_new_project,
            self.ui.lbl_use_new_project,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        self.ui.txt_use_project_name.setEnabled(False)
        gui_elements_to_hide = [
            self.ui.btn_use_next,
            self.ui.list_use_existing_projects,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.disable_text_box(self.ui.txt_use_project_name, self.ui.lbl_use_project_name)
        self.ui.btn_use_add_available_protein_structures.setEnabled(False)
        self.ui.btn_use_remove_selected_protein_structures.setEnabled(False)

    def hide_protein_selection_for_use(self) -> None:
        """Hides the two lists for the protein selection."""
        gui_elements_to_show = [
            self.ui.btn_use_next,
            self.ui.list_use_existing_projects,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        self.ui.txt_use_project_name.setEnabled(True)

        gui_elements_to_hide = [
            self.ui.lbl_use_search,
            self.ui.lbl_use_status_search,
            self.ui.txt_use_search,
            self.ui.btn_use_add_available_protein_structures,
            self.ui.lbl_use_available_protein_structures,
            self.ui.list_use_available_protein_structures,
            self.ui.btn_use_remove_selected_protein_structures,
            self.ui.lbl_use_selected_protein_structures,
            self.ui.list_use_selected_protein_structures,
            self.ui.btn_use_back,
            self.ui.btn_use_create_new_project,
            self.ui.lbl_use_new_project,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_use_project_name, self.ui.lbl_use_project_name)

    def use_enable_add(self) -> None:
        """Enables the add button."""
        self.ui.btn_use_add_available_protein_structures.setEnabled(True)

    def use_enable_remove(self) -> None:
        """Enables the remove button."""
        self.ui.btn_use_remove_selected_protein_structures.setEnabled(True)

    def pre_create_use_project(self) -> None:
        """Sets up the worker for the create_use_project task."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        # copy proteins in new project
        proteins_to_copy = []
        for i in range(self.ui.list_use_selected_protein_structures.count()):
            self.ui.list_use_selected_protein_structures.setCurrentRow(i)
            proteins_to_copy.append(self.ui.list_use_selected_protein_structures.currentItem().text())

        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.CreateUseProjectWorker(self.workspace_path, proteins_to_copy)
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.create_use_project)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        self.ui.lbl_current_project_name.setText(self.ui.txt_use_project_name.text())
        self.status_bar.showMessage(f"Creating new project: {self.ui.txt_use_project_name.text()} ...")
        # save project folder in current workspace
        new_project = project.Project(self.ui.txt_use_project_name.text(), self.workspace_path)
        # new_project.create_project_tree()
        self.app_project = new_project

    def create_use_project(self, proteins_for_new_project) -> None:
        """Post thread method.

        Args:
            proteins_for_new_project (list): the proteins which are in the new project
        """
        # protein_infos = (tools.scan_workspace_for_non_duplicate_proteins(self.workspace_path))[0]
        # for tmp_protein in proteins_to_copy:
        #     for tmp_protein_info in protein_infos:
        #         if tmp_protein_info.name == tmp_protein:
        #             """Var: project_proteins is a list which contains all proteins from a single project"""
        #             xml_deserializer = filesystem_io.XmlDeserializer(
        #                 pathlib.Path(f"{self.workspace_path}/{tmp_protein_info.project_name}.xml"))
        #             for xml_protein in xml_deserializer.xml_root.iter(element_names.PROTEIN):
        #                 if xml_protein.attrib[attribute_names.ID] == tmp_protein_info.id:
        #                     basic_information = xml_protein.attrib
        #                     pdb_lines = []
        #                     session_data_base64 = ""
        #                     for tmp_data in xml_protein:
        #                         if tmp_data.tag == "pdb_data":
        #                             for tmp_atom in tmp_data.findall("atom"):
        #                                 pdb_lines.append(tmp_atom.text)
        #                         elif tmp_data.tag == "session_data":
        #                             session_data_base64 = tmp_data.attrib[attribute_names.PROTEIN_SESSION]
        #                         else:
        #                             raise ValueError
        #                     tmp_protein_obj = protein.Protein(
        #                         molecule_object=basic_information[attribute_names.PROTEIN_MOLECULE_OBJECT],
        #                         pdb_xml_string=xml_protein)
        #                     tmp_protein_obj.set_all_attributes(basic_information, pdb_lines, session_data_base64)
        for tmp_protein_obj in proteins_for_new_project:
            self.app_project.add_existing_protein(tmp_protein_obj)
        self.app_project.serialize_project(
            pathlib.Path(f"{self.workspace_path}/{self.app_project.get_project_name()}.xml"),
        )
        # shows options which can be done with the data in the project folder
        self._project_watcher.current_project = self.app_project
        self._project_watcher.on_home_page = False
        self._project_watcher.show_valid_options(self.ui)
        self.project_scanner.project = self.app_project
        self._init_use_page()
        constants.PYSSA_LOGGER.info(
            f"The project {self.app_project.get_project_name()} was successfully created through a use.",
        )
        self.display_view_page()
        QtWidgets.QApplication.restoreOverrideCursor()

    # </editor-fold>

    # <editor-fold desc="Import, Export functions">
    def import_project(self) -> None:
        """Imports a project.xml into the current workspace."""
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_import_project)
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getOpenFileName(self, "Select a project file to import", "", "XML Files (*.xml)")
        if file_path:
            tmp_project = project.Project("", self.workspace_path)
            tmp_project = tmp_project.deserialize_project(pathlib.Path(file_path), self.app_settings)
            tmp_project.set_workspace_path(self.workspace_path)
            if len(tmp_project.proteins) <= 1:
                if self.app_settings.wsl_install == 0:
                    basic_boxes.ok(
                        "Create new project",
                        "Please install local colabfold to import this project!",
                        QtWidgets.QMessageBox.Warning,
                    )
                    return
                elif self.app_settings.local_colabfold == 0:
                    basic_boxes.ok(
                        "Create new project",
                        "Please install local colabfold to import this project!",
                        QtWidgets.QMessageBox.Warning,
                    )
                    return
            new_filepath = pathlib.Path(f"{self.workspace_path}/{tmp_project.get_project_name()}.xml")
            tmp_project.serialize_project(new_filepath)
            self.app_project = self.app_project.deserialize_project(new_filepath, self.app_settings)
            constants.PYSSA_LOGGER.info(f"Opening the project {self.app_project.get_project_name()}.")
            self._project_watcher.current_project = self.app_project
            self.project_scanner.project = self.app_project
            constants.PYSSA_LOGGER.info(
                f"{self._project_watcher.current_project.get_project_name()} is the current project.",
            )
            self.ui.lbl_current_project_name.setText(self.app_project.get_project_name())
            self._project_watcher.on_home_page = False
            self._project_watcher.show_valid_options(self.ui)
            self.ui.btn_manage_session.show()
            self.display_view_page()
            basic_boxes.ok(
                "Import Project", "The project was successfully imported.", QtWidgets.QMessageBox.Information,
            )

    def export_current_project(self) -> None:
        """Exports the current project to an importable format."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_export_project)
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(self, "Save current project", "", "XML Files (*.xml)")
        if file_path:
            self.app_project.serialize_project(pathlib.Path(file_path))
            basic_boxes.ok(
                "Export Project", "The project was successfully exported.", QtWidgets.QMessageBox.Information,
            )

    # </editor-fold>

    # <editor-fold desc="Close project functions">
    def close_project(self) -> None:
        """Closes the current project."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        cmd.reinitialize()
        self.ui.list_hotspots_choose_protein.clear()
        self._project_watcher.on_home_page = True
        self._project_watcher.current_project = project.Project("", pathlib.Path(""))
        self._project_watcher.show_valid_options(self.ui)
        self.ui.lbl_current_project_name.setText("")
        self._init_all_pages()
        self.results_name = ""
        constants.PYSSA_LOGGER.info(f"The project {self.app_project.get_project_name()} was closed")
        self.display_home_page()

    # </editor-fold>

    # <editor-fold desc="ESMFold Monomer functions">
    def _init_esm_pred_mono_page(self) -> None:
        """Clears all text boxes and sets up the default values for the page."""
        # clears everything
        self.ui.txt_esm_prot_name.clear()
        self.ui.txt_esm_prot_seq.clear()
        for i in range(self.ui.table_esm_prot_to_predict.rowCount()):
            self.ui.table_esm_prot_to_predict.removeRow(i)
        # sets up defaults: Prediction
        self.ui.btn_esm_next.setEnabled(False)
        self.ui.btn_esm_next_2.setEnabled(False)
        self.ui.lbl_esm_prot_name_status.setText("")
        self.ui.lbl_esm_prot_seq_status.setText("")

    def display_esm_pred_mono(self) -> None:
        """Displays the esm_fold monomer page."""
        self._init_esm_pred_mono_page()
        gui_elements_to_show = [
            self.ui.lbl_esm_prot_to_predict,
            self.ui.table_esm_prot_to_predict,
            self.ui.btn_esm_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self.ui.btn_esm_seq_to_predict_remove,
            self.ui.lbl_esm_prot_name,
            self.ui.txt_esm_prot_name,
            self.ui.lbl_esm_prot_name_status,
            self.ui.btn_esm_back,
            self.ui.btn_esm_next,
            self.ui.lbl_esm_prot_seq,
            self.ui.txt_esm_prot_seq,
            self.ui.lbl_esm_prot_seq_status,
            self.ui.btn_esm_back_2,
            self.ui.btn_esm_next_2,
            self.ui.btn_esm_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        styles.color_button_not_ready(self.ui.btn_esm_next)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 1, "ESMFold Monomer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_pred_cloud_monomer_page,
        )

    def cloud_esm_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self.ui.txt_esm_prot_name.text(), self.ui.table_esm_prot_to_predict,
        ):
            self.ui.lbl_esm_prot_name_status.setText("Protein name already used.")
            self.ui.btn_esm_next.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_esm_next)
        else:
            self.ui.btn_esm_next.setEnabled(True)
            tools.validate_protein_name(
                self.ui.txt_esm_prot_name, self.ui.lbl_esm_prot_name_status, self.ui.btn_esm_next,
            )

    def cloud_esm_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self.ui.txt_esm_prot_seq, self.ui.lbl_esm_prot_seq_status, self.ui.btn_esm_next_2,
        )

    def setup_defaults_esm_monomer_prediction(self) -> None:
        """Sets up the default values for the page."""
        # clears everything
        self.ui.txt_esm_prot_name.clear()
        self.ui.txt_esm_prot_seq.clear()
        # sets up defaults: Prediction
        self.ui.btn_esm_next.setEnabled(False)
        self.ui.btn_esm_next_2.setEnabled(False)
        self.ui.lbl_esm_prot_name_status.setText("")
        self.ui.lbl_esm_prot_seq_status.setText("")

    def cloud_esm_add_seq_to_predict(self) -> None:
        """Shows the gui elements to add a sequence to the protein to predict."""
        gui_elements_to_show = [
            self.ui.lbl_esm_prot_to_predict,
            self.ui.table_esm_prot_to_predict,
            self.ui.lbl_esm_prot_name,
            self.ui.txt_esm_prot_name,
            self.ui.lbl_esm_prot_name_status,
            self.ui.btn_esm_back,
            self.ui.btn_esm_next,
        ]
        gui_utils.enable_text_box(self.ui.txt_esm_prot_name, self.ui.lbl_esm_prot_name)
        gui_elements_to_hide = [
            self.ui.btn_esm_seq_to_predict_remove,
            self.ui.btn_esm_seq_to_predict,
            self.ui.lbl_esm_prot_seq,
            self.ui.txt_esm_prot_seq,
            self.ui.lbl_esm_prot_seq_status,
            self.ui.btn_esm_back_2,
            self.ui.btn_esm_next_2,
            self.ui.btn_esm_predict,
        ]
        gui_utils.disable_text_box(self.ui.txt_esm_prot_seq, self.ui.lbl_esm_prot_seq)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_esm_next.setEnabled(False)
        self.ui.txt_esm_prot_name.clear()
        styles.color_button_not_ready(self.ui.btn_esm_next)
        if self.ui.table_esm_prot_to_predict.rowCount() > 0:
            try:
                self.ui.table_esm_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Monomer Prediction in overview table.")

    def cloud_esm_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_esm_prot_to_predict,
            self.ui.table_esm_prot_to_predict,
            self.ui.btn_esm_seq_to_predict_remove,
            self.ui.btn_esm_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self.ui.lbl_esm_prot_name,
            self.ui.txt_esm_prot_name,
            self.ui.lbl_esm_prot_name_status,
            self.ui.btn_esm_back,
            self.ui.btn_esm_next,
            self.ui.lbl_esm_prot_seq,
            self.ui.txt_esm_prot_seq,
            self.ui.lbl_esm_prot_seq_status,
            self.ui.btn_esm_back_2,
            self.ui.btn_esm_next_2,
            self.ui.btn_esm_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.cloud_esm_check_if_table_is_empty()
        self.ui.btn_esm_seq_to_predict_remove.setEnabled(False)

    def cloud_esm_next(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_esm_prot_to_predict,
            self.ui.table_esm_prot_to_predict,
            self.ui.lbl_esm_prot_name,
            self.ui.txt_esm_prot_name,
            self.ui.lbl_esm_prot_seq,
            self.ui.txt_esm_prot_seq,
            self.ui.lbl_esm_prot_seq_status,
            self.ui.btn_esm_back_2,
            self.ui.btn_esm_next_2,
        ]
        gui_utils.enable_text_box(self.ui.txt_esm_prot_seq, self.ui.lbl_esm_prot_seq)
        gui_elements_to_hide = [
            self.ui.btn_esm_seq_to_predict_remove,
            self.ui.btn_esm_seq_to_predict,
            self.ui.lbl_esm_prot_name_status,
            self.ui.btn_esm_back,
            self.ui.btn_esm_next,
            self.ui.btn_esm_predict,
        ]
        gui_utils.disable_text_box(self.ui.txt_esm_prot_name, self.ui.lbl_esm_prot_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.txt_esm_prot_seq.clear()

    def cloud_esm_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_esm_prot_to_predict,
            self.ui.table_esm_prot_to_predict,
            self.ui.lbl_esm_prot_name,
            self.ui.txt_esm_prot_name,
            self.ui.lbl_esm_prot_name_status,
            self.ui.btn_esm_back,
            self.ui.btn_esm_next,
        ]
        gui_elements_to_hide = [
            self.ui.btn_esm_seq_to_predict_remove,
            self.ui.btn_esm_seq_to_predict,
            self.ui.lbl_esm_prot_seq,
            self.ui.txt_esm_prot_seq,
            self.ui.lbl_esm_prot_seq_status,
            self.ui.btn_esm_back_2,
            self.ui.btn_esm_next_2,
            self.ui.btn_esm_predict,
        ]
        gui_utils.enable_text_box(self.ui.txt_esm_prot_name, self.ui.lbl_esm_prot_name)
        gui_utils.disable_text_box(self.ui.txt_esm_prot_seq, self.ui.lbl_esm_prot_seq)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def cloud_esm_add_protein(self) -> None:
        """Adds protein to the list of proteins to predict."""
        self.ui.table_esm_prot_to_predict.setRowCount(self.ui.table_esm_prot_to_predict.rowCount() + 1)
        self.ui.table_esm_prot_to_predict.insertRow(self.ui.table_esm_prot_to_predict.rowCount() + 1)
        self.ui.table_esm_prot_to_predict.setItem(
            self.ui.table_esm_prot_to_predict.rowCount() - 1, 0, QtWidgets.QTableWidgetItem("A"),
        )
        self.ui.table_esm_prot_to_predict.setItem(
            self.ui.table_esm_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self.ui.txt_esm_prot_seq.toPlainText()),
        )
        self.ui.table_esm_prot_to_predict.setVerticalHeaderItem(
            self.ui.table_esm_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self.ui.txt_esm_prot_name.text()),
        )
        self.ui.table_esm_prot_to_predict.resizeColumnsToContents()
        self.cloud_esm_check_if_table_is_empty()
        gui_elements_to_show = [
            self.ui.lbl_esm_prot_to_predict,
            self.ui.table_esm_prot_to_predict,
            self.ui.btn_esm_seq_to_predict_remove,
            self.ui.btn_esm_seq_to_predict,
            self.ui.btn_esm_predict,
        ]
        gui_utils.enable_text_box(self.ui.txt_esm_prot_name, self.ui.lbl_esm_prot_name)
        gui_elements_to_hide = [
            self.ui.lbl_esm_prot_name,
            self.ui.txt_esm_prot_name,
            self.ui.lbl_esm_prot_name_status,
            self.ui.btn_esm_back,
            self.ui.btn_esm_next,
            self.ui.lbl_esm_prot_seq,
            self.ui.txt_esm_prot_seq,
            self.ui.lbl_esm_prot_seq_status,
            self.ui.btn_esm_back_2,
            self.ui.btn_esm_next_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_esm_predict.setEnabled(True)
        self.ui.btn_esm_seq_to_predict_remove.setEnabled(False)
        styles.color_button_ready(self.ui.btn_esm_predict)
        self.setup_defaults_esm_monomer_prediction()

    def cloud_esm_remove(self) -> None:
        """Removes a protein from the list of proteins to predict."""
        self.ui.table_esm_prot_to_predict.removeRow(self.ui.table_esm_prot_to_predict.currentRow())
        gui_elements_to_show = [
            self.ui.lbl_esm_prot_to_predict,
            self.ui.table_esm_prot_to_predict,
            self.ui.btn_esm_seq_to_predict_remove,
            self.ui.btn_esm_seq_to_predict,
            self.ui.btn_esm_predict,
        ]
        gui_utils.enable_text_box(self.ui.txt_esm_prot_name, self.ui.lbl_esm_prot_name)
        gui_elements_to_hide = [
            self.ui.lbl_esm_prot_name,
            self.ui.txt_esm_prot_name,
            self.ui.lbl_esm_prot_name_status,
            self.ui.btn_esm_back,
            self.ui.btn_esm_next,
            self.ui.lbl_esm_prot_seq,
            self.ui.txt_esm_prot_seq,
            self.ui.lbl_esm_prot_seq_status,
            self.ui.btn_esm_back_2,
            self.ui.btn_esm_next_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_esm_seq_to_predict_remove.setEnabled(False)
        self.cloud_esm_check_if_table_is_empty()

    def cloud_esm_item_changed(self) -> None:
        """Enables the remove button."""
        self.ui.btn_esm_seq_to_predict_remove.setEnabled(True)

    def cloud_esm_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self.ui.table_esm_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_esm_predict)
            self.ui.btn_esm_predict.setEnabled(False)
            gui_elements_to_show = [
                self.ui.lbl_esm_prot_to_predict,
                self.ui.table_esm_prot_to_predict,
                self.ui.btn_esm_seq_to_predict,
            ]
            gui_utils.enable_text_box(self.ui.txt_esm_prot_name, self.ui.lbl_esm_prot_name)
            gui_elements_to_hide = [
                self.ui.btn_esm_seq_to_predict_remove,
                self.ui.lbl_esm_prot_name,
                self.ui.txt_esm_prot_name,
                self.ui.lbl_esm_prot_name_status,
                self.ui.btn_esm_back,
                self.ui.btn_esm_next,
                self.ui.lbl_esm_prot_seq,
                self.ui.txt_esm_prot_seq,
                self.ui.lbl_esm_prot_seq_status,
                self.ui.btn_esm_back_2,
                self.ui.btn_esm_next_2,
                self.ui.btn_esm_predict,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            styles.color_button_ready(self.ui.btn_esm_predict)
            self.ui.btn_esm_predict.setEnabled(True)

    def post_predict_esm_monomer(self, output) -> None:
        """Post thread method, for the prediction process.

        Args:
            output (list): the proteins which got predicted
        """
        self.block_box_prediction.destroy(True)
        for tmp_filename in os.listdir(constants.ESMFOLD_PDB_DIR):
            constants.PYSSA_LOGGER.info(
                f"Add protein {tmp_filename} to the current project {self.app_project.get_project_name()}",
            )
            self.app_project.add_existing_protein(
                protein.Protein(
                    tmp_filename.replace(".pdb", ""),
                    path_util.FilePath(pathlib.Path(f"{constants.ESMFOLD_PDB_DIR}/{tmp_filename}")),
                ),
            )
        if len(output) > 0:
            formatted_output = ", ".join(output)
            basic_boxes.ok(
                "ESMFold Prediction",
                f"These protein prediction failed: {formatted_output}.",
                QtWidgets.QMessageBox.Critical,
            )
        else:
            basic_boxes.ok("ESMFold Prediction", "The prediction was successful.", QtWidgets.QMessageBox.Information)
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        self.display_view_page()
        self._project_watcher.show_valid_options(self.ui)

    def predict_esm_monomer(self) -> None:
        """Sets up the worker to predict the proteins with the ESM-Fold."""
        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.EsmFoldWorker(self.ui.table_esm_prot_to_predict)
        self.tmp_thread = task_workers.setup_worker_for_work(
            self.tmp_thread, self.tmp_worker, self.post_predict_esm_monomer,
        )
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        self.block_box_prediction = QtWidgets.QMessageBox()
        self.block_box_prediction = gui_utils.setup_standard_block_box(
            self.block_box_prediction,
            "Structure Prediction",
            "A prediction is currently running.",
        )
        # self.block_box_prediction.setStandardButtons(QtWidgets.QMessageBox.NoButton)
        # self.block_box_prediction.setIcon(QtWidgets.QMessageBox.Information)
        # self.block_box_prediction.setWindowIcon(
        #     PyQt5.constants.PLUGIN_LOGO_ICON_OBJ)
        # styles.set_stylesheet(self.block_box_prediction)
        # self.block_box_prediction.setWindowTitle("Structure Prediction")
        # self.block_box_prediction.setText("A prediction is currently running.")
        self.block_box_prediction.exec_()

    # </editor-fold>

    # <editor-fold desc="Monomer Local Prediction functions">
    def _init_local_pred_mono_page(self) -> None:
        """Clears all text boxes and sets default values for the gui elements."""
        # clears everything
        self.ui.txt_pred_mono_prot_name.clear()
        self.ui.txt_pred_mono_seq_name.clear()
        for i in range(self.ui.table_pred_mono_prot_to_predict.rowCount()):
            self.ui.table_pred_mono_prot_to_predict.removeRow(i)
        # sets up defaults: Prediction
        self.ui.btn_pred_mono_next.setEnabled(False)
        self.ui.btn_pred_mono_add_protein.setEnabled(False)
        self.ui.lbl_pred_mono_prot_name_status.setText("")
        self.ui.lbl_pred_mono_seq_name_status.setText("")

    def display_local_pred_mono(self) -> None:
        """Displays the local prediction monomer page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return
        self._init_local_pred_mono_page()
        gui_elements_to_show = [
            self.ui.lbl_pred_mono_prot_to_predict,
            self.ui.table_pred_mono_prot_to_predict,
            self.ui.btn_pred_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.lbl_pred_mono_prot_name,
            self.ui.txt_pred_mono_prot_name,
            self.ui.lbl_pred_mono_prot_name_status,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
            self.ui.lbl_pred_mono_seq_name,
            self.ui.txt_pred_mono_seq_name,
            self.ui.lbl_pred_mono_seq_name_status,
            self.ui.btn_pred_mono_back_2,
            self.ui.btn_pred_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config,
            self.ui.btn_pred_mono_advanced_config,
            self.ui.btn_pred_mono_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        styles.color_button_not_ready(self.ui.btn_pred_mono_next)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 19, "Local Monomer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_pred_local_monomer_page,
        )

    def local_pred_mono_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self.ui.txt_pred_mono_prot_name.text(), self.ui.table_pred_mono_prot_to_predict,
        ):
            self.ui.lbl_pred_mono_prot_name_status.setText("Protein name already used.")
            self.ui.btn_pred_mono_next.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_pred_mono_next)
        else:
            self.ui.btn_pred_mono_next.setEnabled(True)
            tools.validate_protein_name(
                self.ui.txt_pred_mono_prot_name, self.ui.lbl_pred_mono_prot_name_status, self.ui.btn_pred_mono_next,
            )

    def local_pred_mono_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self.ui.txt_pred_mono_seq_name, self.ui.lbl_pred_mono_seq_name_status, self.ui.btn_pred_mono_add_protein,
        )

    def show_prediction_configuration(self) -> None:
        """Opens the prediction configuration dialog window."""
        config = dialog_advanced_prediction_configurations.DialogAdvancedPredictionConfigurations(
            self.prediction_configuration,
        )
        config.exec_()
        self.prediction_configuration.amber_force_field = config.prediction_config.amber_force_field
        self.prediction_configuration.templates = config.prediction_config.templates

    def setup_defaults_monomer_prediction(self) -> None:
        """Sets up the defaults for the page."""
        # clears everything
        self.ui.txt_pred_mono_prot_name.clear()
        self.ui.txt_pred_mono_seq_name.clear()
        # sets up defaults: Prediction
        self.ui.btn_pred_mono_next.setEnabled(False)
        self.ui.btn_pred_mono_add_protein.setEnabled(False)
        self.ui.lbl_pred_mono_prot_name_status.setText("")
        self.ui.lbl_pred_mono_seq_name_status.setText("")

    def local_pred_mono_add_seq_to_predict(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_mono_prot_to_predict,
            self.ui.table_pred_mono_prot_to_predict,
            self.ui.lbl_pred_mono_prot_name,
            self.ui.txt_pred_mono_prot_name,
            self.ui.lbl_pred_mono_prot_name_status,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_mono_prot_name, self.ui.lbl_pred_mono_prot_name)
        gui_elements_to_hide = [
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.btn_pred_mono_seq_to_predict,
            self.ui.lbl_pred_mono_seq_name,
            self.ui.txt_pred_mono_seq_name,
            self.ui.lbl_pred_mono_seq_name_status,
            self.ui.btn_pred_mono_back_2,
            self.ui.btn_pred_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config,
            self.ui.btn_pred_mono_advanced_config,
            self.ui.btn_pred_mono_predict,
        ]
        gui_utils.disable_text_box(self.ui.txt_pred_mono_seq_name, self.ui.lbl_pred_mono_seq_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_pred_mono_next.setEnabled(False)
        self.ui.txt_pred_mono_prot_name.clear()
        styles.color_button_not_ready(self.ui.btn_pred_mono_next)
        if self.ui.table_pred_mono_prot_to_predict.rowCount() > 0:
            try:
                self.ui.table_pred_mono_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Monomer Prediction in overview table.")

    def local_pred_mono_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_mono_prot_to_predict,
            self.ui.table_pred_mono_prot_to_predict,
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.btn_pred_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self.ui.lbl_pred_mono_prot_name,
            self.ui.txt_pred_mono_prot_name,
            self.ui.lbl_pred_mono_prot_name_status,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
            self.ui.lbl_pred_mono_seq_name,
            self.ui.txt_pred_mono_seq_name,
            self.ui.lbl_pred_mono_seq_name_status,
            self.ui.btn_pred_mono_back_2,
            self.ui.btn_pred_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config,
            self.ui.btn_pred_mono_advanced_config,
            self.ui.btn_pred_mono_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.local_pred_mono_check_if_table_is_empty()
        self.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(False)

    def local_pred_mono_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_mono_prot_to_predict,
            self.ui.table_pred_mono_prot_to_predict,
            self.ui.lbl_pred_mono_prot_name,
            self.ui.txt_pred_mono_prot_name,
            self.ui.lbl_pred_mono_seq_name,
            self.ui.txt_pred_mono_seq_name,
            self.ui.lbl_pred_mono_seq_name_status,
            self.ui.btn_pred_mono_back_2,
            self.ui.btn_pred_mono_add_protein,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_mono_seq_name, self.ui.lbl_pred_mono_seq_name)
        gui_elements_to_hide = [
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.btn_pred_mono_seq_to_predict,
            self.ui.lbl_pred_mono_prot_name_status,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
            self.ui.lbl_pred_mono_advanced_config,
            self.ui.btn_pred_mono_advanced_config,
            self.ui.btn_pred_mono_predict,
        ]
        gui_utils.disable_text_box(self.ui.txt_pred_mono_prot_name, self.ui.lbl_pred_mono_prot_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.txt_pred_mono_seq_name.clear()

    def local_pred_mono_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_mono_prot_to_predict,
            self.ui.table_pred_mono_prot_to_predict,
            self.ui.lbl_pred_mono_prot_name,
            self.ui.txt_pred_mono_prot_name,
            self.ui.lbl_pred_mono_prot_name_status,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.btn_pred_mono_seq_to_predict,
            self.ui.lbl_pred_mono_seq_name,
            self.ui.txt_pred_mono_seq_name,
            self.ui.lbl_pred_mono_seq_name_status,
            self.ui.btn_pred_mono_back_2,
            self.ui.btn_pred_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config,
            self.ui.btn_pred_mono_advanced_config,
            self.ui.btn_pred_mono_predict,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_mono_prot_name, self.ui.lbl_pred_mono_prot_name)
        gui_utils.disable_text_box(self.ui.txt_pred_mono_seq_name, self.ui.lbl_pred_mono_seq_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def local_pred_mono_add_protein(self) -> None:
        """Adds protein to the list of proteins to predict."""
        self.ui.table_pred_mono_prot_to_predict.setRowCount(self.ui.table_pred_mono_prot_to_predict.rowCount() + 1)
        self.ui.table_pred_mono_prot_to_predict.insertRow(self.ui.table_pred_mono_prot_to_predict.rowCount() + 1)
        self.ui.table_pred_mono_prot_to_predict.setItem(
            self.ui.table_pred_mono_prot_to_predict.rowCount() - 1, 0, QtWidgets.QTableWidgetItem("A"),
        )
        self.ui.table_pred_mono_prot_to_predict.setItem(
            self.ui.table_pred_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self.ui.txt_pred_mono_seq_name.toPlainText()),
        )
        self.ui.table_pred_mono_prot_to_predict.setVerticalHeaderItem(
            self.ui.table_pred_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self.ui.txt_pred_mono_prot_name.text()),
        )
        self.ui.table_pred_mono_prot_to_predict.resizeColumnsToContents()
        self.local_pred_mono_check_if_table_is_empty()
        gui_elements_to_show = [
            self.ui.lbl_pred_mono_prot_to_predict,
            self.ui.table_pred_mono_prot_to_predict,
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.btn_pred_mono_seq_to_predict,
            self.ui.lbl_pred_mono_advanced_config,
            self.ui.btn_pred_mono_advanced_config,
            self.ui.btn_pred_mono_predict,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_mono_prot_name, self.ui.lbl_pred_mono_prot_name)
        gui_elements_to_hide = [
            self.ui.lbl_pred_mono_prot_name,
            self.ui.txt_pred_mono_prot_name,
            self.ui.lbl_pred_mono_prot_name_status,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
            self.ui.lbl_pred_mono_seq_name,
            self.ui.txt_pred_mono_seq_name,
            self.ui.lbl_pred_mono_seq_name_status,
            self.ui.btn_pred_mono_back_2,
            self.ui.btn_pred_mono_add_protein,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_pred_mono_predict.setEnabled(True)
        self.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(False)
        styles.color_button_ready(self.ui.btn_pred_mono_predict)
        self.setup_defaults_monomer_prediction()

    def local_pred_mono_remove(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        self.ui.table_pred_mono_prot_to_predict.removeRow(self.ui.table_pred_mono_prot_to_predict.currentRow())
        gui_elements_to_show = [
            self.ui.lbl_pred_mono_prot_to_predict,
            self.ui.table_pred_mono_prot_to_predict,
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.btn_pred_mono_seq_to_predict,
            self.ui.lbl_pred_mono_advanced_config,
            self.ui.btn_pred_mono_advanced_config,
            self.ui.btn_pred_mono_predict,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_mono_prot_name, self.ui.lbl_pred_mono_prot_name)
        gui_elements_to_hide = [
            self.ui.lbl_pred_mono_prot_name,
            self.ui.txt_pred_mono_prot_name,
            self.ui.lbl_pred_mono_prot_name_status,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
            self.ui.lbl_pred_mono_seq_name,
            self.ui.txt_pred_mono_seq_name,
            self.ui.lbl_pred_mono_seq_name_status,
            self.ui.btn_pred_mono_back_2,
            self.ui.btn_pred_mono_add_protein,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(False)
        self.local_pred_mono_check_if_table_is_empty()

    def local_pred_mono_item_changed(self) -> None:
        """Enables the remove button."""
        self.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(True)

    def local_pred_mono_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self.ui.table_pred_mono_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_mono_predict)
            self.ui.btn_pred_mono_predict.setEnabled(False)
            gui_elements_to_show = [
                self.ui.lbl_pred_mono_prot_to_predict,
                self.ui.table_pred_mono_prot_to_predict,
                self.ui.btn_pred_mono_seq_to_predict,
            ]
            gui_utils.enable_text_box(self.ui.txt_pred_mono_prot_name, self.ui.lbl_pred_mono_prot_name)
            gui_elements_to_hide = [
                self.ui.btn_pred_mono_seq_to_predict_remove,
                self.ui.lbl_pred_mono_prot_name,
                self.ui.txt_pred_mono_prot_name,
                self.ui.lbl_pred_mono_prot_name_status,
                self.ui.btn_pred_mono_back,
                self.ui.btn_pred_mono_next,
                self.ui.lbl_pred_mono_seq_name,
                self.ui.txt_pred_mono_seq_name,
                self.ui.lbl_pred_mono_seq_name_status,
                self.ui.btn_pred_mono_back_2,
                self.ui.btn_pred_mono_add_protein,
                self.ui.lbl_pred_mono_advanced_config,
                self.ui.btn_pred_mono_advanced_config,
                self.ui.btn_pred_mono_predict,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            styles.color_button_ready(self.ui.btn_pred_mono_predict)
            self.ui.btn_pred_mono_predict.setEnabled(True)

    def post_prediction_process(self, an_exit_code: int, an_exit_code_description: str) -> None:
        """Process which runs after each prediction job."""
        if an_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction", "Prediction failed because there was an error writing the fasta file(s)!", QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction", "Prediction failed because the fasta file(s) could not be found!", QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction", "Prediction failed because a subprocess failed!", QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction", "Prediction failed because of an unknown error!", QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            if self.prediction_type == constants.PREDICTION_TYPE_PRED:
                self.app_project.serialize_project(self.app_project.get_project_xml_path())
                constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
                self.block_box_prediction.destroy(True)
                basic_boxes.ok(
                    "Structure prediction",
                    "All structure predictions are done. Go to View to check the new proteins.",
                    QtWidgets.QMessageBox.Information,
                )
                constants.PYSSA_LOGGER.info("All structure predictions are done.")
                self._project_watcher.show_valid_options(self.ui)
                self._init_local_pred_mono_page()
                self._init_local_pred_multi_page()
                self.display_view_page()
            elif self.prediction_type == constants.PREDICTION_TYPE_PRED_MONO_ANALYSIS:
                # executes if monomers were successfully predicted
                self.app_project.serialize_project(self.app_project.get_project_xml_path())
                constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
                constants.PYSSA_LOGGER.info("All structure predictions are done.")
                constants.PYSSA_LOGGER.info("Begin analysis process.")
                constants.PYSSA_LOGGER.debug(
                    f"Thread count before analysis worker: {self.threadpool.activeThreadCount()}",
                )

                # self.worker_analysis = workers.AnalysisWorkerPool(
                #    self.ui.list_pred_analysis_mono_overview, self.ui.cb_pred_analysis_mono_images,
                #    self.status_bar, self.app_project, self.app_settings, self._init_mono_pred_analysis_page)
                constants.PYSSA_LOGGER.info("Thread started for analysis process.")
                # self.threadpool.start(self.worker_analysis)
                constants.PYSSA_LOGGER.debug(
                    f"Thread count after analysis worker: {self.threadpool.activeThreadCount()}",
                )

                # <editor-fold desc="Worker setup">
                # TODO: test code below
                # --Begin: worker setup
                self.tmp_thread = QtCore.QThread()
                self.tmp_worker = task_workers.DistanceAnalysisWorker(
                    self.ui.list_pred_analysis_mono_overview,
                    self.ui.cb_pred_analysis_mono_images,
                    self.status_bar,
                    self.app_project,
                    self.app_settings,
                    self._init_mono_pred_analysis_page,
                )
                self.tmp_thread = task_workers.setup_worker_for_work(
                    self.tmp_thread, self.tmp_worker, self.display_view_page,
                )
                self.tmp_worker.finished.connect(self.post_analysis_process)
                self.tmp_thread.start()
                # --End: worker setup

                # </editor-fold>

                if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
                    os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
                self.block_box_prediction.destroy(True)
                self.block_box_analysis.exec_()
                self.display_view_page()
                self._project_watcher.show_valid_options(self.ui)
            elif self.prediction_type == constants.PREDICTION_TYPE_PRED_MULTI_ANALYSIS:
                self.app_project.serialize_project(self.app_project.get_project_xml_path())
                constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
                constants.PYSSA_LOGGER.info("All structure predictions are done.")
                constants.PYSSA_LOGGER.info("Begin analysis process.")
                constants.PYSSA_LOGGER.debug(
                    f"Thread count before analysis worker: {self.threadpool.activeThreadCount()}",
                )

                # self.worker_analysis = workers.AnalysisWorkerPool(
                #    self.ui.list_pred_analysis_multi_overview, self.ui.cb_pred_analysis_multi_images,
                #    self.status_bar, self.app_project, self.app_settings, self._init_multi_pred_analysis_page)
                constants.PYSSA_LOGGER.info("Thread started for analysis process.")
                # self.threadpool.start(self.worker_analysis)
                constants.PYSSA_LOGGER.debug(
                    f"Thread count after analysis worker: {self.threadpool.activeThreadCount()}",
                )

                # <editor-fold desc="Worker setup">
                # TODO: test code below
                # --Begin: worker setup
                self.tmp_thread = QtCore.QThread()
                self.tmp_worker = task_workers.DistanceAnalysisWorker(
                    self.ui.list_pred_analysis_multi_overview,
                    self.ui.cb_pred_analysis_multi_images,
                    self.status_bar,
                    self.app_project,
                    self.app_settings,
                    self._init_multi_pred_analysis_page,
                )
                self.tmp_thread = task_workers.setup_worker_for_work(
                    self.tmp_thread, self.tmp_worker, self.display_view_page,
                )
                self.tmp_worker.finished.connect(self.post_analysis_process)
                self.tmp_thread.start()
                # --End: worker setup

                # </editor-fold>

                if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
                    os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
                self.block_box_prediction.destroy(True)
                self.block_box_analysis.exec_()
                self.display_view_page()
                self._project_watcher.show_valid_options(self.ui)

    def predict_local_monomer(self) -> None:
        """Sets tup the worker for the prediction with the colabfold."""
        self.prediction_type = constants.PREDICTION_TYPE_PRED
        constants.PYSSA_LOGGER.info("Begin prediction process.")

        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.ColabfoldWorker(
            self.ui.table_pred_mono_prot_to_predict, self.prediction_configuration, self.app_project,
        )
        self.tmp_worker.finished.connect(self.post_prediction_process)
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_view_page)
        self.tmp_thread.start()
        constants.PYSSA_LOGGER.info("Thread started for prediction process.")
        # --End: worker setup

        # </editor-fold>

        gui_elements_to_show = []
        gui_elements_to_hide = [
            self.ui.btn_prediction_abort,
            self.ui.btn_use_page,
            self.ui.btn_close_project,
            self.ui.btn_pred_local_monomer_page,
            self.ui.btn_pred_local_multimer_page,
        ]
        gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

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
            return
        else:
            print("Unexpected Error.")
            self.block_box_prediction.close()

    def abort_prediction(self) -> None:
        """Aborts the running prediction."""
        constants.PYSSA_LOGGER.info("Structure prediction process was aborted manually.")
        subprocess.run(["wsl", "--shutdown"])
        constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
        filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
        constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
        basic_boxes.ok("Abort prediction", "The structure prediction was aborted.", QtWidgets.QMessageBox.Information)
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_prediction_abort)
        self._project_watcher.show_valid_options(self.ui)

    # </editor-fold>

    # <editor-fold desc="Multimer Local Prediction functions">
    def _init_local_pred_multi_page(self) -> None:
        """Clears all text boxes and sets default values for the gui elements."""
        # clears everything
        self.ui.txt_pred_multi_prot_name.clear()
        self.ui.txt_pred_multi_prot_seq.clear()
        self.ui.list_pred_multi_prot_seq_overview.clear()
        # sets up defaults: Prediction
        self.ui.btn_pred_multi_next.setEnabled(False)
        self.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        self.ui.lbl_pred_multi_prot_name_status.setText("")
        self.ui.lbl_pred_multi_prot_seq_status.setText("")

    def display_local_pred_multi(self) -> None:
        """Displays the local prediction multimer page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return
        
        gui_elements_to_show = [
            self.ui.lbl_pred_multi_prot_to_predict,
            self.ui.table_pred_multi_prot_to_predict,
            self.ui.btn_pred_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.lbl_pred_multi_prot_name_status,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
            self.ui.lbl_pred_multi_prot_name,
            self.ui.txt_pred_multi_prot_name,
            self.ui.lbl_pred_multi_prot_seq,
            self.ui.txt_pred_multi_prot_seq,
            self.ui.lbl_pred_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add,
            self.ui.btn_pred_multi_prot_seq_add,
            self.ui.lbl_pred_multi_prot_seq_overview,
            self.ui.list_pred_multi_prot_seq_overview,
            self.ui.btn_pred_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_multi_prot_to_predict_2,
            self.ui.btn_pred_multi_back_2,
            self.ui.btn_pred_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_multi_advanced_config,
            self.ui.btn_pred_multi_advanced_config,
            self.ui.btn_pred_multi_predict,
        ]
        for i in range(self.ui.table_pred_multi_prot_to_predict.rowCount() - 1, -1, -1):
            self.ui.table_pred_multi_prot_to_predict.removeRow(i)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 20, "Local Multimer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_pred_local_multimer_page,
        )

    def local_pred_multi_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self.ui.txt_pred_multi_prot_name.text(), self.ui.table_pred_multi_prot_to_predict,
        ):
            self.ui.lbl_pred_multi_prot_name_status.setText("Protein name already used.")
            self.ui.btn_pred_multi_next.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_pred_multi_next)
        else:
            self.ui.btn_pred_multi_next.setEnabled(True)
            tools.validate_protein_name(
                self.ui.txt_pred_multi_prot_name, self.ui.lbl_pred_multi_prot_name_status, self.ui.btn_pred_multi_next,
            )

    def local_pred_multi_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self.ui.txt_pred_multi_prot_seq, self.ui.lbl_pred_multi_prot_seq_status, self.ui.btn_pred_multi_prot_seq_add,
        )

    def local_pred_multi_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self.ui.table_pred_multi_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_multi_predict)
            self.ui.btn_pred_multi_predict.setEnabled(False)
            gui_elements_to_show = [
                self.ui.lbl_pred_multi_prot_to_predict,
                self.ui.table_pred_multi_prot_to_predict,
                self.ui.btn_pred_multi_prot_to_predict_add,
            ]
            gui_elements_to_hide = [
                self.ui.btn_pred_multi_prot_to_predict_remove,
                self.ui.lbl_pred_multi_prot_name_status,
                self.ui.btn_pred_multi_back,
                self.ui.btn_pred_multi_next,
                self.ui.lbl_pred_multi_prot_name,
                self.ui.txt_pred_multi_prot_name,
                self.ui.lbl_pred_multi_prot_seq,
                self.ui.txt_pred_multi_prot_seq,
                self.ui.lbl_pred_multi_prot_seq_status,
                self.ui.lbl_pred_multi_prot_seq_add,
                self.ui.btn_pred_multi_prot_seq_add,
                self.ui.lbl_pred_multi_prot_seq_overview,
                self.ui.list_pred_multi_prot_seq_overview,
                self.ui.btn_pred_multi_prot_seq_overview_remove,
                self.ui.lbl_pred_multi_prot_to_predict_2,
                self.ui.btn_pred_multi_back_2,
                self.ui.btn_pred_multi_prot_to_predict_add_2,
                self.ui.lbl_pred_multi_advanced_config,
                self.ui.btn_pred_multi_advanced_config,
                self.ui.btn_pred_multi_predict,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_multi_predict)
            self.ui.btn_pred_multi_predict.setEnabled(True)
            gui_elements_to_show = [
                self.ui.lbl_pred_multi_prot_to_predict,
                self.ui.table_pred_multi_prot_to_predict,
                self.ui.btn_pred_multi_prot_to_predict_add,
                self.ui.btn_pred_multi_prot_to_predict_remove,
                self.ui.lbl_pred_multi_advanced_config,
                self.ui.btn_pred_multi_advanced_config,
                self.ui.btn_pred_multi_predict,
            ]
            gui_elements_to_hide = [
                self.ui.lbl_pred_multi_prot_name_status,
                self.ui.btn_pred_multi_back,
                self.ui.btn_pred_multi_next,
                self.ui.lbl_pred_multi_prot_name,
                self.ui.txt_pred_multi_prot_name,
                self.ui.lbl_pred_multi_prot_seq,
                self.ui.txt_pred_multi_prot_seq,
                self.ui.lbl_pred_multi_prot_seq_status,
                self.ui.lbl_pred_multi_prot_seq_add,
                self.ui.btn_pred_multi_prot_seq_add,
                self.ui.lbl_pred_multi_prot_seq_overview,
                self.ui.list_pred_multi_prot_seq_overview,
                self.ui.btn_pred_multi_prot_seq_overview_remove,
                self.ui.lbl_pred_multi_prot_to_predict_2,
                self.ui.btn_pred_multi_back_2,
                self.ui.btn_pred_multi_prot_to_predict_add_2,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)

    def local_pred_multi_add_sequence_to_list(self) -> None:
        """Adds the entered sequence to the list of sequences of the protein."""
        self.ui.list_pred_multi_prot_seq_overview.addItem(
            QtWidgets.QListWidgetItem(self.ui.txt_pred_multi_prot_seq.toPlainText()),
        )
        self.local_pred_multi_check_if_list_is_empty()

    def local_pred_multi_remove_sequence_to_list(self) -> None:
        """Removes the entered sequence to the list of sequences of the protein."""
        self.ui.list_pred_multi_prot_seq_overview.takeItem(self.ui.list_pred_multi_prot_seq_overview.currentRow())
        self.local_pred_multi_check_if_list_is_empty()
        self.ui.btn_pred_multi_prot_seq_overview_remove.setEnabled(False)

    def local_pred_multi_check_if_list_is_empty(self) -> None:
        """Checks if the list of sequences of the protein is empty."""
        if self.ui.list_pred_multi_prot_seq_overview.count() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_multi_prot_to_predict_add_2)
            self.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_multi_prot_to_predict_add_2)
            self.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(True)

    def local_pred_multi_add(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_multi_prot_to_predict,
            self.ui.table_pred_multi_prot_to_predict,
            self.ui.lbl_pred_multi_prot_name,
            self.ui.txt_pred_multi_prot_name,
            self.ui.lbl_pred_multi_prot_name_status,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.lbl_pred_multi_prot_seq,
            self.ui.txt_pred_multi_prot_seq,
            self.ui.lbl_pred_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add,
            self.ui.btn_pred_multi_prot_seq_add,
            self.ui.lbl_pred_multi_prot_seq_overview,
            self.ui.list_pred_multi_prot_seq_overview,
            self.ui.btn_pred_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_multi_prot_to_predict_2,
            self.ui.btn_pred_multi_back_2,
            self.ui.btn_pred_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_multi_advanced_config,
            self.ui.btn_pred_multi_advanced_config,
            self.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_multi_prot_name, self.ui.lbl_pred_multi_prot_name)
        gui_utils.disable_text_box(self.ui.txt_pred_multi_prot_seq, self.ui.lbl_pred_multi_prot_seq)
        self.ui.btn_pred_multi_next.setEnabled(False)
        self.ui.txt_pred_multi_prot_name.clear()
        styles.color_button_not_ready(self.ui.btn_pred_multi_next)
        if self.ui.table_pred_multi_prot_to_predict.rowCount() > 0:
            try:
                self.ui.table_pred_multi_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Multimer Prediction in overview table.")

    def local_pred_multi_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_multi_prot_to_predict,
            self.ui.table_pred_multi_prot_to_predict,
            self.ui.btn_pred_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.lbl_pred_multi_prot_name,
            self.ui.txt_pred_multi_prot_name,
            self.ui.lbl_pred_multi_prot_name_status,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
            self.ui.lbl_pred_multi_prot_seq,
            self.ui.txt_pred_multi_prot_seq,
            self.ui.lbl_pred_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add,
            self.ui.btn_pred_multi_prot_seq_add,
            self.ui.lbl_pred_multi_prot_seq_overview,
            self.ui.list_pred_multi_prot_seq_overview,
            self.ui.btn_pred_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_multi_prot_to_predict_2,
            self.ui.btn_pred_multi_back_2,
            self.ui.btn_pred_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_multi_advanced_config,
            self.ui.btn_pred_multi_advanced_config,
            self.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.local_pred_multi_check_if_table_is_empty()

    def local_pred_multi_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_multi_prot_to_predict,
            self.ui.table_pred_multi_prot_to_predict,
            self.ui.lbl_pred_multi_prot_name,
            self.ui.txt_pred_multi_prot_name,
            self.ui.lbl_pred_multi_prot_seq,
            self.ui.txt_pred_multi_prot_seq,
            self.ui.lbl_pred_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add,
            self.ui.btn_pred_multi_prot_seq_add,
            self.ui.lbl_pred_multi_prot_seq_overview,
            self.ui.list_pred_multi_prot_seq_overview,
            self.ui.btn_pred_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_multi_prot_to_predict_2,
            self.ui.btn_pred_multi_back_2,
            self.ui.btn_pred_multi_prot_to_predict_add_2,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.lbl_pred_multi_prot_name_status,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
            self.ui.lbl_pred_multi_advanced_config,
            self.ui.btn_pred_multi_advanced_config,
            self.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_multi_prot_seq, self.ui.lbl_pred_multi_prot_seq)
        gui_utils.disable_text_box(self.ui.txt_pred_multi_prot_name, self.ui.lbl_pred_multi_prot_name)
        self.ui.txt_pred_multi_prot_seq.clear()
        self.ui.list_pred_multi_prot_seq_overview.clear()
        self.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        self.ui.btn_pred_multi_prot_seq_overview_remove.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_pred_multi_prot_to_predict_add_2)

    def local_pred_multi_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_multi_prot_to_predict,
            self.ui.table_pred_multi_prot_to_predict,
            self.ui.lbl_pred_multi_prot_name_status,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
            self.ui.lbl_pred_multi_prot_name,
            self.ui.txt_pred_multi_prot_name,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.lbl_pred_multi_prot_seq,
            self.ui.txt_pred_multi_prot_seq,
            self.ui.lbl_pred_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add,
            self.ui.btn_pred_multi_prot_seq_add,
            self.ui.lbl_pred_multi_prot_seq_overview,
            self.ui.list_pred_multi_prot_seq_overview,
            self.ui.btn_pred_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_multi_prot_to_predict_2,
            self.ui.btn_pred_multi_back_2,
            self.ui.btn_pred_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_multi_advanced_config,
            self.ui.btn_pred_multi_advanced_config,
            self.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_multi_prot_name, self.ui.lbl_pred_multi_prot_name)
        gui_utils.disable_text_box(self.ui.txt_pred_multi_prot_seq, self.ui.lbl_pred_multi_prot_seq)

    def local_pred_multi_prot_seq_overview_item_changed(self) -> None:
        """Enables the remove button of the list sequences of the protein."""
        self.ui.btn_pred_multi_prot_seq_overview_remove.setEnabled(True)

    def local_pred_multi_prot_to_predict_item_changed(self) -> None:
        """Enables the remove button of the table proteins to predict."""
        self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(True)

    def local_pred_multi_prot_to_predict_add_2(self) -> None:
        """Adds the protein to the list of proteins to predict."""
        for i in range(self.ui.list_pred_multi_prot_seq_overview.count()):
            self.ui.table_pred_multi_prot_to_predict.setRowCount(
                self.ui.table_pred_multi_prot_to_predict.rowCount() + 1,
            )
            self.ui.table_pred_multi_prot_to_predict.insertRow(self.ui.table_pred_multi_prot_to_predict.rowCount() + 1)
            tmp_chain_seq = (constants.chain_dict.get(i), self.ui.list_pred_multi_prot_seq_overview.item(i).text())
            self.ui.table_pred_multi_prot_to_predict.setItem(
                self.ui.table_pred_multi_prot_to_predict.rowCount() - 1, 0, QtWidgets.QTableWidgetItem(tmp_chain_seq[0]),
            )
            self.ui.table_pred_multi_prot_to_predict.setItem(
                self.ui.table_pred_multi_prot_to_predict.rowCount() - 1, 1, QtWidgets.QTableWidgetItem(tmp_chain_seq[1]),
            )
            name_item = QtWidgets.QTableWidgetItem(self.ui.txt_pred_multi_prot_name.text())
            self.ui.table_pred_multi_prot_to_predict.setVerticalHeaderItem(
                self.ui.table_pred_multi_prot_to_predict.rowCount() - 1, name_item,
            )
        self.ui.table_pred_multi_prot_to_predict.resizeColumnsToContents()
        self.local_pred_multi_check_if_table_is_empty()
        gui_elements_to_show = [
            self.ui.lbl_pred_multi_prot_to_predict,
            self.ui.table_pred_multi_prot_to_predict,
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.lbl_pred_multi_advanced_config,
            self.ui.btn_pred_multi_advanced_config,
            self.ui.btn_pred_multi_predict,
        ]
        gui_elements_to_hide = [
            self.ui.lbl_pred_multi_prot_name_status,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
            self.ui.lbl_pred_multi_prot_name,
            self.ui.txt_pred_multi_prot_name,
            self.ui.lbl_pred_multi_prot_seq,
            self.ui.txt_pred_multi_prot_seq,
            self.ui.lbl_pred_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add,
            self.ui.btn_pred_multi_prot_seq_add,
            self.ui.lbl_pred_multi_prot_seq_overview,
            self.ui.list_pred_multi_prot_seq_overview,
            self.ui.btn_pred_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_multi_prot_to_predict_2,
            self.ui.btn_pred_multi_back_2,
            self.ui.btn_pred_multi_prot_to_predict_add_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._init_local_pred_multi_page()
        self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)

    def local_pred_multi_remove(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        self.ui.table_pred_multi_prot_to_predict.removeRow(self.ui.table_pred_multi_prot_to_predict.currentRow())
        if self.ui.table_pred_multi_prot_to_predict.rowCount() > 0:
            prot_name = self.ui.table_pred_multi_prot_to_predict.verticalHeaderItem(
                self.ui.table_pred_multi_prot_to_predict.currentRow(),
            ).text()
            for i in range(self.ui.table_pred_multi_prot_to_predict.rowCount()):
                if self.ui.table_pred_multi_prot_to_predict.verticalHeaderItem(i).text() == prot_name:
                    self.ui.table_pred_multi_prot_to_predict.setItem(
                        i, 0, QtWidgets.QTableWidgetItem(constants.chain_dict.get(i)),
                    )
        self.local_pred_multi_check_if_table_is_empty()
        gui_elements_to_show = [
            self.ui.lbl_pred_multi_prot_to_predict,
            self.ui.table_pred_multi_prot_to_predict,
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.lbl_pred_multi_advanced_config,
            self.ui.btn_pred_multi_advanced_config,
            self.ui.btn_pred_multi_predict,
        ]
        gui_elements_to_hide = [
            self.ui.lbl_pred_multi_prot_name_status,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
            self.ui.lbl_pred_multi_prot_name,
            self.ui.txt_pred_multi_prot_name,
            self.ui.lbl_pred_multi_prot_seq,
            self.ui.txt_pred_multi_prot_seq,
            self.ui.lbl_pred_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add,
            self.ui.btn_pred_multi_prot_seq_add,
            self.ui.lbl_pred_multi_prot_seq_overview,
            self.ui.list_pred_multi_prot_seq_overview,
            self.ui.btn_pred_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_multi_prot_to_predict_2,
            self.ui.btn_pred_multi_back_2,
            self.ui.btn_pred_multi_prot_to_predict_add_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        self.local_pred_multi_check_if_table_is_empty()

    # def local_pred_multi_show_protein_overview(self):
    #     if self.ui.table_pred_multi_prot_to_predict.rowCount() == 0:
    #         self.local_pred_multimer_management.show_stage_x(0)
    #     else:
    #         gui_elements_to_show = [
    #             self.ui.btn_pred_multi_prot_to_predict_add,
    #             self.ui.btn_pred_multi_prot_to_predict_remove,
    #         ]
    #         self.local_pred_multimer_management.show_gui_elements_stage_x(
    #             [0, 3], [1, 2], show_specific_elements=gui_elements_to_show
    #         )
    #
    # def local_pred_multi_show_protein_name(self):
    #     gui_elements_to_hide = [
    #         self.ui.btn_pred_multi_prot_to_predict_add,
    #         self.ui.btn_pred_multi_prot_to_predict_remove,
    #     ]
    #     self.local_pred_multimer_management.show_gui_elements_stage_x(
    #         [0, 1], [2, 3], hide_specific_elements=gui_elements_to_hide)
    #     gui_utils.enable_text_box(self.ui.txt_pred_multi_prot_name,
    #                               self.ui.lbl_pred_multi_prot_name)
    #
    # def local_pred_multi_show_protein_sequence(self):
    #     gui_elements_to_hide = [
    #         self.ui.btn_pred_multi_prot_to_predict_add,
    #         self.ui.btn_pred_multi_prot_to_predict_remove,
    #         self.ui.btn_pred_multi_back,
    #         self.ui.btn_pred_multi_next,
    #     ]
    #     self.local_pred_multimer_management.show_gui_elements_stage_x(
    #         [0, 1, 2], [3], hide_specific_elements=gui_elements_to_hide)
    #     gui_utils.disable_text_box(self.ui.txt_pred_multi_prot_name,
    #                                self.ui.lbl_pred_multi_prot_name)

    def predict_local_multimer(self) -> None:
        """Sets up the worker for the prediction with the colabfold."""
        self.prediction_type = constants.PREDICTION_TYPE_PRED
        constants.PYSSA_LOGGER.info("Begin multimer prediction process.")
        # worker = workers.PredictionWorkerPool(self.ui.table_pred_multi_prot_to_predict,
        #                                      self.prediction_configuration, self.app_project)
        # worker.signals.finished.connect(self.post_prediction_process)
        constants.PYSSA_LOGGER.info("Thread started for prediction process.")
        # self.threadpool.start(worker)

        # <editor-fold desc="Worker setup">
        # TODO: test code below
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.ColabfoldWorker(
            self.ui.table_pred_multi_prot_to_predict, self.prediction_configuration, self.app_project,
        )
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_view_page)
        self.tmp_worker.finished.connect(self.post_prediction_process)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        gui_elements_to_show = [
            self.ui.btn_prediction_abort,
        ]
        gui_elements_to_hide = [
            self.ui.btn_use_page,
            self.ui.btn_close_project,
            self.ui.btn_pred_local_monomer_page,
            self.ui.btn_pred_local_multimer_page,
        ]
        gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

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
            return
        else:
            print("Unexpected Error.")
            self.block_box_prediction.close()

    # def validate_local_pred_multi(self):
    #     self.local_pred_multimer_management.create_validation()
    #
    # def show_local_pred_multi_stage_prediction_mode(self):
    #     self.local_pred_multimer_management.show_stage_x(0)
    #
    # # --- single prediction
    # def show_local_pred_multi_stage_protein_name(self):
    #     if self.ui.table_local_pred_multi_prot_overview.isVisible():
    #         self.local_pred_multimer_management.show_stage_x(0)
    #     else:
    #         self.local_pred_multimer_management.show_stage_x(1)
    #
    # def show_local_pred_multi_stage_protein_sequence_single(self):
    #     gui_elements_to_show = [
    #         self.ui.btn_local_pred_multi_add_seq_single
    #     ]
    #     gui_elements_to_hide = [
    #         self.ui.btn_local_pred_multi_next,
    #         self.ui.btn_local_pred_multi_back_prediction_mode,
    #     ]
    #     self.local_pred_multimer_management.show_gui_elements_stage_x([1, 2], [0, 3, 4, 5],
    #                                                                   hide_specific_elements=gui_elements_to_hide,
    #                                                                   show_specific_elements=gui_elements_to_show)
    #     self.local_pred_multimer_management.disable_text_boxes_stage_x([1])
    #     self.local_pred_multimer_management.enable_text_boxes_stage_x([2])
    #
    # def show_local_pred_multi_stage_prediction_single(self):
    #     gui_elements_to_hide = [
    #         self.ui.btn_local_pred_multi_add_seq_single
    #     ]
    #     self.local_pred_multimer_management.show_gui_elements_stage_x([1, 2, 4], [0, 3, 5], hide_specific_elements=gui_elements_to_hide)
    #     self.local_pred_multimer_management.disable_text_boxes_stage_x([1, 2])
    #     self.local_pred_multimer_management.activate_specific_button(self.ui.btn_local_pred_multi_predict)
    #
    # # --- batch prediction
    # def show_local_pred_multi_stage_protein_sequence_batch(self):
    #     self.local_pred_multimer_management.show_gui_elements_stage_x([3], [0, 1, 2, 4, 5])

    # </editor-fold>

    # <editor-fold desc="Monomer Prediction + Analysis functions">
    def _init_mono_pred_analysis_page(self) -> None:
        """Clears all text boxes and sets default values for the gui elements."""
        # <editor-fold desc="Prediction section">
        self.ui.txt_pred_analysis_mono_prot_name.clear()
        self.ui.txt_pred_analysis_mono_seq_name.clear()
        for i in range(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            self.ui.table_pred_analysis_mono_prot_to_predict.removeRow(i)
        # sets up defaults: Prediction
        self.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self.ui.btn_pred_analysis_mono_add_protein.setEnabled(False)
        self.ui.lbl_pred_analysis_mono_prot_name_status.setText("")
        self.ui.lbl_pred_analysis_mono_seq_name_status.setText("")

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self.ui.list_pred_analysis_mono_overview.clear()
        self.ui.btn_pred_analysis_mono_remove.hide()

        # </editor-fold>

    def display_monomer_pred_analysis(self) -> None:
        """Displays the monomer prediction + analysis page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return
        
        self._init_mono_pred_analysis_page()
        self.ui.table_pred_analysis_mono_prot_to_predict.clear()
        self.ui.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(0, QtWidgets.QTableWidgetItem("Chain"))
        self.ui.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(
            1, QtWidgets.QTableWidgetItem("Sequence"),
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_prot_to_predict,
            self.ui.table_pred_analysis_mono_prot_to_predict,
            self.ui.btn_pred_analysis_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self.ui.lbl_pred_analysis_mono_prot_name,
            self.ui.txt_pred_analysis_mono_prot_name,
            self.ui.lbl_pred_analysis_mono_prot_name_status,
            self.ui.btn_pred_analysis_mono_back,
            self.ui.btn_pred_analysis_mono_next,
            self.ui.lbl_pred_analysis_mono_seq_name,
            self.ui.txt_pred_analysis_mono_seq_name,
            self.ui.lbl_pred_analysis_mono_seq_name_status,
            self.ui.btn_pred_analysis_mono_back_2,
            self.ui.btn_pred_analysis_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config_2,
            self.ui.btn_pred_mono_advanced_config_2,
            self.ui.btn_pred_analysis_mono_go_analysis_setup,
            self.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self.ui.tabWidget.currentIndex() == 1:
            self.ui.tabWidget.setCurrentIndex(0)
        self.ui.tabWidget.setTabEnabled(1, False)
        self.ui.tabWidget.setTabEnabled(0, True)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 21, "Monomer Prediction + Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_pred_analysis_monomer_page,
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setEnabled(True)

    # <editor-fold desc="Sections">

    # <editor-fold desc="Prediction section">
    def mono_pred_analysis_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self.ui.txt_pred_analysis_mono_prot_name.text(), self.ui.table_pred_analysis_mono_prot_to_predict,
        ):
            self.ui.lbl_pred_analysis_mono_prot_name_status.setText("Protein name already used.")
            self.ui.btn_pred_analysis_mono_next.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_pred_analysis_mono_next)
        else:
            self.ui.btn_pred_analysis_mono_next.setEnabled(True)
            tools.validate_protein_name(
                self.ui.txt_pred_analysis_mono_prot_name,
                self.ui.lbl_pred_analysis_mono_prot_name_status,
                self.ui.btn_pred_analysis_mono_next,
            )

    def mono_pred_analysis_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self.ui.txt_pred_analysis_mono_seq_name,
            self.ui.lbl_pred_analysis_mono_seq_name_status,
            self.ui.btn_pred_analysis_mono_add_protein,
        )

    def mono_pred_analysis_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_analysis_mono_go_analysis_setup)
            self.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(False)
            gui_elements_to_show = [
                self.ui.lbl_pred_analysis_mono_prot_to_predict,
                self.ui.table_pred_analysis_mono_prot_to_predict,
                self.ui.btn_pred_analysis_mono_seq_to_predict,
            ]
            gui_elements_to_hide = [
                self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                self.ui.lbl_pred_analysis_mono_prot_name,
                self.ui.txt_pred_analysis_mono_prot_name,
                self.ui.lbl_pred_analysis_mono_prot_name_status,
                self.ui.btn_pred_analysis_mono_back,
                self.ui.btn_pred_analysis_mono_next,
                self.ui.lbl_pred_analysis_mono_seq_name,
                self.ui.txt_pred_analysis_mono_seq_name,
                self.ui.lbl_pred_analysis_mono_seq_name_status,
                self.ui.btn_pred_analysis_mono_back_2,
                self.ui.btn_pred_analysis_mono_add_protein,
                self.ui.lbl_pred_mono_advanced_config_2,
                self.ui.btn_pred_mono_advanced_config_2,
                self.ui.btn_pred_analysis_mono_go_analysis_setup,
                self.ui.lbl_pred_analysis_mono_to_analysis_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            gui_elements_to_show = [
                self.ui.lbl_pred_analysis_mono_prot_to_predict,
                self.ui.table_pred_analysis_mono_prot_to_predict,
                self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                self.ui.btn_pred_analysis_mono_seq_to_predict,
                self.ui.lbl_pred_mono_advanced_config_2,
                self.ui.btn_pred_mono_advanced_config_2,
                self.ui.btn_pred_analysis_mono_go_analysis_setup,
                self.ui.lbl_pred_analysis_mono_to_analysis_setup,
            ]
            gui_elements_to_hide = [
                self.ui.lbl_pred_analysis_mono_prot_name,
                self.ui.txt_pred_analysis_mono_prot_name,
                self.ui.lbl_pred_analysis_mono_prot_name_status,
                self.ui.btn_pred_analysis_mono_back,
                self.ui.btn_pred_analysis_mono_next,
                self.ui.lbl_pred_analysis_mono_seq_name,
                self.ui.txt_pred_analysis_mono_seq_name,
                self.ui.lbl_pred_analysis_mono_seq_name_status,
                self.ui.btn_pred_analysis_mono_back_2,
                self.ui.btn_pred_analysis_mono_add_protein,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            styles.color_button_ready(self.ui.btn_pred_analysis_mono_go_analysis_setup)
            self.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(True)

    def setup_defaults_monomer_prediction_analysis(self) -> None:
        """Sets up default values for the prediction tab."""
        # clears everything
        self.ui.txt_pred_analysis_mono_prot_name.clear()
        self.ui.txt_pred_analysis_mono_seq_name.clear()
        # sets up defaults: Prediction
        self.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self.ui.btn_pred_analysis_mono_add_protein.setEnabled(False)
        self.ui.lbl_pred_analysis_mono_prot_name_status.setText("")
        self.ui.lbl_pred_analysis_mono_seq_name_status.setText("")

    def mono_pred_analysis_add_seq_to_predict(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_prot_to_predict,
            self.ui.table_pred_analysis_mono_prot_to_predict,
            self.ui.lbl_pred_analysis_mono_prot_name,
            self.ui.txt_pred_analysis_mono_prot_name,
            self.ui.lbl_pred_analysis_mono_prot_name_status,
            self.ui.btn_pred_analysis_mono_back,
            self.ui.btn_pred_analysis_mono_next,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_mono_prot_name, self.ui.lbl_pred_analysis_mono_prot_name)
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self.ui.btn_pred_analysis_mono_seq_to_predict,
            self.ui.lbl_pred_analysis_mono_seq_name,
            self.ui.txt_pred_analysis_mono_seq_name,
            self.ui.lbl_pred_analysis_mono_seq_name_status,
            self.ui.btn_pred_analysis_mono_back_2,
            self.ui.btn_pred_analysis_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config_2,
            self.ui.btn_pred_mono_advanced_config_2,
            self.ui.btn_pred_analysis_mono_go_analysis_setup,
            self.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_mono_seq_name, self.ui.lbl_pred_analysis_mono_seq_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self.ui.txt_pred_analysis_mono_prot_name.clear()
        styles.color_button_not_ready(self.ui.btn_pred_analysis_mono_next)
        if self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() > 0:
            try:
                self.ui.table_pred_analysis_mono_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Monomer Prediction in overview table.")

    def mono_pred_analysis_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_prot_to_predict,
            self.ui.table_pred_analysis_mono_prot_to_predict,
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self.ui.btn_pred_analysis_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self.ui.lbl_pred_analysis_mono_prot_name,
            self.ui.txt_pred_analysis_mono_prot_name,
            self.ui.lbl_pred_analysis_mono_prot_name_status,
            self.ui.btn_pred_analysis_mono_back,
            self.ui.btn_pred_analysis_mono_next,
            self.ui.lbl_pred_analysis_mono_seq_name,
            self.ui.txt_pred_analysis_mono_seq_name,
            self.ui.lbl_pred_analysis_mono_seq_name_status,
            self.ui.btn_pred_analysis_mono_back_2,
            self.ui.btn_pred_analysis_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config_2,
            self.ui.btn_pred_mono_advanced_config_2,
            self.ui.btn_pred_analysis_mono_go_analysis_setup,
            self.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.mono_pred_analysis_check_if_table_is_empty()
        self.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)

    def mono_pred_analysis_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_prot_to_predict,
            self.ui.table_pred_analysis_mono_prot_to_predict,
            self.ui.lbl_pred_analysis_mono_prot_name,
            self.ui.txt_pred_analysis_mono_prot_name,
            self.ui.lbl_pred_analysis_mono_seq_name,
            self.ui.txt_pred_analysis_mono_seq_name,
            self.ui.lbl_pred_analysis_mono_seq_name_status,
            self.ui.btn_pred_analysis_mono_back_2,
            self.ui.btn_pred_analysis_mono_add_protein,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_mono_seq_name, self.ui.lbl_pred_analysis_mono_seq_name)
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self.ui.btn_pred_analysis_mono_seq_to_predict,
            self.ui.lbl_pred_analysis_mono_prot_name_status,
            self.ui.btn_pred_analysis_mono_back,
            self.ui.btn_pred_analysis_mono_next,
            self.ui.lbl_pred_mono_advanced_config_2,
            self.ui.btn_pred_mono_advanced_config_2,
            self.ui.btn_pred_analysis_mono_go_analysis_setup,
            self.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_mono_prot_name, self.ui.lbl_pred_analysis_mono_prot_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.txt_pred_analysis_mono_seq_name.clear()

    def mono_pred_analysis_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_prot_to_predict,
            self.ui.table_pred_analysis_mono_prot_to_predict,
            self.ui.lbl_pred_analysis_mono_prot_name,
            self.ui.txt_pred_analysis_mono_prot_name,
            self.ui.lbl_pred_analysis_mono_prot_name_status,
            self.ui.btn_pred_analysis_mono_back,
            self.ui.btn_pred_analysis_mono_next,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self.ui.btn_pred_analysis_mono_seq_to_predict,
            self.ui.lbl_pred_analysis_mono_seq_name,
            self.ui.txt_pred_analysis_mono_seq_name,
            self.ui.lbl_pred_analysis_mono_seq_name_status,
            self.ui.btn_pred_analysis_mono_back_2,
            self.ui.btn_pred_analysis_mono_add_protein,
            self.ui.lbl_pred_mono_advanced_config_2,
            self.ui.btn_pred_mono_advanced_config_2,
            self.ui.btn_pred_analysis_mono_go_analysis_setup,
            self.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_mono_prot_name, self.ui.lbl_pred_analysis_mono_prot_name)
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_mono_seq_name, self.ui.lbl_pred_analysis_mono_seq_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def mono_pred_analysis_add_protein(self) -> None:
        """Adds the protein to the list of proteins to predict."""
        self.ui.table_pred_analysis_mono_prot_to_predict.setRowCount(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.insertRow(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1, 0, QtWidgets.QTableWidgetItem("A"),
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self.ui.txt_pred_analysis_mono_seq_name.toPlainText()),
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setVerticalHeaderItem(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self.ui.txt_pred_analysis_mono_prot_name.text()),
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        self.mono_pred_analysis_check_if_table_is_empty()
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_prot_to_predict,
            self.ui.table_pred_analysis_mono_prot_to_predict,
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self.ui.btn_pred_analysis_mono_seq_to_predict,
            self.ui.lbl_pred_mono_advanced_config_2,
            self.ui.btn_pred_mono_advanced_config_2,
            self.ui.btn_pred_analysis_mono_go_analysis_setup,
            self.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_mono_prot_name, self.ui.lbl_pred_analysis_mono_prot_name)
        gui_elements_to_hide = [
            self.ui.lbl_pred_analysis_mono_prot_name,
            self.ui.txt_pred_analysis_mono_prot_name,
            self.ui.lbl_pred_analysis_mono_prot_name_status,
            self.ui.btn_pred_analysis_mono_back,
            self.ui.btn_pred_analysis_mono_next,
            self.ui.lbl_pred_analysis_mono_seq_name,
            self.ui.txt_pred_analysis_mono_seq_name,
            self.ui.lbl_pred_analysis_mono_seq_name_status,
            self.ui.btn_pred_analysis_mono_back_2,
            self.ui.btn_pred_analysis_mono_add_protein,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(True)
        self.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)
        styles.color_button_ready(self.ui.btn_pred_analysis_mono_go_analysis_setup)
        self.setup_defaults_monomer_prediction()

    def mono_pred_analysis_prediction_overview_item_clicked(self) -> None:
        """Enables the remove button."""
        self.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(True)

    def mono_pred_analysis_add_protein_to_predict(self) -> None:
        """Needs to be removed."""
        self.ui.table_pred_analysis_mono_prot_to_predict.setRowCount(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.insertRow(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1, 0, QtWidgets.QTableWidgetItem("A"),
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self.ui.txt_pred_analysis_mono_seq_name.toPlainText()),
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.setVerticalHeaderItem(
            self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self.ui.txt_pred_analysis_mono_prot_name.text()),
        )
        self.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        self.mono_pred_analysis_check_if_table_is_empty()
        self.setup_defaults_monomer_prediction_analysis()

    def mono_pred_analysis_remove_protein_to_predict(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        self.ui.table_pred_analysis_mono_prot_to_predict.removeRow(
            self.ui.table_pred_analysis_mono_prot_to_predict.currentRow(),
        )
        self.mono_pred_analysis_check_if_table_is_empty()
        self.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)

    # </editor-fold>

    # <editor-fold desc="Analysis section">
    def mono_pred_analysis_structure_analysis_add(self) -> None:
        """Shows the gui elements for the selection of the two proteins."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_pred_analysis_mono_prot_struct_1.clear()
        self.ui.lbl_pred_analysis_mono_prot_struct_2.clear()
        self.ui.lbl_pred_analysis_mono_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_pred_analysis_mono_prot_struct_2.setText("Protein structure 2")
        self.fill_mono_pred_analysis_protein_boxes()
        if self.ui.list_pred_analysis_mono_overview.count() > 0:
            try:
                self.ui.list_pred_analysis_mono_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def mono_pred_analysis_structure_analysis_next_2(self) -> None:
        """Shows the gui elements for the chain selection of protein 1."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_pred_analysis_mono_prot_struct_1.setText(self.ui.box_pred_analysis_mono_prot_struct_1.currentText())
        self.ui.lbl_pred_analysis_mono_prot_struct_2.setText(self.ui.box_pred_analysis_mono_prot_struct_2.currentText())
        self.ui.list_pred_analysis_mono_ref_chains.clear()
        self.ui.btn_pred_analysis_mono_next_3.setEnabled(False)
        self.ui.list_pred_analysis_mono_ref_chains.setEnabled(True)

        for i in range(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            if (
                self.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text()
                == self.ui.box_pred_analysis_mono_prot_struct_1.currentText()
            ):
                self.ui.list_pred_analysis_mono_ref_chains.addItem(
                    self.ui.table_pred_analysis_mono_prot_to_predict.item(i, 0).text(),
                )
        if self.ui.list_pred_analysis_mono_ref_chains.count() == 0:
            tmp_protein = self.app_project.search_protein(self.ui.box_pred_analysis_mono_prot_struct_1.currentText())
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self.ui.list_pred_analysis_mono_ref_chains.addItem(tmp_chain.chain_letter)
        if self.ui.list_pred_analysis_mono_ref_chains.count() == 1:
            self.ui.lbl_pred_analysis_mono_ref_chains.setText(
                f"Select chain in protein structure {self.ui.lbl_pred_analysis_mono_prot_struct_1.text()}.",
            )
        else:
            self.ui.lbl_pred_analysis_mono_ref_chains.setText(
                f"Select chains in protein structure {self.ui.lbl_pred_analysis_mono_prot_struct_1.text()}.",
            )

    def mono_pred_analysis_structure_analysis_back_3(self) -> None:
        """Hides the gui elements to select the two proteins."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self.ui.list_pred_analysis_mono_overview.count() > 0:
            self.ui.btn_pred_analysis_mono_remove.show()
            self.ui.btn_pred_analysis_mono_remove.setEnabled(False)
            self.ui.btn_pred_analysis_mono_start.show()
            self.ui.lbl_pred_analysis_mono_images.show()
            self.ui.cb_pred_analysis_mono_images.show()
            styles.color_button_ready(self.ui.btn_pred_analysis_mono_start)

    def mono_pred_analysis_structure_analysis_next_3(self) -> None:
        """Shows the gui elements to select the chains of protein 2."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.list_pred_analysis_mono_model_chains.clear()
        self.ui.list_pred_analysis_mono_ref_chains.setEnabled(False)
        self.ui.btn_pred_analysis_mono_next_4.setEnabled(False)

        for i in range(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            if (
                self.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text()
                == self.ui.box_pred_analysis_mono_prot_struct_2.currentText()
            ):
                self.ui.list_pred_analysis_mono_model_chains.addItem(
                    self.ui.table_pred_analysis_mono_prot_to_predict.item(i, 0).text(),
                )
        if self.ui.list_pred_analysis_mono_model_chains.count() == 0:
            tmp_protein = self.app_project.search_protein(self.ui.box_pred_analysis_mono_prot_struct_2.currentText())
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self.ui.list_pred_analysis_mono_model_chains.addItem(tmp_chain.chain_letter)
        if self.ui.list_pred_analysis_mono_model_chains.count() == 1:
            self.ui.lbl_pred_analysis_mono_model_chains.setText(
                f"Select chain in protein structure {self.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.",
            )
        else:
            self.ui.lbl_pred_analysis_mono_model_chains.setText(
                f"Select {len(self.ui.list_pred_analysis_mono_model_chains.selectedItems())} chains in protein structure {self.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.",
            )

        # tmp_protein = self.app_project.search_protein(self.ui.box_pred_analysis_mono_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self.ui.list_pred_analysis_mono_model_chains.addItem(tmp_chain.chain_letter)
        # if len(self.ui.list_pred_analysis_mono_ref_chains.selectedItems()) == 1:
        #     self.ui.lbl_pred_analysis_mono_model_chains.setText(
        #         f"Select 1 chain in protein structure {self.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.")
        # else:
        #     self.ui.lbl_pred_analysis_mono_model_chains.setText(
        #         f"Select {len(self.ui.list_pred_analysis_mono_ref_chains.selectedItems())} chains in protein structure {self.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.")

    def mono_pred_analysis_structure_analysis_back_4(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_pred_analysis_mono_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_pred_analysis_mono_prot_struct_2.setText("Protein structure 2")

    def mono_pred_analysis_structure_analysis_next_4(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self.ui.lbl_pred_analysis_mono_prot_struct_1.text()
        prot_1_chains = []
        for chain in self.ui.list_pred_analysis_mono_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self.ui.lbl_pred_analysis_mono_prot_struct_2.text()
        prot_2_chains = []
        for chain in self.ui.list_pred_analysis_mono_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self.ui.list_pred_analysis_mono_overview.addItem(item)
        self.ui.btn_pred_analysis_mono_remove.setEnabled(False)
        styles.color_button_ready(self.ui.btn_pred_analysis_mono_start)

    def mono_pred_analysis_structure_analysis_back_5(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.list_pred_analysis_mono_ref_chains.setEnabled(True)

        # tmp_protein = self.app_project.search_protein(self.ui.box_pred_analysis_mono_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self.ui.list_pred_analysis_mono_ref_chains.addItem(tmp_chain.chain_letter)

    def mono_pred_analysis_structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self.ui.btn_pred_analysis_mono_remove.setEnabled(True)

    def fill_mono_pred_analysis_protein_boxes(self) -> None:
        """Fills the combo box of the protein structures."""
        protein_names = []
        for i in range(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            protein_names.append(self.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text())
        for tmp_protein in self.app_project.proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        self.ui.box_pred_analysis_mono_prot_struct_1.clear()
        self.ui.box_pred_analysis_mono_prot_struct_2.clear()
        gui_utils.fill_combo_box(self.ui.box_pred_analysis_mono_prot_struct_1, protein_names)
        gui_utils.fill_combo_box(self.ui.box_pred_analysis_mono_prot_struct_2, protein_names)

    def remove_mono_pred_analysis_analysis_run(self) -> None:
        """Removes a selected protein pair form the list of protein pairs to analyze."""
        self.ui.list_pred_analysis_mono_overview.takeItem(self.ui.list_pred_analysis_mono_overview.currentRow())
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_mono_overview,
            self.ui.list_pred_analysis_mono_overview,
            self.ui.btn_pred_analysis_mono_add,
            self.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_remove,
            self.ui.lbl_pred_analysis_mono_prot_struct_1,
            self.ui.lbl_pred_analysis_mono_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_2,
            self.ui.lbl_pred_analysis_mono_ref_chains,
            self.ui.list_pred_analysis_mono_ref_chains,
            self.ui.btn_pred_analysis_mono_back_4,
            self.ui.btn_pred_analysis_mono_next_3,
            self.ui.box_pred_analysis_mono_prot_struct_1,
            self.ui.box_pred_analysis_mono_prot_struct_2,
            self.ui.btn_pred_analysis_mono_back_3,
            self.ui.btn_pred_analysis_mono_next_2,
            self.ui.lbl_pred_analysis_mono_model_chains,
            self.ui.list_pred_analysis_mono_model_chains,
            self.ui.btn_pred_analysis_mono_back_5,
            self.ui.btn_pred_analysis_mono_next_4,
            self.ui.lbl_pred_analysis_mono_images,
            self.ui.cb_pred_analysis_mono_images,
            self.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self.ui.list_pred_analysis_mono_overview.count() > 0:
            self.ui.btn_pred_analysis_mono_remove.show()
            self.ui.btn_pred_analysis_mono_remove.setEnabled(False)
            self.ui.btn_pred_analysis_mono_start.show()
            self.ui.lbl_pred_analysis_mono_images.show()
            self.ui.cb_pred_analysis_mono_images.show()
            styles.color_button_ready(self.ui.btn_pred_analysis_mono_start)
        # if self.ui.list_pred_analysis_mono_overview.count() == 0:
        #
        #     self.ui.btn_pred_analysis_mono_back_pred_setup.show()
        #     self.ui.btn_pred_analysis_mono_remove.hide()

    def check_mono_pred_analysis_if_same_no_of_chains_selected(self) -> None:
        """Checks if the same number of chains were selected."""
        self.ui.btn_pred_analysis_mono_next_4.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_pred_analysis_mono_next_4)
        if self.no_of_selected_chains == len(self.ui.list_pred_analysis_mono_model_chains.selectedItems()):
            styles.color_button_ready(self.ui.btn_pred_analysis_mono_next_4)
            self.ui.btn_pred_analysis_mono_next_4.setEnabled(True)

        prot_1_name = self.ui.lbl_pred_analysis_mono_prot_struct_1.text()
        prot_1_chains = []
        for chain in self.ui.list_pred_analysis_mono_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self.ui.lbl_pred_analysis_mono_prot_struct_2.text()
        prot_2_chains = []
        for chain in self.ui.list_pred_analysis_mono_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self.ui.list_pred_analysis_mono_overview.count()):
            if analysis_name == self.ui.list_pred_analysis_mono_overview.item(tmp_row).text():
                self.ui.btn_pred_analysis_mono_next_4.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_pred_analysis_mono_next_4)
                return

    def check_mono_pred_analysis_if_prot_structs_are_filled(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self.ui.box_pred_analysis_mono_prot_struct_1.itemText(
            self.ui.box_pred_analysis_mono_prot_struct_1.currentIndex(),
        )
        prot_2 = self.ui.box_pred_analysis_mono_prot_struct_2.itemText(
            self.ui.box_pred_analysis_mono_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self.ui.btn_pred_analysis_mono_next_2.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_mono_next_2.setEnabled(False)

    def count_mono_pred_analysis_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains selected in protein 1."""
        self.no_of_selected_chains = len(self.ui.list_pred_analysis_mono_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self.ui.btn_pred_analysis_mono_next_3.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_mono_next_3.setEnabled(False)

    # </editor-fold>

    # </editor-fold>

    def switch_monomer_pred_analysis_tab(self) -> None:
        """Switches the tabs from prediction to analysis and vice versa."""
        if self.ui.tabWidget.currentIndex() == 0:
            # goes from prediction to analysis
            self.ui.tabWidget.setTabEnabled(1, True)
            self.ui.tabWidget.setTabEnabled(0, False)
            self.ui.tabWidget.setCurrentIndex(1)
            gui_elements_to_show = [
                self.ui.lbl_pred_analysis_mono_overview,
                self.ui.list_pred_analysis_mono_overview,
                self.ui.btn_pred_analysis_mono_add,
                self.ui.btn_pred_analysis_mono_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self.ui.btn_pred_analysis_mono_remove,
                self.ui.lbl_pred_analysis_mono_prot_struct_1,
                self.ui.lbl_pred_analysis_mono_prot_struct_2,
                self.ui.lbl_analysis_batch_vs_2,
                self.ui.lbl_pred_analysis_mono_ref_chains,
                self.ui.list_pred_analysis_mono_ref_chains,
                self.ui.btn_pred_analysis_mono_back_4,
                self.ui.btn_pred_analysis_mono_next_3,
                self.ui.box_pred_analysis_mono_prot_struct_1,
                self.ui.box_pred_analysis_mono_prot_struct_2,
                self.ui.btn_pred_analysis_mono_back_3,
                self.ui.btn_pred_analysis_mono_next_2,
                self.ui.lbl_pred_analysis_mono_model_chains,
                self.ui.list_pred_analysis_mono_model_chains,
                self.ui.btn_pred_analysis_mono_back_5,
                self.ui.btn_pred_analysis_mono_next_4,
                self.ui.lbl_pred_analysis_mono_images,
                self.ui.cb_pred_analysis_mono_images,
                self.ui.btn_pred_analysis_mono_start,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            if self.ui.list_pred_analysis_mono_overview.count() > 0:
                self.ui.btn_pred_analysis_mono_remove.show()
                self.ui.btn_pred_analysis_mono_remove.setEnabled(False)
                self.ui.btn_pred_analysis_mono_start.show()
                self.ui.lbl_pred_analysis_mono_images.show()
                self.ui.cb_pred_analysis_mono_images.show()
                styles.color_button_ready(self.ui.btn_pred_analysis_mono_start)
        else:
            # goes from analysis to prediction
            if self.ui.list_pred_analysis_mono_overview.count() > 0:
                gui_elements_to_show = [
                    self.ui.lbl_pred_analysis_mono_prot_to_predict,
                    self.ui.table_pred_analysis_mono_prot_to_predict,
                    self.ui.btn_pred_analysis_mono_go_analysis_setup,
                    self.ui.lbl_pred_analysis_mono_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                    self.ui.btn_pred_analysis_mono_seq_to_predict,
                    self.ui.lbl_pred_analysis_mono_prot_name,
                    self.ui.txt_pred_analysis_mono_prot_name,
                    self.ui.lbl_pred_analysis_mono_prot_name_status,
                    self.ui.btn_pred_analysis_mono_back,
                    self.ui.btn_pred_analysis_mono_next,
                    self.ui.lbl_pred_analysis_mono_seq_name,
                    self.ui.txt_pred_analysis_mono_seq_name,
                    self.ui.lbl_pred_analysis_mono_seq_name_status,
                    self.ui.btn_pred_analysis_mono_back_2,
                    self.ui.btn_pred_analysis_mono_add_protein,
                    self.ui.lbl_pred_mono_advanced_config_2,
                    self.ui.btn_pred_mono_advanced_config_2,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
            else:
                gui_elements_to_show = [
                    self.ui.lbl_pred_analysis_mono_prot_to_predict,
                    self.ui.table_pred_analysis_mono_prot_to_predict,
                    self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                    self.ui.btn_pred_analysis_mono_seq_to_predict,
                    self.ui.lbl_pred_mono_advanced_config_2,
                    self.ui.btn_pred_mono_advanced_config_2,
                    self.ui.btn_pred_analysis_mono_go_analysis_setup,
                    self.ui.lbl_pred_analysis_mono_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self.ui.lbl_pred_analysis_mono_prot_name,
                    self.ui.txt_pred_analysis_mono_prot_name,
                    self.ui.lbl_pred_analysis_mono_prot_name_status,
                    self.ui.btn_pred_analysis_mono_back,
                    self.ui.btn_pred_analysis_mono_next,
                    self.ui.lbl_pred_analysis_mono_seq_name,
                    self.ui.txt_pred_analysis_mono_seq_name,
                    self.ui.lbl_pred_analysis_mono_seq_name_status,
                    self.ui.btn_pred_analysis_mono_back_2,
                    self.ui.btn_pred_analysis_mono_add_protein,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
                self.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)
            self.ui.tabWidget.setTabEnabled(0, True)
            self.ui.tabWidget.setTabEnabled(1, False)
            self.ui.tabWidget.setCurrentIndex(0)

    def start_monomer_prediction_analysis(self) -> None:
        """Sets up the worker for the prediction of the proteins."""
        # # creating tmp directories in scratch folder to organize prediction inputs and outputs
        # # TODO: is there a more elegant way to do it?
        # if not os.path.exists(pathlib.Path(f"{self.scratch_path}/local_predictions")):
        #     os.mkdir(pathlib.Path(f"{self.scratch_path}/local_predictions"))
        # if not os.path.exists(constants.PREDICTION_FASTA_DIR):
        #     os.mkdir(constants.PREDICTION_FASTA_DIR)
        # if not os.path.exists(constants.PREDICTION_PDB_DIR):
        #     os.mkdir(constants.PREDICTION_PDB_DIR)
        # # creating fasta file
        # predictions: list[tuple[str, str]] = gui_utils.get_prediction_name_and_seq_from_table(self.ui.table_pred_mono_prot_to_predict)
        # prot_entries = []
        # last_header = predictions[0][0]
        # pred_list = prediction_list.PredictionList("", [])
        # for tmp_prediction in predictions:
        #     current_header = tmp_prediction[0]
        #     if last_header == current_header:
        #         pred_list.protein_name = tmp_prediction[0]
        #         pred_list.protein_sequence.append(tmp_prediction[1])
        #     else:
        #         prot_entries.append(pred_list)
        #         pred_list = prediction_list.PredictionList("", [])
        #         last_header = current_header
        # for tmp_prot_to_predict in prot_entries:
        #     tmp_prot_to_predict.write_fasta_file()
        # user_name = os.getlogin()
        # fasta_path = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/fasta"
        # pdb_path = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/pdb"
        # # running prediction script
        # if self.prediction_configuration.templates == "none":
        #     try:
        #         subprocess.run([constants.POWERSHELL_EXE, constants.CONVERT_DOS_TO_UNIX])
        #         subprocess.run(["wsl", constants.COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT,
        #                         fasta_path, pdb_path])
        #         subprocess.run(["wsl", "--shutdown"])
        #     except OSError:
        #         shutil.rmtree(pathlib.Path(f"{self.scratch_path}/local_predictions"))
        #         return
        # else:
        #     try:
        #         subprocess.run([constants.POWERSHELL_EXE, constants.CONVERT_DOS_TO_UNIX])
        #         subprocess.run(["wsl", constants.COLABFOLD_PREDICT_SCRIPT,
        #                         fasta_path, pdb_path])
        #         subprocess.run(["wsl", "--shutdown"])
        #     except OSError:
        #         shutil.rmtree(pathlib.Path(f"{self.scratch_path}/local_predictions"))
        #         return
        # # moving best prediction model
        # prediction_results: list[str] = os.listdir(pathlib.Path(constants.PREDICTION_PDB_DIR))
        # for tmp_prediction in predictions:
        #     for filename in prediction_results:
        #         check = filename.find(f"{tmp_prediction[0]}_relaxed_rank_1")
        #         if check != -1:
        #             src = pathlib.Path(f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{filename}")
        #             dest = pathlib.Path(
        #                 f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}")
        #             shutil.copy(src, dest)
        #             os.rename(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}",
        #                       f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{tmp_prediction[0]}.pdb")
        #             tmp_protein = protein.Protein(tmp_prediction[0], pathlib.Path(self.app_project.get_pdb_path()))
        #             self.app_project.add_existing_protein(tmp_protein)
        #             break
        # shutil.rmtree(pathlib.Path(f"{self.scratch_path}/local_predictions"))
        # try:
        #     tmp_first_prediction = predictions[0]
        #     cmd.load(
        #         f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{tmp_first_prediction[0]}.pdb")
        # except pymol.CmdException:
        #     print("Loading the model failed.")
        #     return
        self.prediction_type = constants.PREDICTION_TYPE_PRED_MONO_ANALYSIS
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        # self.worker_prediction_analysis = workers.PredictionWorkerPool(self.ui.table_pred_analysis_mono_prot_to_predict,
        #                                                               self.prediction_configuration, self.app_project)
        constants.PYSSA_LOGGER.info("Thread started for prediction process.")
        # self.threadpool.start(self.worker_prediction_analysis)

        # <editor-fold desc="Worker setup">
        # TODO: test code below
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.ColabfoldWorker(
            self.ui.table_pred_analysis_mono_prot_to_predict, self.prediction_configuration, self.app_project,
        )
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_view_page)
        self.tmp_worker.finished.connect(self.post_prediction_process)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        gui_elements_to_show = [
            self.ui.btn_prediction_abort,
        ]
        gui_elements_to_hide = [
            self.ui.btn_use_page,
            self.ui.btn_close_project,
            self.ui.btn_pred_local_monomer_page,
            self.ui.btn_pred_local_multimer_page,
        ]
        gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        # constants.PYSSA_LOGGER.info("Begin prediction process.")
        # self.prediction_thread = QThread()
        # constants.PYSSA_LOGGER.info("Created a new prediction thread.")
        # self.prediction_worker = workers.PredictionWorker(self.ui.table_pred_mono_prot_to_predict,
        #                                                   self.prediction_configuration, self.app_project)
        # constants.PYSSA_LOGGER.info("Created a new prediction worker.")
        # self._thread_controller.thread_worker_pairs.get(constants.PREDICTION_TASK).setup_and_run_thread(display_self.block_box_prediction_box2)
        # gui_elements_to_show = [
        #     self.ui.btn_prediction_abort,
        # ]
        # gui_elements_to_hide = [
        #     self.ui.btn_use_page,
        #     self.ui.btn_close_project,
        # ]
        # gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        # #self._project_watcher.show_valid_options(self.ui)
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
            return
        else:
            print("Unexpected Error.")
            self.block_box_prediction.close()

    # </editor-fold>

    # <editor-fold desc="Multimer Prediction + Analysis functions">
    def _init_multi_pred_analysis_page(self) -> None:
        """Clears the text boxes and sets the default values for the gui elements."""
        # <editor-fold desc="Prediction section">
        # clears everything
        self.ui.txt_pred_analysis_multi_prot_name.clear()
        self.ui.txt_pred_analysis_multi_prot_seq.clear()
        self.ui.list_pred_analysis_multi_prot_seq_overview.clear()
        for i in range(self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1, -1, -1):
            self.ui.table_pred_analysis_multi_prot_to_predict.removeRow(i)

        # sets up defaults: Prediction
        self.ui.btn_pred_analysis_multi_next.setEnabled(False)
        self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        self.ui.lbl_pred_analysis_multi_prot_name_status.setText("")
        self.ui.lbl_pred_analysis_multi_prot_seq_status.setText("")

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self.ui.list_pred_analysis_multi_overview.clear()
        self.ui.btn_pred_analysis_multi_remove.hide()

        # </editor-fold>

        # self.multi_pred_analysis_show_protein_overview()

    def display_multimer_pred_analysis(self) -> None:
        """Displays the multimer prediction + analysis page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return
        
        self._init_multi_pred_analysis_page()
        self.ui.table_pred_analysis_multi_prot_to_predict.clear()
        self.ui.table_pred_analysis_multi_prot_to_predict.setHorizontalHeaderItem(
            0,
            QtWidgets.QTableWidgetItem("Chain"),
        )
        self.ui.table_pred_analysis_multi_prot_to_predict.setHorizontalHeaderItem(
            1,
            QtWidgets.QTableWidgetItem("Sequence"),
        )
        self.ui.table_pred_analysis_multi_prot_to_predict.resizeColumnsToContents()
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_prot_to_predict,
            self.ui.table_pred_analysis_multi_prot_to_predict,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self.ui.lbl_pred_analysis_multi_prot_name,
            self.ui.txt_pred_analysis_multi_prot_name,
            self.ui.lbl_pred_analysis_multi_prot_name_status,
            self.ui.btn_pred_analysis_multi_back,
            self.ui.btn_pred_analysis_multi_next,
            self.ui.lbl_pred_analysis_multi_prot_seq,
            self.ui.txt_pred_analysis_multi_prot_seq,
            self.ui.lbl_pred_analysis_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add_2,
            self.ui.btn_pred_analysis_multi_prot_seq_add,
            self.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self.ui.list_pred_analysis_multi_prot_seq_overview,
            self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self.ui.btn_pred_analysis_multi_back_2,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_go_analysis_setup,
            self.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self.ui.tabWidget_2.currentIndex() == 1:
            self.ui.tabWidget_2.setCurrentIndex(0)
        self.ui.tabWidget_2.setTabEnabled(1, False)
        self.ui.tabWidget_2.setTabEnabled(0, True)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 22, "Multimer Prediction + Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_pred_analysis_multimer_page,
        )
        self.ui.table_pred_analysis_multi_prot_to_predict.setEnabled(True)

    # <editor-fold desc="Prediction section">
    def multi_pred_analysis_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self.ui.txt_pred_analysis_multi_prot_name.text(), self.ui.table_pred_analysis_multi_prot_to_predict,
        ):
            self.ui.lbl_pred_analysis_multi_prot_name_status.setText("Protein name already used.")
            self.ui.btn_pred_analysis_multi_next.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_next)
        else:
            self.ui.btn_pred_analysis_multi_next.setEnabled(True)
            tools.validate_protein_name(
                self.ui.txt_pred_analysis_multi_prot_name,
                self.ui.lbl_pred_analysis_multi_prot_name_status,
                self.ui.btn_pred_analysis_multi_next,
            )

    def multi_pred_analysis_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self.ui.txt_pred_analysis_multi_prot_seq,
            self.ui.lbl_pred_analysis_multi_prot_seq_status,
            self.ui.btn_pred_analysis_multi_prot_seq_add,
        )

    def multi_pred_analysis_check_if_list_is_empty(self) -> None:
        """Checks if the list of sequences of the protein is empty."""
        if self.ui.list_pred_analysis_multi_prot_seq_overview.count() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_prot_to_predict_add_2)
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_analysis_multi_prot_to_predict_add_2)
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(True)

    def multi_pred_analysis_add_sequence_to_list(self) -> None:
        """Adds the entered sequence to the sequences of the protein."""
        self.ui.list_pred_analysis_multi_prot_seq_overview.addItem(
            QtWidgets.QListWidgetItem(self.ui.txt_pred_analysis_multi_prot_seq.toPlainText()),
        )
        self.multi_pred_analysis_check_if_list_is_empty()

    def multi_pred_analysis_remove_sequence_to_list(self) -> None:
        """Removes the entered sequence from the sequences of the protein."""
        self.ui.list_pred_analysis_multi_prot_seq_overview.takeItem(
            self.ui.list_pred_analysis_multi_prot_seq_overview.currentRow(),
        )
        self.multi_pred_analysis_check_if_list_is_empty()
        self.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(False)

    def multi_pred_analysis_check_if_table_is_empty(self) -> None:
        """Checks if the list of proteins to predict is empty."""
        if self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_multi_predict)
            gui_elements_to_show = [
                self.ui.lbl_pred_analysis_multi_prot_to_predict,
                self.ui.table_pred_analysis_multi_prot_to_predict,
                self.ui.btn_pred_analysis_multi_prot_to_predict_add,
            ]
            gui_elements_to_hide = [
                self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                self.ui.lbl_pred_analysis_multi_prot_name,
                self.ui.txt_pred_analysis_multi_prot_name,
                self.ui.lbl_pred_analysis_multi_prot_name_status,
                self.ui.btn_pred_analysis_multi_back,
                self.ui.btn_pred_analysis_multi_next,
                self.ui.lbl_pred_analysis_multi_prot_seq,
                self.ui.txt_pred_analysis_multi_prot_seq,
                self.ui.lbl_pred_analysis_multi_prot_seq_status,
                self.ui.lbl_pred_multi_prot_seq_add_2,
                self.ui.btn_pred_analysis_multi_prot_seq_add,
                self.ui.lbl_pred_analysis_multi_prot_seq_overview,
                self.ui.list_pred_analysis_multi_prot_seq_overview,
                self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                self.ui.btn_pred_analysis_multi_back_2,
                self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
                self.ui.lbl_pred_analysis_multi_advanced_config,
                self.ui.btn_pred_analysis_multi_advanced_config,
                self.ui.btn_pred_analysis_multi_go_analysis_setup,
                self.ui.lbl_pred_analysis_multi_to_analysis_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.btn_pred_multi_predict.setEnabled(False)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_analysis_multi_go_analysis_setup)
            self.ui.btn_pred_multi_predict.setEnabled(True)
            gui_elements_to_show = [
                self.ui.lbl_pred_analysis_multi_prot_to_predict,
                self.ui.table_pred_analysis_multi_prot_to_predict,
                self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                self.ui.btn_pred_analysis_multi_prot_to_predict_add,
                self.ui.lbl_pred_analysis_multi_advanced_config,
                self.ui.btn_pred_analysis_multi_advanced_config,
                self.ui.btn_pred_analysis_multi_go_analysis_setup,
                self.ui.lbl_pred_analysis_multi_to_analysis_setup,
            ]
            gui_elements_to_hide = [
                self.ui.lbl_pred_analysis_multi_prot_name,
                self.ui.txt_pred_analysis_multi_prot_name,
                self.ui.lbl_pred_analysis_multi_prot_name_status,
                self.ui.btn_pred_analysis_multi_back,
                self.ui.btn_pred_analysis_multi_next,
                self.ui.lbl_pred_analysis_multi_prot_seq,
                self.ui.txt_pred_analysis_multi_prot_seq,
                self.ui.lbl_pred_analysis_multi_prot_seq_status,
                self.ui.lbl_pred_multi_prot_seq_add_2,
                self.ui.btn_pred_analysis_multi_prot_seq_add,
                self.ui.lbl_pred_analysis_multi_prot_seq_overview,
                self.ui.list_pred_analysis_multi_prot_seq_overview,
                self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                self.ui.btn_pred_analysis_multi_back_2,
                self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_add_protein_to_predict(self) -> None:
        """Adds the proteins to the list of proteins to predict."""
        for i in range(self.ui.list_pred_analysis_multi_prot_seq_overview.count()):
            self.ui.table_pred_analysis_multi_prot_to_predict.setRowCount(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1,
            )
            self.ui.table_pred_analysis_multi_prot_to_predict.insertRow(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1,
            )
            tmp_chain_seq = (
                constants.chain_dict.get(i),
                self.ui.list_pred_analysis_multi_prot_seq_overview.item(i).text(),
            )
            self.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                0,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[0]),
            )
            self.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                1,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[1]),
            )
            name_item = QtWidgets.QTableWidgetItem(self.ui.txt_pred_analysis_multi_prot_name.text())
            self.ui.table_pred_analysis_multi_prot_to_predict.setVerticalHeaderItem(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1, name_item,
            )
        self.ui.table_pred_analysis_multi_prot_to_predict.resizeColumnsToContents()
        self.multi_pred_analysis_check_if_table_is_empty()
        self.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_remove_protein_to_predict(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        if self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 1:
            self.ui.table_pred_analysis_multi_prot_to_predict.removeRow(0)
        else:
            self.ui.table_pred_analysis_multi_prot_to_predict.removeRow(
                self.ui.table_pred_analysis_multi_prot_to_predict.currentRow(),
            )
            prot_name = self.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(
                self.ui.table_pred_analysis_multi_prot_to_predict.currentRow(),
            ).text()
            for i in range(self.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
                if self.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text() == prot_name:
                    self.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                        i, 0, QtWidgets.QTableWidgetItem(constants.chain_dict.get(i)),
                    )
        self.multi_pred_analysis_check_if_table_is_empty()
        self.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_add(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_prot_to_predict,
            self.ui.table_pred_analysis_multi_prot_to_predict,
            self.ui.lbl_pred_analysis_multi_prot_name,
            self.ui.txt_pred_analysis_multi_prot_name,
            self.ui.lbl_pred_analysis_multi_prot_name_status,
            self.ui.btn_pred_analysis_multi_back,
            self.ui.btn_pred_analysis_multi_next,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add,
            self.ui.lbl_pred_analysis_multi_prot_seq,
            self.ui.txt_pred_analysis_multi_prot_seq,
            self.ui.lbl_pred_analysis_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add_2,
            self.ui.btn_pred_analysis_multi_prot_seq_add,
            self.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self.ui.list_pred_analysis_multi_prot_seq_overview,
            self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self.ui.btn_pred_analysis_multi_back_2,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_go_analysis_setup,
            self.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_multi_prot_name, self.ui.lbl_pred_analysis_multi_prot_name)
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_multi_prot_seq, self.ui.lbl_pred_analysis_multi_prot_seq)
        self.ui.btn_pred_analysis_multi_next.setEnabled(False)
        self.ui.txt_pred_analysis_multi_prot_name.clear()
        styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_next)
        if self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() > 0:
            try:
                self.ui.table_pred_analysis_multi_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Multimer Prediction in overview table.")

    def multi_pred_analysis_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_prot_to_predict,
            self.ui.table_pred_analysis_multi_prot_to_predict,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self.ui.lbl_pred_analysis_multi_prot_name,
            self.ui.txt_pred_analysis_multi_prot_name,
            self.ui.lbl_pred_analysis_multi_prot_name_status,
            self.ui.btn_pred_analysis_multi_back,
            self.ui.btn_pred_analysis_multi_next,
            self.ui.lbl_pred_analysis_multi_prot_seq,
            self.ui.txt_pred_analysis_multi_prot_seq,
            self.ui.lbl_pred_analysis_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add_2,
            self.ui.btn_pred_analysis_multi_prot_seq_add,
            self.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self.ui.list_pred_analysis_multi_prot_seq_overview,
            self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self.ui.btn_pred_analysis_multi_back_2,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_go_analysis_setup,
            self.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.multi_pred_analysis_check_if_table_is_empty()

    def multi_pred_analysis_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_prot_to_predict,
            self.ui.table_pred_analysis_multi_prot_to_predict,
            self.ui.lbl_pred_analysis_multi_prot_name,
            self.ui.txt_pred_analysis_multi_prot_name,
            self.ui.lbl_pred_analysis_multi_prot_seq,
            self.ui.txt_pred_analysis_multi_prot_seq,
            self.ui.lbl_pred_analysis_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add_2,
            self.ui.btn_pred_analysis_multi_prot_seq_add,
            self.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self.ui.list_pred_analysis_multi_prot_seq_overview,
            self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self.ui.btn_pred_analysis_multi_back_2,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add,
            self.ui.lbl_pred_analysis_multi_prot_name_status,
            self.ui.btn_pred_analysis_multi_back,
            self.ui.btn_pred_analysis_multi_next,
            self.ui.lbl_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_go_analysis_setup,
            self.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_multi_prot_seq, self.ui.lbl_pred_analysis_multi_prot_seq)
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_multi_prot_name, self.ui.lbl_pred_analysis_multi_prot_name)
        self.ui.txt_pred_analysis_multi_prot_seq.clear()
        self.ui.list_pred_analysis_multi_prot_seq_overview.clear()
        self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        self.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_prot_to_predict_add_2)

    def multi_pred_analysis_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_prot_to_predict,
            self.ui.table_pred_analysis_multi_prot_to_predict,
            self.ui.lbl_pred_analysis_multi_prot_name_status,
            self.ui.btn_pred_analysis_multi_back,
            self.ui.btn_pred_analysis_multi_next,
            self.ui.lbl_pred_analysis_multi_prot_name,
            self.ui.txt_pred_analysis_multi_prot_name,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add,
            self.ui.lbl_pred_analysis_multi_prot_seq,
            self.ui.txt_pred_analysis_multi_prot_seq,
            self.ui.lbl_pred_analysis_multi_prot_seq_status,
            self.ui.lbl_pred_multi_prot_seq_add_2,
            self.ui.btn_pred_analysis_multi_prot_seq_add,
            self.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self.ui.list_pred_analysis_multi_prot_seq_overview,
            self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self.ui.btn_pred_analysis_multi_back_2,
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self.ui.lbl_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_advanced_config,
            self.ui.btn_pred_analysis_multi_go_analysis_setup,
            self.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_multi_prot_name, self.ui.lbl_pred_analysis_multi_prot_name)
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_multi_prot_seq, self.ui.lbl_pred_analysis_multi_prot_seq)

    def multi_pred_analysis_prot_seq_overview_item_changed(self) -> None:
        """Enables the remove button of the list of sequences of the protein."""
        self.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(True)

    def multi_pred_analysis_prot_to_predict_item_changed(self) -> None:
        """Enables the remove button of the list of proteins to predict."""
        self.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(True)

    # <editor-fold desc="old">
    # def multi_pred_analysis_show_protein_overview(self):
    #     if self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 0:
    #         self.multimer_prediction_analysis_management.show_stage_x(0)
    #     else:
    #         gui_elements_to_show = [
    #             self.ui.btn_pred_analysis_multi_prot_to_predict_add,
    #             self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
    #         ]
    #         self.multimer_prediction_analysis_management.show_gui_elements_stage_x(
    #             [0, 3], [1, 2], show_specific_elements=gui_elements_to_show,
    #         )
    #
    # def multi_pred_analysis_show_protein_name(self):
    #     gui_elements_to_hide = [
    #         self.ui.btn_pred_multi_prot_to_predict_add,
    #         self.ui.btn_pred_multi_prot_to_predict_remove,
    #     ]
    #     self.multimer_prediction_analysis_management.show_gui_elements_stage_x(
    #         [0, 1], [2, 3], hide_specific_elements=gui_elements_to_hide)
    #     gui_utils.enable_text_box(self.ui.txt_pred_analysis_multi_prot_name,
    #                               self.ui.lbl_pred_analysis_multi_prot_name)
    #
    # def multi_pred_analysis_show_protein_sequence(self):
    #     gui_elements_to_hide = [
    #         self.ui.btn_pred_multi_prot_to_predict_add,
    #         self.ui.btn_pred_multi_prot_to_predict_remove,
    #         self.ui.btn_pred_analysis_multi_back,
    #         self.ui.btn_pred_analysis_multi_next,
    #     ]
    #     self.multimer_prediction_analysis_management.show_gui_elements_stage_x(
    #         [0, 1, 2], [3], hide_specific_elements=gui_elements_to_hide)
    #     gui_utils.disable_text_box(self.ui.txt_pred_analysis_multi_prot_name,
    #                                self.ui.lbl_pred_analysis_multi_prot_name)
    # </editor-fold>

    # </editor-fold>

    def switch_multimer_pred_analysis_tab(self) -> None:
        """Switches the tabs from prediction to analysis and vice versa."""
        if self.ui.tabWidget_2.currentIndex() == 0:
            # goes from prediction to analysis
            self.ui.tabWidget_2.setCurrentIndex(1)
            gui_elements_to_show = [
                self.ui.lbl_pred_analysis_multi_overview,
                self.ui.list_pred_analysis_multi_overview,
                self.ui.btn_pred_analysis_multi_add,
                self.ui.btn_pred_analysis_multi_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self.ui.btn_pred_analysis_multi_remove,
                self.ui.lbl_pred_analysis_multi_prot_struct_1,
                self.ui.lbl_pred_analysis_multi_prot_struct_2,
                self.ui.lbl_analysis_batch_vs_3,
                self.ui.lbl_pred_analysis_multi_ref_chains,
                self.ui.list_pred_analysis_multi_ref_chains,
                self.ui.btn_pred_analysis_multi_back_4,
                self.ui.btn_pred_analysis_multi_next_3,
                self.ui.box_pred_analysis_multi_prot_struct_1,
                self.ui.box_pred_analysis_multi_prot_struct_2,
                self.ui.btn_pred_analysis_multi_back_3,
                self.ui.btn_pred_analysis_multi_next_2,
                self.ui.lbl_pred_analysis_multi_model_chains,
                self.ui.list_pred_analysis_multi_model_chains,
                self.ui.btn_pred_analysis_multi_back_5,
                self.ui.btn_pred_analysis_multi_next_4,
                self.ui.lbl_pred_analysis_multi_images,
                self.ui.cb_pred_analysis_multi_images,
                self.ui.btn_pred_analysis_multi_start,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.tabWidget_2.setTabEnabled(1, True)
            self.ui.tabWidget_2.setTabEnabled(0, False)
            if self.ui.list_pred_analysis_multi_overview.count() > 0:
                self.ui.btn_pred_analysis_multi_remove.show()
                self.ui.btn_pred_analysis_multi_remove.setEnabled(False)
                self.ui.btn_pred_analysis_multi_start.show()
                self.ui.lbl_pred_analysis_multi_images.show()
                self.ui.cb_pred_analysis_multi_images.show()
                styles.color_button_ready(self.ui.btn_pred_analysis_multi_start)
        else:
            # goes from analysis to prediction
            self.ui.tabWidget_2.setCurrentIndex(0)
            if self.ui.list_pred_analysis_multi_overview.count() > 0:
                gui_elements_to_show = [
                    self.ui.lbl_pred_analysis_multi_prot_to_predict,
                    self.ui.table_pred_analysis_multi_prot_to_predict,
                    self.ui.btn_pred_analysis_multi_go_analysis_setup,
                    self.ui.lbl_pred_analysis_multi_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                    self.ui.btn_pred_analysis_multi_prot_to_predict_add,
                    self.ui.lbl_pred_analysis_multi_advanced_config,
                    self.ui.btn_pred_analysis_multi_advanced_config,
                    self.ui.lbl_pred_analysis_multi_prot_name,
                    self.ui.txt_pred_analysis_multi_prot_name,
                    self.ui.lbl_pred_analysis_multi_prot_name_status,
                    self.ui.btn_pred_analysis_multi_back,
                    self.ui.btn_pred_analysis_multi_next,
                    self.ui.lbl_pred_analysis_multi_prot_seq,
                    self.ui.txt_pred_analysis_multi_prot_seq,
                    self.ui.lbl_pred_analysis_multi_prot_seq_status,
                    self.ui.lbl_pred_multi_prot_seq_add_2,
                    self.ui.btn_pred_analysis_multi_prot_seq_add,
                    self.ui.lbl_pred_analysis_multi_prot_seq_overview,
                    self.ui.list_pred_analysis_multi_prot_seq_overview,
                    self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                    self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                    self.ui.btn_pred_analysis_multi_back_2,
                    self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
            else:
                gui_elements_to_show = [
                    self.ui.lbl_pred_analysis_multi_prot_to_predict,
                    self.ui.table_pred_analysis_multi_prot_to_predict,
                    self.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                    self.ui.btn_pred_analysis_multi_prot_to_predict_add,
                    self.ui.lbl_pred_analysis_multi_advanced_config,
                    self.ui.btn_pred_analysis_multi_advanced_config,
                    self.ui.btn_pred_analysis_multi_go_analysis_setup,
                    self.ui.lbl_pred_analysis_multi_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self.ui.lbl_pred_analysis_multi_prot_name,
                    self.ui.txt_pred_analysis_multi_prot_name,
                    self.ui.lbl_pred_analysis_multi_prot_name_status,
                    self.ui.btn_pred_analysis_multi_back,
                    self.ui.btn_pred_analysis_multi_next,
                    self.ui.lbl_pred_analysis_multi_prot_seq,
                    self.ui.txt_pred_analysis_multi_prot_seq,
                    self.ui.lbl_pred_analysis_multi_prot_seq_status,
                    self.ui.lbl_pred_multi_prot_seq_add_2,
                    self.ui.btn_pred_analysis_multi_prot_seq_add,
                    self.ui.lbl_pred_analysis_multi_prot_seq_overview,
                    self.ui.list_pred_analysis_multi_prot_seq_overview,
                    self.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                    self.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                    self.ui.btn_pred_analysis_multi_back_2,
                    self.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.tabWidget_2.setTabEnabled(0, True)
            self.ui.tabWidget_2.setTabEnabled(1, False)

    # <editor-fold desc="Analysis section">
    def multi_pred_analysis_structure_analysis_add(self) -> None:
        """Shows the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_pred_analysis_multi_prot_struct_1.clear()
        self.ui.lbl_pred_analysis_multi_prot_struct_2.clear()
        self.ui.lbl_pred_analysis_multi_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_pred_analysis_multi_prot_struct_2.setText("Protein structure 2")
        self.fill_multi_pred_analysis_protein_boxes()
        if self.ui.list_pred_analysis_multi_overview.count() > 0:
            try:
                self.ui.list_pred_analysis_multi_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def multi_pred_analysis_structure_analysis_next_2(self) -> None:
        """Shows the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_pred_analysis_multi_prot_struct_1.setText(
            self.ui.box_pred_analysis_multi_prot_struct_1.currentText(),
        )
        self.ui.lbl_pred_analysis_multi_prot_struct_2.setText(
            self.ui.box_pred_analysis_multi_prot_struct_2.currentText(),
        )
        self.ui.list_pred_analysis_multi_ref_chains.clear()
        self.ui.btn_pred_analysis_multi_next_3.setEnabled(False)
        self.ui.list_pred_analysis_multi_ref_chains.setEnabled(True)

        for i in range(self.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            if (
                self.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text()
                == self.ui.box_pred_analysis_multi_prot_struct_1.currentText()
            ):
                self.ui.list_pred_analysis_multi_ref_chains.addItem(
                    self.ui.table_pred_analysis_multi_prot_to_predict.item(i, 0).text(),
                )
        if self.ui.list_pred_analysis_multi_ref_chains.count() == 0:
            tmp_protein = self.app_project.search_protein(self.ui.box_pred_analysis_multi_prot_struct_1.currentText())
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self.ui.list_pred_analysis_multi_ref_chains.addItem(tmp_chain.chain_letter)
        if self.ui.list_pred_analysis_multi_ref_chains.count() == 1:
            self.ui.lbl_pred_analysis_multi_ref_chains.setText(
                f"Select chain in protein structure {self.ui.lbl_pred_analysis_multi_prot_struct_1.text()}.",
            )
        else:
            self.ui.lbl_pred_analysis_multi_ref_chains.setText(
                f"Select chains in protein structure {self.ui.lbl_pred_analysis_multi_prot_struct_1.text()}.",
            )

    def multi_pred_analysis_structure_analysis_back_3(self) -> None:
        """Hides the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self.ui.list_pred_analysis_multi_overview.count() > 0:
            self.ui.btn_pred_analysis_multi_remove.show()
            self.ui.btn_pred_analysis_multi_remove.setEnabled(False)
            self.ui.btn_pred_analysis_multi_start.show()
            self.ui.lbl_pred_analysis_multi_images.show()
            self.ui.cb_pred_analysis_multi_images.show()
            styles.color_button_ready(self.ui.btn_pred_analysis_multi_start)

    def multi_pred_analysis_structure_analysis_next_3(self) -> None:
        """Shows the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.list_pred_analysis_multi_model_chains.clear()
        self.ui.list_pred_analysis_multi_ref_chains.setEnabled(False)
        self.ui.btn_pred_analysis_multi_next_4.setEnabled(False)

        for i in range(self.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            if (
                self.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text()
                == self.ui.box_pred_analysis_multi_prot_struct_2.currentText()
            ):
                self.ui.list_pred_analysis_multi_model_chains.addItem(
                    self.ui.table_pred_analysis_multi_prot_to_predict.item(i, 0).text(),
                )
        if self.ui.list_pred_analysis_multi_model_chains.count() == 0:
            tmp_protein = self.app_project.search_protein(self.ui.box_pred_analysis_multi_prot_struct_2.currentText())
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self.ui.list_pred_analysis_multi_model_chains.addItem(tmp_chain.chain_letter)
        if len(self.ui.list_pred_analysis_multi_ref_chains.selectedItems()) == 1:
            self.ui.lbl_pred_analysis_multi_model_chains.setText(
                f"Select 1 chain in protein structure {self.ui.lbl_pred_analysis_multi_prot_struct_2.text()}.",
            )
        else:
            self.ui.lbl_pred_analysis_multi_model_chains.setText(
                f"Select {len(self.ui.list_pred_analysis_multi_ref_chains.selectedItems())} chains in "
                f"protein structure {self.ui.lbl_pred_analysis_multi_prot_struct_2.text()}.",
            )

    def multi_pred_analysis_structure_analysis_back_4(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_pred_analysis_multi_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_pred_analysis_multi_prot_struct_2.setText("Protein structure 2")

    def multi_pred_analysis_structure_analysis_next_4(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self.ui.lbl_pred_analysis_multi_prot_struct_1.text()
        prot_1_chains = []
        for chain in self.ui.list_pred_analysis_multi_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self.ui.lbl_pred_analysis_multi_prot_struct_2.text()
        prot_2_chains = []
        for chain in self.ui.list_pred_analysis_multi_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self.ui.list_pred_analysis_multi_overview.addItem(item)
        self.ui.btn_pred_analysis_multi_remove.setEnabled(False)
        styles.color_button_ready(self.ui.btn_pred_analysis_multi_start)

    def multi_pred_analysis_structure_analysis_back_5(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.list_pred_analysis_multi_ref_chains.setEnabled(True)

        # tmp_protein = self.app_project.search_protein(self.ui.box_pred_analysis_multi_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self.ui.list_pred_analysis_multi_ref_chains.addItem(tmp_chain.chain_letter)

    def multi_pred_analysis_structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self.ui.btn_pred_analysis_multi_remove.setEnabled(True)

    def fill_multi_pred_analysis_protein_boxes(self) -> None:
        """Fills the combo boxes with the protein names."""
        protein_names = []
        for i in range(self.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            protein_names.append(self.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text())
        for tmp_protein in self.app_project.proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        protein_names = list(set(protein_names))
        self.ui.box_pred_analysis_multi_prot_struct_1.clear()
        self.ui.box_pred_analysis_multi_prot_struct_2.clear()
        gui_utils.fill_combo_box(self.ui.box_pred_analysis_multi_prot_struct_1, protein_names)
        gui_utils.fill_combo_box(self.ui.box_pred_analysis_multi_prot_struct_2, protein_names)

    def remove_multi_pred_analysis_analysis_run(self) -> None:
        """Removes the selected protein pair from the list of protein pairs to analyze."""
        self.ui.list_pred_analysis_multi_overview.takeItem(self.ui.list_pred_analysis_multi_overview.currentRow())
        gui_elements_to_show = [
            self.ui.lbl_pred_analysis_multi_overview,
            self.ui.list_pred_analysis_multi_overview,
            self.ui.btn_pred_analysis_multi_add,
            self.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_multi_remove,
            self.ui.lbl_pred_analysis_multi_prot_struct_1,
            self.ui.lbl_pred_analysis_multi_prot_struct_2,
            self.ui.lbl_analysis_batch_vs_3,
            self.ui.lbl_pred_analysis_multi_ref_chains,
            self.ui.list_pred_analysis_multi_ref_chains,
            self.ui.btn_pred_analysis_multi_back_4,
            self.ui.btn_pred_analysis_multi_next_3,
            self.ui.box_pred_analysis_multi_prot_struct_1,
            self.ui.box_pred_analysis_multi_prot_struct_2,
            self.ui.btn_pred_analysis_multi_back_3,
            self.ui.btn_pred_analysis_multi_next_2,
            self.ui.lbl_pred_analysis_multi_model_chains,
            self.ui.list_pred_analysis_multi_model_chains,
            self.ui.btn_pred_analysis_multi_back_5,
            self.ui.btn_pred_analysis_multi_next_4,
            self.ui.lbl_pred_analysis_multi_images,
            self.ui.cb_pred_analysis_multi_images,
            self.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self.ui.list_pred_analysis_multi_overview.count() > 0:
            self.ui.btn_pred_analysis_multi_remove.show()
            self.ui.btn_pred_analysis_multi_remove.setEnabled(False)
            self.ui.btn_pred_analysis_multi_start.show()
            self.ui.lbl_pred_analysis_multi_images.show()
            self.ui.cb_pred_analysis_multi_images.show()
            styles.color_button_ready(self.ui.btn_pred_analysis_multi_start)
        # if self.ui.list_pred_analysis_multi_overview.count() == 0:
        #
        #     self.ui.btn_pred_analysis_multi_back_pred_setup.show()
        #     self.ui.btn_pred_analysis_multi_remove.hide()

    def check_multi_pred_analysis_if_same_no_of_chains_selected(self) -> None:
        """Checks if the same number of chains were selected."""
        self.ui.btn_pred_analysis_multi_next_4.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_next_4)
        if self.no_of_selected_chains == len(self.ui.list_pred_analysis_multi_model_chains.selectedItems()):
            styles.color_button_ready(self.ui.btn_pred_analysis_multi_next_4)
            self.ui.btn_pred_analysis_multi_next_4.setEnabled(True)

        prot_1_name = self.ui.lbl_pred_analysis_multi_prot_struct_1.text()
        prot_1_chains = []
        for chain in self.ui.list_pred_analysis_multi_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self.ui.lbl_pred_analysis_multi_prot_struct_2.text()
        prot_2_chains = []
        for chain in self.ui.list_pred_analysis_multi_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self.ui.list_pred_analysis_multi_overview.count()):
            if analysis_name == self.ui.list_pred_analysis_multi_overview.item(tmp_row).text():
                self.ui.btn_pred_analysis_multi_next_4.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_next_4)
                return

    def check_multi_pred_analysis_if_prot_structs_are_filled(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self.ui.box_pred_analysis_multi_prot_struct_1.itemText(
            self.ui.box_pred_analysis_multi_prot_struct_1.currentIndex(),
        )
        prot_2 = self.ui.box_pred_analysis_multi_prot_struct_2.itemText(
            self.ui.box_pred_analysis_multi_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self.ui.btn_pred_analysis_multi_next_2.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_multi_next_2.setEnabled(False)

    def count_multi_pred_analysis_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains in protein 1."""
        self.no_of_selected_chains = len(self.ui.list_pred_analysis_multi_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self.ui.btn_pred_analysis_multi_next_3.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_multi_next_3.setEnabled(False)

    def start_multimer_prediction_analysis(self) -> None:
        """Sets up the prediction process."""
        self.prediction_type = constants.PREDICTION_TYPE_PRED_MULTI_ANALYSIS
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        # self.worker_prediction = workers.PredictionWorkerPool(self.ui.table_pred_analysis_multi_prot_to_predict,
        #                                                      self.prediction_configuration, self.app_project)
        constants.PYSSA_LOGGER.info("Thread started for prediction process.")
        # self.threadpool.start(self.worker_prediction)
        gui_elements_to_show = [
            self.ui.btn_prediction_abort,
        ]
        gui_elements_to_hide = [
            self.ui.btn_use_page,
            self.ui.btn_close_project,
            self.ui.btn_pred_local_monomer_page,
            self.ui.btn_pred_local_multimer_page,
        ]

        # <editor-fold desc="Worker setup">
        # TODO: test code below
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.ColabfoldWorker(
            self.ui.table_pred_analysis_multi_prot_to_predict, self.prediction_configuration, self.app_project,
        )
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_view_page)
        self.tmp_worker.finished.connect(self.post_prediction_process)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

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
            return
        else:
            print("Unexpected Error.")
            self.block_box_prediction.close()

    # </editor-fold>

    # <editor-fold desc="Structure Analysis functions">
    def structure_analysis_add(self) -> None:
        """Shows the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self.ui.lbl_analysis_batch_overview,
            self.ui.list_analysis_batch_overview,
            self.ui.lbl_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.lbl_analysis_batch_vs,
            self.ui.lbl_analysis_batch_prot_struct_2,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next,
        ]
        gui_elements_to_hide = [
            self.ui.btn_analysis_batch_remove,
            self.ui.btn_analysis_batch_add,
            self.ui.lbl_analysis_batch_ref_chains,
            self.ui.list_analysis_batch_ref_chains,
            self.ui.btn_analysis_batch_back_2,
            self.ui.btn_analysis_batch_next_2,
            self.ui.lbl_analysis_batch_model_chains,
            self.ui.list_analysis_batch_model_chains,
            self.ui.btn_analysis_batch_back_3,
            self.ui.btn_analysis_batch_next_3,
            self.ui.lbl_analysis_batch_images,
            self.ui.cb_analysis_batch_images,
            self.ui.btn_analysis_batch_start,
        ]

        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_analysis_batch_prot_struct_1.clear()
        self.ui.lbl_analysis_batch_prot_struct_2.clear()
        self.ui.lbl_analysis_batch_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_analysis_batch_prot_struct_2.setText("Protein structure 2")
        self.fill_protein_boxes_batch()
        if self.ui.list_analysis_batch_overview.count() > 0:
            try:
                self.ui.list_analysis_batch_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def structure_analysis_next(self) -> None:
        """Shows the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self.ui.lbl_analysis_batch_overview,
            self.ui.list_analysis_batch_overview,
            self.ui.lbl_analysis_batch_prot_struct_1,
            self.ui.lbl_analysis_batch_prot_struct_2,
            self.ui.lbl_analysis_batch_vs,
            self.ui.lbl_analysis_batch_ref_chains,
            self.ui.list_analysis_batch_ref_chains,
            self.ui.btn_analysis_batch_back_2,
            self.ui.btn_analysis_batch_next_2,
        ]
        gui_elements_to_hide = [
            self.ui.btn_analysis_batch_remove,
            self.ui.btn_analysis_batch_add,
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next,
            self.ui.lbl_analysis_batch_model_chains,
            self.ui.list_analysis_batch_model_chains,
            self.ui.btn_analysis_batch_back_3,
            self.ui.btn_analysis_batch_next_3,
            self.ui.lbl_analysis_batch_images,
            self.ui.cb_analysis_batch_images,
            self.ui.btn_analysis_batch_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_analysis_batch_prot_struct_1.setText(self.ui.box_analysis_batch_prot_struct_1.currentText())
        self.ui.lbl_analysis_batch_prot_struct_2.setText(self.ui.box_analysis_batch_prot_struct_2.currentText())
        self.ui.list_analysis_batch_ref_chains.clear()
        self.ui.btn_analysis_batch_next_2.setEnabled(False)
        self.ui.list_analysis_batch_ref_chains.setEnabled(True)

        tmp_protein = self.app_project.search_protein(self.ui.box_analysis_batch_prot_struct_1.currentText())
        for tmp_chain in tmp_protein.chains:
            if tmp_chain.chain_type == "protein_chain":
                self.ui.list_analysis_batch_ref_chains.addItem(tmp_chain.chain_letter)
        if self.ui.list_analysis_batch_ref_chains.count() == 1:
            self.ui.lbl_analysis_batch_ref_chains.setText(
                f"Select chain in protein structure {self.ui.lbl_analysis_batch_prot_struct_1.text()}.",
            )
        else:
            self.ui.lbl_analysis_batch_ref_chains.setText(
                f"Select chains in protein structure {self.ui.lbl_analysis_batch_prot_struct_1.text()}.",
            )

    def structure_analysis_back(self) -> None:
        """Hides the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self.ui.lbl_analysis_batch_overview,
            self.ui.list_analysis_batch_overview,
            self.ui.btn_analysis_batch_add,
        ]
        gui_elements_to_hide = [
            self.ui.btn_analysis_batch_remove,
            self.ui.lbl_analysis_batch_prot_struct_1,
            self.ui.lbl_analysis_batch_prot_struct_2,
            self.ui.lbl_analysis_batch_vs,
            self.ui.lbl_analysis_batch_ref_chains,
            self.ui.list_analysis_batch_ref_chains,
            self.ui.btn_analysis_batch_back_2,
            self.ui.btn_analysis_batch_next_2,
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next,
            self.ui.lbl_analysis_batch_model_chains,
            self.ui.list_analysis_batch_model_chains,
            self.ui.btn_analysis_batch_back_3,
            self.ui.btn_analysis_batch_next_3,
            self.ui.lbl_analysis_batch_images,
            self.ui.cb_analysis_batch_images,
            self.ui.btn_analysis_batch_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self.ui.list_analysis_batch_overview.count() > 0:
            self.ui.btn_analysis_batch_remove.show()
            self.ui.btn_analysis_batch_remove.setEnabled(False)
            self.ui.btn_analysis_batch_start.show()
            self.ui.lbl_analysis_batch_images.show()
            self.ui.cb_analysis_batch_images.show()
            styles.color_button_ready(self.ui.btn_analysis_batch_start)

    def structure_analysis_next_2(self) -> None:
        """Shows the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self.ui.lbl_analysis_batch_overview,
            self.ui.list_analysis_batch_overview,
            self.ui.lbl_analysis_batch_prot_struct_1,
            self.ui.lbl_analysis_batch_prot_struct_2,
            self.ui.lbl_analysis_batch_vs,
            self.ui.lbl_analysis_batch_ref_chains,
            self.ui.list_analysis_batch_ref_chains,
            self.ui.lbl_analysis_batch_model_chains,
            self.ui.list_analysis_batch_model_chains,
            self.ui.btn_analysis_batch_back_3,
            self.ui.btn_analysis_batch_next_3,
        ]
        gui_elements_to_hide = [
            self.ui.btn_analysis_batch_remove,
            self.ui.btn_analysis_batch_add,
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next,
            self.ui.btn_analysis_batch_back_2,
            self.ui.btn_analysis_batch_next_2,
            self.ui.lbl_analysis_batch_images,
            self.ui.cb_analysis_batch_images,
            self.ui.btn_analysis_batch_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.list_analysis_batch_model_chains.clear()
        self.ui.list_analysis_batch_ref_chains.setEnabled(False)
        self.ui.btn_analysis_batch_next_3.setEnabled(False)

        tmp_protein = self.app_project.search_protein(self.ui.box_analysis_batch_prot_struct_2.currentText())
        for tmp_chain in tmp_protein.chains:
            if tmp_chain.chain_type == "protein_chain":
                self.ui.list_analysis_batch_model_chains.addItem(tmp_chain.chain_letter)
        if len(self.ui.list_analysis_batch_ref_chains.selectedItems()) == 1:
            self.ui.lbl_analysis_batch_model_chains.setText(
                f"Select 1 chain in protein structure {self.ui.lbl_analysis_batch_prot_struct_2.text()}.",
            )
        else:
            self.ui.lbl_analysis_batch_model_chains.setText(
                f"Select {len(self.ui.list_analysis_batch_ref_chains.selectedItems())} chains in "
                f"protein structure {self.ui.lbl_analysis_batch_prot_struct_2.text()}.",
            )

    def structure_analysis_back_2(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self.ui.lbl_analysis_batch_overview,
            self.ui.list_analysis_batch_overview,
            self.ui.lbl_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.lbl_analysis_batch_vs,
            self.ui.lbl_analysis_batch_prot_struct_2,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next,
        ]
        gui_elements_to_hide = [
            self.ui.btn_analysis_batch_remove,
            self.ui.btn_analysis_batch_add,
            self.ui.lbl_analysis_batch_ref_chains,
            self.ui.list_analysis_batch_ref_chains,
            self.ui.btn_analysis_batch_back_2,
            self.ui.btn_analysis_batch_next_2,
            self.ui.lbl_analysis_batch_model_chains,
            self.ui.list_analysis_batch_model_chains,
            self.ui.btn_analysis_batch_back_3,
            self.ui.btn_analysis_batch_next_3,
            self.ui.lbl_analysis_batch_images,
            self.ui.cb_analysis_batch_images,
            self.ui.btn_analysis_batch_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.lbl_analysis_batch_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_analysis_batch_prot_struct_2.setText("Protein structure 2")

    def structure_analysis_next_3(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self.ui.btn_analysis_batch_remove,
            self.ui.btn_analysis_batch_add,
            self.ui.lbl_analysis_batch_overview,
            self.ui.list_analysis_batch_overview,
            self.ui.lbl_analysis_batch_images,
            self.ui.cb_analysis_batch_images,
            self.ui.btn_analysis_batch_start,
        ]
        gui_elements_to_hide = [
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.lbl_analysis_batch_prot_struct_1,
            self.ui.lbl_analysis_batch_prot_struct_2,
            self.ui.lbl_analysis_batch_vs,
            self.ui.lbl_analysis_batch_ref_chains,
            self.ui.list_analysis_batch_ref_chains,
            self.ui.lbl_analysis_batch_model_chains,
            self.ui.list_analysis_batch_model_chains,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next,
            self.ui.btn_analysis_batch_back_2,
            self.ui.btn_analysis_batch_next_2,
            self.ui.btn_analysis_batch_back_3,
            self.ui.btn_analysis_batch_next_3,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self.ui.lbl_analysis_batch_prot_struct_1.text()
        prot_1_chains = []
        for chain in self.ui.list_analysis_batch_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self.ui.lbl_analysis_batch_prot_struct_2.text()
        prot_2_chains = []
        for chain in self.ui.list_analysis_batch_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self.ui.list_analysis_batch_overview.addItem(item)
        self.ui.btn_analysis_batch_remove.setEnabled(False)
        styles.color_button_ready(self.ui.btn_analysis_batch_start)

    def structure_analysis_back_3(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self.ui.lbl_analysis_batch_overview,
            self.ui.list_analysis_batch_overview,
            self.ui.lbl_analysis_batch_prot_struct_1,
            self.ui.lbl_analysis_batch_prot_struct_2,
            self.ui.lbl_analysis_batch_vs,
            self.ui.lbl_analysis_batch_ref_chains,
            self.ui.list_analysis_batch_ref_chains,
            self.ui.btn_analysis_batch_back_2,
            self.ui.btn_analysis_batch_next_2,
        ]
        gui_elements_to_hide = [
            self.ui.btn_analysis_batch_remove,
            self.ui.btn_analysis_batch_add,
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next,
            self.ui.btn_analysis_batch_back_3,
            self.ui.btn_analysis_batch_next_3,
            self.ui.lbl_analysis_batch_images,
            self.ui.cb_analysis_batch_images,
            self.ui.btn_analysis_batch_start,
            self.ui.lbl_analysis_batch_model_chains,
            self.ui.list_analysis_batch_model_chains,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.ui.list_analysis_batch_ref_chains.setEnabled(True)

        # tmp_protein = self.app_project.search_protein(self.ui.box_analysis_batch_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self.ui.list_analysis_batch_ref_chains.addItem(tmp_chain.chain_letter)

    def structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self.ui.btn_analysis_batch_remove.setEnabled(True)

    # <editor-fold desc="Old">
    # def show_batch_analysis_stage_0(self):
    #     gui_page_management.show_analysis_page_stage_0(self.batch_analysis_management,
    #                                                    self.ui.list_analysis_batch_ref_chains,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_1,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_2,
    #                                                    self.ui.list_analysis_batch_model_chains,
    #                                                    self.ui.list_analysis_batch_overview,
    #                                                    self.ui.btn_analysis_batch_remove,
    #                                                    self.ui.btn_analysis_batch_add,
    #                                                    0)
    #
    # def show_batch_analysis_stage_1(self):
    #     gui_page_management.show_analysis_page_stage_1(self.batch_analysis_management,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_1,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_2,
    #                                                    self.ui.btn_analysis_batch_remove,
    #                                                    self.ui.btn_analysis_batch_add,
    #                                                    0,
    #                                                    self.fill_protein_boxes_batch)
    #
    # def show_batch_analysis_stage_2(self):
    #     gui_page_management.show_analysis_page_stage_2(self.app_project,
    #                                                    self.batch_analysis_management,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_1,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_2,
    #                                                    self.ui.box_analysis_batch_prot_struct_1,
    #                                                    self.ui.box_analysis_batch_prot_struct_2,
    #                                                    self.ui.lbl_analysis_batch_ref_chains,
    #                                                    self.ui.list_analysis_batch_ref_chains,
    #                                                    self.ui.btn_analysis_batch_next,
    #                                                    self.ui.btn_analysis_batch_next_2,
    #                                                    self.ui.btn_analysis_batch_back,
    #                                                    0)
    #
    # def show_batch_analysis_stage_3(self):
    #     gui_page_management.show_analysis_page_stage_3(self.app_project,
    #                                                    self.batch_analysis_management,
    #                                                    self.ui.list_analysis_batch_overview,
    #                                                    self.ui.btn_analysis_batch_add,
    #                                                    self.ui.btn_analysis_batch_remove,
    #                                                    self.ui.btn_analysis_batch_next,
    #                                                    self.ui.btn_analysis_batch_back,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_1,
    #                                                    self.ui.lbl_analysis_batch_prot_struct_2,
    #                                                    self.ui.box_analysis_batch_prot_struct_1,
    #                                                    self.ui.box_analysis_batch_prot_struct_2,
    #                                                    self.ui.btn_analysis_batch_next_2,
    #                                                    self.ui.btn_analysis_batch_back_2,
    #                                                    self.ui.lbl_analysis_batch_model_chains,
    #                                                    self.ui.list_analysis_batch_model_chains,
    #                                                    self.ui.btn_analysis_batch_next_3,
    #                                                    self.no_of_selected_chains,
    #                                                    0)
    # </editor-fold>

    def fill_protein_boxes_batch(self) -> None:
        """Fills the combo boxes with the protein names."""
        proteins = []
        for tmp_protein in self.app_project.proteins:
            proteins.append(tmp_protein.get_molecule_object())
        proteins.insert(0, "")
        self.ui.box_analysis_batch_prot_struct_1.clear()
        self.ui.box_analysis_batch_prot_struct_2.clear()
        gui_utils.fill_combo_box(self.ui.box_analysis_batch_prot_struct_1, proteins)
        gui_utils.fill_combo_box(self.ui.box_analysis_batch_prot_struct_2, proteins)

    def remove_analysis_run(self) -> None:
        """Removes the selected protein pair from the list of protein pairs to analyze."""
        self.ui.list_analysis_batch_overview.takeItem(self.ui.list_analysis_batch_overview.currentRow())
        if self.ui.list_analysis_batch_overview.count() == 0:
            gui_elements_to_show = [
                self.ui.lbl_analysis_batch_overview,
                self.ui.list_analysis_batch_overview,
                self.ui.btn_analysis_batch_add,
            ]

            gui_elements_to_hide = [
                self.ui.btn_analysis_batch_remove,
                self.ui.lbl_analysis_batch_prot_struct_1,
                self.ui.lbl_analysis_batch_prot_struct_2,
                self.ui.lbl_analysis_batch_vs,
                self.ui.lbl_analysis_batch_ref_chains,
                self.ui.list_analysis_batch_ref_chains,
                self.ui.btn_analysis_batch_back_2,
                self.ui.btn_analysis_batch_next_2,
                self.ui.box_analysis_batch_prot_struct_1,
                self.ui.box_analysis_batch_prot_struct_2,
                self.ui.btn_analysis_batch_back,
                self.ui.btn_analysis_batch_next,
                self.ui.lbl_analysis_batch_model_chains,
                self.ui.list_analysis_batch_model_chains,
                self.ui.btn_analysis_batch_back_3,
                self.ui.btn_analysis_batch_next_3,
                self.ui.lbl_analysis_batch_images,
                self.ui.cb_analysis_batch_images,
                self.ui.btn_analysis_batch_start,
            ]

            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.btn_analysis_batch_remove.hide()
        else:
            if self.ui.list_analysis_batch_overview.count() > 0:
                try:
                    self.ui.list_analysis_batch_overview.currentItem().setSelected(False)
                except AttributeError:
                    constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")
        self.ui.btn_analysis_batch_remove.setEnabled(False)

    def check_if_same_no_of_chains_selected_batch(self) -> None:
        """Checks if the same number of proteins were selected."""
        self.ui.btn_analysis_batch_next_3.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_analysis_batch_next_3)

        if self.no_of_selected_chains == len(self.ui.list_analysis_batch_model_chains.selectedItems()):
            styles.color_button_ready(self.ui.btn_analysis_batch_next_3)
            self.ui.btn_analysis_batch_next_3.setEnabled(True)

        prot_1_name = self.ui.lbl_analysis_batch_prot_struct_1.text()
        prot_1_chains = []
        for chain in self.ui.list_analysis_batch_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self.ui.lbl_analysis_batch_prot_struct_2.text()
        prot_2_chains = []
        for chain in self.ui.list_analysis_batch_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self.ui.list_analysis_batch_overview.count()):
            if analysis_name == self.ui.list_analysis_batch_overview.item(tmp_row).text():
                self.ui.btn_analysis_batch_next_3.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_analysis_batch_next_3)
                return

    def check_if_prot_structs_are_filled_batch(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self.ui.box_analysis_batch_prot_struct_1.itemText(
            self.ui.box_analysis_batch_prot_struct_1.currentIndex(),
        )
        prot_2 = self.ui.box_analysis_batch_prot_struct_2.itemText(
            self.ui.box_analysis_batch_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self.ui.btn_analysis_batch_next.setEnabled(True)
        else:
            self.ui.btn_analysis_batch_next.setEnabled(False)

    def count_batch_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains of protein 1."""
        self.no_of_selected_chains = len(self.ui.list_analysis_batch_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self.ui.btn_analysis_batch_next_2.setEnabled(True)
        else:
            self.ui.btn_analysis_batch_next_2.setEnabled(False)

    def post_analysis_process(self) -> None:
        """Post process after the analysis thread finished."""
        constants.PYSSA_LOGGER.debug("post_analysis_process() started ...")
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
        self.block_box_analysis.destroy(True)
        basic_boxes.ok(
            "Structure analysis",
            "All structure analysis' are done. Go to results to check the new results.",
            QtWidgets.QMessageBox.Information,
        )
        constants.PYSSA_LOGGER.info("All structure analysis' are done.")
        self._project_watcher.show_valid_options(self.ui)
        self.display_view_page()
        self._init_batch_analysis_page()

    # def thread_func_run_analysis(self, callback):
    #     worker = workers.AnalysisWorkerPool(
    #         self.ui.list_analysis_batch_overview, self.ui.cb_analysis_images,
    #         self.status_bar, self.app_project, self.app_settings, self._init_batch_analysis_page)
    #     worker.run()
    #     callback()
    #
    # def analysis_callback(self):
    #     self.post_analysis_process()

    def start_process_batch(self) -> None:
        """Sets up the worker for the analysis task."""
        # batch_analysis = []
        # for row_no in range(self.ui.list_analysis_batch_overview.count()):
        #     tmp_batch_analysis = self.ui.list_analysis_batch_overview.item(row_no).text()
        #     separator_index = tmp_batch_analysis.find("_vs_")
        #     prot_1 = tmp_batch_analysis[:separator_index]
        #     if prot_1.find(";") != -1:
        #         prot_1_name = prot_1[:prot_1.find(";")]
        #         prot_1_chains = prot_1[prot_1.find(";") + 1:].split(",")
        #     else:
        #         prot_1_name = prot_1
        #         prot_1_chains = None
        #     prot_2 = tmp_batch_analysis[separator_index + 4:]
        #     if prot_2.find(";") != -1:
        #         prot_2_name = prot_2[:prot_2.find(";")]
        #         prot_2_chains = prot_2[prot_2.find(";") + 1:].split(",")
        #     else:
        #         prot_2_name = prot_2
        #         prot_2_chains = None
        #
        #     tmp_prot_1 = protein_analysis_info.ProteinAnalysisInfo(prot_1_name, prot_1_chains, tmp_batch_analysis)
        #     tmp_prot_2 = protein_analysis_info.ProteinAnalysisInfo(prot_2_name, prot_2_chains, tmp_batch_analysis)
        #     batch_analysis.append((tmp_prot_1, tmp_prot_2))
        #
        # transformer = data_transformer.DataTransformer(self.ui)
        # # contains analysis-ready data format: list(tuple(prot_1, prot_2, export_dir), ...)
        # batch_data = transformer.transform_data_for_analysis(self.app_project, batch_analysis)
        #
        # for analysis_data in batch_data:
        #     if not os.path.exists(analysis_data[2]):
        #         os.mkdir(analysis_data[2])
        #     else:
        #         basic_boxes.ok("Single Analysis", f"The structure analysis: {analysis_data[3]} already exists!", QtWidgets.QMessageBox.Critical)
        #         self._init_batch_analysis_page()
        #         return
        #
        #     cmd.reinitialize()
        #     structure_analysis_obj = structure_analysis.StructureAnalysis(
        #         reference_protein=[analysis_data[0]], model_proteins=[analysis_data[1]],
        #         ref_chains=analysis_data[0].chains, model_chains=analysis_data[1].chains,
        #         export_dir=analysis_data[2], cycles=self.app_settings.get_cycles(),
        #         cutoff=self.app_settings.get_cutoff(),
        #     )
        #     if self.ui.cb_analysis_images.isChecked():
        #         structure_analysis_obj.response_create_images = True
        #     structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.ref_chains,
        #                                                          structure_analysis_obj.reference_protein)
        #     structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.model_chains,
        #                                                          structure_analysis_obj.model_proteins)
        #     protein_pairs = structure_analysis_obj.create_protein_pairs()
        #     structure_analysis_obj.do_analysis_in_pymol(protein_pairs, self.status_bar)
        #     protein_pairs[0].name = analysis_data[3]
        #     protein_pairs[0].cutoff = self.app_settings.cutoff
        #     self.app_project.add_protein_pair(protein_pairs[0])
        #     protein_pairs[0].serialize_protein_pair(self.app_project.get_objects_protein_pairs_path())
        #     self.app_project.serialize_project(self.app_project.folder_paths["project"], "project")
        constants.PYSSA_LOGGER.info("Begin analysis process.")
        # analysis_thread = threading.Thread(target=self.thread_func_run_analysis, args=(self.analysis_callback,))
        # analysis_thread.start()
        # self.block_box_analysis.exec_()
        # self.worker_analysis = workers.AnalysisWorkerPool(
        #    self.ui.list_analysis_batch_overview, self.ui.cb_analysis_images,
        #    self.status_bar, self.app_project, self.app_settings, self._init_batch_analysis_page)
        constants.PYSSA_LOGGER.info("Thread started for analysis process.")
        # self.threadpool.start(self.worker_analysis)

        # <editor-fold desc="Worker setup">
        # TODO: test code below
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.DistanceAnalysisWorker(
            self.ui.list_analysis_batch_overview,
            self.ui.cb_analysis_images,
            self.status_bar,
            self.app_project,
            self.app_settings,
            self._init_batch_analysis_page,
        )
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_use_page)
        self.tmp_worker.finished.connect(self.post_analysis_process)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
            os.mkdir(constants.SCRATCH_DIR_ANALYSIS)

        # gui_elements_to_show = [
        #     self.ui.btn_analysis_abort,
        # ]
        # gui_elements_to_hide = [
        #     self.ui.btn_save_project,
        #     self.ui.btn_edit_page,
        #     self.ui.btn_view_page,
        #     self.ui.btn_use_page,
        #     self.ui.btn_export_project,
        #     self.ui.btn_close_project,
        #     self.ui.btn_batch_analysis_page,
        #     self.ui.btn_image_analysis_page,
        #     self.ui.btn_results_page,
        #     self.ui.lbl_handle_pymol_session,
        #     self.ui.btn_image_page,
        #     self.ui.btn_hotspots_page,
        # ]
        # gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

        # constants.PYSSA_LOGGER.info("Begin analysis process.")
        # analysis_thread = QThread()
        # constants.PYSSA_LOGGER.info("Created a new analysis thread.")
        # analysis_worker = workers.AnalysisWorker(self.ui.list_analysis_batch_overview,
        #                                               self.ui.cb_analysis_images,
        #                                               self.status_bar,
        #                                               self.app_project,
        #                                               self.app_settings,
        #                                               self._init_batch_analysis_page)
        # constants.PYSSA_LOGGER.info("Created a new analysis worker.")
        # try:
        #     self._thread_controller.create_and_add_new_thread_worker_pair(constants.ANALYSIS_TASK,
        #                                                                   analysis_thread,
        #                                                                   analysis_worker)
        # except ValueError:
        #     constants.PYSSA_LOGGER.info("A valid thread_worker_pair already exists.")
        # for key in self._thread_controller.thread_worker_pairs:
        #     constants.PYSSA_LOGGER.debug(f"The main task {key} is currently controlled by the _thread_controller.")

        # self._thread_controller.thread_worker_pairs.get(constants.ANALYSIS_TASK).worker.update_attributes(
        #     self.ui.list_analysis_batch_overview, self.ui.cb_analysis_images,
        #     self.status_bar, self.app_project, self.app_settings,
        # )
        # self._thread_controller.thread_worker_pairs.get(constants.ANALYSIS_TASK).setup_and_run_thread(self.post_analysis_process)
        self.block_box_analysis.exec_()
        self.display_view_page()

        # self.ui.btn_analysis_start.setEnabled(False)
        # self.status_bar.showMessage("Protein structure analysis started ...")
        # cmd.reinitialize()
        # data_transformer_analysis = data_transformer.DataTransformer(self.ui)
        # transformed_analysis_data = data_transformer_analysis.transform_to_analysis(self.app_project)
        #

        #
        #
        # self.status_bar.showMessage("Checking user input ...")
        # job = job_utils.Job(self.ui.txt_batch_job_name.text(), self.workspace_path)
        # raw_models_input = self.ui.txt_batch_load_model.toPlainText()
        # models = raw_models_input.split("\n")
        #
        #
        #
        # # runs analysis with project creation
        # for model in models:
        #     cmd.reinitialize()
        #     model_file_info = Qt.QtCore.QFileInfo(model)
        #     MODEL_OBJ_NAME = model_file_info.baseName()
        #     MODEL_DIR = model_file_info.canonicalPath()
        #
        #     project = project.Project(f"project_of_{MODEL_OBJ_NAME}",
        #                                           f"{self.workspace_path}/{job.get_job_name()}")
        #     project.create_project_tree()
        #     project.set_job_name(job.get_job_name())
        #     project.set_pdb_file(self.ui.txt_batch_load_reference.text())
        #     project.set_pdb_id(self.ui.txt_batch_load_reference.text())
        #     project.set_pdb_model(MODEL_OBJ_NAME)
        #     project.set_ref_chains(self.ui.txt_batch_chain_ref.text())
        #     project.set_model_chains((self.ui.txt_batch_chain_model.text()))
        #     project.create_xml_file()
        #     job.add_project_to_job(project)
        #
        #     # gets reference filename and filepath
        #     if len(self.ui.txt_batch_load_reference.text()) == 4:
        #         tmp_protein = core.Protein(self.ui.txt_batch_load_reference.text(),
        #                                    export_data_dir=project.get_pdb_path())
        #         tmp_protein.clean_pdb_file()
        #         REFERENCE_OBJ_NAME = self.ui.txt_batch_load_reference.text()
        #         REFERENCE_DIR = project.get_pdb_path()
        #     else:
        #         ref_file_info = Qt.QtCore.QFileInfo(self.ui.txt_batch_load_reference.text())
        #         REFERENCE_OBJ_NAME = ref_file_info.baseName()
        #         REFERENCE_DIR = ref_file_info.canonicalPath()
        #     reference_protein: list[core.Protein] = [core.Protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]
        #     model_proteins: list[core.Protein] = [core.Protein(MODEL_OBJ_NAME, MODEL_DIR)]
        #     export_dir = project.get_results_path()
        #     structure_analysis = structure_analysis.StructureAnalysis(
        #         reference_protein, model_proteins,
        #         project.get_ref_chains().split(","), project.get_model_chains().split(","),
        #         export_dir, cycles=global_variables.global_var_settings_obj.get_cycles(),
        #         cutoff=global_variables.global_var_settings_obj.get_cutoff(),
        #     )
        #     structure_analysis.create_selection_for_proteins(structure_analysis.ref_chains,
        #                                                      structure_analysis.reference_protein)
        #     structure_analysis.create_selection_for_proteins(structure_analysis.model_chains,
        #                                                      structure_analysis.model_proteins)
        #     structure_analysis.do_analysis_in_pymol(structure_analysis.create_protein_pairs(),
        #                                             self.status_bar, self.ui.progress_bar_batch)
        # job.create_xml_file()

    # </editor-fold>

    # <editor-fold desc="Analysis Images">
    def display_image_analysis_page(self) -> None:
        """Displays the analysis image work area."""
        # get all protein pairs without images
        self.ui.list_analysis_images_struct_analysis.clear()
        self.ui.list_analysis_images_creation_struct_analysis.clear()
        for tmp_protein_pair in self.app_project.protein_pairs:
            if len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) == 0:
                self.ui.list_analysis_images_struct_analysis.addItem(tmp_protein_pair.name)
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button, self.ui.btn_image_analysis_page,
        )
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 23, "Analysis Images")
        self.ui.btn_add_analysis_images_struct_analysis.setEnabled(False)
        self.ui.btn_remove_analysis_images_creation_struct_analysis.setEnabled(False)
        self.ui.btn_start_automatic_image_creation.setEnabled(False)

    def analysis_images_enable_add(self) -> None:
        """Enables the add button."""
        self.ui.btn_add_analysis_images_struct_analysis.setEnabled(True)

    def analysis_images_enable_remove(self) -> None:
        """Enables the remove button."""
        self.ui.btn_remove_analysis_images_creation_struct_analysis.setEnabled(True)

    def add_protein_pair_to_image_creation_queue(self) -> None:
        """Adds a protein pair from the list of protein pairs to make images of."""
        protein_pair_to_add = self.ui.list_analysis_images_struct_analysis.currentItem().text()
        self.ui.list_analysis_images_creation_struct_analysis.addItem(protein_pair_to_add)
        self.ui.list_analysis_images_struct_analysis.takeItem(self.ui.list_analysis_images_struct_analysis.currentRow())
        self.ui.btn_add_analysis_images_struct_analysis.setEnabled(False)
        self.analysis_images_check_if_creation_can_start()

    def remove_protein_pair_from_image_creation_queue(self) -> None:
        """Removes a protein pair from the list of protein pairs to make images of."""
        protein_pair_to_remove = self.ui.list_analysis_images_creation_struct_analysis.currentItem()
        self.ui.list_analysis_images_creation_struct_analysis.takeItem(
            self.ui.list_analysis_images_creation_struct_analysis.currentRow(),
        )
        self.ui.list_analysis_images_struct_analysis.addItem(protein_pair_to_remove)
        self.ui.btn_remove_analysis_images_creation_struct_analysis.setEnabled(False)
        self.analysis_images_check_if_creation_can_start()

    def analysis_images_check_if_creation_can_start(self) -> None:
        """Checks if the list of protein pairs which get images are empty."""
        if self.ui.list_analysis_images_creation_struct_analysis.count() > 0:
            self.ui.btn_start_automatic_image_creation.setEnabled(True)
            styles.color_button_ready(self.ui.btn_start_automatic_image_creation)
        else:
            self.ui.btn_start_automatic_image_creation.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_start_automatic_image_creation)

    def post_image_creation_process(self) -> None:
        """Post method after the image creation task is finished."""
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
        self.block_box_images.destroy(True)
        basic_boxes.ok(
            "Analysis Images",
            "All images of all analysis' have been created. Go to results to check the new results.",
            QtWidgets.QMessageBox.Information,
        )
        constants.PYSSA_LOGGER.info("All images of all analysis' have been created.")
        self._init_analysis_image_page()
        self.display_view_page()
        self._project_watcher.show_valid_options(self.ui)

    def start_automatic_image_creation(self) -> None:
        """Sets up the worker for the image creation process."""
        constants.PYSSA_LOGGER.info("Begin image creation process.")
        # self.worker_image_creation = workers.BatchImageWorkerPool(
        #    self.ui.list_analysis_images_struct_analysis, self.ui.list_analysis_images_creation_struct_analysis,
        #    self.status_bar, self.app_project)
        constants.PYSSA_LOGGER.info("Thread started for image creation process.")
        # self.threadpool.start(self.worker_image_creation)

        # <editor-fold desc="Worker setup">
        # TODO: test code below
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.BatchImageWorker(
            self.ui.list_analysis_images_struct_analysis,
            self.ui.list_analysis_images_creation_struct_analysis,
            self.status_bar,
            self.app_project,
        )
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_view_page)
        self.tmp_worker.finished.connect(self.post_image_creation_process)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
            os.mkdir(constants.SCRATCH_DIR_ANALYSIS)

        self.ui.list_analysis_images_struct_analysis.setEnabled(False)
        self.ui.list_analysis_images_creation_struct_analysis.setEnabled(False)

        # gui_elements_to_show = [
        #     self.ui.btn_analysis_abort,
        # ]
        # gui_elements_to_hide = [
        #     self.ui.btn_save_project,
        #     self.ui.btn_edit_page,
        #     self.ui.btn_view_page,
        #     self.ui.btn_use_page,
        #     self.ui.btn_export_project,
        #     self.ui.btn_close_project,
        #     self.ui.btn_batch_analysis_page,
        #     self.ui.btn_image_analysis_page,
        #     self.ui.btn_results_page,
        #     self.ui.lbl_handle_pymol_session,
        #     self.ui.btn_image_page,
        #     self.ui.btn_hotspots_page,
        # ]
        # gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

        self.block_box_images.exec_()

    # </editor-fold>

    # <editor-fold desc="Results page functions">
    def show_analysis_results_options(self) -> None:
        """Shows the combo box of protein pairs."""
        self.results_management.show_stage_x(0)

    def show_results_interactions(self, gui_elements_to_show=None, gui_elements_to_hide=None) -> None:
        """Shows the gui elements of the results page."""
        if gui_elements_to_hide is not None:
            self.results_management.show_gui_elements_stage_x(
                [0, 1],
                [],
                show_specific_elements=[self.ui.lbl_results_analysis_options, self.ui.cb_results_analysis_options],
                hide_specific_elements=gui_elements_to_hide,
            )
        else:
            self.results_management.show_gui_elements_stage_x(
                [0, 1],
                [],
                show_specific_elements=[self.ui.lbl_results_analysis_options, self.ui.cb_results_analysis_options],
            )

    @staticmethod
    def post_load_results() -> None:
        """Prints that the results were loaded."""
        print("Results were loaded.")

    def pre_load_results(self) -> None:
        """Sets up the worker for the result loading process."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        shutil.rmtree(constants.CACHE_DIR)
        os.mkdir(constants.CACHE_DIR)
        os.mkdir(constants.CACHE_IMAGES)
        self.results_name = self.ui.cb_results_analysis_options.currentText()
        if self.results_name == "":
            self.show_analysis_results_options()
            return
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        tmp_protein_pair = self.app_project.search_protein_pair(self.results_name)

        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.LoadResultsWorker(tmp_protein_pair, self.app_project.get_project_xml_path())
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.load_results)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        # distance_data: dict[str, np.ndarray] = tmp_protein_pair.distance_analysis.analysis_results.distance_data
        # distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
        self.ui.list_results_interest_regions.clear()
        self.status_bar.showMessage(f"Loading results of {self.results_name} ...")
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)

    def load_results(self, images_type) -> None:
        """Loads the results of the selected protein pair.

        Args:
            images_type: the type of images which were made during analysis
        """
        if images_type == constants.IMAGES_ALL:
            gui_elements_to_show = [
                self.ui.lbl_results_analysis_options,
                self.ui.cb_results_analysis_options,
                self.ui.lbl_results_rmsd,
                self.ui.txt_results_rmsd,
                self.ui.lbl_color_rmsd,
                self.ui.btn_color_rmsd,
                self.ui.lbl_results_aligned_residues,
                self.ui.txt_results_aligned_residues,
                self.ui.lbl_results_distance_plot,
                self.ui.btn_view_distance_plot,
                self.ui.lbl_results_distance_histogram,
                self.ui.btn_view_distance_histogram,
                self.ui.lbl_results_distance_table,
                self.ui.btn_view_distance_table,
                self.ui.lbl_results_structure_alignment,
                self.ui.btn_view_struct_alignment,
                self.ui.lbl_results_interest_regions,
                self.ui.list_results_interest_regions,
                self.ui.btn_view_interesting_region,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            for tmp_filename in os.listdir(constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                self.ui.list_results_interest_regions.addItem(tmp_filename)
        elif images_type == constants.IMAGES_STRUCT_ALN_ONLY:
            gui_elements_to_show = [
                self.ui.lbl_results_analysis_options,
                self.ui.cb_results_analysis_options,
                self.ui.lbl_results_rmsd,
                self.ui.txt_results_rmsd,
                self.ui.lbl_color_rmsd,
                self.ui.btn_color_rmsd,
                self.ui.lbl_results_aligned_residues,
                self.ui.txt_results_aligned_residues,
                self.ui.lbl_results_distance_plot,
                self.ui.btn_view_distance_plot,
                self.ui.lbl_results_distance_histogram,
                self.ui.btn_view_distance_histogram,
                self.ui.lbl_results_distance_table,
                self.ui.btn_view_distance_table,
                self.ui.lbl_results_structure_alignment,
                self.ui.btn_view_struct_alignment,
            ]
            gui_elements_to_hide = [
                self.ui.lbl_results_interest_regions,
                self.ui.list_results_interest_regions,
                self.ui.btn_view_interesting_region,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        elif images_type == constants.IMAGES_NONE:
            gui_elements_to_show = [
                self.ui.lbl_results_analysis_options,
                self.ui.cb_results_analysis_options,
                self.ui.lbl_results_rmsd,
                self.ui.txt_results_rmsd,
                self.ui.lbl_color_rmsd,
                self.ui.btn_color_rmsd,
                self.ui.lbl_results_aligned_residues,
                self.ui.txt_results_aligned_residues,
                self.ui.lbl_results_distance_plot,
                self.ui.btn_view_distance_plot,
                self.ui.lbl_results_distance_histogram,
                self.ui.btn_view_distance_histogram,
                self.ui.lbl_results_distance_table,
                self.ui.btn_view_distance_table,
            ]
            gui_elements_to_hide = [
                self.ui.lbl_results_structure_alignment,
                self.ui.btn_view_struct_alignment,
                self.ui.lbl_results_interest_regions,
                self.ui.list_results_interest_regions,
                self.ui.btn_view_interesting_region,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            raise ValueError("Illegal argument.")

        tmp_protein_pair = self.app_project.search_protein_pair(self.results_name)
        self.ui.list_results_interest_regions.sortItems()
        self.ui.txt_results_rmsd.setText(str(tmp_protein_pair.distance_analysis.analysis_results.rmsd))
        

        # # <editor-fold desc="Histogram check">
        # # check if histogram can be created
        # distance_list.sort()
        # x, y = np.histogram(distance_list, bins=np.arange(0, distance_list[len(distance_list) - 1], 0.25))
        # if x.size != y.size:
        #     x = np.resize(x, (1, y.size))
        # # this conversion is needed for the pyqtgraph library!
        # x = x.tolist()
        # try:
        #     x = x[0]
        # except IndexError:
        #     # histogram could not be created
        #     gui_elements_to_hide.append(self.ui.lbl_results_distance_histogram)
        #     gui_elements_to_hide.append(self.ui.btn_view_distance_histogram)
        #     gui_utils.hide_gui_elements(gui_elements_to_hide)
        #
        # # </editor-fold>

        cmd.reinitialize()
        tmp_protein_pair.load_pymol_session()
        self.current_session = current_session.CurrentSession(
            "protein_pair", tmp_protein_pair.name, tmp_protein_pair.pymol_session,
        )
        
        self.ui.txt_results_aligned_residues.setText(
            str(tmp_protein_pair.distance_analysis.analysis_results.aligned_aa),
        )
        cmd.scene(f"{tmp_protein_pair.protein_1.get_molecule_object()}-{tmp_protein_pair.protein_2.get_molecule_object()}",
                  action="recall")
        QtWidgets.QApplication.restoreOverrideCursor()
        self.status_bar.showMessage(f"Current workspace: {str(self.workspace_path)}")

    def color_protein_pair_by_rmsd(self) -> None:
        """Colors the residues in 5 colors depending on their distance to the reference."""
        tmp_protein_pair: "protein_pair.ProteinPair" = self.app_project.search_protein_pair(self.results_name)

        cutoff_1 = 0.5
        cutoff_2 = 1.0
        cutoff_3 = 2
        cutoff_4 = 4
        cutoff_5 = 6

        color_1 = "br0"
        color_2 = "br2"
        color_3 = "br4"
        color_4 = "br6"
        color_5 = "br8"
        color_6 = "red"

        cmd.color("hydrogen", tmp_protein_pair.protein_2.get_molecule_object())

        i: int = 0
        for distance_value in tmp_protein_pair.distance_analysis.analysis_results.distance_data.get("distance"):
            if distance_value <= cutoff_1:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                cmd.color(color_1, atom1)
                cmd.color(color_1, atom2)
                # cmd.do(f"color {color_1}, {atom1}")
                # cmd.do(f"color {color_1}, {atom2}")
                i += 1

            elif distance_value <= cutoff_2:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                cmd.color(color_2, atom1)
                cmd.color(color_2, atom2)
                i += 1

            elif distance_value <= cutoff_3:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}/CA"
                )
                atom2: str = (
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}/CA"
                )
                # coloring
                cmd.color(color_3, atom1)
                cmd.color(color_3, atom2)
                i += 1

            elif distance_value <= cutoff_4:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                cmd.color(color_4, atom1)
                cmd.color(color_4, atom2)
                i += 1

            elif distance_value <= cutoff_5:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                cmd.color(color_5, atom1)
                cmd.color(color_5, atom2)
                i += 1

            elif distance_value > cutoff_5:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                cmd.color(color_6, f"({atom1})")
                cmd.color(color_6, f"({atom2})")
                i += 1

        # hide unnecessary representations
        # fixme: it might be a problem to hide any representation at this point
        # cmd.hide("cartoon", tmp_protein_pair.protein_1.get_molecule_object())
        # cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}")
        # cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}")

    def display_structure_alignment(self) -> None:
        """Opens a window which displays the image of the structure alignment."""
        png_dialog = QtWidgets.QDialog(self)
        label = QtWidgets.QLabel(self)
        pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/{self.results_name}")
        self.ui.cb_results_analysis_options.currentText()
        pixmap = QtGui.QPixmap(
            f"{constants.CACHE_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{self.ui.cb_results_analysis_options.currentText()}",
        )
        # TO-DO: Create setting for min. image size
        pixmap = pixmap.scaled(450, 450, transformMode=QtCore.Qt.SmoothTransformation)
        label.setPixmap(pixmap)
        label.setScaledContents(True)
        png_dialog_layout = QtWidgets.QHBoxLayout()
        png_dialog_layout.addWidget(label)
        png_dialog.setLayout(png_dialog_layout)
        png_dialog.setWindowTitle("Image of: structure alignment")
        png_dialog.show()

    def display_distance_plot(self) -> None:
        """Opens a window which displays the distance plot."""
        protein_pair_of_analysis = self.app_project.search_protein_pair(
            self.ui.cb_results_analysis_options.currentText(),
        )
        self.distance_plot_dialog = dialog_distance_plot.DialogDistancePlot(protein_pair_of_analysis)
        self.distance_plot_dialog.setWindowModality(Qt.WindowModal)
        self.is_distance_plot_open = True
        self.distance_plot_dialog.exec_()
        # try:
        #     protein_pair_of_analysis = self.app_project.search_protein_pair(self.ui.cb_results_analysis_options.currentText())
        #     dialog = dialog_distance_plot.DialogDistancePlot(protein_pair_of_analysis)
        #     dialog.exec_()
        # except:
        #     constants.PYSSA_LOGGER.error("The distance plot could not be created, due to an known bug.")
        #     basic_boxes.ok("Display distance plot", "There was a problem with the display of the distance plot.\n"
        #                                             "Try closing the project, or restarting the application.", QtWidgets.QMessageBox.Error)

    def display_distance_histogram(self) -> None:
        """Opens a window which displays the distance histogram."""
        # item = self.ui.project_list.selectedItems()
        # if item is None:
        #     raise ValueError
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        protein_pair_of_analysis = self.app_project.search_protein_pair(
            self.ui.cb_results_analysis_options.currentText(),
        )
        dialog = dialog_distance_histogram.DialogDistanceHistogram(protein_pair_of_analysis)
        dialog.exec_()

    def display_interesting_region(self) -> None:
        """Displays an image of an interesting region."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        png_dialog = QtWidgets.QDialog(self)
        label = QtWidgets.QLabel(self)
        file_name = self.ui.list_results_interest_regions.currentItem().text()
        pixmap = QtGui.QPixmap(f"{constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{file_name}")
        # TO-DO: Create setting for min. image size
        pixmap = pixmap.scaled(450, 450, transformMode=QtCore.Qt.SmoothTransformation)
        label.setPixmap(pixmap)
        label.setScaledContents(True)
        png_dialog_layout = QtWidgets.QHBoxLayout()
        png_dialog_layout.addWidget(label)
        png_dialog.setLayout(png_dialog_layout)
        png_dialog.setWindowTitle(f"Image of: {file_name}")
        png_dialog.show()

    def display_distance_table(self) -> None:
        """Displays the distances in a table."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        csv_model = QtGui.QStandardItemModel()
        csv_model.setColumnCount(7)
        labels = [
            "Residue pair no.",
            "Protein 1 Chain",
            "Protein 1 Position",
            "Protein 1 Residue",
            "Protein 2 Chain",
            "Protein 2 Position",
            "Protein 2 Residue",
            "Distance in Å",
        ]
        csv_model.setHorizontalHeaderLabels(labels)
        table_dialog = QtWidgets.QDialog(self)
        table_view = QtWidgets.QTableView()
        table_view.setModel(csv_model)

        tmp_protein_pair = self.app_project.search_protein_pair(self.ui.cb_results_analysis_options.currentText())
        csv_filepath = pathlib.Path(f"{constants.CACHE_CSV_DIR}/{tmp_protein_pair.name}.csv")
        if not os.path.exists(constants.CACHE_CSV_DIR):
            os.mkdir(constants.CACHE_CSV_DIR)
        tmp_protein_pair = self.app_project.search_protein_pair(self.results_name)

        distance_data = tmp_protein_pair.distance_analysis.analysis_results.distance_data
        distance_data_array = np.array(
            [
                distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI],
                distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES],
            ],
        )
        distance_data_array_transpose = distance_data_array.transpose()
        with open(csv_filepath, mode="w", newline="") as file:
            writer = csv.writer(file, delimiter=",")
            writer.writerows(distance_data_array_transpose)

        file.close()

        with open(csv_filepath, "r", encoding="utf-8") as csv_file:
            i = 0
            for line in csv_file:
                tmp_list = line.split(",")
                # tmp_list.pop(0)
                standard_item_list = []
                pair_no_item = QtGui.QStandardItem()
                pair_no_item.setData(int(tmp_list[0]), role=QtCore.Qt.DisplayRole)
                ref_chain_item = QtGui.QStandardItem()
                ref_chain_item.setData(str(tmp_list[1]), role=QtCore.Qt.DisplayRole)
                ref_pos_item = QtGui.QStandardItem()
                ref_pos_item.setData(int(tmp_list[2]), role=QtCore.Qt.DisplayRole)
                ref_resi_item = QtGui.QStandardItem()
                ref_resi_item.setData(str(tmp_list[3]), role=QtCore.Qt.DisplayRole)
                model_chain_item = QtGui.QStandardItem()
                model_chain_item.setData(str(tmp_list[4]), role=QtCore.Qt.DisplayRole)
                model_pos_item = QtGui.QStandardItem()
                model_pos_item.setData(int(tmp_list[5]), role=QtCore.Qt.DisplayRole)
                model_resi_item = QtGui.QStandardItem()
                model_resi_item.setData(str(tmp_list[6]), role=QtCore.Qt.DisplayRole)
                distance_item = QtGui.QStandardItem()
                distance_item.setData(float(tmp_list[7]), role=QtCore.Qt.DisplayRole)
                standard_item_list.append(pair_no_item)
                standard_item_list.append(ref_chain_item)
                standard_item_list.append(ref_pos_item)
                standard_item_list.append(ref_resi_item)
                standard_item_list.append(model_chain_item)
                standard_item_list.append(model_pos_item)
                standard_item_list.append(model_resi_item)
                standard_item_list.append(distance_item)
                csv_model.insertRow(i, standard_item_list)
            i += 1
        csv_file.close()
        csv_model.removeRow(0)
        table_view.setAlternatingRowColors(True)
        table_view.resizeColumnsToContents()
        table_view.verticalHeader().setVisible(False)
        table_view.setSortingEnabled(True)
        table_view.sortByColumn(0, QtCore.Qt.AscendingOrder)
        table_dialog_layout = QtWidgets.QHBoxLayout()
        table_dialog_layout.addWidget(table_view)
        table_dialog.setLayout(table_dialog_layout)
        # styles
        stylesheet = """
        QDialog {background-color: #F6F4F8;}
        QTableView {background-color: white;}
        """
        table_dialog.setStyleSheet(stylesheet)
        table_dialog.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        table_dialog.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        table_dialog.setModal(True)
        table_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)  # disables editing of cells
        table_dialog.setWindowTitle("Distances of Structure Alignment")
        table_dialog.show()

    # </editor-fold>

    # <editor-fold desc="Manage page functions">
    def choose_manage_open_protein(self) -> None:
        """Shows the different configuration gui elements."""
        if self.ui.box_manage_choose_protein.currentText() != "":
            gui_elements_to_show = [
                self.ui.lbl_manage_choose_color,
                self.ui.box_manage_choose_color,
                self.ui.lbl_manage_choose_representation,
                self.ui.box_manage_choose_representation,
                self.ui.lbl_manage_choose_bg_color,
                self.ui.box_manage_choose_bg_color,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            if cmd.select(name="disulfides", selection=f"{self.ui.box_manage_choose_protein.currentText()} & byres (resn CYS and name SG) within 2 of (resn CYS and name SG)") > 0:
                gui_elements_to_show = [
                    self.ui.lbl_disulfid_bond_1,
                    self.ui.lbl_disulfid_bond_2,
                    self.ui.btn_disulfid_bond_show,
                    self.ui.btn_disulfid_bond_hide,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
        else:
            gui_elements_to_hide = [
                self.ui.lbl_manage_choose_color,
                self.ui.box_manage_choose_color,
                self.ui.lbl_manage_choose_representation,
                self.ui.box_manage_choose_representation,
                self.ui.lbl_manage_choose_bg_color,
                self.ui.box_manage_choose_bg_color,
                self.ui.lbl_disulfid_bond_1,
                self.ui.lbl_disulfid_bond_2,
                self.ui.btn_disulfid_bond_show,
                self.ui.btn_disulfid_bond_hide,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        
    def choose_manage_color_selected_protein(self) -> None:
        """Sets the protein color."""
        input = self.ui.box_manage_choose_protein.currentText()
        tmp_protein = self.app_project.search_protein(input)
        tmp_protein.color_protein_in_pymol(
            self.ui.box_manage_choose_color.currentText(), f"/{tmp_protein.get_molecule_object()}",
        )

    def choose_manage_representation(self) -> None:
        """Sets the representation."""
        input = self.ui.box_manage_choose_protein.currentText()
        tmp_protein = self.app_project.search_protein(input)
        tmp_selection = f"/{tmp_protein.get_molecule_object()}"
        if self.ui.box_manage_choose_representation.currentIndex() == 0:
            print("Please select a representation.")
            self.status_bar.showMessage("Please select a representation.")
        elif self.ui.box_manage_choose_representation.currentIndex() == 1:
            cmd.show("cartoon", tmp_selection)
            cmd.hide("ribbon", tmp_selection)
        elif self.ui.box_manage_choose_representation.currentIndex() == 2:
            cmd.show("ribbon", tmp_selection)
            cmd.hide("cartoon", tmp_selection)
        else:
            print("Missing implementation!")

    def choose_manage_bg_color(self) -> None:
        """Sets the background color."""
        if self.ui.box_manage_choose_bg_color.currentIndex() == 0:
            print("Please select a background color.")
            self.status_bar.showMessage("Please select a background color.")
        elif self.ui.box_manage_choose_bg_color.currentIndex() == 1:
            cmd.bg_color("black")
        elif self.ui.box_manage_choose_bg_color.currentIndex() == 2:
            cmd.bg_color("white")
        else:
            print("Missing implementation!")
    
    def show_disulfid_bonds_as_sticks(self):
        cmd.select(name="disulfides", selection=f"{self.ui.box_manage_choose_protein.currentText()} & byres (resn CYS and name SG) within 2 of (resn CYS and name SG)")
        cmd.color(color="atomic", selection="disulfides and not elem C")
        cmd.set("valence", 0)  # this needs to be better implemented
        cmd.show("sticks", "disulfides")
        cmd.hide("sticks", "elem H")
    
    def hide_disulfid_bonds_as_sticks(self):
        cmd.select(name="disulfides",
                   selection=f"{self.ui.box_manage_choose_protein.currentText()} & byres (resn CYS and name SG) within 2 of (resn CYS and name SG)")
        cmd.hide("sticks", "disulfides")
    
    # </editor-fold>

    # <editor-fold desc="Image page functions">
    def show_representation(self) -> None:
        """Sets the representation."""
        if self.ui.box_representation.currentIndex() == 0:
            print("Please select a representation.")
            self.status_bar.showMessage("Please select a representation.")
        elif self.ui.box_representation.currentIndex() == 1:
            cmd.show("cartoon", "all")
            cmd.hide("ribbon", "all")
        elif self.ui.box_representation.currentIndex() == 2:
            cmd.show("ribbon", "all")
            cmd.hide("cartoon", "all")
        else:
            print("Missing implementation!")

    def choose_bg_color(self) -> None:
        """Sets the background color."""
        if self.ui.box_bg_color.currentIndex() == 0:
            print("Please select a background color.")
            self.status_bar.showMessage("Please select a background color.")
        elif self.ui.box_bg_color.currentIndex() == 1:
            cmd.bg_color("black")
        elif self.ui.box_bg_color.currentIndex() == 2:
            cmd.bg_color("white")
        else:
            print("Missing implementation!")

    def choose_renderer(self) -> None:
        """Sets the renderer."""
        if self.ui.box_renderer.currentIndex() == 0:
            self.status_bar.showMessage("Please select a renderer.")
            self.ui.cb_ray_tracing.hide()
            self.ui.label_26.hide()
        elif self.ui.box_renderer.currentIndex() == 1:
            self.renderer = "-1"
            self.ui.cb_ray_tracing.show()
            self.ui.label_26.show()
        elif self.ui.box_renderer.currentIndex() == 2:
            self.renderer = "0"
            self.ui.cb_ray_tracing.show()
            self.ui.label_26.show()
        else:
            print("Missing implementation!")

    def choose_ray_trace_mode(self) -> None:
        """Sets the ray-trace mode."""
        if self.ui.box_ray_trace_mode.currentIndex() == 0:
            print("Please select a Ray-Trace-Mode.")
            self.status_bar.showMessage("Please select a Ray-Trace-Mode.")
        elif self.ui.box_ray_trace_mode.currentIndex() == 1:
            cmd.set("ray_trace_mode", 0)
        elif self.ui.box_ray_trace_mode.currentIndex() == 2:
            cmd.set("ray_trace_mode", 1)
        elif self.ui.box_ray_trace_mode.currentIndex() == 3:
            cmd.set("ray_trace_mode", 2)
        elif self.ui.box_ray_trace_mode.currentIndex() == 4:
            cmd.set("ray_trace_mode", 3)
        else:
            print("Missing implementation!")

    def choose_ray_texture(self) -> None:
        """Sets the ray texture."""
        if self.ui.box_ray_texture.currentIndex() == 0:
            print("Please select a Ray Texture.")
            self.status_bar.showMessage("Please select a Ray Texture.")
        elif self.ui.box_ray_texture.currentIndex() == 1:
            cmd.set("ray_texture", 0)
        elif self.ui.box_ray_texture.currentIndex() == 2:
            cmd.set("ray_texture", 1)
        elif self.ui.box_ray_texture.currentIndex() == 3:
            cmd.set("ray_texture", 2)
        elif self.ui.box_ray_texture.currentIndex() == 4:
            cmd.set("ray_texture", 3)
        elif self.ui.box_ray_texture.currentIndex() == 5:
            cmd.set("ray_texture", 4)
        elif self.ui.box_ray_texture.currentIndex() == 6:
            cmd.set("ray_texture", 5)
        else:
            print("Missing implementation!")
    
    def decide_ray_tracing(self) -> None:
        """Sets the ray-tracing options."""
        if self.ui.cb_ray_tracing.isChecked():
            self.ui.cb_transparent_bg.hide()
            self.ui.label_23.hide()
            self.ui.box_renderer.setEnabled(False)
            self.ui.label_10.show()
            self.ui.box_ray_trace_mode.show()
            self.ui.label_14.show()
            self.ui.box_ray_texture.show()
            cmd.set("ray_opaque_background", "off")
        else:
            self.ui.cb_transparent_bg.show()
            self.ui.label_23.show()
            self.ui.box_renderer.setEnabled(True)
            self.ui.label_10.hide()
            self.ui.box_ray_trace_mode.hide()
            self.ui.label_14.hide()
            self.ui.box_ray_texture.hide()
    
    def decide_transparent_bg(self) -> None:
        """Sets the transparent background."""
        if self.ui.cb_transparent_bg.isChecked():
            cmd.set("opaque_background", "off")
        else:
            cmd.set("opaque_background", "on")

    @staticmethod
    def update_scene() -> None:
        """Updates the current selected PyMOL scene."""
        cmd.scene(key="auto", action="update")

    def save_scene(self) -> None:
        """Saves the current view as a new PyMOL scene."""
        # returns tuple with (name, bool)
        scene_name = QtWidgets.QInputDialog.getText(self, "Save Scene", "Enter scene name:")
        if scene_name[1]:
            cmd.scene(key=scene_name[0], action="append")

    def post_preview_image(self) -> None:
        """Hides the block box of the preview process."""
        self.block_box_uni.hide()
        self.block_box_uni.destroy(True)
        self.status_bar.showMessage("Finished preview of ray-traced image.")
        QtWidgets.QApplication.restoreOverrideCursor()

    def preview_image(self) -> None:
        """Previews the image."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        if self.ui.cb_ray_tracing.isChecked():
            self.status_bar.showMessage("Preview ray-traced image ...")
            # <editor-fold desc="Worker setup">
            # --Begin: worker setup
            self.tmp_thread = QtCore.QThread()
            self.tmp_worker = task_workers.PreviewRayImageWorker(self.renderer)
            self.tmp_thread = task_workers.setup_worker_for_work(
                self.tmp_thread, self.tmp_worker, self.display_view_page,
            )
            self.tmp_worker.finished.connect(self.post_preview_image)
            self.tmp_thread.start()
            # --End: worker setup

            # </editor-fold>
            gui_utils.setup_standard_block_box(
                self.block_box_uni, "Preview ray-trace image", "Creating preview for the ray-traced image ...",
            )
            self.block_box_uni.exec_()
        else:
            self.status_bar.showMessage("Preview draw image ...")
            cmd.draw(2400, 2400)
            self.status_bar.showMessage("Finished preview of drawn image.")
            QtWidgets.QApplication.restoreOverrideCursor()

    def post_save_image(self) -> None:
        """Displays a message box which informs that the process has finished."""
        self.block_box_uni.hide()
        self.block_box_uni.destroy(True)
        self.status_bar.showMessage("Finished image creation.")
        QtWidgets.QApplication.restoreOverrideCursor()
        basic_boxes.ok("Finished image creation", "The image has been created.", QtWidgets.QMessageBox.Information)

    def save_image(self) -> None:
        """Saves the image as a png file."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        if self.ui.cb_ray_tracing.isChecked():
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display(
                        "info", "No file has been selected.", self.status_bar, "No file has been selected.",
                    )
                    return
                self.status_bar.showMessage("Creating ray-traced image ...")

                # <editor-fold desc="Worker setup">
                # --Begin: worker setup
                self.tmp_thread = QtCore.QThread()
                self.tmp_worker = task_workers.SaveRayImageWorker(self.renderer, full_file_name[0])
                self.tmp_thread = task_workers.setup_worker_for_work(
                    self.tmp_thread, self.tmp_worker, self.display_view_page,
                )
                self.tmp_worker.finished.connect(self.post_save_image)
                self.tmp_thread.start()
                # --End: worker setup

                # </editor-fold>
                gui_utils.setup_standard_block_box(
                    self.block_box_uni, "Save ray-trace image", "Creating the ray-traced image ...",
                )
                self.block_box_uni.exec_()

                # cmd.ray(2400, 2400, renderer=int(self.renderer))
                # cmd.png(full_file_name[0], dpi=300)

            except FileExistsError:
                tools.quick_log_and_display("error", "File exists already.", self.status_bar, "File exists already.")
            except pymol.CmdException:
                tools.quick_log_and_display(
                    "error",
                    "Unexpected Error from PyMOL while saving the " "an image",
                    self.status_bar,
                    "Unexpected Error from PyMOL",
                )
        else:
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display(
                        "info", "No file has been selected.", self.status_bar, "No file has been selected.",
                    )
                    return
                self.status_bar.showMessage("Creating draw image ...")
                cmd.draw(2400, 2400)
                cmd.png(full_file_name[0], dpi=300)
                self.status_bar.showMessage("Finished image creation.")
                basic_boxes.ok(
                    "Finished image creation", "The image has been created.", QtWidgets.QMessageBox.Information,
                )
            except FileExistsError:
                tools.quick_log_and_display("error", "File exists already.", self.status_bar, "File exists already.")
            except pymol.CmdException:
                tools.quick_log_and_display(
                    "error",
                    "Unexpected Error from PyMOL while saving the " "an image",
                    self.status_bar,
                    "Unexpected Error from PyMOL",
                )
            finally:
                QtWidgets.QApplication.restoreOverrideCursor()

    # </editor-fold>

    # <editor-fold desc="Hotspots page functions">
    def open_protein(self) -> None:
        """Loads the selected protein's pymol session."""
        try:
            test = self.ui.list_hotspots_choose_protein.currentItem().text()
        except AttributeError:
            return
        if self.ui.list_hotspots_choose_protein.currentItem().text() != "":
            tools.ask_to_save_pymol_session(self.app_project, self.current_session)

            input = self.ui.list_hotspots_choose_protein.currentItem().text()
            if input.find("_vs_") == -1:
                # one protein is selected
                tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
                tmp_protein.load_protein_pymol_session()
                self.current_session.name = tmp_protein.get_molecule_object()
                self.current_session.type = "protein"
                self.current_session.session = tmp_protein.pymol_session
                cmd.set("seq_view", 1)
                # tmp_protein.pymol_selection.set_selections_without_chains_ca()
                # tmp_model = cmd.get_model(tmp_protein.pymol_selection.selection_string)
                # first_amoino_acid_no = tmp_model.atom[0].resi
                # self.ui.sp_hotspots_resi_no.setMinimum(int(first_amoino_acid_no))
                # tmp_sequence = cmd.get_fastastr('all')
                # self.ui.sp_hotspots_resi_no.setMaximum(len(tmp_sequence))
                gui_elements_to_show = [
                    self.ui.lbl_hotspots_resi_show,
                    self.ui.btn_hotspots_resi_show,
                    self.ui.lbl_hotspots_resi_hide,
                    self.ui.btn_hotspots_resi_hide,
                    self.ui.lbl_hotspots_resi_zoom,
                    self.ui.btn_hotspots_resi_zoom,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
            else:
                # protein pair is selected
                tmp_protein_pair = self.app_project.search_protein_pair(input)
                gui_elements_to_show = [
                    self.ui.lbl_hotspots_resi_show,
                    self.ui.btn_hotspots_resi_show,
                    self.ui.lbl_hotspots_resi_hide,
                    self.ui.btn_hotspots_resi_hide,
                    self.ui.lbl_hotspots_resi_zoom,
                    self.ui.btn_hotspots_resi_zoom,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                tmp_protein_pair.load_pymol_session()
                self.current_session.name = tmp_protein_pair.name
                self.current_session.type = "protein_pair"
                self.current_session.session = tmp_protein_pair.pymol_session
                cmd.set("seq_view", 1)
                # prot1_no = int(tools.check_first_seq_no(tmp_protein_pair.protein_1))
                # prot2_no = int(tools.check_first_seq_no(tmp_protein_pair.protein_2))
                # if prot1_no > prot2_no:
                #     self.ui.sp_hotspots_resi_no.setMinimum(prot2_no)
                #     self.used_proteins.prot_1 = tmp_protein_pair.protein_1.get_molecule_object()
                #     self.used_proteins.prot_1_seq_start = prot1_no
                #     self.used_proteins.prot_1_seq_end = tools.check_seq_length(tmp_protein_pair.protein_1)
                #     self.used_proteins.prot_2 = tmp_protein_pair.protein_2.get_molecule_object()
                #     self.used_proteins.prot_2_seq_start = prot2_no
                #     self.used_proteins.prot_2_seq_end = tools.check_seq_length(tmp_protein_pair.protein_2)
                # else:
                #     self.ui.sp_hotspots_resi_no.setMinimum(prot1_no)
                #     self.ui.sp_hotspots_resi_no.setMinimum(prot2_no)
                #     self.used_proteins.prot_1 = tmp_protein_pair.protein_1.get_molecule_object()
                #     self.used_proteins.prot_1_seq_start = prot1_no
                #     self.used_proteins.prot_1_seq_end = tools.check_seq_length(tmp_protein_pair.protein_1)
                #     self.used_proteins.prot_2 = tmp_protein_pair.protein_2.get_molecule_object()
                #     self.used_proteins.prot_2_seq_start = prot2_no
                #     self.used_proteins.prot_2_seq_end = tools.check_seq_length(tmp_protein_pair.protein_2)
                #
                # prot1_seq_len = int(tools.check_seq_length(tmp_protein_pair.protein_1))
                # prot2_seq_len = int(tools.check_seq_length(tmp_protein_pair.protein_2))
                # if prot1_seq_len > prot2_seq_len:
                #     self.ui.sp_hotspots_resi_no.setMaximum(prot2_seq_len)
                #     self.ui.sp_hotspots_resi_no.setMinimum(prot2_no)
                #     self.used_proteins.prot_1 = tmp_protein_pair.protein_1.get_molecule_object()
                #     self.used_proteins.prot_1_seq_start = prot1_no
                #     self.used_proteins.prot_1_seq_end = tools.check_seq_length(tmp_protein_pair.protein_1)
                #     self.used_proteins.prot_2 = tmp_protein_pair.protein_2.get_molecule_object()
                #     self.used_proteins.prot_2_seq_start = prot2_no
                #     self.used_proteins.prot_2_seq_end = tools.check_seq_length(tmp_protein_pair.protein_2)
                # else:
                #     self.ui.sp_hotspots_resi_no.setMaximum(prot1_seq_len)
                #     self.ui.sp_hotspots_resi_no.setMinimum(prot2_no)
                #     self.used_proteins.prot_1 = tmp_protein_pair.protein_1.get_molecule_object()
                #     self.used_proteins.prot_1_seq_start = prot1_no
                #     self.used_proteins.prot_1_seq_end = tools.check_seq_length(tmp_protein_pair.protein_1)
                #     self.used_proteins.prot_2 = tmp_protein_pair.protein_2.get_molecule_object()
                #     self.used_proteins.prot_2_seq_start = prot2_no
                #     self.used_proteins.prot_2_seq_end = tools.check_seq_length(tmp_protein_pair.protein_2)
        self.ui.btn_manage_session.show()

    @staticmethod
    def hide_resi_sticks() -> None:
        """Hides the balls and sticks representation of the pymol selection."""
        session_util.check_if_sele_is_empty()
        cmd.hide(representation="sticks", selection="sele")

    @staticmethod
    def show_resi_sticks() -> None:
        """Shows the pymol selection as sticks."""
        session_util.check_if_sele_is_empty()
        cmd.show(representation="sticks", selection="sele and not hydrogens")
        cmd.select(name="sele", selection="sele and not hydrogens")
        cmd.color(color="atomic", selection="sele and not elem C")
        cmd.set("valence", 0)  # this needs to be better implemented

    @staticmethod
    def zoom_resi_position() -> None:
        """Zooms to the pymol selection."""
        session_util.check_if_sele_is_empty()
        cmd.zoom(selection="sele", buffer=8.0, state=0, complete=0)

    # </editor-fold>


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    styles.set_stylesheet(app)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())
