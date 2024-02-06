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
"""Module for the main view presenter of the pyssa plugin."""
import logging
import os
import shutil
import subprocess
import sys
import pathlib
import csv

import numpy as np
import pymol
from typing import TYPE_CHECKING
from pymol import cmd
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.gui.ui.dialogs import dialog_settings_global
from pyssa.gui.ui.dialogs import dialog_help
from pyssa.gui.ui.dialogs import dialog_startup
from pyssa.gui.ui.dialogs import dialog_distance_histogram
from pyssa.gui.ui.dialogs import dialog_about
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.dialogs import dialog_tutorial_videos
from pyssa.gui.ui.dialogs import dialog_rename_protein
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import plot_view, add_protein_view

from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures import project_watcher
from pyssa.internal.data_structures import settings
from pyssa.internal.data_structures.data_classes import prediction_configuration, prediction_protein_info
from pyssa.internal.data_structures.data_classes import image_state
from pyssa.internal.data_structures.data_classes import stage
from pyssa.internal.data_structures.data_classes import current_session
from pyssa.internal.data_structures.data_classes import main_window_state
from pyssa.internal.data_structures.data_classes import results_state
from pyssa.internal.thread import task_workers, tasks

from pyssa.io_pyssa import safeguard, filesystem_helpers
from pyssa.io_pyssa import filesystem_io
from pyssa.io_pyssa import path_util

from pyssa.util import pyssa_keys, prediction_util
from pyssa.util import main_window_util
from pyssa.util import exit_codes
from pyssa.util import globals
from pyssa.util import session_util
from pyssa.util import exception
from pyssa.util import constants
from pyssa.util import input_validator
from pyssa.util import gui_page_management
from pyssa.util import tools
from pyssa.util import gui_utils
from pyssa.logging_pyssa import log_handlers
from pyssa.presenter import main_presenter_async

if TYPE_CHECKING:
    from pyssa.gui.ui.views import main_view

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class MainPresenter:
    """Class for main presenter of the pyssa plugin."""

    """
    The main view of the pyssa plugin.
    """
    _view: "main_view.MainView"

    """
    The path to the active workspace.
    """
    _workspace_path: pathlib.Path

    """
    The application settings object.
    """
    _application_settings: "settings.Settings"

    """
    The currently active project.
    """
    _current_project: "project.Project"

    """
    A watcher for the project state.
    """
    _project_watcher: "project_watcher.ProjectWatcher"

    """
    The active task of the application.
    """
    _active_task: tasks.Task

    def __init__(self, a_view: "main_view.MainView") -> None:
        """Constructor.

        Args:
            a_view: the main view of the pyssa plugin.

        Raises:
            IllegalArgumentException: If an argument is illegal.
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_view, logger)

        # </editor-fold>

        # <editor-fold desc="Initialize class attributes">
        self._view = a_view
        self._workspace_path = constants.DEFAULT_WORKSPACE_PATH
        # TODO: settings and project objects get not initialized properly
        self._application_settings = settings.Settings("", "")
        self._current_project: "project.Project" = project.Project()
        self._project_watcher: "project_watcher.ProjectWatcher" = project_watcher.ProjectWatcher(self._current_project)

        # </editor-fold>

        # Check OS
        main_window_util.check_operating_system()
        # Create log dir
        # <editor-fold desc="Program directory check">
        filesystem_helpers.create_directory(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa"))
        filesystem_helpers.create_directory(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/logs"))
        constants.PYSSA_LOGGER.info("Checked program and logs directory.")

        # </editor-fold>
        # Configure settings
        # <editor-fold desc="Setup App Settings">
        self._application_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        if not os.path.exists(constants.SETTINGS_FULL_FILEPATH):
            constants.PYSSA_LOGGER.info("Settings file not found, open configuration dialog.")
            # Configuration dialog to setup setting file
            dialog = dialog_startup.DialogStartup()
            dialog.exec_()

            # checks if the cancel button was pressed
            if dialog_startup.global_var_terminate_app == 1:
                os.remove(constants.SETTINGS_FULL_FILEPATH)
                constants.PYSSA_LOGGER.info("Configuration dialog closed, and removed new settings file.")
                sys.exit()

            self._application_settings.app_launch = 1
            self._application_settings.workspace_path = pathlib.Path(dialog_startup.global_var_startup_workspace)

            constants.PYSSA_LOGGER.info("Demo projects are getting downloaded and extracted ...")
            import zipfile

            with zipfile.ZipFile(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"), "r") as zip_ref:
                zip_ref.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
            constants.PYSSA_LOGGER.info(
                "Demo projects are downloaded and extracted.\n Import of demo projects started ...",
            )

            path_of_demo_projects = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
            tmp_project = project.Project("", dialog_startup.global_var_startup_workspace)
            for tmp_filename in os.listdir(path_of_demo_projects):
                try:
                    tmp_project = tmp_project.deserialize_project(
                        pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}"),
                        self._application_settings,
                    )
                except exception.IllegalArgumentError:
                    constants.PYSSA_LOGGER.warning(
                        "The workspace path does not exist on this system, " "but this is due to the demo projects.",
                    )
                tmp_project.set_workspace_path(dialog_startup.global_var_startup_workspace)
                new_filepath = pathlib.Path(f"{dialog_startup.global_var_startup_workspace}/{tmp_filename}")
                tmp_project.serialize_project(new_filepath)
            constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
            try:
                os.remove(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"))
            except FileNotFoundError:
                constants.PYSSA_LOGGER.warning("Zip archive of demo projects could not be found!")
            constants.PYSSA_LOGGER.info("Serialize settings ...")
            self._application_settings.serialize_settings()
            constants.PYSSA_LOGGER.info("Serialize settings finished.")

            QtWidgets.QApplication.restoreOverrideCursor()

        globals.g_settings = main_window_util.setup_app_settings(self._application_settings)
        self._application_settings = globals.g_settings
        # </editor-fold>
        # Check version number
        main_window_util.check_version_number()

        # <editor-fold desc="Class attributes">
        self._workspace_path = self._application_settings.workspace_path
        self._workspace_status = f"Current workspace: {str(self._workspace_path)}"
        self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")
        self.prediction_configuration = prediction_configuration.PredictionConfiguration(True, "pdb70")
        self.results_name = ""
        self.no_of_selected_chains = 0
        self.plot_dialog = QtWidgets.QDialog(self._view)
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

        self.main_window_state = main_window_state.MainWindowState(
            results_state.ResultsState(),
            image_state.ImageState(),
        )

        constants.PYSSA_LOGGER.info("Setup class attributes finished.")
        # </editor-fold>

        # sets up the status bar
        self._setup_statusbar()
        tools.create_directory(constants.SETTINGS_DIR, "scratch")
        self._setup_default_configuration()

        # if len(os.listdir(constants.LOG_PATH)) > 0:
        #     self.open_change_log()

        constants.PYSSA_LOGGER.info("Setup rest of GUI related elements ...")

        # <editor-fold desc="GUI page management">
        # -- Gui page management vars
        self.batch_analysis_management: gui_page_management.GuiPageManagement
        self.results_management: gui_page_management.GuiPageManagement

        # management functions
        self._create_batch_analysis_management()
        self._create_results_management()

        # </editor-fold>

        # <editor-fold desc="Block box definitions">
        self.block_box_analysis = basic_boxes.no_buttons(
            "Analysis",
            "An analysis is currently running, please wait.",
            QtWidgets.QMessageBox.Information,
        )
        self.block_box_prediction: QtWidgets.QMessageBox = QtWidgets.QMessageBox()
        self.block_box_images = basic_boxes.no_buttons(
            "Analysis Images",
            "Images getting created, please wait.",
            QtWidgets.QMessageBox.Information,
        )
        self.block_box_uni = basic_boxes.no_buttons("Generic", "Generic", QtWidgets.QMessageBox.Information)

        # </editor-fold>

        # configure gui element properties
        self._view.ui.txt_results_aligned_residues.setAlignment(QtCore.Qt.AlignRight)
        self._view.ui.table_pred_mono_prot_to_predict.setSizeAdjustPolicy(
            QtWidgets.QAbstractScrollArea.AdjustToContents,
        )
        self._view.ui.table_pred_mono_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self._view.ui.table_pred_multi_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self._view.ui.list_pred_analysis_multi_ref_chains.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection,
        )
        self._view.ui.list_pred_analysis_multi_model_chains.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection,
        )

        # helper attributes
        self.pymol_session_specs = {
            pyssa_keys.SESSION_SPEC_PROTEIN: [0, ""],
            pyssa_keys.SESSION_SPEC_COLOR: [0, ""],
            pyssa_keys.SESSION_SPEC_REPRESENTATION: [0, ""],
            pyssa_keys.SESSION_SPEC_BG_COLOR: [0, ""],
        }

        # <editor-fold desc="Setup defaults for pages">
        self._init_fill_combo_boxes()
        self._init_new_page()
        self._init_use_page()
        self._init_local_pred_mono_page()
        self._init_local_pred_multi_page()
        self._init_batch_analysis_page()
        self._view.ui.action_toggle_notebook_visibility.setVisible(False)
        self._view.ui.action_settings_model_w_off_colab_notebook.setVisible(False)
        self.last_sidebar_button = QtWidgets.QPushButton()
        self._view.ui.table_pred_mono_prot_to_predict.setEditTriggers(
            self._view.ui.table_pred_mono_prot_to_predict.NoEditTriggers,
        )
        self._view.ui.table_pred_multi_prot_to_predict.setEditTriggers(
            self._view.ui.table_pred_multi_prot_to_predict.NoEditTriggers,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setEditTriggers(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.NoEditTriggers,
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.setEditTriggers(
            self._view.ui.table_pred_analysis_multi_prot_to_predict.NoEditTriggers,
        )

        # TODO: temp gui changes (need to be changed in the designer)
        self._view.ui.lbl_hotspots_resi_show.setText("Residue(s) as sticks")
        self._view.ui.lbl_hotspots_resi_hide.setText("Residue(s) as sticks")

        # fixme: should the pdf documentation be accessible through the pyssa gui?
        self._view.ui.action_help_docs_pdf.setText("Documentation")
        self._view.ui.action_help_docs_pdf.setVisible(True)
        self._view.ui.action_help_docs.setText("Tutorials")
        self._view.ui.action_help_docs.setVisible(True)

        if self._application_settings.wsl_install == 1 and self._application_settings.local_colabfold == 0:
            self._view.ui.action_install_from_file.setVisible(True)
        else:
            self._view.ui.action_install_from_file.setVisible(False)
        # </editor-fold>

        self._connect_all_gui_elements()

        self._project_watcher.show_valid_options(self._view.ui)
        self._project_watcher.check_workspace_for_projects(self._workspace_path, self._view.ui)

        self.threadpool = QtCore.QThreadPool()
        # create scratch and cache dir
        try:
            filesystem_helpers.delete_directory(constants.SCRATCH_DIR)
        except Exception as e:
            constants.PYSSA_LOGGER.warning(f"Scratch path could not be deleted. {e}")
        filesystem_helpers.create_directory(constants.SCRATCH_DIR)
        filesystem_helpers.create_directory(constants.CACHE_DIR)

    def start_worker_thread(self, worker_obj: object, post_process_func) -> None:  # noqa: ANN001
        """Sets up the worker, moves the worker to a thread and starts the thread.

        Args:
            worker_obj: an object of type <work>Worker (QObject).
            post_process_func: a function which sould be executed if the worker is finished.
        """
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = worker_obj
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, post_process_func)
        self.tmp_thread.start()

    def _setup_statusbar(self) -> None:
        """Sets up the status bar and fills it with the current workspace."""
        self._view.setStatusBar(self._view.status_bar)
        self._view.status_bar.showMessage(str(self._workspace_label.text()))

    def _setup_default_configuration(self) -> None:
        """Sets up the default values for specific gui elements."""
        self._view.ui.lbl_current_project_name.setText("")
        # menu
        # side menu

        # new project page
        self._view.ui.btn_new_create_project.setEnabled(False)
        self._view.ui.cb_new_add_reference.setCheckable(False)
        self._view.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
        # open project page
        self._view.ui.lbl_open_status_search.setText("")
        self._view.ui.btn_open_open_project.setEnabled(False)
        # delete project page
        self._view.ui.lbl_delete_status_search.setText("")
        self._view.ui.btn_delete_delete_project.setEnabled(False)
        # edit project page

        # view project page

        # use project page

        # new sequence page
        # sequence vs .pdb page
        self._view.ui.btn_s_v_p_start.setEnabled(False)
        self._view.ui.list_s_v_p_ref_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        # single analysis page
        self._view.ui.lbl_analysis_model_chains.hide()
        self._view.ui.list_analysis_model_chains.hide()
        self._view.ui.btn_analysis_back.hide()
        self._view.ui.btn_analysis_start.hide()
        # batch analysis page

        # results page

        # image page

    def _connect_all_gui_elements(self) -> None:
        """Connects all gui elements with their corresponding slots."""
        self._view.ui.btn_info.clicked.connect(self.open_page_information)

        # <editor-fold desc="Menu">
        self._view.ui.action_file_quit.triggered.connect(self.quit_app)
        self._view.ui.action_file_restore_settings.triggered.connect(self.restore_settings)
        self._view.ui.action_settings_edit_all.triggered.connect(self.open_settings_global)
        self._view.ui.action_help_docs.triggered.connect(self.open_tutorial)
        self._view.ui.action_help_docs_pdf.triggered.connect(self.open_documentation)
        self._view.ui.action_help_about.triggered.connect(self.open_about)
        self._view.ui.action_settings_open_logs.triggered.connect(self.open_logs)
        self._view.ui.action_settings_clear_logs.triggered.connect(self.clear_all_log_files)
        self._view.ui.action_help_changelog.triggered.connect(self.open_release_notes_in_standard_application)
        # </editor-fold>

        # <editor-fold desc="Side Menu">
        self._view.ui.btn_new_page.clicked.connect(self.display_new_page)
        self._view.ui.btn_open_page.clicked.connect(self.display_open_page)
        self._view.ui.btn_delete_page.clicked.connect(self.display_delete_page)
        self._view.ui.btn_save_project.clicked.connect(self.save_project)
        self._view.ui.btn_edit_page.clicked.connect(self.display_edit_page)
        self._view.ui.btn_view_page.clicked.connect(self.display_view_page)
        self._view.ui.btn_use_page.clicked.connect(self.display_use_page)
        self._view.ui.btn_import_project.clicked.connect(self.import_project)
        self._view.ui.btn_export_project.clicked.connect(self.export_current_project)
        self._view.ui.btn_close_project.clicked.connect(self.close_project)
        self._view.ui.btn_pred_cloud_monomer_page.clicked.connect(self.display_esm_pred_mono)
        self._view.ui.btn_pred_local_monomer_page.clicked.connect(self.display_local_pred_mono)
        self._view.ui.btn_pred_local_multimer_page.clicked.connect(self.display_local_pred_multi)
        self._view.ui.btn_prediction_abort.clicked.connect(self.abort_prediction)
        self._view.ui.btn_pred_analysis_monomer_page.clicked.connect(self.display_monomer_pred_analysis)
        self._view.ui.btn_pred_analysis_multimer_page.clicked.connect(self.display_multimer_pred_analysis)
        self._view.ui.btn_batch_analysis_page.clicked.connect(self.display_job_analysis_page)
        self._view.ui.btn_image_analysis_page.clicked.connect(self.display_image_analysis_page)
        self._view.ui.btn_results_page.clicked.connect(self.display_results_page)
        # self._view.ui.btn_analysis_abort.clicked.connect(self.abort_analysis)
        self._view.ui.btn_manage_session.clicked.connect(self.display_manage_pymol_session)
        self._view.ui.btn_image_page.clicked.connect(self.display_image_page)
        self._view.ui.btn_hotspots_page.clicked.connect(self.display_hotspots_page)

        # </editor-fold>

        # <editor-fold desc="New project page">
        self._view.ui.btn_new_choose_reference.clicked.connect(self.load_reference_in_project)
        self._view.ui.txt_new_project_name.textChanged.connect(self.validate_project_name)
        self._view.ui.txt_new_choose_reference.textChanged.connect(self.validate_reference_in_project)
        self._view.ui.cb_new_add_reference.stateChanged.connect(self.show_add_reference)
        self._view.ui.btn_new_create_project.clicked.connect(self.create_new_project)

        # </editor-fold>

        # <editor-fold desc="Open project page">
        self._view.ui.btn_open_open_project.clicked.connect(self.open_project)
        self._view.ui.list_open_projects.doubleClicked.connect(self.open_project)
        self._view.ui.txt_open_search.textChanged.connect(self.validate_open_search)
        self._view.ui.txt_open_selected_project.textChanged.connect(self.activate_open_button)
        self._view.ui.list_open_projects.currentItemChanged.connect(self.select_project_from_open_list)

        # </editor-fold>

        # <editor-fold desc="Delete project page">
        self._view.ui.btn_delete_delete_project.clicked.connect(self.delete_project)
        self._view.ui.txt_delete_search.textChanged.connect(self.validate_delete_search)
        self._view.ui.txt_delete_selected_projects.textChanged.connect(self.activate_delete_button)
        self._view.ui.list_delete_projects.currentItemChanged.connect(self.select_project_from_delete_list)

        # </editor-fold>

        # <editor-fold desc="Edit project page">
        self._view.ui.btn_edit_page.clicked.connect(self.display_edit_page)
        self._view.ui.list_edit_project_proteins.currentItemChanged.connect(self.check_for_cleaning)
        self._view.ui.btn_edit_clean_new_prot.clicked.connect(self.clean_protein_new)
        self._view.ui.btn_edit_clean_update_prot.clicked.connect(self.clean_protein_update)
        self._view.ui.btn_edit_project_delete.clicked.connect(self.delete_protein)
        self._view.ui.btn_edit_existing_protein_struct.clicked.connect(self.add_existing_protein)
        self._view.ui.btn_edit_project_save.clicked.connect(self.save_selected_protein_structure_as_pdb_file)
        self._view.ui.btn_edit_protein_rename.clicked.connect(self.rename_selected_protein_structure)
        # </editor-fold>

        # <editor-fold desc="View project page">
        self._view.ui.btn_view_project_show.clicked.connect(self.view_sequence)
        self._view.ui.btn_view_project_show_structure.clicked.connect(self.view_structure)
        self._view.ui.list_view_project_proteins.doubleClicked.connect(self.view_sequence)
        self._view.ui.list_view_project_proteins.itemClicked.connect(self.view_show_options)
        # </editor-fold>

        # <editor-fold desc="Use project page">
        self._view.ui.txt_use_project_name.textChanged.connect(self.validate_use_project_name)
        self._view.ui.btn_use_next.clicked.connect(self.show_protein_selection_for_use)
        self._view.ui.txt_use_search.textChanged.connect(self.validate_use_search)
        self._view.ui.btn_use_add_available_protein_structures.clicked.connect(
            self.add_protein_structure_to_new_project,
        )
        self._view.ui.list_use_available_protein_structures.doubleClicked.connect(
            self.add_protein_structure_to_new_project,
        )
        self._view.ui.list_use_available_protein_structures.itemClicked.connect(self.use_enable_add)
        self._view.ui.btn_use_remove_selected_protein_structures.clicked.connect(
            self.remove_protein_structure_to_new_project,
        )
        self._view.ui.list_use_selected_protein_structures.doubleClicked.connect(
            self.remove_protein_structure_to_new_project,
        )
        self._view.ui.list_use_selected_protein_structures.itemClicked.connect(self.use_enable_remove)
        self._view.ui.btn_use_back.clicked.connect(self.hide_protein_selection_for_use)
        self._view.ui.btn_use_create_new_project.clicked.connect(self.pre_create_use_project)

        # </editor-fold>

        # <editor-fold desc="ESMFold Monomer Prediction page">
        self._view.ui.btn_esm_seq_to_predict.clicked.connect(self.cloud_esm_add_seq_to_predict)
        self._view.ui.btn_esm_seq_to_predict_remove.clicked.connect(self.cloud_esm_remove)
        self._view.ui.btn_esm_next.clicked.connect(self.cloud_esm_next)
        self._view.ui.btn_esm_back.clicked.connect(self.cloud_esm_back)
        self._view.ui.btn_esm_next_2.clicked.connect(self.cloud_esm_add_protein)
        self._view.ui.btn_esm_back_2.clicked.connect(self.cloud_esm_back_2)
        self._view.ui.txt_esm_prot_name.textChanged.connect(self.cloud_esm_validate_protein_name)
        self._view.ui.txt_esm_prot_seq.textChanged.connect(self.cloud_esm_validate_protein_sequence)
        self._view.ui.btn_esm_predict.clicked.connect(self.predict_esm_monomer)

        self._view.ui.table_esm_prot_to_predict.itemSelectionChanged.connect(self.cloud_esm_item_changed)
        # </editor-fold>

        # <editor-fold desc="Monomer local prediction page">
        self._view.ui.btn_pred_mono_seq_to_predict.clicked.connect(self.local_pred_mono_add_seq_to_predict)
        self._view.ui.btn_pred_mono_seq_to_predict_remove.clicked.connect(self.local_pred_mono_remove)
        self._view.ui.btn_pred_mono_next.clicked.connect(self.local_pred_mono_next)
        self._view.ui.btn_pred_mono_back.clicked.connect(self.local_pred_mono_back)
        self._view.ui.btn_pred_mono_add_protein.clicked.connect(self.local_pred_mono_add_protein)
        self._view.ui.btn_pred_mono_back_2.clicked.connect(self.local_pred_mono_back_2)
        self._view.ui.txt_pred_mono_prot_name.textChanged.connect(self.local_pred_mono_validate_protein_name)
        self._view.ui.txt_pred_mono_seq_name.textChanged.connect(self.local_pred_mono_validate_protein_sequence)
        self._view.ui.btn_pred_mono_advanced_config.clicked.connect(self.show_prediction_configuration)
        self._view.ui.btn_pred_mono_predict.clicked.connect(self.predict_local_monomer)

        self._view.ui.table_pred_mono_prot_to_predict.itemSelectionChanged.connect(self.local_pred_mono_item_changed)
        # </editor-fold>

        # <editor-fold desc="Multimer prediction page">
        self._view.ui.btn_pred_multi_prot_to_predict_add.clicked.connect(self.local_pred_multi_add)
        self._view.ui.btn_pred_multi_prot_to_predict_remove.clicked.connect(self.local_pred_multi_remove)
        self._view.ui.btn_pred_multi_next.clicked.connect(self.local_pred_multi_next)
        self._view.ui.btn_pred_multi_back.clicked.connect(self.local_pred_multi_back)
        self._view.ui.btn_pred_multi_prot_to_predict_add_2.clicked.connect(self.local_pred_multi_prot_to_predict_add_2)
        self._view.ui.btn_pred_multi_back_2.clicked.connect(self.local_pred_multi_back_2)
        self._view.ui.txt_pred_multi_prot_name.textChanged.connect(self.local_pred_multi_validate_protein_name)
        self._view.ui.txt_pred_multi_prot_seq.textChanged.connect(self.local_pred_multi_validate_protein_sequence)
        self._view.ui.btn_pred_multi_prot_seq_add.clicked.connect(self.local_pred_multi_add_sequence_to_list)
        self._view.ui.btn_pred_multi_prot_seq_overview_remove.clicked.connect(
            self.local_pred_multi_remove_sequence_to_list,
        )
        self._view.ui.btn_pred_multi_advanced_config.clicked.connect(self.show_prediction_configuration)
        self._view.ui.btn_pred_multi_predict.clicked.connect(self.predict_local_multimer)

        self._view.ui.list_pred_multi_prot_seq_overview.itemClicked.connect(
            self.local_pred_multi_prot_seq_overview_item_changed,
        )
        self._view.ui.table_pred_multi_prot_to_predict.itemSelectionChanged.connect(
            self.local_pred_multi_prot_to_predict_item_changed,
        )
        # </editor-fold>

        # <editor-fold desc="Monomer Prediction + Analysis page">
        # <editor-fold desc="Prediction section">
        self._view.ui.btn_pred_analysis_mono_seq_to_predict.clicked.connect(self.mono_pred_analysis_add_seq_to_predict)
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.clicked.connect(
            self.mono_pred_analysis_remove_protein_to_predict,
        )
        self._view.ui.btn_pred_analysis_mono_next.clicked.connect(self.mono_pred_analysis_next)
        self._view.ui.btn_pred_analysis_mono_back.clicked.connect(self.mono_pred_analysis_back)
        self._view.ui.btn_pred_analysis_mono_add_protein.clicked.connect(self.mono_pred_analysis_add_protein)
        self._view.ui.btn_pred_analysis_mono_back_2.clicked.connect(self.mono_pred_analysis_back_2)
        self._view.ui.txt_pred_analysis_mono_prot_name.textChanged.connect(
            self.mono_pred_analysis_validate_protein_name,
        )
        self._view.ui.txt_pred_analysis_mono_seq_name.textChanged.connect(
            self.mono_pred_analysis_validate_protein_sequence,
        )
        self._view.ui.btn_pred_mono_advanced_config_2.clicked.connect(self.show_prediction_configuration)
        self._view.ui.btn_pred_analysis_mono_go_analysis_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self._view.ui.table_pred_analysis_mono_prot_to_predict.itemSelectionChanged.connect(
            self.mono_pred_analysis_prediction_overview_item_clicked,
        )

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.btn_pred_analysis_mono_add.clicked.connect(self.mono_pred_analysis_structure_analysis_add)
        self._view.ui.btn_pred_analysis_mono_remove.clicked.connect(self.remove_mono_pred_analysis_analysis_run)
        self._view.ui.btn_pred_analysis_mono_back_3.clicked.connect(self.mono_pred_analysis_structure_analysis_back_3)
        self._view.ui.btn_pred_analysis_mono_next_2.clicked.connect(self.mono_pred_analysis_structure_analysis_next_2)
        self._view.ui.btn_pred_analysis_mono_back_4.clicked.connect(self.mono_pred_analysis_structure_analysis_back_4)
        self._view.ui.btn_pred_analysis_mono_next_3.clicked.connect(self.mono_pred_analysis_structure_analysis_next_3)
        self._view.ui.btn_pred_analysis_mono_back_5.clicked.connect(self.mono_pred_analysis_structure_analysis_back_5)
        self._view.ui.btn_pred_analysis_mono_next_4.clicked.connect(self.mono_pred_analysis_structure_analysis_next_4)
        self._view.ui.box_pred_analysis_mono_prot_struct_1.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.box_pred_analysis_mono_prot_struct_2.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.list_pred_analysis_mono_ref_chains.itemSelectionChanged.connect(
            self.count_mono_pred_analysis_selected_chains_for_prot_struct_1,
        )
        self._view.ui.list_pred_analysis_mono_model_chains.itemSelectionChanged.connect(
            self.check_mono_pred_analysis_if_same_no_of_chains_selected,
        )
        self._view.ui.btn_pred_analysis_mono_back_pred_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self._view.ui.btn_pred_analysis_mono_start.clicked.connect(self.start_monomer_prediction_analysis)
        self._view.ui.list_pred_analysis_mono_overview.clicked.connect(
            self.mono_pred_analysis_structure_analysis_overview_clicked,
        )

        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Multimer Prediction + Analysis page">
        # <editor-fold desc="Prediction section">
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add.clicked.connect(self.multi_pred_analysis_add)
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.clicked.connect(
            self.multi_pred_analysis_remove_protein_to_predict,
        )
        self._view.ui.btn_pred_analysis_multi_next.clicked.connect(self.multi_pred_analysis_next)
        self._view.ui.btn_pred_analysis_multi_back.clicked.connect(self.multi_pred_analysis_back)
        self._view.ui.btn_pred_analysis_multi_prot_seq_add.clicked.connect(
            self.multi_pred_analysis_add_sequence_to_list,
        )
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.clicked.connect(
            self.multi_pred_analysis_remove_sequence_to_list,
        )
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.clicked.connect(
            self.multi_pred_analysis_add_protein_to_predict,
        )
        self._view.ui.btn_pred_analysis_multi_back_2.clicked.connect(self.multi_pred_analysis_back_2)

        self._view.ui.txt_pred_analysis_multi_prot_name.textChanged.connect(
            self.multi_pred_analysis_validate_protein_name,
        )
        self._view.ui.txt_pred_analysis_multi_prot_seq.textChanged.connect(
            self.multi_pred_analysis_validate_protein_sequence,
        )
        self._view.ui.btn_pred_analysis_multi_advanced_config.clicked.connect(self.show_prediction_configuration)
        self._view.ui.btn_pred_analysis_multi_go_analysis_setup.clicked.connect(self.switch_multimer_pred_analysis_tab)

        self._view.ui.list_pred_analysis_multi_prot_seq_overview.clicked.connect(
            self.multi_pred_analysis_prot_seq_overview_item_changed,
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.itemSelectionChanged.connect(
            self.multi_pred_analysis_prot_to_predict_item_changed,
        )

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.btn_pred_analysis_multi_add.clicked.connect(self.multi_pred_analysis_structure_analysis_add)
        self._view.ui.btn_pred_analysis_multi_remove.clicked.connect(self.remove_multi_pred_analysis_analysis_run)
        self._view.ui.btn_pred_analysis_multi_back_3.clicked.connect(self.multi_pred_analysis_structure_analysis_back_3)
        self._view.ui.btn_pred_analysis_multi_next_2.clicked.connect(self.multi_pred_analysis_structure_analysis_next_2)
        self._view.ui.btn_pred_analysis_multi_back_4.clicked.connect(self.multi_pred_analysis_structure_analysis_back_4)
        self._view.ui.btn_pred_analysis_multi_next_3.clicked.connect(self.multi_pred_analysis_structure_analysis_next_3)
        self._view.ui.btn_pred_analysis_multi_back_5.clicked.connect(self.multi_pred_analysis_structure_analysis_back_5)
        self._view.ui.btn_pred_analysis_multi_next_4.clicked.connect(self.multi_pred_analysis_structure_analysis_next_4)
        self._view.ui.box_pred_analysis_multi_prot_struct_1.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.box_pred_analysis_multi_prot_struct_2.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.list_pred_analysis_multi_ref_chains.itemSelectionChanged.connect(
            self.count_multi_pred_analysis_selected_chains_for_prot_struct_1,
        )
        self._view.ui.list_pred_analysis_multi_model_chains.itemSelectionChanged.connect(
            self.check_multi_pred_analysis_if_same_no_of_chains_selected,
        )
        self._view.ui.btn_pred_analysis_multi_back_pred_setup.clicked.connect(self.switch_multimer_pred_analysis_tab)
        self._view.ui.btn_pred_analysis_multi_start.clicked.connect(self.start_multimer_prediction_analysis)
        self._view.ui.list_pred_analysis_multi_overview.clicked.connect(
            self.multi_pred_analysis_structure_analysis_overview_clicked,
        )
        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Structure analysis page">
        self._view.ui.btn_analysis_batch_add.clicked.connect(self.structure_analysis_add)
        self._view.ui.btn_analysis_batch_remove.clicked.connect(self.remove_analysis_run)
        self._view.ui.btn_analysis_batch_back.clicked.connect(self.structure_analysis_back)
        self._view.ui.btn_analysis_batch_next.clicked.connect(self.structure_analysis_next)
        self._view.ui.btn_analysis_batch_back_2.clicked.connect(self.structure_analysis_back_2)
        self._view.ui.btn_analysis_batch_next_2.clicked.connect(self.structure_analysis_next_2)
        self._view.ui.btn_analysis_batch_back_3.clicked.connect(self.structure_analysis_back_3)
        self._view.ui.btn_analysis_batch_next_3.clicked.connect(self.structure_analysis_next_3)
        self._view.ui.box_analysis_batch_prot_struct_1.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch,
        )
        self._view.ui.box_analysis_batch_prot_struct_2.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch,
        )
        self._view.ui.list_analysis_batch_ref_chains.itemSelectionChanged.connect(
            self.count_batch_selected_chains_for_prot_struct_1,
        )
        self._view.ui.list_analysis_batch_model_chains.itemSelectionChanged.connect(
            self.check_if_same_no_of_chains_selected_batch,
        )
        self._view.ui.btn_analysis_batch_start.clicked.connect(self.start_process_batch)
        self._view.ui.list_analysis_batch_overview.clicked.connect(self.structure_analysis_overview_clicked)

        # </editor-fold>

        # <editor-fold desc="Analysis images page">
        self._view.ui.btn_add_analysis_images_struct_analysis.clicked.connect(
            self.add_protein_pair_to_image_creation_queue,
        )
        self._view.ui.list_analysis_images_struct_analysis.doubleClicked.connect(
            self.add_protein_pair_to_image_creation_queue,
        )
        self._view.ui.list_analysis_images_struct_analysis.clicked.connect(self.analysis_images_enable_add)
        self._view.ui.btn_remove_analysis_images_creation_struct_analysis.clicked.connect(
            self.remove_protein_pair_from_image_creation_queue,
        )
        self._view.ui.list_analysis_images_creation_struct_analysis.doubleClicked.connect(
            self.remove_protein_pair_from_image_creation_queue,
        )
        self._view.ui.list_analysis_images_creation_struct_analysis.clicked.connect(self.analysis_images_enable_remove)
        self._view.ui.btn_start_automatic_image_creation.clicked.connect(self.start_automatic_image_creation)

        # </editor-fold>

        # <editor-fold desc="Results page">
        self._view.ui.cb_results_analysis_options.currentIndexChanged.connect(self.load_results)
        self._view.ui.btn_color_rmsd.clicked.connect(self.color_protein_pair_by_rmsd)
        self._view.ui.btn_view_struct_alignment.clicked.connect(self.display_structure_alignment)
        self._view.ui.btn_view_distance_plot.clicked.connect(self.display_distance_plot)
        self._view.ui.btn_view_distance_histogram.clicked.connect(self.display_distance_histogram)
        self._view.ui.btn_view_interesting_region.clicked.connect(self.display_interesting_region)
        self._view.ui.btn_view_distance_table.clicked.connect(self.display_distance_table)

        # </editor-fold>

        # <editor-fold desc="Manage">
        self._view.ui.box_manage_choose_protein.activated.connect(self.choose_manage_open_protein)
        self._view.ui.box_manage_choose_color.activated.connect(self.choose_manage_color_selected_protein)
        self._view.ui.box_manage_choose_representation.activated.connect(self.choose_manage_representation)
        self._view.ui.box_manage_choose_bg_color.activated.connect(self.choose_manage_bg_color)
        self._view.ui.btn_disulfid_bond_show.clicked.connect(self.show_disulfid_bonds_as_sticks)
        self._view.ui.btn_disulfid_bond_hide.clicked.connect(self.hide_disulfid_bonds_as_sticks)

        # </editor-fold>

        # <editor-fold desc="Image page">
        self._view.ui.btn_update_scene.clicked.connect(self.update_scene)
        self._view.ui.btn_save_scene.clicked.connect(self.save_scene)
        self._view.ui.btn_save_image.clicked.connect(self.save_image)
        self._view.ui.btn_preview_image.clicked.connect(self.preview_image)
        self._view.ui.box_representation.activated.connect(self.show_representation)
        self._view.ui.box_bg_color.activated.connect(self.choose_bg_color)
        self._view.ui.box_renderer.activated.connect(self.choose_renderer)
        self._view.ui.box_ray_trace_mode.activated.connect(self.choose_ray_trace_mode)
        self._view.ui.box_ray_texture.activated.connect(self.choose_ray_texture)
        self._view.ui.cb_transparent_bg.stateChanged.connect(self.decide_transparent_bg)
        self._view.ui.cb_ray_tracing.stateChanged.connect(self.decide_ray_tracing)

        # </editor-fold>

        # <editor-fold desc="Hotspots page">
        self._view.ui.list_hotspots_choose_protein.currentItemChanged.connect(self.open_protein_for_hotspots)
        self._view.ui.btn_hotspots_resi_show.clicked.connect(self.show_resi_sticks)
        self._view.ui.btn_hotspots_resi_hide.clicked.connect(self.hide_resi_sticks)
        self._view.ui.btn_hotspots_resi_zoom.clicked.connect(self.zoom_resi_position)

        # </editor-fold>

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

    def quit_app(self) -> None:
        """Closes the entire plugin."""
        self._view.quit_app()

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

    # <editor-fold desc="Open extra views">
    def open_settings_global(self) -> None:
        """Opens the dialog for the global settings."""
        dialog = dialog_settings_global.DialogSettingsGlobal()
        dialog.exec_()
        self._application_settings = self._application_settings.deserialize_settings()
        globals.g_settings = self._application_settings
        self._workspace_path = globals.g_settings.workspace_path
        self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")
        self._setup_statusbar()

    def open_logs(self) -> None:
        """Opens a file explorer with all log files and can open a log file in the default application."""
        file_dialog = QtWidgets.QFileDialog()
        log_path = str(constants.LOG_PATH)
        file_dialog.setDirectory(log_path)
        file_path, _ = file_dialog.getOpenFileName(self._view, "Select a log file to open", "", "LOG File (*.log)")
        if file_path:
            os.startfile(file_path)

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

    def open_page_information(self) -> None:
        """Opens the message box, to display extra information based on the page."""
        with open(
            f"{constants.PAGE_HELP_PATHS_DICT[self._view.ui.lbl_page_title.text()]}",
            "r",
            encoding="utf-8",
        ) as file:
            html_content = file.read()
            file.close()
        tmp_dialog = dialog_help.DialogHelp(html_content)
        tmp_dialog.exec_()

    @staticmethod
    def open_release_notes_in_standard_application() -> None:
        """Opens the release notes in the default app."""
        os.startfile(constants.CHANGELOG_HTML_PATH)

    # </editor-fold>

    # def open_change_log(self) -> None:
    #     """Opens change log based on the last run pyssa version.
    #
    #     Notes:
    #         Change log opens only if the last pyssa version is older than the current one.
    #     """
    #     last_version = f"v{self.get_version_from_latest_log_file(self.get_filepath_of_latest_log_file())}"
    #     if last_version != constants.VERSION_NUMBER:
    #         self.open_release_notes_in_standard_application()
    #
    # def get_version_from_latest_log_file(self, a_filepath: pathlib.Path) -> str:
    #     """Gets the pyssa version of the latest log file.
    #
    #     Args:
    #         a_filepath: the filepath to the latest log file.
    #     """
    #     with open(str(a_filepath), "r", encoding="utf-8") as file:
    #         file_content = file.read()
    #         file.close()
    #     # Define the regex pattern to extract the version number
    #     pattern = r"PySSA started with version v(\d+\.\d+\.\d+)"
    #     # Search for the pattern in the text
    #     match = re.search(pattern, file_content)
    #     # Check if a match is found
    #     if match:
    #         version_number = match.group(1)
    #         print("Extracted version number:", version_number)
    #     else:
    #         print("Version number not found in the text.")
    #         version_number = None
    #     return version_number
    #
    # def get_filepath_of_latest_log_file(self) -> pathlib.Path:
    #     """Gets the filepath of the latest log file.
    #
    #     Raises:
    #         FileNotFoundError: If no log files can be found in the log directory.
    #     """
    #     # Get a list of files in the directory
    #     files = [f for f in os.listdir(constants.LOG_PATH) if os.path.isfile(os.path.join(constants.LOG_PATH, f))]
    #     # Filter files to include only log files (adjust the extension accordingly)
    #     log_files = [f for f in files if f.endswith(".log")]
    #     # Check if there are any log files
    #     if not log_files:
    #         raise FileNotFoundError("No log files found in the directory.")
    #     # Get the full path of each log file
    #     log_files_paths = [os.path.join(constants.LOG_PATH, f) for f in log_files]
    #     # Get the latest log file based on modification time
    #     latest_log_file_index = log_files_paths.index(max(log_files_paths, key=os.path.getmtime)) - 1
    #     return pathlib.Path(log_files_paths[latest_log_file_index])

    # <editor-fold desc="Logic related methods">

    # <editor-fold desc="New project page functions">
    def create_new_project(self) -> None:
        """Creates a new project with the content of the new page."""
        # <editor-fold desc="Checks">
        if self._application_settings.wsl_install == 0:
            basic_boxes.ok(
                "Create new project",
                "Please install local colabfold to create a project!",
                QtWidgets.QMessageBox.Warning,
            )
            return
        if self._application_settings.local_colabfold == 0:
            basic_boxes.ok(
                "Create new project",
                "Please install local colabfold to create a project!",
                QtWidgets.QMessageBox.Warning,
            )
            return

        # </editor-fold>

        self._view.wait_spinner.start()
        an_add_protein_flag: bool = False
        if self._view.ui.cb_new_add_reference.checkState() == 2:
            an_add_protein_flag = True

        self._active_task = tasks.Task(
            target=main_presenter_async.create_new_project,
            args=(
                self._view.ui.txt_new_project_name.text(),
                self._workspace_path,
                an_add_protein_flag,
                self._view.ui.txt_new_choose_reference.text(),
            ),
            post_func=self.__await_create_new_project,
        )
        self._active_task.start()
        self.update_status("Creating new project ...")

    def __await_create_new_project(self, a_result: tuple) -> None:
        self._current_project: "project.Project" = a_result[1]
        constants.PYSSA_LOGGER.info(f"Created the project {self._current_project.get_project_name()}.")
        self._view.ui.cb_new_add_reference.setCheckState(0)
        self._project_watcher.current_project = self._current_project
        constants.PYSSA_LOGGER.info(
            f"{self._project_watcher.current_project.get_project_name()} is the current project.",
        )
        self._project_watcher.on_home_page = False
        # update gui
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.ui.lbl_current_project_name.setText(self._current_project.get_project_name())
        self.main_window_state = main_window_state.MainWindowState(
            results_state.ResultsState(),
            image_state.ImageState(),
        )
        self.display_view_page()
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)

    # </editor-fold>

    # <editor-fold desc="Open project page functions">
    def open_project(self) -> None:
        """Initiates the task to open an existing project."""
        self._view.wait_spinner.start()
        self._active_task = tasks.Task(
            target=main_presenter_async.open_project,
            args=(
                self._workspace_path,
                self._view.ui.txt_open_selected_project.text(),
                self._application_settings,
            ),
            post_func=self.__await_open_project,
        )
        self._active_task.start()
        self.update_status("Opening existing project ...")

    def __await_open_project(self, a_result: tuple) -> None:
        self._current_project = a_result[1]
        self._project_watcher.current_project = self._current_project
        constants.PYSSA_LOGGER.info(
            f"{self._project_watcher.current_project.get_project_name()} is the current project.",
        )
        self._view.ui.lbl_current_project_name.setText(self._current_project.get_project_name())
        self._project_watcher.on_home_page = False
        self._project_watcher.show_valid_options(self._view.ui)
        cmd.reinitialize()
        self._view.ui.btn_manage_session.hide()
        self.display_view_page()
        self.main_window_state = main_window_state.MainWindowState(
            results_state.ResultsState(),
            image_state.ImageState(),
        )
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)

    # </editor-fold>

    # <editor-fold desc="Delete project page functions">
    def delete_project(self) -> None:
        """Deletes an existing project."""
        # popup message which warns the user that the selected project gets deleted
        response: bool = gui_utils.warning_message_project_gets_deleted()
        tmp_project_name = self._view.ui.txt_delete_selected_projects.text()
        if response is True:
            os.remove(pathlib.Path(f"{self._workspace_path}/{self._view.ui.txt_delete_selected_projects.text()}"))
            if self._view.ui.txt_delete_selected_projects.text() == self._view.ui.lbl_current_project_name.text():
                self._view.ui.lbl_current_project_name.clear()
            self._view.ui.txt_delete_selected_projects.clear()
            # update list
            self._view.ui.list_delete_projects.clear()
            # pre-process
            self._view.status_bar.showMessage(self._workspace_label.text())
            self._view.ui.list_delete_projects.clear()
            self._view.status_bar.showMessage(self._workspace_label.text())
            tools.scan_workspace_for_valid_projects(self._workspace_path, self._view.ui.list_delete_projects)
            constants.PYSSA_LOGGER.info(f"The project {tmp_project_name} was successfully deleted.")
            if self._view.ui.list_delete_projects.count() == 0:
                self.display_home_page()
                self._project_watcher.check_workspace_for_projects(self._workspace_path, self._view.ui)
        else:
            constants.PYSSA_LOGGER.info("No project has been deleted. No changes were made.")

    # </editor-fold>

    # <editor-fold desc="Save project functions">
    def save_project(self) -> None:
        """Saves the project.xml."""
        self._view.wait_spinner.start()
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_save_project,
        )
        tools.ask_to_save_pymol_session(self._current_project, self.current_session, self._application_settings)
        self._active_task = tasks.Task(
            target=main_presenter_async.save_project,
            args=(self._current_project, 0),
            post_func=self.__await_save_project,
        )
        self._active_task.start()
        self.update_status("Saving current project ...")

    def __await_save_project(self, result: tuple) -> None:
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)
        basic_boxes.ok("Save Project", "The project was successfully saved.", QtWidgets.QMessageBox.Information)

    # </editor-fold>

    # <editor-fold desc="Edit project page functions">
    def check_for_cleaning(self) -> None:
        """Checks if the selected protein can be cleaned."""
        self._view.wait_spinner.start()
        tools.ask_to_save_pymol_session(self._current_project, self.current_session, self._application_settings)
        try:
            tmp_protein_name: str = self._view.ui.list_edit_project_proteins.currentItem().text()
        except AttributeError:
            self._view.wait_spinner.stop()
            return
        self._active_task = tasks.Task(
            target=main_presenter_async.check_for_cleaning,
            args=(
                tmp_protein_name,
                self._current_project,
            ),
            post_func=self.__await_check_for_cleaning,
        )
        self._active_task.start()
        self.update_status("Checking protein properties ...")

    def __await_check_for_cleaning(self, result: tuple) -> None:
        # check if selected protein contains any organic or solvent molecules which can be removed
        tmp_is_cleanable: bool = result[1]
        tmp_is_in_protein_pair: bool = result[2]
        if tmp_is_cleanable:
            gui_elements_to_show = [
                self._view.ui.lbl_edit_clean_new_prot,
                self._view.ui.btn_edit_clean_new_prot,
                self._view.ui.lbl_edit_clean_update_prot,
                self._view.ui.btn_edit_clean_update_prot,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
        else:
            gui_elements_to_hide = [
                self._view.ui.lbl_edit_clean_new_prot,
                self._view.ui.btn_edit_clean_new_prot,
                self._view.ui.lbl_edit_clean_update_prot,
                self._view.ui.btn_edit_clean_update_prot,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        # check if selected protein is in any existing protein pair
        if tmp_is_in_protein_pair:
            gui_elements_to_hide = [
                self._view.ui.label_12,
                self._view.ui.btn_edit_project_delete,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            gui_elements_to_show = [
                self._view.ui.label_12,
                self._view.ui.btn_edit_project_delete,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
        self._view.ui.btn_edit_project_save.show()
        self._view.ui.label_13.show()
        self._view.ui.label_15.show()
        self._view.ui.btn_edit_protein_rename.show()
        self._view.ui.btn_manage_session.show()
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)

    def clean_protein_new(self) -> None:
        """Cleans the selected protein structure and creates a new cleaned structure."""
        self._view.wait_spinner.start()
        self._active_task = tasks.Task(
            target=main_presenter_async.clean_protein_new,
            args=(
                self._view.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""),
                self._current_project,
            ),
            post_func=self.__await_clean_protein_new,
        )
        self._active_task.start()
        self.update_status("Duplicating and cleaning protein ...")

    def __await_clean_protein_new(self, result: tuple) -> None:
        self._init_edit_page()
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)

    def clean_protein_update(self) -> None:
        """Cleans the selected protein structure."""
        self._view.wait_spinner.start()
        if basic_boxes.yes_or_no(
            "Clean protein",
            "Are you sure you want to clean this protein?\n" "This will remove all organic and solvent components!",
            QtWidgets.QMessageBox.Information,
        ):
            self._active_task = tasks.Task(
                target=main_presenter_async.clean_protein_update,
                args=(
                    self._view.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""),
                    self._current_project,
                ),
                post_func=self.__await_clean_protein_update,
            )
            self._active_task.start()
            self.update_status("Cleaning protein ...")
        else:
            constants.PYSSA_LOGGER.info("No protein has been cleaned.")
            self._view.wait_spinner.stop()

    def __await_clean_protein_update(self) -> None:
        self._init_edit_page()
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)

    def delete_protein(self) -> None:
        """Deletes the selected protein structure."""
        self._view.wait_spinner.start()
        if gui_utils.warning_message_protein_gets_deleted():
            self._active_task = tasks.Task(
                target=main_presenter_async.delete_protein,
                args=(
                    self._view.ui.list_edit_project_proteins.currentItem().text(),
                    self._current_project,
                ),
                post_func=self.__await_delete_protein,
            )
            self._active_task.start()
            self.update_status("Deleting protein ...")
        else:
            constants.PYSSA_LOGGER.info("No protein was deleted.")
            self._view.wait_spinner.stop()

    def __await_delete_protein(self) -> None:
        self._init_edit_page()
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)

    def add_existing_protein(self) -> None:
        """Opens a dialog to adds an existing protein structure to the project."""
        self.tmp_dialog = view_add_protein.AddProteinView()
        self.tmp_dialog.return_value.connect(self.post_add_existing_protein)
        self.tmp_dialog.show()

    def post_add_existing_protein(self, return_value: tuple) -> None:
        """Adds an existing protein structure to the current project either by id or from the filesystem.

        Args:
            return_value: a tuple consisting of the filepath or PDB id and the length of the first one.
        """
        self._view.wait_spinner.start()
        if return_value[1] > 0:
            self._active_task = tasks.Task(
                target=main_presenter_async.add_existing_protein_to_project,
                args=(return_value, self._current_project),
                post_func=self.__await_post_add_existing_protein,
            )
            self._active_task.start()
            self.update_status("Adding protein to current project ...")
        else:
            self._view.wait_spinner.stop()
            self.display_edit_page()

    def __await_post_add_existing_protein(self, result: tuple) -> None:
        self._current_project = result[1]
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)
        self.display_edit_page()

    def save_selected_protein_structure_as_pdb_file(self) -> None:
        """Saves selected protein as pdb file."""
        self._view.wait_spinner.start()
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(
            self._view,
            "Save protein structure",
            "",
            "Protein Data Bank File (*.pdb)",
        )
        if file_path:
            self._active_task = tasks.Task(
                target=main_presenter_async.save_selected_protein_structure_as_pdb_file,
                args=(self._view.ui.list_edit_project_proteins.currentItem().text(), self._current_project, file_path),
                post_func=self.__await_save_selected_protein_structure_as_pdb_file,
            )
            self._active_task.start()
        else:
            self._view.wait_spinner.stop()

    def __await_save_selected_protein_structure_as_pdb_file(self, result: tuple) -> None:
        self._view.wait_spinner.stop()
        if result[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            basic_boxes.ok(
                "Save protein structure",
                "Saving the protein as .pdb file failed!",
                QtWidgets.QMessageBox.Error,
            )
        elif result[0] == exit_codes.EXIT_CODE_ZERO[0]:
            basic_boxes.ok(
                "Save protein structure",
                "The protein was successfully saved as .pdb file.",
                QtWidgets.QMessageBox.Information,
            )
        else:
            basic_boxes.ok(
                "Save protein structure",
                "Saving the protein as .pdb file failed with an unexpected error!",
                QtWidgets.QMessageBox.Error,
            )
        self._project_watcher.show_valid_options(self._view.ui)
        self.update_status(self._workspace_status)

    def rename_selected_protein_structure(self) -> None:
        """Opens a new view to rename the selected protein."""
        self._view.wait_spinner.start()
        self.tmp_dialog = dialog_rename_protein.DialogRenameProtein(self._workspace_path)
        self.tmp_dialog.return_value.connect(self.post_rename_selected_protein_structure)
        self.tmp_dialog.show()

    def post_rename_selected_protein_structure(self, return_value: tuple) -> None:
        """Renames a selected protein structure."""
        if return_value[1] is True:
            self._active_task = tasks.Task(
                target=main_presenter_async.rename_selected_protein_structure,
                args=(
                    self._view.ui.list_edit_project_proteins.currentItem().text(),
                    return_value[0],
                    self._current_project,
                ),
                post_func=self.__await_post_rename_selected_protein_structure,
            )
            self._active_task.start()
        else:
            self._view.wait_spinner.stop()

    def __await_post_rename_selected_protein_structure(self, result: tuple) -> None:
        self._init_edit_page()
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.wait_spinner.stop()

    # </editor-fold>

    # <editor-fold desc="View project page functions">
    def view_structure(self) -> None:
        """Displays the structure of the selected protein in pymol."""
        protein_name = self._view.ui.list_view_project_proteins.currentItem().text()
        tools.ask_to_save_pymol_session(self._current_project, self.current_session, self._application_settings)
        cmd.reinitialize()
        self._view.ui.btn_manage_session.show()
        try:
            self._current_project.search_protein(protein_name).load_protein_pymol_session()
            constants.PYSSA_LOGGER.info("Loaded PyMOL session of protein %s", protein_name)
        except pymol.CmdException:
            constants.PYSSA_LOGGER.error("Error while loading protein in PyMOL!")
            return
        self.current_session = current_session.CurrentSession(
            "protein",
            protein_name,
            self._current_project.search_protein(protein_name).pymol_session,
        )
        print(self.current_session)

    # </editor-fold>

    # <editor-fold desc="Use project page functions">
    def display_use_page(self) -> None:
        """Displays the use project page."""
        self._view.wait_spinner.start()
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False

        self.start_worker_thread(
            task_workers.LoadUsePageWorker(
                self._workspace_path,
                self._current_project.convert_list_of_proteins_to_list_of_protein_infos(),
            ),
            self.post_display_use_page,
        )
        self._init_use_page()

    def post_display_use_page(self, return_value) -> None:  # noqa: ANN001
        """Displays the use project page, after cpu intense task (post thread method).

        Args:
            return_value: the value which gets returned from the thread process
        """
        # this for-loop is necessary for eliminating all proteins which are in the current project from the ones which
        # are available
        for tmp_item in return_value[0]:
            self._view.ui.list_use_available_protein_structures.addItem(tmp_item)
        for i in range(self._view.ui.list_use_selected_protein_structures.count()):
            self._view.ui.list_use_selected_protein_structures.setCurrentRow(i)
            tmp_prot_name = self._view.ui.list_use_selected_protein_structures.currentItem().text()
            if tmp_prot_name in return_value[1]:
                return_value[1].remove(tmp_prot_name)

        for tmp_item in return_value[1]:
            self._view.ui.list_use_available_protein_structures.addItem(tmp_item)
        for tmp_project_name in return_value[2]:
            self._view.ui.list_use_existing_projects.addItem(QtWidgets.QListWidgetItem(tmp_project_name))

        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 14, "Use existing project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self._view.ui.btn_use_page)
        self._view.wait_spinner.stop()

    def pre_create_use_project(self) -> None:
        """Sets up the worker for the create_use_project task."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        # copy proteins in new project
        proteins_to_copy = []
        for i in range(self._view.ui.list_use_selected_protein_structures.count()):
            self._view.ui.list_use_selected_protein_structures.setCurrentRow(i)
            proteins_to_copy.append(self._view.ui.list_use_selected_protein_structures.currentItem().text())

        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.CreateUseProjectWorker(self._workspace_path, proteins_to_copy)
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.create_use_project)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        self._view.ui.lbl_current_project_name.setText(self._view.ui.txt_use_project_name.text())
        self._view.status_bar.showMessage(f"Creating new project: {self._view.ui.txt_use_project_name.text()} ...")
        # save project folder in current workspace
        new_project = project.Project(self._view.ui.txt_use_project_name.text(), self._workspace_path)
        # new_project.create_project_tree()
        self._current_project = new_project

    def create_use_project(self, proteins_for_new_project) -> None:  # noqa: ANN001
        """Post thread method.

        Args:
            proteins_for_new_project (list): the proteins which are in the new project
        """
        for tmp_protein_obj in proteins_for_new_project:
            self._current_project.add_existing_protein(tmp_protein_obj)
        self._current_project.serialize_project(
            pathlib.Path(f"{self._workspace_path}/{self._current_project.get_project_name()}.xml"),
        )
        # shows options which can be done with the data in the project folder
        self._project_watcher.current_project = self._current_project
        self._project_watcher.on_home_page = False
        self._project_watcher.show_valid_options(self._view.ui)
        self.project_scanner.project = self._current_project
        self._init_use_page()
        constants.PYSSA_LOGGER.info(
            f"The project {self._current_project.get_project_name()} was successfully created through a use.",
        )
        self.display_view_page()
        QtWidgets.QApplication.restoreOverrideCursor()

    # </editor-fold>

    # <editor-fold desc="Import, Export functions">
    def import_project(self) -> None:
        """Imports a project.xml into the current workspace."""
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_import_project,
        )
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
            tmp_project = project.Project("", self._workspace_path)
            tmp_project = tmp_project.deserialize_project(pathlib.Path(file_path), self._application_settings)
            tmp_project.set_workspace_path(self._workspace_path)
            if len(tmp_project.proteins) <= 1:
                if self._application_settings.wsl_install == 0:
                    basic_boxes.ok(
                        "Create new project",
                        "Please install local colabfold to import this project!",
                        QtWidgets.QMessageBox.Warning,
                    )
                    return
                elif self._application_settings.local_colabfold == 0:  # noqa: RET505
                    basic_boxes.ok(
                        "Create new project",
                        "Please install local colabfold to import this project!",
                        QtWidgets.QMessageBox.Warning,
                    )
                    return
            new_filepath = pathlib.Path(f"{self._workspace_path}/{tmp_project.get_project_name()}.xml")
            tmp_project.serialize_project(new_filepath)
            self._current_project = self._current_project.deserialize_project(new_filepath, self._application_settings)
            constants.PYSSA_LOGGER.info(f"Opening the project {self._current_project.get_project_name()}.")
            self._project_watcher.current_project = self._current_project
            self.project_scanner.project = self._current_project
            constants.PYSSA_LOGGER.info(
                f"{self._project_watcher.current_project.get_project_name()} is the current project.",
            )
            self._view.ui.lbl_current_project_name.setText(self._current_project.get_project_name())
            self._project_watcher.on_home_page = False
            self._project_watcher.show_valid_options(self._view.ui)
            self._view.ui.btn_manage_session.show()
            self.display_view_page()
            basic_boxes.ok(
                "Import Project",
                "The project was successfully imported.",
                QtWidgets.QMessageBox.Information,
            )

    def export_current_project(self) -> None:
        """Exports the current project to an importable format."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_export_project,
        )
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(self._view, "Save current project", "", "XML Files (*.xml)")
        if file_path:
            self._current_project.serialize_project(pathlib.Path(file_path))
            basic_boxes.ok(
                "Export Project",
                "The project was successfully exported.",
                QtWidgets.QMessageBox.Information,
            )

    # </editor-fold>

    # <editor-fold desc="Close project functions">
    def close_project(self) -> None:
        """Closes the current project."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        tools.ask_to_save_pymol_session(self._current_project, self.current_session, self._application_settings)
        cmd.reinitialize()
        self._view.ui.list_hotspots_choose_protein.clear()
        self._project_watcher.on_home_page = True
        self._project_watcher.current_project = project.Project("", pathlib.Path(""))
        self._project_watcher.show_valid_options(self._view.ui)
        self._view.ui.lbl_current_project_name.setText("")
        self._init_all_pages()
        self.results_name = ""
        constants.PYSSA_LOGGER.info(f"The project {self._current_project.get_project_name()} was closed")
        self.display_home_page()

    # </editor-fold>

    # <editor-fold desc="ESMFold Monomer functions">
    def post_predict_esm_monomer(self, output) -> None:  # noqa: ANN001
        """Post thread method, for the prediction process.

        Args:
            output (list): the proteins which got predicted
        """
        self.block_box_prediction.destroy(True)
        for tmp_filename in os.listdir(constants.ESMFOLD_PDB_DIR):
            constants.PYSSA_LOGGER.info(
                f"Add protein {tmp_filename} to the current project {self._current_project.get_project_name()}",
            )
            self._current_project.add_existing_protein(
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
        self._current_project.serialize_project(self._current_project.get_project_xml_path())
        self.display_view_page()
        self._project_watcher.show_valid_options(self._view.ui)

    def predict_esm_monomer(self) -> None:
        """Sets up the worker to predict the proteins with the ESM-Fold."""
        # <editor-fold desc="Worker setup">
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.EsmFoldWorker(self._view.ui.table_esm_prot_to_predict)
        self.tmp_thread = task_workers.setup_worker_for_work(
            self.tmp_thread,
            self.tmp_worker,
            self.post_predict_esm_monomer,
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
    def post_prediction_process(self, an_exit_code: int, an_exit_code_description: str) -> None:
        """Process which runs after each prediction job."""
        if an_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because there was an error writing the fasta file(s)!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because the fasta file(s) could not be found!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because a subprocess failed!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because of an unknown error!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(f"Prediction ended with exit code {an_exit_code}: {an_exit_code_description}")
        elif an_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            if self.prediction_type == constants.PREDICTION_TYPE_PRED:
                self._current_project.serialize_project(self._current_project.get_project_xml_path())
                constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
                self.block_box_prediction.destroy(True)
                basic_boxes.ok(
                    "Structure prediction",
                    "All structure predictions are done. Go to View to check the new proteins.",
                    QtWidgets.QMessageBox.Information,
                )
                constants.PYSSA_LOGGER.info("All structure predictions are done.")
                self._project_watcher.show_valid_options(self._view.ui)
                self._init_local_pred_mono_page()
                self._init_local_pred_multi_page()
                self.display_view_page()
            elif self.prediction_type == constants.PREDICTION_TYPE_PRED_MONO_ANALYSIS:
                # executes if monomers were successfully predicted
                self._current_project.serialize_project(self._current_project.get_project_xml_path())
                constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
                constants.PYSSA_LOGGER.info("All structure predictions are done.")
                constants.PYSSA_LOGGER.info("Begin analysis process.")
                constants.PYSSA_LOGGER.debug(
                    f"Thread count before analysis worker: {self.threadpool.activeThreadCount()}",
                )

                # self.worker_analysis = workers.AnalysisWorkerPool(
                #    self._view.ui.list_pred_analysis_mono_overview, self._view.ui.cb_pred_analysis_mono_images,
                #    self._view.status_bar, self._current_project, self._application_settings, self._init_mono_pred_analysis_page)
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
                    self._view.ui.list_pred_analysis_mono_overview,
                    self._view.ui.cb_pred_analysis_mono_images,
                    self._view.status_bar,
                    self._current_project,
                    self._application_settings,
                    self._init_mono_pred_analysis_page,
                )
                self.tmp_thread = task_workers.setup_worker_for_work(
                    self.tmp_thread,
                    self.tmp_worker,
                    self.display_view_page,
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
                self._project_watcher.show_valid_options(self._view.ui)
            elif self.prediction_type == constants.PREDICTION_TYPE_PRED_MULTI_ANALYSIS:
                self._current_project.serialize_project(self._current_project.get_project_xml_path())
                constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
                constants.PYSSA_LOGGER.info("All structure predictions are done.")
                constants.PYSSA_LOGGER.info("Begin analysis process.")
                constants.PYSSA_LOGGER.debug(
                    f"Thread count before analysis worker: {self.threadpool.activeThreadCount()}",
                )

                # self.worker_analysis = workers.AnalysisWorkerPool(
                #    self._view.ui.list_pred_analysis_multi_overview, self._view.ui.cb_pred_analysis_multi_images,
                #    self._view.status_bar, self._current_project, self._application_settings, self._init_multi_pred_analysis_page)
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
                    self._view.ui.list_pred_analysis_multi_overview,
                    self._view.ui.cb_pred_analysis_multi_images,
                    self._view.status_bar,
                    self._current_project,
                    self._application_settings,
                    self._init_multi_pred_analysis_page,
                )
                self.tmp_thread = task_workers.setup_worker_for_work(
                    self.tmp_thread,
                    self.tmp_worker,
                    self.display_view_page,
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
                self._project_watcher.show_valid_options(self._view.ui)
        try:
            shutil.rmtree(pathlib.Path(f"{constants.SCRATCH_DIR}/local_predictions"))
        except Exception as e:
            constants.PYSSA_LOGGER.warning(f"Local predictions scratch path could not be deleted. {e}")

    def predict_local_monomer(self) -> None:
        """Sets tup the worker for the prediction with the colabfold."""
        self._view.wait_spinner.start()
        self.prediction_type = constants.PREDICTION_TYPE_PRED
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        self.update_status("Begin prediction process ...")
        predictions: list[
            prediction_protein_info.PredictionProteinInfo
        ] = prediction_util.get_prediction_name_and_seq_from_table(self._view.ui.table_pred_mono_prot_to_predict)

        self._active_task = tasks.Task(
            target=main_presenter_async.predict_protein_with_colabfold,
            args=(
                predictions,
                self.prediction_configuration,
                self._current_project,
            ),
            post_func=self.__await_predict_protein_with_colabfold,
        )
        self._active_task.start()

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
            return
        else:  # noqa: RET505
            self.block_box_prediction.close()
            self._view.wait_spinner.stop()

    def __await_predict_protein_with_colabfold(self, result: tuple) -> None:
        """Process which runs after each prediction job."""
        tmp_exit_code = result[0]
        tmp_exit_code_description = [1]
        if tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because there was an error writing the fasta file(s)!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
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
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
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
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
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
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            if self.prediction_type == constants.PREDICTION_TYPE_PRED:
                self._current_project.serialize_project(self._current_project.get_project_xml_path())
                constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
                self.block_box_prediction.destroy(True)
                basic_boxes.ok(
                    "Structure prediction",
                    "All structure predictions are done. Go to View to check the new proteins.",
                    QtWidgets.QMessageBox.Information,
                )
                constants.PYSSA_LOGGER.info("All structure predictions are done.")
                self._project_watcher.show_valid_options(self._view.ui)
                self._init_local_pred_mono_page()
                self._init_local_pred_multi_page()
                self.display_view_page()
        self._view.wait_spinner.stop()

    def abort_prediction(self) -> None:
        """Aborts the running prediction."""
        constants.PYSSA_LOGGER.info("Structure prediction process was aborted manually.")
        subprocess.run(["wsl", "--shutdown"])
        constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
        filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
        constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
        basic_boxes.ok("Abort prediction", "The structure prediction was aborted.", QtWidgets.QMessageBox.Information)
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_prediction_abort,
        )
        self._project_watcher.show_valid_options(self._view.ui)

    # </editor-fold>

    # <editor-fold desc="Multimer Local Prediction functions">
    def predict_local_multimer(self) -> None:
        """Sets up the worker for the prediction with the colabfold."""
        self._view.wait_spinner.start()
        self.prediction_type = constants.PREDICTION_TYPE_PRED
        constants.PYSSA_LOGGER.info("Begin multimer prediction process.")
        self.update_status("Begin prediction process ...")
        predictions: list[
            prediction_protein_info.PredictionProteinInfo
        ] = prediction_util.get_prediction_name_and_seq_from_table(self._view.ui.table_pred_multi_prot_to_predict)

        self._active_task = tasks.Task(
            target=main_presenter_async.predict_protein_with_colabfold,
            args=(
                predictions,
                self.prediction_configuration,
                self._current_project,
            ),
            post_func=self.__await_predict_protein_with_colabfold,
        )
        self._active_task.start()

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

    # </editor-fold>

    # <editor-fold desc="Monomer Prediction + Analysis functions">
    def start_monomer_prediction_analysis(self) -> None:
        """Sets up the worker for the prediction of the proteins."""
        self._view.wait_spinner.start()
        self.prediction_type = constants.PREDICTION_TYPE_PRED_MONO_ANALYSIS
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        self.update_status("Begin prediction process ...")
        predictions: list[
            prediction_protein_info.PredictionProteinInfo
        ] = prediction_util.get_prediction_name_and_seq_from_table(
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
        )

        self._active_task = tasks.Task(
            target=main_presenter_async.predict_protein_with_colabfold,
            args=(
                predictions,
                self.prediction_configuration,
                self._current_project,
            ),
            post_func=self.__await_monomer_prediction_for_subsequent_analysis,
        )
        self._active_task.start()

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

    def __await_monomer_prediction_for_subsequent_analysis(self, result: tuple) -> None:
        tmp_exit_code = result[0]
        tmp_exit_code_description = [1]
        if tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            self.block_box_prediction.destroy(True)
            constants.PYSSA_LOGGER.info("All structure predictions are done.")
            self.update_status("All structure predictions are done.")
            constants.PYSSA_LOGGER.info("Begin analysis process.")
            self.update_status("Begin analysis process ...")

            tmp_raw_analysis_run_names: list = []
            for row_no in range(self._view.ui.list_pred_analysis_mono_overview.count()):
                tmp_raw_analysis_run_names.append(self._view.ui.list_pred_analysis_mono_overview.item(row_no).text())

            self._active_task = tasks.Task(
                target=main_presenter_async.run_distance_analysis,
                args=(
                    tmp_raw_analysis_run_names,
                    self._current_project,
                    self._application_settings,
                    self._view.ui.cb_pred_analysis_mono_images.isChecked(),
                ),
                post_func=self.post_analysis_process,
            )
            self._active_task.start()

            if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
                os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
            self.block_box_analysis.exec_()
            self.display_view_page()

        elif tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because there was an error writing the fasta file(s)!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because the fasta file(s) could not be found!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because a subprocess failed!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
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
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()

    def mono_pred_analysis_structure_analysis_next_2(self) -> None:
        """Shows the gui elements for the chain selection of protein 1."""
        self._view.wait_spinner.start()
        tmp_proteins_to_predict: list[str] = []
        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            tmp_proteins_to_predict.append(
                self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text(),
            )
        self._active_task = tasks.Task(
            target=main_presenter_async.check_chains_for_subsequent_analysis,
            args=(
                self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText(),
                self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText(),
                self._current_project,
                tmp_proteins_to_predict,
            ),
            post_func=self.__await_mono_pred_analysis_structure_analysis_next_2,
        )
        self._active_task.start()

    def __await_mono_pred_analysis_structure_analysis_next_2(self, result: tuple) -> None:
        _, tmp_analysis_name = result
        if tmp_analysis_name != "":
            gui_elements_to_show = [
                self._view.ui.btn_pred_analysis_mono_remove,
                self._view.ui.btn_pred_analysis_mono_add,
                self._view.ui.lbl_pred_analysis_mono_overview,
                self._view.ui.list_pred_analysis_mono_overview,
                self._view.ui.lbl_pred_analysis_mono_images,
                self._view.ui.cb_pred_analysis_mono_images,
                self._view.ui.btn_pred_analysis_mono_start,
                self._view.ui.btn_pred_analysis_mono_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.box_pred_analysis_mono_prot_struct_1,
                self._view.ui.box_pred_analysis_mono_prot_struct_2,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_2,
                self._view.ui.lbl_pred_analysis_mono_ref_chains,
                self._view.ui.list_pred_analysis_mono_ref_chains,
                self._view.ui.lbl_pred_analysis_mono_model_chains,
                self._view.ui.list_pred_analysis_mono_model_chains,
                self._view.ui.btn_pred_analysis_mono_back_3,
                self._view.ui.btn_pred_analysis_mono_next_2,
                self._view.ui.btn_pred_analysis_mono_back_4,
                self._view.ui.btn_pred_analysis_mono_next_3,
                self._view.ui.btn_pred_analysis_mono_back_5,
                self._view.ui.btn_pred_analysis_mono_next_4,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            item = QtWidgets.QListWidgetItem(tmp_analysis_name)
            self._view.ui.list_pred_analysis_mono_overview.addItem(item)
            self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_mono_overview,
                self._view.ui.list_pred_analysis_mono_overview,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_2,
                self._view.ui.lbl_pred_analysis_mono_ref_chains,
                self._view.ui.list_pred_analysis_mono_ref_chains,
                self._view.ui.btn_pred_analysis_mono_back_4,
                self._view.ui.btn_pred_analysis_mono_next_3,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_mono_remove,
                self._view.ui.btn_pred_analysis_mono_add,
                self._view.ui.box_pred_analysis_mono_prot_struct_1,
                self._view.ui.box_pred_analysis_mono_prot_struct_2,
                self._view.ui.btn_pred_analysis_mono_back_3,
                self._view.ui.btn_pred_analysis_mono_next_2,
                self._view.ui.lbl_pred_analysis_mono_model_chains,
                self._view.ui.list_pred_analysis_mono_model_chains,
                self._view.ui.btn_pred_analysis_mono_back_5,
                self._view.ui.btn_pred_analysis_mono_next_4,
                self._view.ui.lbl_pred_analysis_mono_images,
                self._view.ui.cb_pred_analysis_mono_images,
                self._view.ui.btn_pred_analysis_mono_start,
                self._view.ui.btn_pred_analysis_mono_back_pred_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1.setText(
                self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText(),
            )
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2.setText(
                self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText(),
            )
            self._view.ui.list_pred_analysis_mono_ref_chains.clear()
            self._view.ui.btn_pred_analysis_mono_next_3.setEnabled(False)
            self._view.ui.list_pred_analysis_mono_ref_chains.setEnabled(True)

            for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
                if (
                    self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text()
                    == self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText()
                ):
                    self._view.ui.list_pred_analysis_mono_ref_chains.addItem(
                        self._view.ui.table_pred_analysis_mono_prot_to_predict.item(i, 0).text(),
                    )
            if self._view.ui.list_pred_analysis_mono_ref_chains.count() == 0:
                tmp_protein = self._current_project.search_protein(
                    self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText(),
                )
                for tmp_chain in tmp_protein.chains:
                    if tmp_chain.chain_type == "protein_chain":
                        self._view.ui.list_pred_analysis_mono_ref_chains.addItem(tmp_chain.chain_letter)
            if self._view.ui.list_pred_analysis_mono_ref_chains.count() == 1:
                self._view.ui.lbl_pred_analysis_mono_ref_chains.setText(
                    f"Select chain in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()}.",
                )
            else:
                self._view.ui.lbl_pred_analysis_mono_ref_chains.setText(
                    f"Select chains in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()}.",
                )
        self._view.wait_spinner.stop()

    # </editor-fold>

    # <editor-fold desc="Multimer Prediction + Analysis functions">
    def start_multimer_prediction_analysis(self) -> None:
        """Sets up the prediction process."""
        self._view.wait_spinner.start()
        self.prediction_type = constants.PREDICTION_TYPE_PRED_MULTI_ANALYSIS
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        self.update_status("Begin prediction process ...")
        predictions: list[
            prediction_protein_info.PredictionProteinInfo
        ] = prediction_util.get_prediction_name_and_seq_from_table(
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
        )

        self._active_task = tasks.Task(
            target=main_presenter_async.predict_protein_with_colabfold,
            args=(
                predictions,
                self.prediction_configuration,
                self._current_project,
            ),
            post_func=self.__await_multimer_prediction_for_subsequent_analysis,
        )
        self._active_task.start()

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

    def __await_multimer_prediction_for_subsequent_analysis(self, result: tuple) -> None:
        tmp_exit_code = result[0]
        tmp_exit_code_description = [1]
        if tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            self.block_box_prediction.destroy(True)
            constants.PYSSA_LOGGER.info("All structure predictions are done.")
            self.update_status("All structure predictions are done.")
            constants.PYSSA_LOGGER.info("Begin analysis process.")
            self.update_status("Begin analysis process ...")
            tmp_raw_analysis_run_names: list = []
            for row_no in range(self._view.ui.list_pred_analysis_multi_overview.count()):
                tmp_raw_analysis_run_names.append(self._view.ui.list_pred_analysis_multi_overview.item(row_no).text())

            self._active_task = tasks.Task(
                target=main_presenter_async.run_distance_analysis,
                args=(
                    tmp_raw_analysis_run_names,
                    self._current_project,
                    self._application_settings,
                    self._view.ui.cb_pred_analysis_multi_images.isChecked(),
                ),
                post_func=self.post_analysis_process,
            )
            self._active_task.start()

            if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
                os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
            self.block_box_analysis.exec_()
            self.display_view_page()

        elif tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because there was an error writing the fasta file(s)!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because the fasta file(s) could not be found!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            basic_boxes.ok(
                "Prediction",
                "Prediction failed because a subprocess failed!",
                QtWidgets.QMessageBox.Critical,
            )
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
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
            self.display_view_page()
            self._project_watcher.show_valid_options(self._view.ui)
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()

    def multi_pred_analysis_structure_analysis_next_2(self) -> None:
        """Shows the gui elements to select the chains in protein 1."""
        self._view.wait_spinner.start()
        tmp_proteins_to_predict: list[str] = []
        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            tmp_proteins_to_predict.append(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text(),
            )
        self._active_task = tasks.Task(
            target=main_presenter_async.check_chains_for_subsequent_analysis,
            args=(
                self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText(),
                self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText(),
                self._current_project,
                tmp_proteins_to_predict,
            ),
            post_func=self.__await_multi_pred_analysis_structure_analysis_next_2,
        )
        self._active_task.start()

    def __await_multi_pred_analysis_structure_analysis_next_2(self, result: tuple) -> None:
        # fixme: something is not quite right with the detection of the chains!
        _, tmp_analysis_name = result
        if tmp_analysis_name != "":
            gui_elements_to_show = [
                self._view.ui.btn_pred_analysis_multi_remove,
                self._view.ui.btn_pred_analysis_multi_add,
                self._view.ui.lbl_pred_analysis_multi_overview,
                self._view.ui.list_pred_analysis_multi_overview,
                self._view.ui.lbl_pred_analysis_multi_images,
                self._view.ui.cb_pred_analysis_multi_images,
                self._view.ui.btn_pred_analysis_multi_start,
                self._view.ui.btn_pred_analysis_multi_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.box_pred_analysis_multi_prot_struct_1,
                self._view.ui.box_pred_analysis_multi_prot_struct_2,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_3,
                self._view.ui.lbl_pred_analysis_multi_ref_chains,
                self._view.ui.list_pred_analysis_multi_ref_chains,
                self._view.ui.lbl_pred_analysis_multi_model_chains,
                self._view.ui.list_pred_analysis_multi_model_chains,
                self._view.ui.btn_pred_analysis_multi_back_3,
                self._view.ui.btn_pred_analysis_multi_next_2,
                self._view.ui.btn_pred_analysis_multi_back_4,
                self._view.ui.btn_pred_analysis_multi_next_3,
                self._view.ui.btn_pred_analysis_multi_back_5,
                self._view.ui.btn_pred_analysis_multi_next_4,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            item = QtWidgets.QListWidgetItem(tmp_analysis_name)
            self._view.ui.list_pred_analysis_multi_overview.addItem(item)
            self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_multi_overview,
                self._view.ui.list_pred_analysis_multi_overview,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_3,
                self._view.ui.lbl_pred_analysis_multi_ref_chains,
                self._view.ui.list_pred_analysis_multi_ref_chains,
                self._view.ui.btn_pred_analysis_multi_back_4,
                self._view.ui.btn_pred_analysis_multi_next_3,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_multi_remove,
                self._view.ui.btn_pred_analysis_multi_add,
                self._view.ui.box_pred_analysis_multi_prot_struct_1,
                self._view.ui.box_pred_analysis_multi_prot_struct_2,
                self._view.ui.btn_pred_analysis_multi_back_3,
                self._view.ui.btn_pred_analysis_multi_next_2,
                self._view.ui.lbl_pred_analysis_multi_model_chains,
                self._view.ui.list_pred_analysis_multi_model_chains,
                self._view.ui.btn_pred_analysis_multi_back_5,
                self._view.ui.btn_pred_analysis_multi_next_4,
                self._view.ui.lbl_pred_analysis_multi_images,
                self._view.ui.cb_pred_analysis_multi_images,
                self._view.ui.btn_pred_analysis_multi_start,
                self._view.ui.btn_pred_analysis_multi_back_pred_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1.setText(
                self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText(),
            )
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2.setText(
                self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText(),
            )
            self._view.ui.list_pred_analysis_multi_ref_chains.clear()
            self._view.ui.btn_pred_analysis_multi_next_3.setEnabled(False)
            self._view.ui.list_pred_analysis_multi_ref_chains.setEnabled(True)

            for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
                if (
                    self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text()
                    == self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText()
                ):
                    self._view.ui.list_pred_analysis_multi_ref_chains.addItem(
                        self._view.ui.table_pred_analysis_multi_prot_to_predict.item(i, 0).text(),
                    )
            if self._view.ui.list_pred_analysis_multi_ref_chains.count() == 0:
                tmp_protein = self._current_project.search_protein(
                    self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText(),
                )
                for tmp_chain in tmp_protein.chains:
                    if tmp_chain.chain_type == "protein_chain":
                        self._view.ui.list_pred_analysis_multi_ref_chains.addItem(tmp_chain.chain_letter)
            if self._view.ui.list_pred_analysis_multi_ref_chains.count() == 1:
                self._view.ui.lbl_pred_analysis_multi_ref_chains.setText(
                    f"Select chain in protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()}.",
                )
            else:
                self._view.ui.lbl_pred_analysis_multi_ref_chains.setText(
                    f"Select chains in protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()}.",
                )
        self._view.wait_spinner.stop()

    # </editor-fold>

    # <editor-fold desc="Structure Analysis functions">
    def post_analysis_process(self, an_exit_code: tuple[int, str]) -> None:
        """Post process after the analysis thread finished."""
        constants.PYSSA_LOGGER.debug("post_analysis_process() started ...")
        if an_exit_code[0] == exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[0]:
            self.block_box_analysis.destroy(True)
            basic_boxes.ok(
                "Distance analysis",
                "Distance analysis failed because there was an error during the analysis!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
            )
        elif an_exit_code[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_analysis.destroy(True)
            basic_boxes.ok(
                "Distance analysis",
                "Distance analysis failed because of an unknown error!",
                QtWidgets.QMessageBox.Critical,
            )
            constants.PYSSA_LOGGER.error(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
            )
        elif an_exit_code[0] == exit_codes.EXIT_CODE_ZERO[0]:
            self._current_project.serialize_project(self._current_project.get_project_xml_path())
            constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
            self.block_box_analysis.destroy(True)
            basic_boxes.ok(
                "Structure analysis",
                "All structure analysis' are done. Go to results to check the new results.",
                QtWidgets.QMessageBox.Information,
            )
            constants.PYSSA_LOGGER.info("All structure analysis' are done.")

        self._project_watcher.show_valid_options(self._view.ui)
        self.display_view_page()
        self._init_batch_analysis_page()

    def start_process_batch(self) -> None:
        """Sets up the worker for the analysis task."""
        constants.PYSSA_LOGGER.info("Begin analysis process.")

        tmp_raw_analysis_run_names: list = []
        for row_no in range(self._view.ui.list_analysis_batch_overview.count()):
            tmp_raw_analysis_run_names.append(self._view.ui.list_analysis_batch_overview.item(row_no).text())
        self._active_task = tasks.Task(
            target=main_presenter_async.run_distance_analysis,
            args=(
                tmp_raw_analysis_run_names,
                self._current_project,
                self._application_settings,
                self._view.ui.cb_analysis_images.isChecked(),
            ),
            post_func=self.post_analysis_process,
        )
        self._active_task.start()

        if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
            os.mkdir(constants.SCRATCH_DIR_ANALYSIS)

        self.block_box_analysis.exec_()
        self.display_view_page()

    def structure_analysis_next(self) -> None:
        """Shows the gui elements to select the chains in protein 1."""
        self._view.wait_spinner.start()
        self._active_task = tasks.Task(
            target=main_presenter_async.check_chains_for_analysis,
            args=(
                self._view.ui.box_analysis_batch_prot_struct_1.currentText(),
                self._view.ui.box_analysis_batch_prot_struct_2.currentText(),
                self._current_project,
            ),
            post_func=self.__await_structure_analysis_next,
        )
        self._active_task.start()

    def __await_structure_analysis_next(self, result: tuple) -> None:
        _, tmp_is_only_one_chain, tmp_analysis_run_name, tmp_protein_1, tmp_protein_2 = result

        if tmp_is_only_one_chain:
            gui_elements_to_show = [
                self._view.ui.btn_analysis_batch_remove,
                self._view.ui.btn_analysis_batch_add,
                self._view.ui.lbl_analysis_batch_overview,
                self._view.ui.list_analysis_batch_overview,
                self._view.ui.lbl_analysis_batch_images,
                self._view.ui.cb_analysis_batch_images,
                self._view.ui.btn_analysis_batch_start,
            ]
            gui_elements_to_hide = [
                self._view.ui.box_analysis_batch_prot_struct_1,
                self._view.ui.box_analysis_batch_prot_struct_2,
                self._view.ui.lbl_analysis_batch_prot_struct_1,
                self._view.ui.lbl_analysis_batch_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs,
                self._view.ui.lbl_analysis_batch_ref_chains,
                self._view.ui.list_analysis_batch_ref_chains,
                self._view.ui.lbl_analysis_batch_model_chains,
                self._view.ui.list_analysis_batch_model_chains,
                self._view.ui.btn_analysis_batch_back,
                self._view.ui.btn_analysis_batch_next,
                self._view.ui.btn_analysis_batch_back_2,
                self._view.ui.btn_analysis_batch_next_2,
                self._view.ui.btn_analysis_batch_back_3,
                self._view.ui.btn_analysis_batch_next_3,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            item = QtWidgets.QListWidgetItem(tmp_analysis_run_name)
            self._view.ui.list_analysis_batch_overview.addItem(item)
            self._view.ui.btn_analysis_batch_remove.setEnabled(False)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_analysis_batch_overview,
                self._view.ui.list_analysis_batch_overview,
                self._view.ui.lbl_analysis_batch_prot_struct_1,
                self._view.ui.lbl_analysis_batch_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs,
                self._view.ui.lbl_analysis_batch_ref_chains,
                self._view.ui.list_analysis_batch_ref_chains,
                self._view.ui.btn_analysis_batch_back_2,
                self._view.ui.btn_analysis_batch_next_2,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_analysis_batch_remove,
                self._view.ui.btn_analysis_batch_add,
                self._view.ui.box_analysis_batch_prot_struct_1,
                self._view.ui.box_analysis_batch_prot_struct_2,
                self._view.ui.btn_analysis_batch_back,
                self._view.ui.btn_analysis_batch_next,
                self._view.ui.lbl_analysis_batch_model_chains,
                self._view.ui.list_analysis_batch_model_chains,
                self._view.ui.btn_analysis_batch_back_3,
                self._view.ui.btn_analysis_batch_next_3,
                self._view.ui.lbl_analysis_batch_images,
                self._view.ui.cb_analysis_batch_images,
                self._view.ui.btn_analysis_batch_start,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.lbl_analysis_batch_prot_struct_1.setText(
                self._view.ui.box_analysis_batch_prot_struct_1.currentText(),
            )
            self._view.ui.lbl_analysis_batch_prot_struct_2.setText(
                self._view.ui.box_analysis_batch_prot_struct_2.currentText(),
            )
            self._view.ui.list_analysis_batch_ref_chains.clear()
            self._view.ui.btn_analysis_batch_next_2.setEnabled(False)
            self._view.ui.list_analysis_batch_ref_chains.setEnabled(True)

            for tmp_chain in tmp_protein_1.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self._view.ui.list_analysis_batch_ref_chains.addItem(tmp_chain.chain_letter)

            if self._view.ui.list_analysis_batch_ref_chains.count() == 1:
                self._view.ui.lbl_analysis_batch_ref_chains.setText(
                    f"Select chain in protein structure {self._view.ui.lbl_analysis_batch_prot_struct_1.text()}.",
                )
            else:
                self._view.ui.lbl_analysis_batch_ref_chains.setText(
                    f"Select chains in protein structure {self._view.ui.lbl_analysis_batch_prot_struct_1.text()}.",
                )
        self._view.wait_spinner.stop()

    # </editor-fold>

    # <editor-fold desc="Analysis Images">
    def post_image_creation_process(self) -> None:
        """Post method after the image creation task is finished."""
        self._current_project.serialize_project(self._current_project.get_project_xml_path())
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
        self._project_watcher.show_valid_options(self._view.ui)

    def start_automatic_image_creation(self) -> None:
        """Sets up the worker for the image creation process."""
        constants.PYSSA_LOGGER.info("Begin image creation process.")
        # self.worker_image_creation = workers.BatchImageWorkerPool(
        #    self._view.ui.list_analysis_images_struct_analysis, self._view.ui.list_analysis_images_creation_struct_analysis,
        #    self._view.status_bar, self._current_project)
        constants.PYSSA_LOGGER.info("Thread started for image creation process.")
        # self.threadpool.start(self.worker_image_creation)

        # <editor-fold desc="Worker setup">
        # TODO: test code below
        # --Begin: worker setup
        self.tmp_thread = QtCore.QThread()
        self.tmp_worker = task_workers.BatchImageWorker(
            self._view.ui.list_analysis_images_struct_analysis,
            self._view.ui.list_analysis_images_creation_struct_analysis,
            self._view.status_bar,
            self._current_project,
        )
        self.tmp_thread = task_workers.setup_worker_for_work(self.tmp_thread, self.tmp_worker, self.display_view_page)
        self.tmp_worker.finished.connect(self.post_image_creation_process)
        self.tmp_thread.start()
        # --End: worker setup

        # </editor-fold>

        if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
            os.mkdir(constants.SCRATCH_DIR_ANALYSIS)

        self._view.ui.list_analysis_images_struct_analysis.setEnabled(False)
        self._view.ui.list_analysis_images_creation_struct_analysis.setEnabled(False)

        # gui_elements_to_show = [
        #     self._view.ui.btn_analysis_abort,
        # ]
        # gui_elements_to_hide = [
        #     self._view.ui.btn_save_project,
        #     self._view.ui.btn_edit_page,
        #     self._view.ui.btn_view_page,
        #     self._view.ui.btn_use_page,
        #     self._view.ui.btn_export_project,
        #     self._view.ui.btn_close_project,
        #     self._view.ui.btn_batch_analysis_page,
        #     self._view.ui.btn_image_analysis_page,
        #     self._view.ui.btn_results_page,
        #     self._view.ui.lbl_handle_pymol_session,
        #     self._view.ui.btn_image_page,
        #     self._view.ui.btn_hotspots_page,
        # ]
        # gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

        self.block_box_images.exec_()

    # </editor-fold>

    # <editor-fold desc="Results page functions">
    def show_analysis_results_options(self) -> None:
        """Shows the combo box of protein pairs."""
        self.results_management.show_stage_x(0)

    def show_results_interactions(self, gui_elements_to_show: list = None, gui_elements_to_hide: list = None) -> None:
        """Shows the gui elements of the results page.

        Args:
            gui_elements_to_show: a list of gui elements which should get displayed.
            gui_elements_to_hide: a list of gui elements which should get hidden.
        """
        if gui_elements_to_hide is not None:
            self.results_management.show_gui_elements_stage_x(
                [0, 1],
                [],
                show_specific_elements=[
                    self._view.ui.lbl_results_analysis_options,
                    self._view.ui.cb_results_analysis_options,
                ],
                hide_specific_elements=gui_elements_to_hide,
            )
        else:
            self.results_management.show_gui_elements_stage_x(
                [0, 1],
                [],
                show_specific_elements=[
                    self._view.ui.lbl_results_analysis_options,
                    self._view.ui.cb_results_analysis_options,
                ],
            )

    def load_results(self) -> None:
        """Sets up the worker for the result loading process."""
        self._view.wait_spinner.start()
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self.results_name = self._view.ui.cb_results_analysis_options.currentText()
        if self.results_name == "":
            self.show_analysis_results_options()
            self._view.wait_spinner.stop()
            return
        tools.ask_to_save_pymol_session(self._current_project, self.current_session, self._application_settings)

        self._active_task = tasks.Task(
            target=main_presenter_async.load_results,
            args=(
                self._current_project,
                self.results_name,
            ),
            post_func=self.__await_load_results,
        )
        self._active_task.start()

        self._view.ui.list_results_interest_regions.clear()
        self._view.status_bar.showMessage(f"Loading results of {self.results_name} ...")

    def __await_load_results(self, result: tuple) -> None:
        _, images_type, self.current_session, tmp_rmsd, tmp_no_aligned_residues = result

        if images_type == constants.IMAGES_ALL:
            gui_elements_to_show = [
                self._view.ui.lbl_results_analysis_options,
                self._view.ui.cb_results_analysis_options,
                self._view.ui.lbl_results_rmsd,
                self._view.ui.txt_results_rmsd,
                self._view.ui.lbl_color_rmsd,
                self._view.ui.btn_color_rmsd,
                self._view.ui.lbl_results_aligned_residues,
                self._view.ui.txt_results_aligned_residues,
                self._view.ui.lbl_results_distance_plot,
                self._view.ui.btn_view_distance_plot,
                self._view.ui.lbl_results_distance_histogram,
                self._view.ui.btn_view_distance_histogram,
                self._view.ui.lbl_results_distance_table,
                self._view.ui.btn_view_distance_table,
                self._view.ui.lbl_results_structure_alignment,
                self._view.ui.btn_view_struct_alignment,
                self._view.ui.lbl_results_interest_regions,
                self._view.ui.list_results_interest_regions,
                self._view.ui.btn_view_interesting_region,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            for tmp_filename in os.listdir(constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                self._view.ui.list_results_interest_regions.addItem(tmp_filename)
            self._view.ui.list_results_interest_regions.sortItems()
        elif images_type == constants.IMAGES_STRUCT_ALN_ONLY:
            gui_elements_to_show = [
                self._view.ui.lbl_results_analysis_options,
                self._view.ui.cb_results_analysis_options,
                self._view.ui.lbl_results_rmsd,
                self._view.ui.txt_results_rmsd,
                self._view.ui.lbl_color_rmsd,
                self._view.ui.btn_color_rmsd,
                self._view.ui.lbl_results_aligned_residues,
                self._view.ui.txt_results_aligned_residues,
                self._view.ui.lbl_results_distance_plot,
                self._view.ui.btn_view_distance_plot,
                self._view.ui.lbl_results_distance_histogram,
                self._view.ui.btn_view_distance_histogram,
                self._view.ui.lbl_results_distance_table,
                self._view.ui.btn_view_distance_table,
                self._view.ui.lbl_results_structure_alignment,
                self._view.ui.btn_view_struct_alignment,
            ]
            gui_elements_to_hide = [
                self._view.ui.lbl_results_interest_regions,
                self._view.ui.list_results_interest_regions,
                self._view.ui.btn_view_interesting_region,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.list_results_interest_regions.sortItems()
        elif images_type == constants.IMAGES_NONE:
            gui_elements_to_show = [
                self._view.ui.lbl_results_analysis_options,
                self._view.ui.cb_results_analysis_options,
                self._view.ui.lbl_results_rmsd,
                self._view.ui.txt_results_rmsd,
                self._view.ui.lbl_color_rmsd,
                self._view.ui.btn_color_rmsd,
                self._view.ui.lbl_results_aligned_residues,
                self._view.ui.txt_results_aligned_residues,
                self._view.ui.lbl_results_distance_plot,
                self._view.ui.btn_view_distance_plot,
                self._view.ui.lbl_results_distance_histogram,
                self._view.ui.btn_view_distance_histogram,
                self._view.ui.lbl_results_distance_table,
                self._view.ui.btn_view_distance_table,
            ]
            gui_elements_to_hide = [
                self._view.ui.lbl_results_structure_alignment,
                self._view.ui.btn_view_struct_alignment,
                self._view.ui.lbl_results_interest_regions,
                self._view.ui.list_results_interest_regions,
                self._view.ui.btn_view_interesting_region,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            raise ValueError("Illegal argument.")

        self._view.ui.txt_results_rmsd.setText(str(tmp_rmsd))
        self._view.ui.txt_results_aligned_residues.setText(str(tmp_no_aligned_residues))
        self._view.status_bar.showMessage(f"Current workspace: {str(self._workspace_path)}")
        self.main_window_state.results_page.results_name = self._view.ui.cb_results_analysis_options.currentText()
        self._view.wait_spinner.stop()

    def color_protein_pair_by_rmsd(self) -> None:
        """Colors the residues in 5 colors depending on their distance to the reference."""
        self._view.wait_spinner.start()
        self._active_task = tasks.Task(
            target=main_presenter_async.color_protein_pair_by_rmsd_value,
            args=(
                self._current_project,
                self.results_name,
            ),
            post_func=self.__await_color_protein_pair_by_rmsd,
        )
        self._active_task.start()
        self._view.status_bar.showMessage("Coloring protein pair by RMSD value ...")
        # hide unnecessary representations
        # fixme: it might be a problem to hide any representation at this point
        # cmd.hide("cartoon", tmp_protein_pair.protein_1.get_molecule_object())
        # cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}")
        # cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}")

    def __await_color_protein_pair_by_rmsd(self, result: tuple) -> None:
        self._view.wait_spinner.stop()

    # </editor-fold>

    # <editor-fold desc="Display page functions">
    def display_structure_alignment(self) -> None:
        """Opens a window which displays the image of the structure alignment."""
        png_dialog = QtWidgets.QDialog(self._view)
        label = QtWidgets.QLabel(self._view)
        pathlib.Path(
            f"{self._workspace_path}/{self._view.ui.lbl_current_project_name.text()}/results/{self.results_name}",
        )
        self._view.ui.cb_results_analysis_options.currentText()
        pixmap = QtGui.QPixmap(
            f"{constants.CACHE_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{self._view.ui.cb_results_analysis_options.currentText()}",
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
        protein_pair_of_analysis = self._current_project.search_protein_pair(
            self._view.ui.cb_results_analysis_options.currentText(),
        )

        tmp_dialog = plot_view.PlotView(protein_pair_of_analysis, self._current_project, self._view.ui.cb_results_analysis_options.currentText())
        tmp_dialog.exec_()

        # self.distance_plot_dialog = dialog_distance_plot.DialogDistancePlot(protein_pair_of_analysis, self._current_project,
        #                                                                     self._view.ui.cb_results_analysis_options.currentText())
        # self.distance_plot_dialog.setWindowModality(Qt.WindowModal)
        # self.is_distance_plot_open = True
        # self.distance_plot_dialog.exec_()

    def display_distance_histogram(self) -> None:
        """Opens a window which displays the distance histogram."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        protein_pair_of_analysis = self._current_project.search_protein_pair(
            self._view.ui.cb_results_analysis_options.currentText(),
        )
        dialog = dialog_distance_histogram.DialogDistanceHistogram(protein_pair_of_analysis)
        dialog.exec_()

    def display_interesting_region(self) -> None:
        """Displays an image of an interesting region."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        png_dialog = QtWidgets.QDialog(self._view)
        label = QtWidgets.QLabel(self._view)
        file_name = self._view.ui.list_results_interest_regions.currentItem().text()
        pixmap = QtGui.QPixmap(f"{constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{file_name}")
        # TODO: Create setting for min. image size
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
        table_dialog = QtWidgets.QDialog(self._view)
        table_view = QtWidgets.QTableView()
        table_view.setModel(csv_model)

        tmp_protein_pair = self._current_project.search_protein_pair(
            self._view.ui.cb_results_analysis_options.currentText(),
        )
        csv_filepath = pathlib.Path(f"{constants.CACHE_CSV_DIR}/{tmp_protein_pair.name}.csv")
        if not os.path.exists(constants.CACHE_CSV_DIR):
            os.mkdir(constants.CACHE_CSV_DIR)
        tmp_protein_pair = self._current_project.search_protein_pair(self.results_name)

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
        if self._view.ui.box_manage_choose_protein.currentText() != "":
            gui_elements_to_show = [
                self._view.ui.lbl_manage_choose_color,
                self._view.ui.box_manage_choose_color,
                self._view.ui.lbl_manage_choose_representation,
                self._view.ui.box_manage_choose_representation,
                self._view.ui.lbl_manage_choose_bg_color,
                self._view.ui.box_manage_choose_bg_color,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
            if (
                cmd.select(
                    name="disulfides",
                    selection=f"{self._view.ui.box_manage_choose_protein.currentText()} & {tmp_pymol_selection_option}",
                )
                > 0
            ):
                gui_elements_to_show = [
                    self._view.ui.lbl_disulfid_bond_1,
                    self._view.ui.lbl_disulfid_bond_2,
                    self._view.ui.btn_disulfid_bond_show,
                    self._view.ui.btn_disulfid_bond_hide,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
        else:
            gui_elements_to_hide = [
                self._view.ui.lbl_manage_choose_color,
                self._view.ui.box_manage_choose_color,
                self._view.ui.lbl_manage_choose_representation,
                self._view.ui.box_manage_choose_representation,
                self._view.ui.lbl_manage_choose_bg_color,
                self._view.ui.box_manage_choose_bg_color,
                self._view.ui.lbl_disulfid_bond_1,
                self._view.ui.lbl_disulfid_bond_2,
                self._view.ui.btn_disulfid_bond_show,
                self._view.ui.btn_disulfid_bond_hide,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)

    def choose_manage_color_selected_protein(self) -> None:
        """Sets the protein color."""
        tmp_input = self._view.ui.box_manage_choose_protein.currentText()
        tmp_protein = self._current_project.search_protein(tmp_input)
        tmp_protein.color_protein_in_pymol(
            self._view.ui.box_manage_choose_color.currentText(),
            f"/{tmp_protein.get_molecule_object()}",
        )

    def choose_manage_representation(self) -> None:
        """Sets the representation."""
        tmp_input = self._view.ui.box_manage_choose_protein.currentText()
        tmp_protein = self._current_project.search_protein(tmp_input)
        tmp_selection = f"/{tmp_protein.get_molecule_object()}"
        if self._view.ui.box_manage_choose_representation.currentIndex() == 0:
            print("Please select a representation.")
            self._view.status_bar.showMessage("Please select a representation.")
        elif self._view.ui.box_manage_choose_representation.currentIndex() == 1:
            cmd.show("cartoon", tmp_selection)
            cmd.hide("ribbon", tmp_selection)
        elif self._view.ui.box_manage_choose_representation.currentIndex() == 2:
            cmd.show("ribbon", tmp_selection)
            cmd.hide("cartoon", tmp_selection)
        else:
            print("Missing implementation!")

    def choose_manage_bg_color(self) -> None:
        """Sets the background color."""
        if self._view.ui.box_manage_choose_bg_color.currentIndex() == 0:
            print("Please select a background color.")
            self._view.status_bar.showMessage("Please select a background color.")
        elif self._view.ui.box_manage_choose_bg_color.currentIndex() == 1:
            cmd.bg_color("black")
        elif self._view.ui.box_manage_choose_bg_color.currentIndex() == 2:
            cmd.bg_color("white")
        else:
            print("Missing implementation!")

    def show_disulfid_bonds_as_sticks(self) -> None:
        """Shows all disulfid bonds within the pymol session."""
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        cmd.select(
            name="disulfides",
            selection=f"{self._view.ui.box_manage_choose_protein.currentText()} & {tmp_pymol_selection_option}",
        )
        cmd.color(color="atomic", selection="disulfides and not elem C")
        cmd.set("valence", 0)  # this needs to be better implemented
        cmd.show("sticks", "disulfides")
        cmd.hide("sticks", "elem H")

    def hide_disulfid_bonds_as_sticks(self) -> None:
        """Hides all disulfid bonds within the pymol session."""
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        cmd.select(
            name="disulfides",
            selection=f"{self._view.ui.box_manage_choose_protein.currentText()} & {tmp_pymol_selection_option}",
        )
        cmd.hide("sticks", "disulfides")

    # </editor-fold>

    # <editor-fold desc="Image page functions">
    def show_representation(self) -> None:
        """Sets the representation."""
        if self._view.ui.box_representation.currentIndex() == 0:
            self.main_window_state.image_page.representation = ""
            self._view.status_bar.showMessage("Please select a representation.")
        elif self._view.ui.box_representation.currentIndex() == 1:
            cmd.show("cartoon", "all")
            cmd.hide("ribbon", "all")
            self.main_window_state.image_page.representation = "cartoon"
        elif self._view.ui.box_representation.currentIndex() == 2:
            cmd.show("ribbon", "all")
            cmd.hide("cartoon", "all")
            self.main_window_state.image_page.representation = "ribbon"
        else:
            print("Missing implementation!")

    def choose_bg_color(self) -> None:
        """Sets the background color."""
        if self._view.ui.box_bg_color.currentIndex() == 0:
            self.main_window_state.image_page.background_color = ""
            self._view.status_bar.showMessage("Please select a background color.")
        elif self._view.ui.box_bg_color.currentIndex() == 1:
            cmd.bg_color("black")
            self.main_window_state.image_page.background_color = "black"
        elif self._view.ui.box_bg_color.currentIndex() == 2:
            cmd.bg_color("white")
            self.main_window_state.image_page.background_color = "white"
        else:
            print("Missing implementation!")

    def choose_renderer(self) -> None:
        """Sets the renderer."""
        if self._view.ui.box_renderer.currentIndex() == 0:
            self._view.status_bar.showMessage("Please select a renderer.")
            self._view.ui.cb_ray_tracing.hide()
            self._view.ui.label_26.hide()
        elif self._view.ui.box_renderer.currentIndex() == 1:
            self.renderer = "-1"
            self._view.ui.cb_ray_tracing.show()
            self._view.ui.label_26.show()
        elif self._view.ui.box_renderer.currentIndex() == 2:
            self.renderer = "0"
            self._view.ui.cb_ray_tracing.show()
            self._view.ui.label_26.show()
        else:
            print("Missing implementation!")
        self.main_window_state.image_page.renderer = self._view.ui.box_renderer.currentIndex()

    def choose_ray_trace_mode(self) -> None:
        """Sets the ray-trace mode."""
        if self._view.ui.box_ray_trace_mode.currentIndex() == 0:
            self._view.status_bar.showMessage("Please select a Ray-Trace-Mode.")
        elif self._view.ui.box_ray_trace_mode.currentIndex() == 1:
            cmd.set("ray_trace_mode", 0)
        elif self._view.ui.box_ray_trace_mode.currentIndex() == 2:
            cmd.set("ray_trace_mode", 1)
        elif self._view.ui.box_ray_trace_mode.currentIndex() == 3:
            cmd.set("ray_trace_mode", 2)
        elif self._view.ui.box_ray_trace_mode.currentIndex() == 4:
            cmd.set("ray_trace_mode", 3)
        else:
            print("Missing implementation!")
        self.main_window_state.image_page.ray_trace_mode = self._view.ui.box_ray_trace_mode.currentIndex()

    def choose_ray_texture(self) -> None:
        """Sets the ray texture."""
        if self._view.ui.box_ray_texture.currentIndex() == 0:
            print("Please select a Ray Texture.")
            self._view.status_bar.showMessage("Please select a Ray Texture.")
        elif self._view.ui.box_ray_texture.currentIndex() == 1:
            cmd.set("ray_texture", 0)
        elif self._view.ui.box_ray_texture.currentIndex() == 2:
            cmd.set("ray_texture", 1)
        elif self._view.ui.box_ray_texture.currentIndex() == 3:
            cmd.set("ray_texture", 2)
        elif self._view.ui.box_ray_texture.currentIndex() == 4:
            cmd.set("ray_texture", 3)
        elif self._view.ui.box_ray_texture.currentIndex() == 5:
            cmd.set("ray_texture", 4)
        elif self._view.ui.box_ray_texture.currentIndex() == 6:
            cmd.set("ray_texture", 5)
        else:
            print("Missing implementation!")
        self.main_window_state.image_page.ray_texture = self._view.ui.box_ray_texture.currentIndex()

    def decide_ray_tracing(self) -> None:
        """Sets the ray-tracing options."""
        print(self._view.ui.cb_ray_tracing.isChecked())
        if self._view.ui.cb_ray_tracing.isChecked():
            self._view.ui.cb_transparent_bg.hide()
            self._view.ui.label_23.hide()
            self._view.ui.box_renderer.setEnabled(False)
            self._view.ui.label_10.show()
            self._view.ui.box_ray_trace_mode.show()
            self._view.ui.label_14.show()
            self._view.ui.box_ray_texture.show()
            cmd.set("ray_opaque_background", "off")
        else:
            self._view.ui.cb_transparent_bg.show()
            self._view.ui.label_23.show()
            self._view.ui.box_renderer.setEnabled(True)
            self._view.ui.label_10.hide()
            self._view.ui.box_ray_trace_mode.hide()
            self._view.ui.label_14.hide()
            self._view.ui.box_ray_texture.hide()
        self.main_window_state.image_page.ray_tracing = self._view.ui.cb_ray_tracing.isChecked()

    def decide_transparent_bg(self) -> None:
        """Sets the transparent background."""
        if self._view.ui.cb_transparent_bg.isChecked():
            cmd.set("opaque_background", "off")
        else:
            cmd.set("opaque_background", "on")
        self.main_window_state.image_page.transparent_background = self._view.ui.cb_transparent_bg.isChecked()

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

    def post_preview_image(self) -> None:
        """Hides the block box of the preview process."""
        self.block_box_uni.hide()
        self.block_box_uni.destroy(True)
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
            self.block_box_uni.exec_()
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

    # <editor-fold desc="Hotspots page functions">
    def open_protein_for_hotspots(self) -> None:
        """Loads the selected protein's pymol session."""
        self._view.wait_spinner.start()
        try:
            tmp_name = self._view.ui.list_hotspots_choose_protein.currentItem().text()
        except AttributeError:
            self._view.wait_spinner.stop()
            return
        if self._view.ui.list_hotspots_choose_protein.currentItem().text() != "":
            tools.ask_to_save_pymol_session(self._current_project, self.current_session, self._application_settings)
            self._active_task = tasks.Task(
                target=main_presenter_async.open_protein_for_hotspots,
                args=(
                    tmp_name,
                    self._current_project,
                    self.current_session,
                ),
                post_func=self.__await_open_protein_for_hotspots,
            )
            self._active_task.start()
            self.update_status("Loading PyMOL session ...")
        else:
            self._view.wait_spinner.stop()
            self.update_status(self._workspace_status)

    def __await_open_protein_for_hotspots(self, result: tuple) -> None:
        _, is_protein, is_protein_pair, self.current_session = result
        if is_protein and not is_protein_pair:
            gui_elements_to_show = [
                self._view.ui.lbl_hotspots_resi_show,
                self._view.ui.btn_hotspots_resi_show,
                self._view.ui.lbl_hotspots_resi_hide,
                self._view.ui.btn_hotspots_resi_hide,
                self._view.ui.lbl_hotspots_resi_zoom,
                self._view.ui.btn_hotspots_resi_zoom,
                self._view.ui.btn_manage_session,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
        elif is_protein_pair and not is_protein:
            gui_elements_to_show = [
                self._view.ui.lbl_hotspots_resi_show,
                self._view.ui.btn_hotspots_resi_show,
                self._view.ui.lbl_hotspots_resi_hide,
                self._view.ui.btn_hotspots_resi_hide,
                self._view.ui.lbl_hotspots_resi_zoom,
                self._view.ui.btn_hotspots_resi_zoom,
                self._view.ui.btn_manage_session,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
        else:
            pass
        self._view.wait_spinner.stop()
        self.update_status(self._workspace_status)

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

    # </editor-fold>

    # <editor-fold desc="GUI related methods">

    # <editor-fold desc="GUI page management functions">
    def _create_batch_analysis_management(self) -> None:
        """Creates a list of gui stages."""
        # gui element management
        tmp_stages = [
            # add a prot analysis: stage 0
            stage.Stage(
                {
                    "label_batch_analysis_overview": self._view.ui.lbl_analysis_batch_overview,
                    "box_protein_structure_1": self._view.ui.list_analysis_batch_overview,
                },
                {
                    "add_button": self._view.ui.btn_analysis_batch_add,
                    "remove_button": self._view.ui.btn_analysis_batch_remove,
                },
            ),
            # choose protein structures: stage 1
            stage.Stage(
                {
                    "label_protein_structure_1": self._view.ui.lbl_analysis_batch_prot_struct_1,
                    "box_protein_structure_1": self._view.ui.box_analysis_batch_prot_struct_1,
                    "label_vs": self._view.ui.lbl_analysis_batch_vs,
                    "label_protein_structure_2": self._view.ui.lbl_analysis_batch_prot_struct_2,
                    "box_protein_structure_2": self._view.ui.box_analysis_batch_prot_struct_2,
                },
                {
                    "next_button": self._view.ui.btn_analysis_batch_next,
                    "back_button": self._view.ui.btn_analysis_batch_back,
                },
            ),
            # choose chains from prot structure 1: stage 2
            stage.Stage(
                {
                    "label_protein_structure_1_chains": self._view.ui.lbl_analysis_batch_ref_chains,
                    "list_protein_structure_1_chains": self._view.ui.list_analysis_batch_ref_chains,
                },
                {
                    "back_button": self._view.ui.btn_analysis_batch_back_2,
                    "next_button": self._view.ui.btn_analysis_batch_next_2,
                },
            ),
            # choose chains from prot structure 2: stage 3
            stage.Stage(
                {
                    "label_protein_structure_2_chains": self._view.ui.lbl_analysis_batch_model_chains,
                    "list_protein_structure_2_chains": self._view.ui.list_analysis_batch_model_chains,
                },
                {
                    "back_button": self._view.ui.btn_analysis_batch_back_3,
                    "next_button": self._view.ui.btn_analysis_batch_next_3,
                },
            ),
            # start batch run: stage 4
            stage.Stage(
                {
                    "label_images": self._view.ui.lbl_analysis_batch_images,
                    "checkbox_images": self._view.ui.cb_analysis_batch_images,
                },
                {
                    "start_button": self._view.ui.btn_analysis_batch_start,
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
                    "label_analysis_options": self._view.ui.lbl_results_analysis_options,
                    "box_results_analysis_options": self._view.ui.cb_results_analysis_options,
                },
                {
                    "": None,
                },
            ),
            # choose chains from prot structure 1: stage 1
            stage.Stage(
                {
                    "label_results_rmsd": self._view.ui.lbl_results_rmsd,
                    "text_results_rmsd": self._view.ui.txt_results_rmsd,
                    "label_color_rmsd": self._view.ui.lbl_color_rmsd,
                    "button_color_rmsd": self._view.ui.btn_color_rmsd,
                    "label_results_aligned_residues": self._view.ui.lbl_results_aligned_residues,
                    "text_results_aligned_residues": self._view.ui.txt_results_aligned_residues,
                    "label_results_distance_plot": self._view.ui.lbl_results_distance_plot,
                    "button_view_distance_plot": self._view.ui.btn_view_distance_plot,
                    "label_results_distance_histogram": self._view.ui.lbl_results_distance_histogram,
                    "button_view_distance_histogram": self._view.ui.btn_view_distance_histogram,
                    "label_results_distance_table": self._view.ui.lbl_results_distance_table,
                    "button_view_distance_table": self._view.ui.btn_view_distance_table,
                    "label_results_structure_alignment": self._view.ui.lbl_results_structure_alignment,
                    "button_view_struct_alignment": self._view.ui.btn_view_struct_alignment,
                    "label_results_interest_regions": self._view.ui.lbl_results_interest_regions,
                    "list_results_interest_regions": self._view.ui.list_results_interest_regions,
                    "button_results_interest_regions": self._view.ui.btn_view_interesting_region,
                },
                {
                    "": None,
                },
            ),
        ]
        self.results_management = gui_page_management.GuiPageManagement(tmp_stages)

    # </editor-fold>

    # <editor-fold desc="Page init functions">
    def _init_fill_combo_boxes(self) -> None:
        """Fills all combo boxes of the plugin."""
        item_list_representation = [
            "",
            "cartoon",
            "ribbon",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_representation, item_list_representation)
        gui_utils.fill_combo_box(self._view.ui.box_manage_choose_representation, item_list_representation)
        # combo box BgColor
        item_list_bg_color = [
            "",
            "black",
            "white",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_bg_color, item_list_bg_color)
        gui_utils.fill_combo_box(self._view.ui.box_manage_choose_bg_color, item_list_bg_color)
        # combo box Renderer
        item_list_renderer = [
            "",
            "default renderer",
            "PyMOL internal renderer",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_renderer, item_list_renderer)
        # combo box RayTraceMode
        item_list_ray_trace_mode = [
            "",
            "normal color",
            "normal color + black outline",
            "black outline only",
            "quantized color + black outline",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_ray_trace_mode, item_list_ray_trace_mode)
        # combo box Ray Texture
        item_list_ray_texture = [
            "",
            "None",
            "Matte 1",
            "Matte 2",
            "Swirl 1",
            "Fiber",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_ray_texture, item_list_ray_texture)
        # combo box protein Colors
        gui_utils.fill_combo_box(self._view.ui.box_manage_choose_color, constants.PYMOL_COLORS)

    def _init_new_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self._view.ui.txt_new_project_name.clear()
        self._view.ui.txt_new_choose_reference.clear()
        self._view.ui.lbl_new_status_project_name.setText("")
        self._view.ui.lbl_new_status_choose_reference.setText("")
        self._view.ui.cb_new_add_reference.setCheckState(0)
        self._view.ui.btn_new_create_project.setEnabled(False)
        styles.color_button_not_ready(self._view.ui.btn_new_create_project)

    def _init_use_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        gui_elements = [
            self._view.ui.lbl_use_search,
            self._view.ui.lbl_use_status_search,
            self._view.ui.txt_use_search,
            self._view.ui.btn_use_add_available_protein_structures,
            self._view.ui.lbl_use_available_protein_structures,
            self._view.ui.list_use_available_protein_structures,
            self._view.ui.btn_use_remove_selected_protein_structures,
            self._view.ui.lbl_use_selected_protein_structures,
            self._view.ui.list_use_selected_protein_structures,
            self._view.ui.btn_use_back,
            self._view.ui.btn_use_create_new_project,
            self._view.ui.lbl_use_new_project,
        ]
        gui_utils.hide_gui_elements(gui_elements)
        self._view.ui.txt_use_project_name.clear()
        self._view.ui.lbl_use_status_project_name.setText("")
        self._view.ui.txt_use_search.clear()
        self._view.ui.lbl_use_status_search.setText("")
        self._view.ui.list_use_available_protein_structures.clear()
        self._view.ui.list_use_selected_protein_structures.clear()
        self._view.ui.list_use_existing_projects.clear()
        self._view.ui.btn_use_next.setEnabled(False)
        self.hide_protein_selection_for_use()

    def _init_edit_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self._view.ui.list_edit_project_proteins.clear()
        gui_elements_to_hide = [
            self._view.ui.lbl_edit_clean_new_prot,
            self._view.ui.btn_edit_clean_new_prot,
            self._view.ui.lbl_edit_clean_update_prot,
            self._view.ui.btn_edit_clean_update_prot,
            self._view.ui.label_12,
            self._view.ui.btn_edit_project_delete,
            self._view.ui.label_15,
            self._view.ui.btn_edit_protein_rename,
            self._view.ui.btn_edit_project_save,
            self._view.ui.label_13,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.fill_list_view_with_protein_names(self._current_project, self._view.ui.list_edit_project_proteins)
        # self.project_scanner.scan_project_for_valid_proteins(self._view.ui.list_edit_project_proteins)

    def _init_sequence_vs_pdb_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self._view.ui.list_s_v_p_ref_chains.clear()
        # # sets up defaults: Prediction + Analysis
        #
        # #self._view.ui.label_18.hide()
        # #self._view.ui.txt_prediction_project_name.hide()
        # #self._view.ui.lbl_prediction_status_project_name.hide()
        # # stage 1
        # self._view.ui.lbl_prediction_status_project_name.setText("")
        #
        # #self._view.ui.list_widget_projects.hide()
        # #self._view.ui.btn_prediction_next_1.hide()
        # # stage 2
        # self._view.ui.lbl_prediction_load_reference.hide()
        # self._view.ui.txt_prediction_load_reference.clear()
        # self._view.ui.txt_prediction_load_reference.hide()
        # self._view.ui.lbl_prediction_status_load_reference.setText("")
        # self._view.ui.btn_prediction_back_2.setEnabled(False)
        # self._view.ui.btn_prediction_back_2.hide()
        # self._view.ui.btn_prediction_next_2.setEnabled(False)
        # self._view.ui.btn_prediction_next_2.hide()
        #
        # # stage 3
        # self._view.ui.lbl_prediction_ref_chains.hide()
        # self._view.ui.list_widget_ref_chains.hide()
        # # TO-DO: stylesheet needs to be integrated into the styles.css
        # self._view.ui.list_widget_ref_chains.setStyleSheet("border-style: solid;"
        #                                              "border-width: 2px;"
        #                                              "border-radius: 8px;"
        #                                              "border-color: #DCDBE3;")
        # self._view.ui.list_widget_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        # self._view.ui.btn_prediction_back_3.setEnabled(False)
        # self._view.ui.btn_prediction_back_3.hide()
        # self._view.ui.btn_prediction_start.setEnabled(False)
        # self._view.ui.btn_prediction_start.hide()
        # # TO-DO: needs to be removed if model chains are implemented at the right spot
        # self._view.ui.lbl_prediction_model_chains.hide()
        # self._view.ui.txt_prediction_chain_model.hide()

    def _init_results_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        # stage 1
        self._view.ui.list_results_interest_regions.clear()
        self._view.ui.txt_results_rmsd.clear()
        self._view.ui.txt_results_aligned_residues.clear()

    def _init_analysis_image_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        self._view.ui.list_analysis_images_struct_analysis.setEnabled(True)
        self._view.ui.list_analysis_images_creation_struct_analysis.setEnabled(True)
        self.display_image_analysis_page()

    def _init_image_page(self) -> None:
        """Hides everything which is needed."""
        self._view.ui.label_10.hide()
        self._view.ui.box_ray_trace_mode.hide()
        self._view.ui.label_14.hide()
        self._view.ui.box_ray_texture.hide()

    def _init_batch_analysis_page(self) -> None:
        """Clears all text fields and hides everything which is needed."""
        # sets up defaults: Batch
        self.batch_analysis_management.show_stage_x(0)
        self._view.ui.list_analysis_batch_overview.clear()
        self._view.ui.btn_analysis_batch_remove.hide()

    def _init_all_pages(self) -> None:
        """Collection of all init methods to reset all page defaults."""
        self._init_local_pred_mono_page()
        self._init_local_pred_multi_page()
        self._init_mono_pred_analysis_page()
        self._init_multi_pred_analysis_page()
        self._init_sequence_vs_pdb_page()
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
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 0, "Home")

    def display_job_analysis_page(self) -> None:
        """Displays the job analysis work area."""
        self._view.ui.list_analysis_batch_ref_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self._view.ui.list_analysis_batch_model_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        # regular work area opening
        self._init_batch_analysis_page()
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 4, "Structure Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_batch_analysis_page,
        )

    def display_results_page(self) -> None:
        """Displays the results work area."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        results = []
        results.insert(0, "")

        for tmp_protein_pair in self._current_project.protein_pairs:
            results.append(tmp_protein_pair.name)
            if tmp_protein_pair.name == self.results_name:
                pass
        self._view.ui.cb_results_analysis_options.clear()
        gui_utils.fill_combo_box(self._view.ui.cb_results_analysis_options, results)

        desired_string: str = self.main_window_state.results_page.results_name
        # Find the index of the string in the combo box
        index = self._view.ui.cb_results_analysis_options.findText(desired_string)
        # Set the current index of the combo box
        if index != -1:
            self._view.ui.cb_results_analysis_options.setCurrentIndex(index)
        else:
            print(f"String '{desired_string}' not found in the combo box.")

        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 5, "Results")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_results_page,
        )

    def display_image_page(self) -> None:
        """Displays the image work area."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self._init_image_page()

        if self.main_window_state.image_page.transparent_background:
            self._view.ui.cb_transparent_bg.setChecked(True)
        else:
            self._view.ui.cb_transparent_bg.setChecked(False)
        if self.main_window_state.image_page.ray_tracing:
            self._view.ui.cb_ray_tracing.setChecked(True)
        else:
            self._view.ui.cb_ray_tracing.setChecked(False)
        try:
            self._view.ui.box_representation.setCurrentIndex(
                self._view.ui.box_representation.findText(self.main_window_state.image_page.representation),
            )
            self._view.ui.box_bg_color.setCurrentIndex(
                self._view.ui.box_bg_color.findText(self.main_window_state.image_page.background_color),
            )
            self._view.ui.box_renderer.setCurrentIndex(self.main_window_state.image_page.renderer)
            self._view.ui.box_ray_trace_mode.setCurrentIndex(self.main_window_state.image_page.ray_trace_mode)
            self._view.ui.box_ray_texture.setCurrentIndex(self.main_window_state.image_page.ray_texture)
        except IndexError:
            constants.PYSSA_LOGGER.error("An option is invalid.")

        if self._view.ui.box_renderer.currentText() == "":
            self._view.ui.cb_ray_tracing.hide()
            self._view.ui.label_26.hide()
        else:
            self._view.ui.cb_ray_tracing.show()
            self._view.ui.label_26.show()

        if self._view.ui.cb_ray_tracing.isChecked():
            self._view.ui.label_10.show()
            self._view.ui.box_ray_trace_mode.show()
            self._view.ui.label_14.show()
            self._view.ui.box_ray_texture.show()
        else:
            self._view.ui.label_10.hide()
            self._view.ui.box_ray_trace_mode.hide()
            self._view.ui.label_14.hide()
            self._view.ui.box_ray_texture.hide()

        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self._view.ui.btn_image_page)
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 6, "Image")

    def display_new_page(self) -> None:
        """Displays the new project work area."""
        self._init_new_page()
        self._view.ui.list_new_projects.clear()
        # pre-process
        gui_elements_to_hide = [
            self._view.ui.lbl_new_choose_reference,
            self._view.ui.txt_new_choose_reference,
            self._view.ui.btn_new_choose_reference,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.status_bar.showMessage(self._workspace_label.text())
        tools.scan_workspace_for_valid_projects(self._workspace_path, self._view.ui.list_new_projects)
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 7, "Create new project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self._view.ui.btn_new_page)

    def display_open_page(self) -> None:
        """Displays the open project work area."""
        self._view.ui.txt_open_search.clear()
        self._view.ui.txt_open_selected_project.clear()
        if safeguard.Safeguard.check_filepath(self._workspace_path):
            self._view.ui.list_open_projects.clear()
            # pre-process
            self._view.status_bar.showMessage(self._workspace_label.text())
            try:
                tools.scan_workspace_for_valid_projects(self._workspace_path, self._view.ui.list_open_projects)
            except PermissionError:
                gui_utils.error_dialog_settings(
                    "The settings file is corrupted. Please restore the settings!",
                    "",
                    self._application_settings,
                )
                self.display_home_page()
                return
            tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 8, "Open existing project")
            self.last_sidebar_button = styles.color_sidebar_buttons(
                self.last_sidebar_button,
                self._view.ui.btn_open_page,
            )
        else:
            gui_utils.error_dialog_settings(
                "The settings file is corrupted. Please restore the settings!",
                "",
                self._application_settings,
            )
            self.display_home_page()

    def display_delete_page(self) -> None:
        """Displays the "delete" project work area."""
        self._view.ui.txt_delete_search.clear()
        self._view.ui.txt_delete_selected_projects.clear()
        self._view.ui.list_delete_projects.clear()
        # pre-process
        self._view.status_bar.showMessage(self._workspace_label.text())
        tools.scan_workspace_for_valid_projects(self._workspace_path, self._view.ui.list_delete_projects)
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 9, "Delete existing project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self._view.ui.btn_delete_page)

    def display_edit_page(self) -> None:
        """Displays the edit project page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        # pre-process
        self._view.status_bar.showMessage(self._workspace_label.text())
        self._init_edit_page()
        tools.switch_page(
            self._view.ui.stackedWidget,
            self._view.ui.lbl_page_title,
            13,
            "Edit proteins of current project",
        )
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self._view.ui.btn_edit_page)
        if len(self._view.ui.list_edit_project_proteins) >= 2:
            self._view.ui.btn_edit_existing_protein_struct.hide()
            self._view.ui.lbl_edit_existing_protein_struct.hide()
        else:
            self._view.ui.btn_edit_existing_protein_struct.show()
            self._view.ui.lbl_edit_existing_protein_struct.show()

    def display_hotspots_page(self) -> None:
        """Displays the hotspots page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        try:
            print(self._view.ui.list_hotspots_choose_protein.currentItem().text())
        except AttributeError:
            self._view.ui.list_hotspots_choose_protein.clear()
            gui_elements_to_hide = [
                self._view.ui.lbl_hotspots_resi_no,
                self._view.ui.sp_hotspots_resi_no,
                self._view.ui.lbl_hotspots_resi_show,
                self._view.ui.btn_hotspots_resi_show,
                self._view.ui.lbl_hotspots_resi_hide,
                self._view.ui.btn_hotspots_resi_hide,
                self._view.ui.lbl_hotspots_resi_zoom,
                self._view.ui.btn_hotspots_resi_zoom,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            gui_utils.fill_list_view_with_protein_names(
                self._current_project,
                self._view.ui.list_hotspots_choose_protein,
            )
            gui_utils.fill_list_view_with_protein_pair_names(
                self._current_project,
                self._view.ui.list_hotspots_choose_protein,
            )
        # self.project_scanner.scan_project_for_valid_proteins(self._view.ui.list_hotspots_choose_protein)
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 18, "Hotspots")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_hotspots_page,
        )

    def display_manage_pymol_session(self) -> None:
        """Displays the manage pymol session page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self._view.ui.box_manage_choose_protein.clear()
        pymol_objs = cmd.get_object_list()
        pymol_objs.insert(0, "")
        for tmp_object in pymol_objs:
            self._view.ui.box_manage_choose_protein.addItem(tmp_object)
        self._view.ui.box_manage_choose_protein.setCurrentIndex(
            self.pymol_session_specs[pyssa_keys.SESSION_SPEC_PROTEIN][0],
        )
        self._view.ui.box_manage_choose_color.setCurrentIndex(
            self.pymol_session_specs[pyssa_keys.SESSION_SPEC_COLOR][0],
        )
        self._view.ui.box_manage_choose_representation.setCurrentIndex(
            self.pymol_session_specs[pyssa_keys.SESSION_SPEC_REPRESENTATION][0],
        )
        self._view.ui.box_manage_choose_bg_color.setCurrentIndex(
            self.pymol_session_specs[pyssa_keys.SESSION_SPEC_BG_COLOR][0],
        )
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 24, "Manage PyMOL session")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_manage_session,
        )
        gui_elements_to_hide = [
            self._view.ui.lbl_manage_choose_color,
            self._view.ui.box_manage_choose_color,
            self._view.ui.lbl_manage_choose_representation,
            self._view.ui.box_manage_choose_representation,
            self._view.ui.lbl_manage_choose_bg_color,
            self._view.ui.box_manage_choose_bg_color,
            self._view.ui.lbl_disulfid_bond_1,
            self._view.ui.lbl_disulfid_bond_2,
            self._view.ui.btn_disulfid_bond_show,
            self._view.ui.btn_disulfid_bond_hide,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    # </editor-fold>

    # <editor-fold desc="New project">
    def show_add_reference(self) -> None:
        """Shows the reference input section."""
        # checkbox is checked
        self._view.ui.cb_new_add_reference.checkState()
        if self._view.ui.cb_new_add_reference.checkState() == 2:
            self._view.ui.txt_new_choose_reference.clear()
            self._view.ui.txt_new_choose_reference.setStyleSheet("background-color: white")
            self._view.ui.lbl_new_choose_reference.show()
            self._view.ui.txt_new_choose_reference.show()
            self._view.ui.btn_new_choose_reference.show()
            self._view.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_new_create_project)
            # check internet connectivity
            if not tools.check_internet_connectivity():
                gui_utils.no_internet_dialog()
                self._view.ui.txt_new_choose_reference.setEnabled(False)
                self._view.ui.lbl_new_status_choose_reference.setText("You cannot enter a PDB ID (no internet).")
                return
            self._view.ui.txt_new_choose_reference.setEnabled(True)
            self._view.ui.lbl_new_status_choose_reference.setText("")
        else:
            self._view.ui.lbl_new_choose_reference.hide()
            self._view.ui.txt_new_choose_reference.hide()
            self._view.ui.btn_new_choose_reference.hide()
            self._view.ui.lbl_new_status_choose_reference.setText("")
            self._view.ui.btn_new_create_project.setEnabled(True)

    def load_reference_in_project(self) -> None:
        """Loads a reference in a new project."""
        try:
            # open file dialog
            file_name = QtWidgets.QFileDialog.getOpenFileName(
                self._view,
                "Open Reference",
                QtCore.QDir.homePath(),
                "PDB Files (*.pdb)",
            )
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self._view.ui.txt_new_choose_reference.setText(str(file_name[0]))
            self._view.ui.txt_new_choose_reference.setEnabled(False)
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #000000")
            self._view.ui.btn_new_create_project.setEnabled(True)
        except ValueError:
            print("No file has been selected.")

    def validate_reference_in_project(self) -> None:
        """Checks if the entered reference protein is valid or not."""
        if len(self._view.ui.txt_new_choose_reference.text()) == 0:
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
            self._view.ui.lbl_new_status_choose_reference.setText("")
            self._view.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_new_create_project)
        elif len(self._view.ui.txt_new_choose_reference.text()) < 4:
            self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
            styles.color_button_not_ready(self._view.ui.btn_new_create_project)
            self._view.ui.btn_new_create_project.setEnabled(False)
            self._view.ui.lbl_new_status_choose_reference.setText("")
        # checks if a pdb id was entered
        elif len(self._view.ui.txt_new_choose_reference.text()) == 4:
            pdb_id = self._view.ui.txt_new_choose_reference.text().upper()
            try:
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                cmd.fetch(pdb_id, type="pdb", path=constants.SCRATCH_DIR)
                os.remove(f"{constants.SCRATCH_DIR}/{pdb_id}.pdb")
                cmd.reinitialize()
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #000000")
                self._view.ui.btn_new_create_project.setEnabled(True)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                return
            except FileNotFoundError:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.lbl_new_status_choose_reference.setText("Invalid PDB ID.")
                self._view.ui.btn_new_create_project.setEnabled(False)
                return
        else:
            if self._view.ui.txt_new_choose_reference.text().find("/") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_new_create_project)

            elif self._view.ui.txt_new_choose_reference.text().find("\\") == -1:
                self._view.ui.txt_new_choose_reference.setStyleSheet("color: #FC5457")
                self._view.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_new_create_project)

    def validate_project_name(self) -> None:
        """Validates the input of the project name in real-time."""
        input_validator.InputValidator.validate_project_name(
            self._view.ui.list_new_projects,
            self._view.ui.txt_new_project_name,
            self._view.ui.lbl_new_status_project_name,
            self._view.ui.btn_new_create_project,
            self._view.ui.cb_new_add_reference,
        )

    # </editor-fold>

    # <editor-fold desc="Open project">
    def validate_open_search(self) -> None:
        """Validates the input of the project name in real-time."""
        if self._view.ui.list_open_projects.currentItem() is not None:
            self._view.ui.list_open_projects.currentItem().setSelected(False)
        # set color for lineEdit
        input_validator.InputValidator.validate_search_input(
            self._view.ui.list_open_projects,
            self._view.ui.txt_open_search,
            self._view.ui.lbl_open_status_search,
            self._view.ui.txt_open_selected_project,
        )

    def select_project_from_open_list(self) -> None:
        """Sets the selected project name in the text box."""
        try:
            self._view.ui.txt_open_selected_project.setText(self._view.ui.list_open_projects.currentItem().text())
        except AttributeError:
            self._view.ui.txt_open_selected_project.setText("")

    def activate_open_button(self) -> None:
        """Activates the open button."""
        if self._view.ui.txt_open_selected_project.text() == "":
            self._view.ui.btn_open_open_project.setEnabled(False)
        else:
            self._view.ui.btn_open_open_project.setEnabled(True)

    # </editor-fold>

    # <editor-fold desc="Delete project">
    def select_project_from_delete_list(self) -> None:
        """Selects a project from the project list on the delete page."""
        try:
            self._view.ui.txt_delete_selected_projects.setText(self._view.ui.list_delete_projects.currentItem().text())
        except AttributeError:
            self._view.ui.txt_delete_selected_projects.setText("")

    def activate_delete_button(self) -> None:
        """Activates the delete button."""
        if self._view.ui.txt_delete_selected_projects.text() == "":
            self._view.ui.btn_delete_delete_project.setEnabled(False)
        else:
            self._view.ui.btn_delete_delete_project.setEnabled(True)

    def validate_delete_search(self) -> None:
        """Validates the input of the project name in real-time."""
        if self._view.ui.list_delete_projects.currentItem() is not None:
            self._view.ui.list_delete_projects.currentItem().setSelected(False)
        # set color for lineEdit
        input_validator.InputValidator.validate_search_input(
            self._view.ui.list_delete_projects,
            self._view.ui.txt_delete_search,
            self._view.ui.lbl_delete_status_search,
            self._view.ui.txt_delete_selected_projects,
        )

    # </editor-fold>

    # <editor-fold desc="View page">
    def display_view_page(self) -> None:
        """Displays the edit project page."""
        if self.is_distance_plot_open:
            self.distance_plot_dialog.close()
            self.is_distance_plot_open = False
        self._view.ui.list_view_project_proteins.clear()
        self._view.ui.txtedit_view_sequence.clear()
        # pre-process
        self._view.status_bar.showMessage(self._workspace_label.text())
        # list all proteins from pdb directory
        gui_utils.fill_list_view_with_protein_names(self._current_project, self._view.ui.list_view_project_proteins)

        tools.switch_page(
            self._view.ui.stackedWidget,
            self._view.ui.lbl_page_title,
            11,
            "View proteins of current project",
        )
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self._view.ui.btn_view_page)
        gui_elements_to_hide = [
            self._view.ui.btn_view_project_show,
            self._view.ui.btn_view_project_show_structure,
            self._view.ui.txtedit_view_sequence,
            self._view.ui.label_9,
            self._view.ui.label_11,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def view_show_options(self) -> None:
        """Controls which gui elements get shown."""
        gui_elements_to_show = [
            # self._view.ui.btn_view_project_show,
            self._view.ui.btn_view_project_show_structure,
            self._view.ui.txtedit_view_sequence,
            # self._view.ui.label_9,
            self._view.ui.label_11,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        self.view_sequence()

    def view_sequence(self) -> None:
        """Displays the sequence of the selected protein in a text box."""
        tmp_protein_basename = self._view.ui.list_view_project_proteins.currentItem().text()
        tmp_protein_sequences = self._current_project.search_protein(tmp_protein_basename).get_protein_sequences()
        self._view.ui.txtedit_view_sequence.clear()
        for tmp_sequence in tmp_protein_sequences:
            self._view.ui.txtedit_view_sequence.append("".join(tmp_sequence.sequence))
        # fixme: experimental sequence viewer gui
        # dialog = dialog_sequence_viewer.SequenceViewer(tmp_protein_sequences, tmp_protein_filename)
        # dialog.exec_()

    # </editor-fold>

    # <editor-fold desc="Use page">
    def validate_use_project_name(self) -> None:
        """Validates the input of the project name in real-time."""
        input_validator.InputValidator.validate_project_name(
            self._view.ui.list_use_existing_projects,
            self._view.ui.txt_use_project_name,
            self._view.ui.lbl_use_status_project_name,
            self._view.ui.btn_use_next,
        )

    def validate_use_search(self) -> None:
        """Validates the input of the protein name in real-time."""
        message = "Protein structure does not exists."
        input_validator.InputValidator.validate_search_input(
            self._view.ui.list_use_available_protein_structures,
            self._view.ui.txt_use_search,
            self._view.ui.lbl_use_status_search,
            status_message=message,
        )

    def add_protein_structure_to_new_project(self) -> None:
        """Adds the selected protein to the list which is used to create the new project."""
        prot_to_add = self._view.ui.list_use_available_protein_structures.currentItem().text()
        self._view.ui.list_use_selected_protein_structures.addItem(prot_to_add)
        self._view.ui.list_use_available_protein_structures.takeItem(
            self._view.ui.list_use_available_protein_structures.currentRow(),
        )
        self._view.ui.btn_use_add_available_protein_structures.setEnabled(False)
        if self._view.ui.list_use_available_protein_structures.count() > 0:
            try:
                self._view.ui.list_use_available_protein_structures.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in use available proteins list on Use page.")

    def remove_protein_structure_to_new_project(self) -> None:
        """Removes the selected protein from the list which is used to create the new project."""
        prot_to_remove = self._view.ui.list_use_selected_protein_structures.currentItem()
        self._view.ui.list_use_selected_protein_structures.takeItem(
            self._view.ui.list_use_selected_protein_structures.currentRow(),
        )
        self._view.ui.list_use_available_protein_structures.addItem(prot_to_remove)
        self._view.ui.btn_use_remove_selected_protein_structures.setEnabled(False)
        if self._view.ui.list_use_selected_protein_structures.count() > 0:
            try:
                self._view.ui.list_use_selected_protein_structures.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in use selected proteins list on Use page.")

    def show_protein_selection_for_use(self) -> None:
        """Shows the two lists for the protein selection."""
        gui_elements_to_show = [
            self._view.ui.lbl_use_search,
            self._view.ui.lbl_use_status_search,
            self._view.ui.txt_use_search,
            self._view.ui.btn_use_add_available_protein_structures,
            self._view.ui.lbl_use_available_protein_structures,
            self._view.ui.list_use_available_protein_structures,
            self._view.ui.btn_use_remove_selected_protein_structures,
            self._view.ui.lbl_use_selected_protein_structures,
            self._view.ui.list_use_selected_protein_structures,
            self._view.ui.btn_use_back,
            self._view.ui.btn_use_create_new_project,
            self._view.ui.lbl_use_new_project,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        self._view.ui.txt_use_project_name.setEnabled(False)
        gui_elements_to_hide = [
            self._view.ui.btn_use_next,
            self._view.ui.list_use_existing_projects,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.disable_text_box(self._view.ui.txt_use_project_name, self._view.ui.lbl_use_project_name)
        self._view.ui.btn_use_add_available_protein_structures.setEnabled(False)
        self._view.ui.btn_use_remove_selected_protein_structures.setEnabled(False)

    def hide_protein_selection_for_use(self) -> None:
        """Hides the two lists for the protein selection."""
        gui_elements_to_show = [
            self._view.ui.btn_use_next,
            self._view.ui.list_use_existing_projects,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        self._view.ui.txt_use_project_name.setEnabled(True)

        gui_elements_to_hide = [
            self._view.ui.lbl_use_search,
            self._view.ui.lbl_use_status_search,
            self._view.ui.txt_use_search,
            self._view.ui.btn_use_add_available_protein_structures,
            self._view.ui.lbl_use_available_protein_structures,
            self._view.ui.list_use_available_protein_structures,
            self._view.ui.btn_use_remove_selected_protein_structures,
            self._view.ui.lbl_use_selected_protein_structures,
            self._view.ui.list_use_selected_protein_structures,
            self._view.ui.btn_use_back,
            self._view.ui.btn_use_create_new_project,
            self._view.ui.lbl_use_new_project,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self._view.ui.txt_use_project_name, self._view.ui.lbl_use_project_name)

    def use_enable_add(self) -> None:
        """Enables the add button."""
        self._view.ui.btn_use_add_available_protein_structures.setEnabled(True)

    def use_enable_remove(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_use_remove_selected_protein_structures.setEnabled(True)

    # </editor-fold>

    # <editor-fold desc="ESM fold">
    def _init_esm_pred_mono_page(self) -> None:
        """Clears all text boxes and sets up the default values for the page."""
        # clears everything
        self._view.ui.txt_esm_prot_name.clear()
        self._view.ui.txt_esm_prot_seq.clear()
        for i in range(self._view.ui.table_esm_prot_to_predict.rowCount()):
            self._view.ui.table_esm_prot_to_predict.removeRow(i)
        # sets up defaults: Prediction
        self._view.ui.btn_esm_next.setEnabled(False)
        self._view.ui.btn_esm_next_2.setEnabled(False)
        self._view.ui.lbl_esm_prot_name_status.setText("")
        self._view.ui.lbl_esm_prot_seq_status.setText("")

    def display_esm_pred_mono(self) -> None:
        """Displays the esm_fold monomer page."""
        self._init_esm_pred_mono_page()
        gui_elements_to_show = [
            self._view.ui.lbl_esm_prot_to_predict,
            self._view.ui.table_esm_prot_to_predict,
            self._view.ui.btn_esm_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_esm_seq_to_predict_remove,
            self._view.ui.lbl_esm_prot_name,
            self._view.ui.txt_esm_prot_name,
            self._view.ui.lbl_esm_prot_name_status,
            self._view.ui.btn_esm_back,
            self._view.ui.btn_esm_next,
            self._view.ui.lbl_esm_prot_seq,
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_back_2,
            self._view.ui.btn_esm_next_2,
            self._view.ui.btn_esm_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        styles.color_button_not_ready(self._view.ui.btn_esm_next)
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 1, "ESMFold Monomer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_pred_cloud_monomer_page,
        )

    def cloud_esm_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self._view.ui.txt_esm_prot_name.text(),
            self._view.ui.table_esm_prot_to_predict,
        ):
            self._view.ui.lbl_esm_prot_name_status.setText("Protein name already used.")
            self._view.ui.btn_esm_next.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_esm_next)
        else:
            self._view.ui.btn_esm_next.setEnabled(True)
            tools.validate_protein_name(
                self._view.ui.txt_esm_prot_name,
                self._view.ui.lbl_esm_prot_name_status,
                self._view.ui.btn_esm_next,
            )

    def cloud_esm_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_next_2,
        )

    def setup_defaults_esm_monomer_prediction(self) -> None:
        """Sets up the default values for the page."""
        # clears everything
        self._view.ui.txt_esm_prot_name.clear()
        self._view.ui.txt_esm_prot_seq.clear()
        # sets up defaults: Prediction
        self._view.ui.btn_esm_next.setEnabled(False)
        self._view.ui.btn_esm_next_2.setEnabled(False)
        self._view.ui.lbl_esm_prot_name_status.setText("")
        self._view.ui.lbl_esm_prot_seq_status.setText("")

    def cloud_esm_add_seq_to_predict(self) -> None:
        """Shows the gui elements to add a sequence to the protein to predict."""
        gui_elements_to_show = [
            self._view.ui.lbl_esm_prot_to_predict,
            self._view.ui.table_esm_prot_to_predict,
            self._view.ui.lbl_esm_prot_name,
            self._view.ui.txt_esm_prot_name,
            self._view.ui.lbl_esm_prot_name_status,
            self._view.ui.btn_esm_back,
            self._view.ui.btn_esm_next,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_esm_prot_name, self._view.ui.lbl_esm_prot_name)
        gui_elements_to_hide = [
            self._view.ui.btn_esm_seq_to_predict_remove,
            self._view.ui.btn_esm_seq_to_predict,
            self._view.ui.lbl_esm_prot_seq,
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_back_2,
            self._view.ui.btn_esm_next_2,
            self._view.ui.btn_esm_predict,
        ]
        gui_utils.disable_text_box(self._view.ui.txt_esm_prot_seq, self._view.ui.lbl_esm_prot_seq)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_esm_next.setEnabled(False)
        self._view.ui.txt_esm_prot_name.clear()
        styles.color_button_not_ready(self._view.ui.btn_esm_next)
        if self._view.ui.table_esm_prot_to_predict.rowCount() > 0:
            try:
                self._view.ui.table_esm_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Monomer Prediction in overview table.")

    def cloud_esm_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_esm_prot_to_predict,
            self._view.ui.table_esm_prot_to_predict,
            self._view.ui.btn_esm_seq_to_predict_remove,
            self._view.ui.btn_esm_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.lbl_esm_prot_name,
            self._view.ui.txt_esm_prot_name,
            self._view.ui.lbl_esm_prot_name_status,
            self._view.ui.btn_esm_back,
            self._view.ui.btn_esm_next,
            self._view.ui.lbl_esm_prot_seq,
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_back_2,
            self._view.ui.btn_esm_next_2,
            self._view.ui.btn_esm_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.cloud_esm_check_if_table_is_empty()
        self._view.ui.btn_esm_seq_to_predict_remove.setEnabled(False)

    def cloud_esm_next(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_esm_prot_to_predict,
            self._view.ui.table_esm_prot_to_predict,
            self._view.ui.lbl_esm_prot_name,
            self._view.ui.txt_esm_prot_name,
            self._view.ui.lbl_esm_prot_seq,
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_back_2,
            self._view.ui.btn_esm_next_2,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_esm_prot_seq, self._view.ui.lbl_esm_prot_seq)
        gui_elements_to_hide = [
            self._view.ui.btn_esm_seq_to_predict_remove,
            self._view.ui.btn_esm_seq_to_predict,
            self._view.ui.lbl_esm_prot_name_status,
            self._view.ui.btn_esm_back,
            self._view.ui.btn_esm_next,
            self._view.ui.btn_esm_predict,
        ]
        gui_utils.disable_text_box(self._view.ui.txt_esm_prot_name, self._view.ui.lbl_esm_prot_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.txt_esm_prot_seq.clear()

    def cloud_esm_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_esm_prot_to_predict,
            self._view.ui.table_esm_prot_to_predict,
            self._view.ui.lbl_esm_prot_name,
            self._view.ui.txt_esm_prot_name,
            self._view.ui.lbl_esm_prot_name_status,
            self._view.ui.btn_esm_back,
            self._view.ui.btn_esm_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_esm_seq_to_predict_remove,
            self._view.ui.btn_esm_seq_to_predict,
            self._view.ui.lbl_esm_prot_seq,
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_back_2,
            self._view.ui.btn_esm_next_2,
            self._view.ui.btn_esm_predict,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_esm_prot_name, self._view.ui.lbl_esm_prot_name)
        gui_utils.disable_text_box(self._view.ui.txt_esm_prot_seq, self._view.ui.lbl_esm_prot_seq)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def cloud_esm_add_protein(self) -> None:
        """Adds protein to the list of proteins to predict."""
        self._view.ui.table_esm_prot_to_predict.setRowCount(self._view.ui.table_esm_prot_to_predict.rowCount() + 1)
        self._view.ui.table_esm_prot_to_predict.insertRow(self._view.ui.table_esm_prot_to_predict.rowCount() + 1)
        self._view.ui.table_esm_prot_to_predict.setItem(
            self._view.ui.table_esm_prot_to_predict.rowCount() - 1,
            0,
            QtWidgets.QTableWidgetItem("A"),
        )
        self._view.ui.table_esm_prot_to_predict.setItem(
            self._view.ui.table_esm_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_esm_prot_seq.toPlainText()),
        )
        self._view.ui.table_esm_prot_to_predict.setVerticalHeaderItem(
            self._view.ui.table_esm_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_esm_prot_name.text()),
        )
        self._view.ui.table_esm_prot_to_predict.resizeColumnsToContents()
        self.cloud_esm_check_if_table_is_empty()
        gui_elements_to_show = [
            self._view.ui.lbl_esm_prot_to_predict,
            self._view.ui.table_esm_prot_to_predict,
            self._view.ui.btn_esm_seq_to_predict_remove,
            self._view.ui.btn_esm_seq_to_predict,
            self._view.ui.btn_esm_predict,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_esm_prot_name, self._view.ui.lbl_esm_prot_name)
        gui_elements_to_hide = [
            self._view.ui.lbl_esm_prot_name,
            self._view.ui.txt_esm_prot_name,
            self._view.ui.lbl_esm_prot_name_status,
            self._view.ui.btn_esm_back,
            self._view.ui.btn_esm_next,
            self._view.ui.lbl_esm_prot_seq,
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_back_2,
            self._view.ui.btn_esm_next_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_esm_predict.setEnabled(True)
        self._view.ui.btn_esm_seq_to_predict_remove.setEnabled(False)
        self.setup_defaults_esm_monomer_prediction()

    def cloud_esm_remove(self) -> None:
        """Removes a protein from the list of proteins to predict."""
        self._view.ui.table_esm_prot_to_predict.removeRow(self._view.ui.table_esm_prot_to_predict.currentRow())
        gui_elements_to_show = [
            self._view.ui.lbl_esm_prot_to_predict,
            self._view.ui.table_esm_prot_to_predict,
            self._view.ui.btn_esm_seq_to_predict_remove,
            self._view.ui.btn_esm_seq_to_predict,
            self._view.ui.btn_esm_predict,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_esm_prot_name, self._view.ui.lbl_esm_prot_name)
        gui_elements_to_hide = [
            self._view.ui.lbl_esm_prot_name,
            self._view.ui.txt_esm_prot_name,
            self._view.ui.lbl_esm_prot_name_status,
            self._view.ui.btn_esm_back,
            self._view.ui.btn_esm_next,
            self._view.ui.lbl_esm_prot_seq,
            self._view.ui.txt_esm_prot_seq,
            self._view.ui.lbl_esm_prot_seq_status,
            self._view.ui.btn_esm_back_2,
            self._view.ui.btn_esm_next_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_esm_seq_to_predict_remove.setEnabled(False)
        self.cloud_esm_check_if_table_is_empty()

    def cloud_esm_item_changed(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_esm_seq_to_predict_remove.setEnabled(True)

    def cloud_esm_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self._view.ui.table_esm_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self._view.ui.btn_esm_predict)
            self._view.ui.btn_esm_predict.setEnabled(False)
            gui_elements_to_show = [
                self._view.ui.lbl_esm_prot_to_predict,
                self._view.ui.table_esm_prot_to_predict,
                self._view.ui.btn_esm_seq_to_predict,
            ]
            gui_utils.enable_text_box(self._view.ui.txt_esm_prot_name, self._view.ui.lbl_esm_prot_name)
            gui_elements_to_hide = [
                self._view.ui.btn_esm_seq_to_predict_remove,
                self._view.ui.lbl_esm_prot_name,
                self._view.ui.txt_esm_prot_name,
                self._view.ui.lbl_esm_prot_name_status,
                self._view.ui.btn_esm_back,
                self._view.ui.btn_esm_next,
                self._view.ui.lbl_esm_prot_seq,
                self._view.ui.txt_esm_prot_seq,
                self._view.ui.lbl_esm_prot_seq_status,
                self._view.ui.btn_esm_back_2,
                self._view.ui.btn_esm_next_2,
                self._view.ui.btn_esm_predict,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            self._view.ui.btn_esm_predict.setEnabled(True)

    # </editor-fold>

    # <editor-fold desc="Monomer local prediction">
    def _init_local_pred_mono_page(self) -> None:
        """Clears all text boxes and sets default values for the gui elements."""
        # clears everything
        self._view.ui.txt_pred_mono_prot_name.clear()
        self._view.ui.txt_pred_mono_seq_name.clear()
        for i in range(self._view.ui.table_pred_mono_prot_to_predict.rowCount()):
            self._view.ui.table_pred_mono_prot_to_predict.removeRow(i)
        # sets up defaults: Prediction
        self._view.ui.btn_pred_mono_next.setEnabled(False)
        self._view.ui.btn_pred_mono_add_protein.setEnabled(False)
        self._view.ui.lbl_pred_mono_prot_name_status.setText("")
        self._view.ui.lbl_pred_mono_seq_name_status.setText("")

    def display_local_pred_mono(self) -> None:
        """Displays the local prediction monomer page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return
        self._init_local_pred_mono_page()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_mono_prot_to_predict,
            self._view.ui.table_pred_mono_prot_to_predict,
            self._view.ui.btn_pred_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_mono_seq_to_predict_remove,
            self._view.ui.lbl_pred_mono_prot_name,
            self._view.ui.txt_pred_mono_prot_name,
            self._view.ui.lbl_pred_mono_prot_name_status,
            self._view.ui.btn_pred_mono_back,
            self._view.ui.btn_pred_mono_next,
            self._view.ui.lbl_pred_mono_seq_name,
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_back_2,
            self._view.ui.btn_pred_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        styles.color_button_not_ready(self._view.ui.btn_pred_mono_next)
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 19, "Local Monomer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_pred_local_monomer_page,
        )

    def local_pred_mono_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self._view.ui.txt_pred_mono_prot_name.text(),
            self._view.ui.table_pred_mono_prot_to_predict,
        ):
            self._view.ui.lbl_pred_mono_prot_name_status.setText("Protein name already used.")
            self._view.ui.btn_pred_mono_next.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_pred_mono_next)
        else:
            self._view.ui.btn_pred_mono_next.setEnabled(True)
            tools.validate_protein_name(
                self._view.ui.txt_pred_mono_prot_name,
                self._view.ui.lbl_pred_mono_prot_name_status,
                self._view.ui.btn_pred_mono_next,
            )

    def local_pred_mono_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_add_protein,
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
        self._view.ui.txt_pred_mono_prot_name.clear()
        self._view.ui.txt_pred_mono_seq_name.clear()
        # sets up defaults: Prediction
        self._view.ui.btn_pred_mono_next.setEnabled(False)
        self._view.ui.btn_pred_mono_add_protein.setEnabled(False)
        self._view.ui.lbl_pred_mono_prot_name_status.setText("")
        self._view.ui.lbl_pred_mono_seq_name_status.setText("")

    def local_pred_mono_add_seq_to_predict(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_mono_prot_to_predict,
            self._view.ui.table_pred_mono_prot_to_predict,
            self._view.ui.lbl_pred_mono_prot_name,
            self._view.ui.txt_pred_mono_prot_name,
            self._view.ui.lbl_pred_mono_prot_name_status,
            self._view.ui.btn_pred_mono_back,
            self._view.ui.btn_pred_mono_next,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_pred_mono_prot_name, self._view.ui.lbl_pred_mono_prot_name)
        gui_elements_to_hide = [
            self._view.ui.btn_pred_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_mono_seq_to_predict,
            self._view.ui.lbl_pred_mono_seq_name,
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_back_2,
            self._view.ui.btn_pred_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_predict,
        ]
        gui_utils.disable_text_box(self._view.ui.txt_pred_mono_seq_name, self._view.ui.lbl_pred_mono_seq_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_mono_next.setEnabled(False)
        self._view.ui.txt_pred_mono_prot_name.clear()
        styles.color_button_not_ready(self._view.ui.btn_pred_mono_next)
        if self._view.ui.table_pred_mono_prot_to_predict.rowCount() > 0:
            try:
                self._view.ui.table_pred_mono_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Monomer Prediction in overview table.")

    def local_pred_mono_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_mono_prot_to_predict,
            self._view.ui.table_pred_mono_prot_to_predict,
            self._view.ui.btn_pred_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_mono_prot_name,
            self._view.ui.txt_pred_mono_prot_name,
            self._view.ui.lbl_pred_mono_prot_name_status,
            self._view.ui.btn_pred_mono_back,
            self._view.ui.btn_pred_mono_next,
            self._view.ui.lbl_pred_mono_seq_name,
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_back_2,
            self._view.ui.btn_pred_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.local_pred_mono_check_if_table_is_empty()
        self._view.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(False)

    def local_pred_mono_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_mono_prot_to_predict,
            self._view.ui.table_pred_mono_prot_to_predict,
            self._view.ui.lbl_pred_mono_prot_name,
            self._view.ui.txt_pred_mono_prot_name,
            self._view.ui.lbl_pred_mono_seq_name,
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_back_2,
            self._view.ui.btn_pred_mono_add_protein,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_pred_mono_seq_name, self._view.ui.lbl_pred_mono_seq_name)
        gui_elements_to_hide = [
            self._view.ui.btn_pred_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_mono_seq_to_predict,
            self._view.ui.lbl_pred_mono_prot_name_status,
            self._view.ui.btn_pred_mono_back,
            self._view.ui.btn_pred_mono_next,
            self._view.ui.lbl_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_predict,
        ]
        gui_utils.disable_text_box(self._view.ui.txt_pred_mono_prot_name, self._view.ui.lbl_pred_mono_prot_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.txt_pred_mono_seq_name.clear()

    def local_pred_mono_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_mono_prot_to_predict,
            self._view.ui.table_pred_mono_prot_to_predict,
            self._view.ui.lbl_pred_mono_prot_name,
            self._view.ui.txt_pred_mono_prot_name,
            self._view.ui.lbl_pred_mono_prot_name_status,
            self._view.ui.btn_pred_mono_back,
            self._view.ui.btn_pred_mono_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_mono_seq_to_predict,
            self._view.ui.lbl_pred_mono_seq_name,
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_back_2,
            self._view.ui.btn_pred_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_predict,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_pred_mono_prot_name, self._view.ui.lbl_pred_mono_prot_name)
        gui_utils.disable_text_box(self._view.ui.txt_pred_mono_seq_name, self._view.ui.lbl_pred_mono_seq_name)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def local_pred_mono_add_protein(self) -> None:
        """Adds protein to the list of proteins to predict."""
        self._view.ui.table_pred_mono_prot_to_predict.setRowCount(
            self._view.ui.table_pred_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_mono_prot_to_predict.insertRow(
            self._view.ui.table_pred_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_mono_prot_to_predict.rowCount() - 1,
            0,
            QtWidgets.QTableWidgetItem("A"),
        )
        self._view.ui.table_pred_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_mono_seq_name.toPlainText()),
        )
        self._view.ui.table_pred_mono_prot_to_predict.setVerticalHeaderItem(
            self._view.ui.table_pred_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_mono_prot_name.text()),
        )
        self._view.ui.table_pred_mono_prot_to_predict.resizeColumnsToContents()
        self.local_pred_mono_check_if_table_is_empty()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_mono_prot_to_predict,
            self._view.ui.table_pred_mono_prot_to_predict,
            self._view.ui.btn_pred_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_mono_seq_to_predict,
            self._view.ui.lbl_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_predict,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_pred_mono_prot_name, self._view.ui.lbl_pred_mono_prot_name)
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_mono_prot_name,
            self._view.ui.txt_pred_mono_prot_name,
            self._view.ui.lbl_pred_mono_prot_name_status,
            self._view.ui.btn_pred_mono_back,
            self._view.ui.btn_pred_mono_next,
            self._view.ui.lbl_pred_mono_seq_name,
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_back_2,
            self._view.ui.btn_pred_mono_add_protein,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_mono_predict.setEnabled(True)
        self._view.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(False)
        self.setup_defaults_monomer_prediction()

    def local_pred_mono_remove(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        self._view.ui.table_pred_mono_prot_to_predict.removeRow(
            self._view.ui.table_pred_mono_prot_to_predict.currentRow(),
        )
        gui_elements_to_show = [
            self._view.ui.lbl_pred_mono_prot_to_predict,
            self._view.ui.table_pred_mono_prot_to_predict,
            self._view.ui.btn_pred_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_mono_seq_to_predict,
            self._view.ui.lbl_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_advanced_config,
            self._view.ui.btn_pred_mono_predict,
        ]
        gui_utils.enable_text_box(self._view.ui.txt_pred_mono_prot_name, self._view.ui.lbl_pred_mono_prot_name)
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_mono_prot_name,
            self._view.ui.txt_pred_mono_prot_name,
            self._view.ui.lbl_pred_mono_prot_name_status,
            self._view.ui.btn_pred_mono_back,
            self._view.ui.btn_pred_mono_next,
            self._view.ui.lbl_pred_mono_seq_name,
            self._view.ui.txt_pred_mono_seq_name,
            self._view.ui.lbl_pred_mono_seq_name_status,
            self._view.ui.btn_pred_mono_back_2,
            self._view.ui.btn_pred_mono_add_protein,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(False)
        self.local_pred_mono_check_if_table_is_empty()

    def local_pred_mono_item_changed(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_pred_mono_seq_to_predict_remove.setEnabled(True)

    def local_pred_mono_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self._view.ui.table_pred_mono_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self._view.ui.btn_pred_mono_predict)
            self._view.ui.btn_pred_mono_predict.setEnabled(False)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_mono_prot_to_predict,
                self._view.ui.table_pred_mono_prot_to_predict,
                self._view.ui.btn_pred_mono_seq_to_predict,
            ]
            gui_utils.enable_text_box(self._view.ui.txt_pred_mono_prot_name, self._view.ui.lbl_pred_mono_prot_name)
            gui_elements_to_hide = [
                self._view.ui.btn_pred_mono_seq_to_predict_remove,
                self._view.ui.lbl_pred_mono_prot_name,
                self._view.ui.txt_pred_mono_prot_name,
                self._view.ui.lbl_pred_mono_prot_name_status,
                self._view.ui.btn_pred_mono_back,
                self._view.ui.btn_pred_mono_next,
                self._view.ui.lbl_pred_mono_seq_name,
                self._view.ui.txt_pred_mono_seq_name,
                self._view.ui.lbl_pred_mono_seq_name_status,
                self._view.ui.btn_pred_mono_back_2,
                self._view.ui.btn_pred_mono_add_protein,
                self._view.ui.lbl_pred_mono_advanced_config,
                self._view.ui.btn_pred_mono_advanced_config,
                self._view.ui.btn_pred_mono_predict,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            self._view.ui.btn_pred_mono_predict.setEnabled(True)

    # </editor-fold>

    # <editor-fold desc="Multimer local prediction">
    def _init_local_pred_multi_page(self) -> None:
        """Clears all text boxes and sets default values for the gui elements."""
        # clears everything
        self._view.ui.txt_pred_multi_prot_name.clear()
        self._view.ui.txt_pred_multi_prot_seq.clear()
        self._view.ui.list_pred_multi_prot_seq_overview.clear()
        # sets up defaults: Prediction
        self._view.ui.btn_pred_multi_next.setEnabled(False)
        self._view.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        self._view.ui.lbl_pred_multi_prot_name_status.setText("")
        self._view.ui.lbl_pred_multi_prot_seq_status.setText("")

    def display_local_pred_multi(self) -> None:
        """Displays the local prediction multimer page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return

        gui_elements_to_show = [
            self._view.ui.lbl_pred_multi_prot_to_predict,
            self._view.ui.table_pred_multi_prot_to_predict,
            self._view.ui.btn_pred_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_multi_prot_to_predict_remove,
            self._view.ui.lbl_pred_multi_prot_name_status,
            self._view.ui.btn_pred_multi_back,
            self._view.ui.btn_pred_multi_next,
            self._view.ui.lbl_pred_multi_prot_name,
            self._view.ui.txt_pred_multi_prot_name,
            self._view.ui.lbl_pred_multi_prot_seq,
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add,
            self._view.ui.btn_pred_multi_prot_seq_add,
            self._view.ui.lbl_pred_multi_prot_seq_overview,
            self._view.ui.list_pred_multi_prot_seq_overview,
            self._view.ui.btn_pred_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_multi_prot_to_predict_2,
            self._view.ui.btn_pred_multi_back_2,
            self._view.ui.btn_pred_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_predict,
        ]
        for i in range(self._view.ui.table_pred_multi_prot_to_predict.rowCount() - 1, -1, -1):
            self._view.ui.table_pred_multi_prot_to_predict.removeRow(i)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 20, "Local Multimer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_pred_local_multimer_page,
        )

    def local_pred_multi_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self._view.ui.txt_pred_multi_prot_name.text(),
            self._view.ui.table_pred_multi_prot_to_predict,
        ):
            self._view.ui.lbl_pred_multi_prot_name_status.setText("Protein name already used.")
            self._view.ui.btn_pred_multi_next.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_pred_multi_next)
        else:
            self._view.ui.btn_pred_multi_next.setEnabled(True)
            tools.validate_protein_name(
                self._view.ui.txt_pred_multi_prot_name,
                self._view.ui.lbl_pred_multi_prot_name_status,
                self._view.ui.btn_pred_multi_next,
            )

    def local_pred_multi_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.btn_pred_multi_prot_seq_add,
        )

    def local_pred_multi_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self._view.ui.table_pred_multi_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self._view.ui.btn_pred_multi_predict)
            self._view.ui.btn_pred_multi_predict.setEnabled(False)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_multi_prot_to_predict,
                self._view.ui.table_pred_multi_prot_to_predict,
                self._view.ui.btn_pred_multi_prot_to_predict_add,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_multi_prot_to_predict_remove,
                self._view.ui.lbl_pred_multi_prot_name_status,
                self._view.ui.btn_pred_multi_back,
                self._view.ui.btn_pred_multi_next,
                self._view.ui.lbl_pred_multi_prot_name,
                self._view.ui.txt_pred_multi_prot_name,
                self._view.ui.lbl_pred_multi_prot_seq,
                self._view.ui.txt_pred_multi_prot_seq,
                self._view.ui.lbl_pred_multi_prot_seq_status,
                self._view.ui.lbl_pred_multi_prot_seq_add,
                self._view.ui.btn_pred_multi_prot_seq_add,
                self._view.ui.lbl_pred_multi_prot_seq_overview,
                self._view.ui.list_pred_multi_prot_seq_overview,
                self._view.ui.btn_pred_multi_prot_seq_overview_remove,
                self._view.ui.lbl_pred_multi_prot_to_predict_2,
                self._view.ui.btn_pred_multi_back_2,
                self._view.ui.btn_pred_multi_prot_to_predict_add_2,
                self._view.ui.lbl_pred_multi_advanced_config,
                self._view.ui.btn_pred_multi_advanced_config,
                self._view.ui.btn_pred_multi_predict,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        else:
            self._view.ui.btn_pred_multi_predict.setEnabled(True)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_multi_prot_to_predict,
                self._view.ui.table_pred_multi_prot_to_predict,
                self._view.ui.btn_pred_multi_prot_to_predict_add,
                self._view.ui.btn_pred_multi_prot_to_predict_remove,
                self._view.ui.lbl_pred_multi_advanced_config,
                self._view.ui.btn_pred_multi_advanced_config,
                self._view.ui.btn_pred_multi_predict,
            ]
            gui_elements_to_hide = [
                self._view.ui.lbl_pred_multi_prot_name_status,
                self._view.ui.btn_pred_multi_back,
                self._view.ui.btn_pred_multi_next,
                self._view.ui.lbl_pred_multi_prot_name,
                self._view.ui.txt_pred_multi_prot_name,
                self._view.ui.lbl_pred_multi_prot_seq,
                self._view.ui.txt_pred_multi_prot_seq,
                self._view.ui.lbl_pred_multi_prot_seq_status,
                self._view.ui.lbl_pred_multi_prot_seq_add,
                self._view.ui.btn_pred_multi_prot_seq_add,
                self._view.ui.lbl_pred_multi_prot_seq_overview,
                self._view.ui.list_pred_multi_prot_seq_overview,
                self._view.ui.btn_pred_multi_prot_seq_overview_remove,
                self._view.ui.lbl_pred_multi_prot_to_predict_2,
                self._view.ui.btn_pred_multi_back_2,
                self._view.ui.btn_pred_multi_prot_to_predict_add_2,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)

    def local_pred_multi_add_sequence_to_list(self) -> None:
        """Adds the entered sequence to the list of sequences of the protein."""
        self._view.ui.list_pred_multi_prot_seq_overview.addItem(
            QtWidgets.QListWidgetItem(self._view.ui.txt_pred_multi_prot_seq.toPlainText()),
        )
        self.local_pred_multi_check_if_list_is_empty()

    def local_pred_multi_remove_sequence_to_list(self) -> None:
        """Removes the entered sequence to the list of sequences of the protein."""
        self._view.ui.list_pred_multi_prot_seq_overview.takeItem(
            self._view.ui.list_pred_multi_prot_seq_overview.currentRow(),
        )
        self.local_pred_multi_check_if_list_is_empty()
        self._view.ui.btn_pred_multi_prot_seq_overview_remove.setEnabled(False)

    def local_pred_multi_check_if_list_is_empty(self) -> None:
        """Checks if the list of sequences of the protein is empty."""
        if self._view.ui.list_pred_multi_prot_seq_overview.count() == 0:
            self._view.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        else:
            self._view.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(True)

    def local_pred_multi_add(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_multi_prot_to_predict,
            self._view.ui.table_pred_multi_prot_to_predict,
            self._view.ui.lbl_pred_multi_prot_name,
            self._view.ui.txt_pred_multi_prot_name,
            self._view.ui.lbl_pred_multi_prot_name_status,
            self._view.ui.btn_pred_multi_back,
            self._view.ui.btn_pred_multi_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_multi_prot_seq,
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add,
            self._view.ui.btn_pred_multi_prot_seq_add,
            self._view.ui.lbl_pred_multi_prot_seq_overview,
            self._view.ui.list_pred_multi_prot_seq_overview,
            self._view.ui.btn_pred_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_multi_prot_to_predict_2,
            self._view.ui.btn_pred_multi_back_2,
            self._view.ui.btn_pred_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self._view.ui.txt_pred_multi_prot_name, self._view.ui.lbl_pred_multi_prot_name)
        gui_utils.disable_text_box(self._view.ui.txt_pred_multi_prot_seq, self._view.ui.lbl_pred_multi_prot_seq)
        self._view.ui.btn_pred_multi_next.setEnabled(False)
        self._view.ui.txt_pred_multi_prot_name.clear()
        styles.color_button_not_ready(self._view.ui.btn_pred_multi_next)
        if self._view.ui.table_pred_multi_prot_to_predict.rowCount() > 0:
            try:
                self._view.ui.table_pred_multi_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Multimer Prediction in overview table.")

    def local_pred_multi_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_multi_prot_to_predict,
            self._view.ui.table_pred_multi_prot_to_predict,
            self._view.ui.btn_pred_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_multi_prot_to_predict_remove,
            self._view.ui.lbl_pred_multi_prot_name,
            self._view.ui.txt_pred_multi_prot_name,
            self._view.ui.lbl_pred_multi_prot_name_status,
            self._view.ui.btn_pred_multi_back,
            self._view.ui.btn_pred_multi_next,
            self._view.ui.lbl_pred_multi_prot_seq,
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add,
            self._view.ui.btn_pred_multi_prot_seq_add,
            self._view.ui.lbl_pred_multi_prot_seq_overview,
            self._view.ui.list_pred_multi_prot_seq_overview,
            self._view.ui.btn_pred_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_multi_prot_to_predict_2,
            self._view.ui.btn_pred_multi_back_2,
            self._view.ui.btn_pred_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.local_pred_multi_check_if_table_is_empty()

    def local_pred_multi_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_multi_prot_to_predict,
            self._view.ui.table_pred_multi_prot_to_predict,
            self._view.ui.lbl_pred_multi_prot_name,
            self._view.ui.txt_pred_multi_prot_name,
            self._view.ui.lbl_pred_multi_prot_seq,
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add,
            self._view.ui.btn_pred_multi_prot_seq_add,
            self._view.ui.lbl_pred_multi_prot_seq_overview,
            self._view.ui.list_pred_multi_prot_seq_overview,
            self._view.ui.btn_pred_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_multi_prot_to_predict_2,
            self._view.ui.btn_pred_multi_back_2,
            self._view.ui.btn_pred_multi_prot_to_predict_add_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_multi_prot_name_status,
            self._view.ui.btn_pred_multi_back,
            self._view.ui.btn_pred_multi_next,
            self._view.ui.lbl_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self._view.ui.txt_pred_multi_prot_seq, self._view.ui.lbl_pred_multi_prot_seq)
        gui_utils.disable_text_box(self._view.ui.txt_pred_multi_prot_name, self._view.ui.lbl_pred_multi_prot_name)
        self._view.ui.txt_pred_multi_prot_seq.clear()
        self._view.ui.list_pred_multi_prot_seq_overview.clear()
        self._view.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        self._view.ui.btn_pred_multi_prot_seq_overview_remove.setEnabled(False)
        styles.color_button_not_ready(self._view.ui.btn_pred_multi_prot_to_predict_add_2)

    def local_pred_multi_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_multi_prot_to_predict,
            self._view.ui.table_pred_multi_prot_to_predict,
            self._view.ui.lbl_pred_multi_prot_name_status,
            self._view.ui.btn_pred_multi_back,
            self._view.ui.btn_pred_multi_next,
            self._view.ui.lbl_pred_multi_prot_name,
            self._view.ui.txt_pred_multi_prot_name,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_multi_prot_seq,
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add,
            self._view.ui.btn_pred_multi_prot_seq_add,
            self._view.ui.lbl_pred_multi_prot_seq_overview,
            self._view.ui.list_pred_multi_prot_seq_overview,
            self._view.ui.btn_pred_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_multi_prot_to_predict_2,
            self._view.ui.btn_pred_multi_back_2,
            self._view.ui.btn_pred_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_predict,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(self._view.ui.txt_pred_multi_prot_name, self._view.ui.lbl_pred_multi_prot_name)
        gui_utils.disable_text_box(self._view.ui.txt_pred_multi_prot_seq, self._view.ui.lbl_pred_multi_prot_seq)

    def local_pred_multi_prot_seq_overview_item_changed(self) -> None:
        """Enables the remove button of the list sequences of the protein."""
        self._view.ui.btn_pred_multi_prot_seq_overview_remove.setEnabled(True)

    def local_pred_multi_prot_to_predict_item_changed(self) -> None:
        """Enables the remove button of the table proteins to predict."""
        self._view.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(True)

    def local_pred_multi_prot_to_predict_add_2(self) -> None:
        """Adds the protein to the list of proteins to predict."""
        for i in range(self._view.ui.list_pred_multi_prot_seq_overview.count()):
            self._view.ui.table_pred_multi_prot_to_predict.setRowCount(
                self._view.ui.table_pred_multi_prot_to_predict.rowCount() + 1,
            )
            self._view.ui.table_pred_multi_prot_to_predict.insertRow(
                self._view.ui.table_pred_multi_prot_to_predict.rowCount() + 1,
            )
            tmp_chain_seq = (
                constants.chain_dict.get(i),
                self._view.ui.list_pred_multi_prot_seq_overview.item(i).text(),
            )
            self._view.ui.table_pred_multi_prot_to_predict.setItem(
                self._view.ui.table_pred_multi_prot_to_predict.rowCount() - 1,
                0,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[0]),
            )
            self._view.ui.table_pred_multi_prot_to_predict.setItem(
                self._view.ui.table_pred_multi_prot_to_predict.rowCount() - 1,
                1,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[1]),
            )
            name_item = QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_multi_prot_name.text())
            self._view.ui.table_pred_multi_prot_to_predict.setVerticalHeaderItem(
                self._view.ui.table_pred_multi_prot_to_predict.rowCount() - 1,
                name_item,
            )
        self._view.ui.table_pred_multi_prot_to_predict.resizeColumnsToContents()
        self.local_pred_multi_check_if_table_is_empty()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_multi_prot_to_predict,
            self._view.ui.table_pred_multi_prot_to_predict,
            self._view.ui.btn_pred_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_multi_prot_name_status,
            self._view.ui.btn_pred_multi_back,
            self._view.ui.btn_pred_multi_next,
            self._view.ui.lbl_pred_multi_prot_name,
            self._view.ui.txt_pred_multi_prot_name,
            self._view.ui.lbl_pred_multi_prot_seq,
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add,
            self._view.ui.btn_pred_multi_prot_seq_add,
            self._view.ui.lbl_pred_multi_prot_seq_overview,
            self._view.ui.list_pred_multi_prot_seq_overview,
            self._view.ui.btn_pred_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_multi_prot_to_predict_2,
            self._view.ui.btn_pred_multi_back_2,
            self._view.ui.btn_pred_multi_prot_to_predict_add_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._init_local_pred_multi_page()
        self._view.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)

    def local_pred_multi_remove(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        self._view.ui.table_pred_multi_prot_to_predict.removeRow(
            self._view.ui.table_pred_multi_prot_to_predict.currentRow(),
        )
        if self._view.ui.table_pred_multi_prot_to_predict.rowCount() > 0:
            prot_name = self._view.ui.table_pred_multi_prot_to_predict.verticalHeaderItem(
                self._view.ui.table_pred_multi_prot_to_predict.currentRow(),
            ).text()
            for i in range(self._view.ui.table_pred_multi_prot_to_predict.rowCount()):
                if self._view.ui.table_pred_multi_prot_to_predict.verticalHeaderItem(i).text() == prot_name:
                    self._view.ui.table_pred_multi_prot_to_predict.setItem(
                        i,
                        0,
                        QtWidgets.QTableWidgetItem(constants.chain_dict.get(i)),
                    )
        self.local_pred_multi_check_if_table_is_empty()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_multi_prot_to_predict,
            self._view.ui.table_pred_multi_prot_to_predict,
            self._view.ui.btn_pred_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_advanced_config,
            self._view.ui.btn_pred_multi_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_multi_prot_name_status,
            self._view.ui.btn_pred_multi_back,
            self._view.ui.btn_pred_multi_next,
            self._view.ui.lbl_pred_multi_prot_name,
            self._view.ui.txt_pred_multi_prot_name,
            self._view.ui.lbl_pred_multi_prot_seq,
            self._view.ui.txt_pred_multi_prot_seq,
            self._view.ui.lbl_pred_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add,
            self._view.ui.btn_pred_multi_prot_seq_add,
            self._view.ui.lbl_pred_multi_prot_seq_overview,
            self._view.ui.list_pred_multi_prot_seq_overview,
            self._view.ui.btn_pred_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_multi_prot_to_predict_2,
            self._view.ui.btn_pred_multi_back_2,
            self._view.ui.btn_pred_multi_prot_to_predict_add_2,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        self.local_pred_multi_check_if_table_is_empty()

    # </editor-fold>

    # <editor-fold desc="Monomer prediction + analysis">
    def _init_mono_pred_analysis_page(self) -> None:
        """Clears all text boxes and sets default values for the gui elements."""
        # <editor-fold desc="Prediction section">
        self._view.ui.txt_pred_analysis_mono_prot_name.clear()
        self._view.ui.txt_pred_analysis_mono_seq_name.clear()
        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            self._view.ui.table_pred_analysis_mono_prot_to_predict.removeRow(i)
        # sets up defaults: Prediction
        self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self._view.ui.btn_pred_analysis_mono_add_protein.setEnabled(False)
        self._view.ui.lbl_pred_analysis_mono_prot_name_status.setText("")
        self._view.ui.lbl_pred_analysis_mono_seq_name_status.setText("")

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.list_pred_analysis_mono_overview.clear()
        self._view.ui.btn_pred_analysis_mono_remove.hide()

        # </editor-fold>

    def display_monomer_pred_analysis(self) -> None:
        """Displays the monomer prediction + analysis page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return

        self._init_mono_pred_analysis_page()
        self._view.ui.table_pred_analysis_mono_prot_to_predict.clear()
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(
            0,
            QtWidgets.QTableWidgetItem("Chain"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(
            1,
            QtWidgets.QTableWidgetItem("Sequence"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.tabWidget.currentIndex() == 1:
            self._view.ui.tabWidget.setCurrentIndex(0)
        self._view.ui.tabWidget.setTabEnabled(1, False)
        self._view.ui.tabWidget.setTabEnabled(0, True)
        tools.switch_page(
            self._view.ui.stackedWidget,
            self._view.ui.lbl_page_title,
            21,
            "Monomer Prediction + Analysis",
        )
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_pred_analysis_monomer_page,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setEnabled(True)

    # <editor-fold desc="Sections">

    # <editor-fold desc="Prediction section">
    def mono_pred_analysis_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self._view.ui.txt_pred_analysis_mono_prot_name.text(),
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
        ):
            self._view.ui.lbl_pred_analysis_mono_prot_name_status.setText("Protein name already used.")
            self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_pred_analysis_mono_next)
        else:
            self._view.ui.btn_pred_analysis_mono_next.setEnabled(True)
            tools.validate_protein_name(
                self._view.ui.txt_pred_analysis_mono_prot_name,
                self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                self._view.ui.btn_pred_analysis_mono_next,
            )

    def mono_pred_analysis_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_add_protein,
        )

    def mono_pred_analysis_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self._view.ui.btn_pred_analysis_mono_go_analysis_setup)
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(False)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                self._view.ui.table_pred_analysis_mono_prot_to_predict,
                self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                self._view.ui.lbl_pred_analysis_mono_prot_name,
                self._view.ui.txt_pred_analysis_mono_prot_name,
                self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                self._view.ui.btn_pred_analysis_mono_back,
                self._view.ui.btn_pred_analysis_mono_next,
                self._view.ui.lbl_pred_analysis_mono_seq_name,
                self._view.ui.txt_pred_analysis_mono_seq_name,
                self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                self._view.ui.btn_pred_analysis_mono_back_2,
                self._view.ui.btn_pred_analysis_mono_add_protein,
                self._view.ui.lbl_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                self._view.ui.table_pred_analysis_mono_prot_to_predict,
                self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                self._view.ui.btn_pred_analysis_mono_seq_to_predict,
                self._view.ui.lbl_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.lbl_pred_analysis_mono_prot_name,
                self._view.ui.txt_pred_analysis_mono_prot_name,
                self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                self._view.ui.btn_pred_analysis_mono_back,
                self._view.ui.btn_pred_analysis_mono_next,
                self._view.ui.lbl_pred_analysis_mono_seq_name,
                self._view.ui.txt_pred_analysis_mono_seq_name,
                self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                self._view.ui.btn_pred_analysis_mono_back_2,
                self._view.ui.btn_pred_analysis_mono_add_protein,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(True)

    def setup_defaults_monomer_prediction_analysis(self) -> None:
        """Sets up default values for the prediction tab."""
        # clears everything
        self._view.ui.txt_pred_analysis_mono_prot_name.clear()
        self._view.ui.txt_pred_analysis_mono_seq_name.clear()
        # sets up defaults: Prediction
        self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self._view.ui.btn_pred_analysis_mono_add_protein.setEnabled(False)
        self._view.ui.lbl_pred_analysis_mono_prot_name_status.setText("")
        self._view.ui.lbl_pred_analysis_mono_seq_name_status.setText("")

    def mono_pred_analysis_add_seq_to_predict(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
        )
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self._view.ui.txt_pred_analysis_mono_prot_name.clear()
        styles.color_button_not_ready(self._view.ui.btn_pred_analysis_mono_next)
        if self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() > 0:
            try:
                self._view.ui.table_pred_analysis_mono_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Monomer Prediction in overview table.")

    def mono_pred_analysis_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.mono_pred_analysis_check_if_table_is_empty()
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)

    def mono_pred_analysis_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
        )
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.txt_pred_analysis_mono_seq_name.clear()

    def mono_pred_analysis_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
        )
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def mono_pred_analysis_add_protein(self) -> None:
        """Adds the protein to the list of proteins to predict."""
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setRowCount(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.insertRow(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            0,
            QtWidgets.QTableWidgetItem("A"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_seq_name.toPlainText()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setVerticalHeaderItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_prot_name.text()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        self.mono_pred_analysis_check_if_table_is_empty()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(True)
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)
        self.setup_defaults_monomer_prediction()

    def mono_pred_analysis_prediction_overview_item_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(True)

    def mono_pred_analysis_add_protein_to_predict(self) -> None:
        """Needs to be removed."""
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setRowCount(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.insertRow(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            0,
            QtWidgets.QTableWidgetItem("A"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_seq_name.toPlainText()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setVerticalHeaderItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_prot_name.text()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        self.mono_pred_analysis_check_if_table_is_empty()
        self.setup_defaults_monomer_prediction_analysis()

    def mono_pred_analysis_remove_protein_to_predict(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        self._view.ui.table_pred_analysis_mono_prot_to_predict.removeRow(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.currentRow(),
        )
        self.mono_pred_analysis_check_if_table_is_empty()
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)

    # </editor-fold>

    # <editor-fold desc="Analysis section">
    def mono_pred_analysis_structure_analysis_add(self) -> None:
        """Shows the gui elements for the selection of the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            self._view.ui.lbl_pred_analysis_mono_images,
            self._view.ui.cb_pred_analysis_mono_images,
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_mono_prot_struct_1.clear()
        self._view.ui.lbl_pred_analysis_mono_prot_struct_2.clear()
        self._view.ui.lbl_pred_analysis_mono_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_mono_prot_struct_2.setText("Protein structure 2")
        self.fill_mono_pred_analysis_protein_boxes()
        if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
            try:
                self._view.ui.list_pred_analysis_mono_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def mono_pred_analysis_structure_analysis_back_3(self) -> None:
        """Hides the gui elements to select the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            self._view.ui.lbl_pred_analysis_mono_images,
            self._view.ui.cb_pred_analysis_mono_images,
            self._view.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
            self._view.ui.btn_pred_analysis_mono_remove.show()
            self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_mono_start.show()
            self._view.ui.lbl_pred_analysis_mono_images.show()
            self._view.ui.cb_pred_analysis_mono_images.show()

    def mono_pred_analysis_structure_analysis_next_3(self) -> None:
        """Shows the gui elements to select the chains of protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.lbl_pred_analysis_mono_images,
            self._view.ui.cb_pred_analysis_mono_images,
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_mono_model_chains.clear()
        self._view.ui.list_pred_analysis_mono_ref_chains.setEnabled(False)
        self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(False)

        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            if (
                self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text()
                == self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText()
            ):
                self._view.ui.list_pred_analysis_mono_model_chains.addItem(
                    self._view.ui.table_pred_analysis_mono_prot_to_predict.item(i, 0).text(),
                )
        if self._view.ui.list_pred_analysis_mono_model_chains.count() == 0:
            tmp_protein = self._current_project.search_protein(
                self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText(),
            )
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self._view.ui.list_pred_analysis_mono_model_chains.addItem(tmp_chain.chain_letter)
        if self._view.ui.list_pred_analysis_mono_model_chains.count() == 1:
            self._view.ui.lbl_pred_analysis_mono_model_chains.setText(
                f"Select chain in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.",
            )
        else:
            self._view.ui.lbl_pred_analysis_mono_model_chains.setText(
                f"Select {len(self._view.ui.list_pred_analysis_mono_model_chains.selectedItems())} chains "
                f"in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.",
            )

    def mono_pred_analysis_structure_analysis_back_4(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            self._view.ui.lbl_pred_analysis_mono_images,
            self._view.ui.cb_pred_analysis_mono_images,
            self._view.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_mono_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_mono_prot_struct_2.setText("Protein structure 2")

    def mono_pred_analysis_structure_analysis_next_4(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_images,
            self._view.ui.cb_pred_analysis_mono_images,
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self._view.ui.list_pred_analysis_mono_overview.addItem(item)
        self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)

    def mono_pred_analysis_structure_analysis_back_5(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            self._view.ui.lbl_pred_analysis_mono_images,
            self._view.ui.cb_pred_analysis_mono_images,
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_mono_ref_chains.setEnabled(True)

        # tmp_protein = self._current_project.search_protein(self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self._view.ui.list_pred_analysis_mono_ref_chains.addItem(tmp_chain.chain_letter)

    def mono_pred_analysis_structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_pred_analysis_mono_remove.setEnabled(True)

    def fill_mono_pred_analysis_protein_boxes(self) -> None:
        """Fills the combo box of the protein structures."""
        protein_names = []
        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            protein_names.append(self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text())
        for tmp_protein in self._current_project.proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        self._view.ui.box_pred_analysis_mono_prot_struct_1.clear()
        self._view.ui.box_pred_analysis_mono_prot_struct_2.clear()
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_mono_prot_struct_1, protein_names)
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_mono_prot_struct_2, protein_names)

    def remove_mono_pred_analysis_analysis_run(self) -> None:
        """Removes a selected protein pair form the list of protein pairs to analyze."""
        self._view.ui.list_pred_analysis_mono_overview.takeItem(
            self._view.ui.list_pred_analysis_mono_overview.currentRow(),
        )
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            self._view.ui.lbl_pred_analysis_mono_images,
            self._view.ui.cb_pred_analysis_mono_images,
            self._view.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
            self._view.ui.btn_pred_analysis_mono_remove.show()
            self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_mono_start.show()
            self._view.ui.lbl_pred_analysis_mono_images.show()
            self._view.ui.cb_pred_analysis_mono_images.show()
        # if self._view.ui.list_pred_analysis_mono_overview.count() == 0:
        #
        #     self._view.ui.btn_pred_analysis_mono_back_pred_setup.show()
        #     self._view.ui.btn_pred_analysis_mono_remove.hide()

    def check_mono_pred_analysis_if_same_no_of_chains_selected(self) -> None:
        """Checks if the same number of chains were selected."""
        self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(False)
        if self.no_of_selected_chains == len(self._view.ui.list_pred_analysis_mono_model_chains.selectedItems()):
            self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(True)

        prot_1_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self._view.ui.list_pred_analysis_mono_overview.count()):
            if analysis_name == self._view.ui.list_pred_analysis_mono_overview.item(tmp_row).text():
                self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_pred_analysis_mono_next_4)
                return

    def check_mono_pred_analysis_if_prot_structs_are_filled(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self._view.ui.box_pred_analysis_mono_prot_struct_1.itemText(
            self._view.ui.box_pred_analysis_mono_prot_struct_1.currentIndex(),
        )
        prot_2 = self._view.ui.box_pred_analysis_mono_prot_struct_2.itemText(
            self._view.ui.box_pred_analysis_mono_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self._view.ui.btn_pred_analysis_mono_next_2.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_mono_next_2.setEnabled(False)

    def count_mono_pred_analysis_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains selected in protein 1."""
        self.no_of_selected_chains = len(self._view.ui.list_pred_analysis_mono_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self._view.ui.btn_pred_analysis_mono_next_3.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_mono_next_3.setEnabled(False)

    # </editor-fold>

    # </editor-fold>

    def switch_monomer_pred_analysis_tab(self) -> None:
        """Switches the tabs from prediction to analysis and vice versa."""
        if self._view.ui.tabWidget.currentIndex() == 0:
            # goes from prediction to analysis
            self._view.ui.tabWidget.setTabEnabled(1, True)
            self._view.ui.tabWidget.setTabEnabled(0, False)
            self._view.ui.tabWidget.setCurrentIndex(1)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_mono_overview,
                self._view.ui.list_pred_analysis_mono_overview,
                self._view.ui.btn_pred_analysis_mono_add,
                self._view.ui.btn_pred_analysis_mono_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_mono_remove,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_2,
                self._view.ui.lbl_pred_analysis_mono_ref_chains,
                self._view.ui.list_pred_analysis_mono_ref_chains,
                self._view.ui.btn_pred_analysis_mono_back_4,
                self._view.ui.btn_pred_analysis_mono_next_3,
                self._view.ui.box_pred_analysis_mono_prot_struct_1,
                self._view.ui.box_pred_analysis_mono_prot_struct_2,
                self._view.ui.btn_pred_analysis_mono_back_3,
                self._view.ui.btn_pred_analysis_mono_next_2,
                self._view.ui.lbl_pred_analysis_mono_model_chains,
                self._view.ui.list_pred_analysis_mono_model_chains,
                self._view.ui.btn_pred_analysis_mono_back_5,
                self._view.ui.btn_pred_analysis_mono_next_4,
                self._view.ui.lbl_pred_analysis_mono_images,
                self._view.ui.cb_pred_analysis_mono_images,
                self._view.ui.btn_pred_analysis_mono_start,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
                self._view.ui.btn_pred_analysis_mono_remove.show()
                self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
                self._view.ui.btn_pred_analysis_mono_start.show()
                self._view.ui.lbl_pred_analysis_mono_images.show()
                self._view.ui.cb_pred_analysis_mono_images.show()
        else:
            # goes from analysis to prediction
            if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
                gui_elements_to_show = [
                    self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                    self._view.ui.table_pred_analysis_mono_prot_to_predict,
                    self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                    self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                    self._view.ui.btn_pred_analysis_mono_seq_to_predict,
                    self._view.ui.lbl_pred_analysis_mono_prot_name,
                    self._view.ui.txt_pred_analysis_mono_prot_name,
                    self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                    self._view.ui.btn_pred_analysis_mono_back,
                    self._view.ui.btn_pred_analysis_mono_next,
                    self._view.ui.lbl_pred_analysis_mono_seq_name,
                    self._view.ui.txt_pred_analysis_mono_seq_name,
                    self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                    self._view.ui.btn_pred_analysis_mono_back_2,
                    self._view.ui.btn_pred_analysis_mono_add_protein,
                    self._view.ui.lbl_pred_mono_advanced_config_2,
                    self._view.ui.btn_pred_mono_advanced_config_2,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
            else:
                gui_elements_to_show = [
                    self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                    self._view.ui.table_pred_analysis_mono_prot_to_predict,
                    self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                    self._view.ui.btn_pred_analysis_mono_seq_to_predict,
                    self._view.ui.lbl_pred_mono_advanced_config_2,
                    self._view.ui.btn_pred_mono_advanced_config_2,
                    self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                    self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self._view.ui.lbl_pred_analysis_mono_prot_name,
                    self._view.ui.txt_pred_analysis_mono_prot_name,
                    self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                    self._view.ui.btn_pred_analysis_mono_back,
                    self._view.ui.btn_pred_analysis_mono_next,
                    self._view.ui.lbl_pred_analysis_mono_seq_name,
                    self._view.ui.txt_pred_analysis_mono_seq_name,
                    self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                    self._view.ui.btn_pred_analysis_mono_back_2,
                    self._view.ui.btn_pred_analysis_mono_add_protein,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
                self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)
            self._view.ui.tabWidget.setTabEnabled(0, True)
            self._view.ui.tabWidget.setTabEnabled(1, False)
            self._view.ui.tabWidget.setCurrentIndex(0)

    # </editor-fold>

    # <editor-fold desc="Multimer prediction + analysis">
    def _init_multi_pred_analysis_page(self) -> None:
        """Clears the text boxes and sets the default values for the gui elements."""
        # <editor-fold desc="Prediction section">
        # clears everything
        self._view.ui.txt_pred_analysis_multi_prot_name.clear()
        self._view.ui.txt_pred_analysis_multi_prot_seq.clear()
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.clear()
        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1, -1, -1):
            self._view.ui.table_pred_analysis_multi_prot_to_predict.removeRow(i)

        # sets up defaults: Prediction
        self._view.ui.btn_pred_analysis_multi_next.setEnabled(False)
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        self._view.ui.lbl_pred_analysis_multi_prot_name_status.setText("")
        self._view.ui.lbl_pred_analysis_multi_prot_seq_status.setText("")

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.list_pred_analysis_multi_overview.clear()
        self._view.ui.btn_pred_analysis_multi_remove.hide()

        # </editor-fold>

        # self.multi_pred_analysis_show_protein_overview()

    def display_multimer_pred_analysis(self) -> None:
        """Displays the multimer prediction + analysis page."""
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return

        self._init_multi_pred_analysis_page()
        self._view.ui.table_pred_analysis_multi_prot_to_predict.clear()
        self._view.ui.table_pred_analysis_multi_prot_to_predict.setHorizontalHeaderItem(
            0,
            QtWidgets.QTableWidgetItem("Chain"),
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.setHorizontalHeaderItem(
            1,
            QtWidgets.QTableWidgetItem("Sequence"),
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.resizeColumnsToContents()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add_2,
            self._view.ui.btn_pred_analysis_multi_prot_seq_add,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self._view.ui.list_pred_analysis_multi_prot_seq_overview,
            self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self._view.ui.btn_pred_analysis_multi_back_2,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.tabWidget_2.currentIndex() == 1:
            self._view.ui.tabWidget_2.setCurrentIndex(0)
        self._view.ui.tabWidget_2.setTabEnabled(1, False)
        self._view.ui.tabWidget_2.setTabEnabled(0, True)
        tools.switch_page(
            self._view.ui.stackedWidget,
            self._view.ui.lbl_page_title,
            22,
            "Multimer Prediction + Analysis",
        )
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_pred_analysis_multimer_page,
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.setEnabled(True)

    # <editor-fold desc="Prediction section">
    def multi_pred_analysis_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
            self._view.ui.txt_pred_analysis_multi_prot_name.text(),
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
        ):
            self._view.ui.lbl_pred_analysis_multi_prot_name_status.setText("Protein name already used.")
            self._view.ui.btn_pred_analysis_multi_next.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_next)
        else:
            self._view.ui.btn_pred_analysis_multi_next.setEnabled(True)
            tools.validate_protein_name(
                self._view.ui.txt_pred_analysis_multi_prot_name,
                self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                self._view.ui.btn_pred_analysis_multi_next,
            )

    def multi_pred_analysis_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
            self._view.ui.btn_pred_analysis_multi_prot_seq_add,
        )

    def multi_pred_analysis_check_if_list_is_empty(self) -> None:
        """Checks if the list of sequences of the protein is empty."""
        if self._view.ui.list_pred_analysis_multi_prot_seq_overview.count() == 0:
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        else:
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(True)

    def multi_pred_analysis_add_sequence_to_list(self) -> None:
        """Adds the entered sequence to the sequences of the protein."""
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.addItem(
            QtWidgets.QListWidgetItem(self._view.ui.txt_pred_analysis_multi_prot_seq.toPlainText()),
        )
        self.multi_pred_analysis_check_if_list_is_empty()

    def multi_pred_analysis_remove_sequence_to_list(self) -> None:
        """Removes the entered sequence from the sequences of the protein."""
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.takeItem(
            self._view.ui.list_pred_analysis_multi_prot_seq_overview.currentRow(),
        )
        self.multi_pred_analysis_check_if_list_is_empty()
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(False)

    def multi_pred_analysis_check_if_table_is_empty(self) -> None:
        """Checks if the list of proteins to predict is empty."""
        if self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self._view.ui.btn_pred_multi_predict)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                self._view.ui.table_pred_analysis_multi_prot_to_predict,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                self._view.ui.lbl_pred_analysis_multi_prot_name,
                self._view.ui.txt_pred_analysis_multi_prot_name,
                self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                self._view.ui.btn_pred_analysis_multi_back,
                self._view.ui.btn_pred_analysis_multi_next,
                self._view.ui.lbl_pred_analysis_multi_prot_seq,
                self._view.ui.txt_pred_analysis_multi_prot_seq,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
                self._view.ui.lbl_pred_multi_prot_seq_add_2,
                self._view.ui.btn_pred_analysis_multi_prot_seq_add,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
                self._view.ui.list_pred_analysis_multi_prot_seq_overview,
                self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                self._view.ui.btn_pred_analysis_multi_back_2,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
                self._view.ui.lbl_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_pred_multi_predict.setEnabled(False)
            self._view.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        else:
            self._view.ui.btn_pred_multi_predict.setEnabled(True)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                self._view.ui.table_pred_analysis_multi_prot_to_predict,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
                self._view.ui.lbl_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.lbl_pred_analysis_multi_prot_name,
                self._view.ui.txt_pred_analysis_multi_prot_name,
                self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                self._view.ui.btn_pred_analysis_multi_back,
                self._view.ui.btn_pred_analysis_multi_next,
                self._view.ui.lbl_pred_analysis_multi_prot_seq,
                self._view.ui.txt_pred_analysis_multi_prot_seq,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
                self._view.ui.lbl_pred_multi_prot_seq_add_2,
                self._view.ui.btn_pred_analysis_multi_prot_seq_add,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
                self._view.ui.list_pred_analysis_multi_prot_seq_overview,
                self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                self._view.ui.btn_pred_analysis_multi_back_2,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_add_protein_to_predict(self) -> None:
        """Adds the proteins to the list of proteins to predict."""
        for i in range(self._view.ui.list_pred_analysis_multi_prot_seq_overview.count()):
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setRowCount(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1,
            )
            self._view.ui.table_pred_analysis_multi_prot_to_predict.insertRow(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1,
            )
            tmp_chain_seq = (
                constants.chain_dict.get(i),
                self._view.ui.list_pred_analysis_multi_prot_seq_overview.item(i).text(),
            )
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                0,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[0]),
            )
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                1,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[1]),
            )
            name_item = QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_multi_prot_name.text())
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setVerticalHeaderItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                name_item,
            )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.resizeColumnsToContents()
        self.multi_pred_analysis_check_if_table_is_empty()
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_remove_protein_to_predict(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        if self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 1:
            self._view.ui.table_pred_analysis_multi_prot_to_predict.removeRow(0)
        else:
            self._view.ui.table_pred_analysis_multi_prot_to_predict.removeRow(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.currentRow(),
            )
            prot_name = self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.currentRow(),
            ).text()
            for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
                if self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text() == prot_name:
                    self._view.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                        i,
                        0,
                        QtWidgets.QTableWidgetItem(constants.chain_dict.get(i)),
                    )
        self.multi_pred_analysis_check_if_table_is_empty()
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_add(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add_2,
            self._view.ui.btn_pred_analysis_multi_prot_seq_add,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self._view.ui.list_pred_analysis_multi_prot_seq_overview,
            self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self._view.ui.btn_pred_analysis_multi_back_2,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
        )
        self._view.ui.btn_pred_analysis_multi_next.setEnabled(False)
        self._view.ui.txt_pred_analysis_multi_prot_name.clear()
        styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_next)
        if self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() > 0:
            try:
                self._view.ui.table_pred_analysis_multi_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Multimer Prediction in overview table.")

    def multi_pred_analysis_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add_2,
            self._view.ui.btn_pred_analysis_multi_prot_seq_add,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self._view.ui.list_pred_analysis_multi_prot_seq_overview,
            self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self._view.ui.btn_pred_analysis_multi_back_2,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.multi_pred_analysis_check_if_table_is_empty()

    def multi_pred_analysis_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add_2,
            self._view.ui.btn_pred_analysis_multi_prot_seq_add,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self._view.ui.list_pred_analysis_multi_prot_seq_overview,
            self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self._view.ui.btn_pred_analysis_multi_back_2,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
        )
        self._view.ui.txt_pred_analysis_multi_prot_seq.clear()
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.clear()
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(False)
        styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2)

    def multi_pred_analysis_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
            self._view.ui.lbl_pred_multi_prot_seq_add_2,
            self._view.ui.btn_pred_analysis_multi_prot_seq_add,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
            self._view.ui.list_pred_analysis_multi_prot_seq_overview,
            self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
            self._view.ui.btn_pred_analysis_multi_back_2,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
        )

    def multi_pred_analysis_prot_seq_overview_item_changed(self) -> None:
        """Enables the remove button of the list of sequences of the protein."""
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(True)

    def multi_pred_analysis_prot_to_predict_item_changed(self) -> None:
        """Enables the remove button of the list of proteins to predict."""
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(True)

    # </editor-fold>

    def switch_multimer_pred_analysis_tab(self) -> None:
        """Switches the tabs from prediction to analysis and vice versa."""
        if self._view.ui.tabWidget_2.currentIndex() == 0:
            # goes from prediction to analysis
            self._view.ui.tabWidget_2.setCurrentIndex(1)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_multi_overview,
                self._view.ui.list_pred_analysis_multi_overview,
                self._view.ui.btn_pred_analysis_multi_add,
                self._view.ui.btn_pred_analysis_multi_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_multi_remove,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_3,
                self._view.ui.lbl_pred_analysis_multi_ref_chains,
                self._view.ui.list_pred_analysis_multi_ref_chains,
                self._view.ui.btn_pred_analysis_multi_back_4,
                self._view.ui.btn_pred_analysis_multi_next_3,
                self._view.ui.box_pred_analysis_multi_prot_struct_1,
                self._view.ui.box_pred_analysis_multi_prot_struct_2,
                self._view.ui.btn_pred_analysis_multi_back_3,
                self._view.ui.btn_pred_analysis_multi_next_2,
                self._view.ui.lbl_pred_analysis_multi_model_chains,
                self._view.ui.list_pred_analysis_multi_model_chains,
                self._view.ui.btn_pred_analysis_multi_back_5,
                self._view.ui.btn_pred_analysis_multi_next_4,
                self._view.ui.lbl_pred_analysis_multi_images,
                self._view.ui.cb_pred_analysis_multi_images,
                self._view.ui.btn_pred_analysis_multi_start,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.tabWidget_2.setTabEnabled(1, True)
            self._view.ui.tabWidget_2.setTabEnabled(0, False)
            if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
                self._view.ui.btn_pred_analysis_multi_remove.show()
                self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
                self._view.ui.btn_pred_analysis_multi_start.show()
                self._view.ui.lbl_pred_analysis_multi_images.show()
                self._view.ui.cb_pred_analysis_multi_images.show()
        else:
            # goes from analysis to prediction
            self._view.ui.tabWidget_2.setCurrentIndex(0)
            if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
                gui_elements_to_show = [
                    self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                    self._view.ui.table_pred_analysis_multi_prot_to_predict,
                    self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                    self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                    self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
                    self._view.ui.lbl_pred_analysis_multi_advanced_config,
                    self._view.ui.btn_pred_analysis_multi_advanced_config,
                    self._view.ui.lbl_pred_analysis_multi_prot_name,
                    self._view.ui.txt_pred_analysis_multi_prot_name,
                    self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                    self._view.ui.btn_pred_analysis_multi_back,
                    self._view.ui.btn_pred_analysis_multi_next,
                    self._view.ui.lbl_pred_analysis_multi_prot_seq,
                    self._view.ui.txt_pred_analysis_multi_prot_seq,
                    self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
                    self._view.ui.lbl_pred_multi_prot_seq_add_2,
                    self._view.ui.btn_pred_analysis_multi_prot_seq_add,
                    self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
                    self._view.ui.list_pred_analysis_multi_prot_seq_overview,
                    self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                    self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                    self._view.ui.btn_pred_analysis_multi_back_2,
                    self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
            else:
                gui_elements_to_show = [
                    self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                    self._view.ui.table_pred_analysis_multi_prot_to_predict,
                    self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                    self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
                    self._view.ui.lbl_pred_analysis_multi_advanced_config,
                    self._view.ui.btn_pred_analysis_multi_advanced_config,
                    self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                    self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
                ]
                gui_elements_to_hide = [
                    self._view.ui.lbl_pred_analysis_multi_prot_name,
                    self._view.ui.txt_pred_analysis_multi_prot_name,
                    self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                    self._view.ui.btn_pred_analysis_multi_back,
                    self._view.ui.btn_pred_analysis_multi_next,
                    self._view.ui.lbl_pred_analysis_multi_prot_seq,
                    self._view.ui.txt_pred_analysis_multi_prot_seq,
                    self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
                    self._view.ui.lbl_pred_multi_prot_seq_add_2,
                    self._view.ui.btn_pred_analysis_multi_prot_seq_add,
                    self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
                    self._view.ui.list_pred_analysis_multi_prot_seq_overview,
                    self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                    self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                    self._view.ui.btn_pred_analysis_multi_back_2,
                    self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.tabWidget_2.setTabEnabled(0, True)
            self._view.ui.tabWidget_2.setTabEnabled(1, False)

    # <editor-fold desc="Analysis section">
    def multi_pred_analysis_structure_analysis_add(self) -> None:
        """Shows the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.lbl_pred_analysis_multi_images,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_multi_prot_struct_1.clear()
        self._view.ui.lbl_pred_analysis_multi_prot_struct_2.clear()
        self._view.ui.lbl_pred_analysis_multi_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_multi_prot_struct_2.setText("Protein structure 2")
        self.fill_multi_pred_analysis_protein_boxes()
        if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
            try:
                self._view.ui.list_pred_analysis_multi_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def multi_pred_analysis_structure_analysis_back_3(self) -> None:
        """Hides the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.lbl_pred_analysis_multi_images,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
            self._view.ui.btn_pred_analysis_multi_remove.show()
            self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_multi_start.show()
            self._view.ui.lbl_pred_analysis_multi_images.show()
            self._view.ui.cb_pred_analysis_multi_images.show()

    def multi_pred_analysis_structure_analysis_next_3(self) -> None:
        """Shows the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.lbl_pred_analysis_multi_images,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_multi_model_chains.clear()
        self._view.ui.list_pred_analysis_multi_ref_chains.setEnabled(False)
        self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(False)

        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            if (
                self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text()
                == self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText()
            ):
                self._view.ui.list_pred_analysis_multi_model_chains.addItem(
                    self._view.ui.table_pred_analysis_multi_prot_to_predict.item(i, 0).text(),
                )
        if self._view.ui.list_pred_analysis_multi_model_chains.count() == 0:
            tmp_protein = self._current_project.search_protein(
                self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText(),
            )
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self._view.ui.list_pred_analysis_multi_model_chains.addItem(tmp_chain.chain_letter)
        if len(self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems()) == 1:
            self._view.ui.lbl_pred_analysis_multi_model_chains.setText(
                f"Select 1 chain in protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()}.",
            )
        else:
            self._view.ui.lbl_pred_analysis_multi_model_chains.setText(
                f"Select {len(self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems())} chains in "
                f"protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()}.",
            )

    def multi_pred_analysis_structure_analysis_back_4(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.lbl_pred_analysis_multi_images,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_multi_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_multi_prot_struct_2.setText("Protein structure 2")

    def multi_pred_analysis_structure_analysis_next_4(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_images,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self._view.ui.list_pred_analysis_multi_overview.addItem(item)
        self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)

    def multi_pred_analysis_structure_analysis_back_5(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.lbl_pred_analysis_multi_images,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_multi_ref_chains.setEnabled(True)

        # tmp_protein = self._current_project.search_protein(self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self._view.ui.list_pred_analysis_multi_ref_chains.addItem(tmp_chain.chain_letter)

    def multi_pred_analysis_structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_pred_analysis_multi_remove.setEnabled(True)

    def fill_multi_pred_analysis_protein_boxes(self) -> None:
        """Fills the combo boxes with the protein names."""
        protein_names = []
        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            protein_names.append(self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text())
        for tmp_protein in self._current_project.proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        protein_names = list(set(protein_names))
        self._view.ui.box_pred_analysis_multi_prot_struct_1.clear()
        self._view.ui.box_pred_analysis_multi_prot_struct_2.clear()
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_multi_prot_struct_1, protein_names)
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_multi_prot_struct_2, protein_names)

    def remove_multi_pred_analysis_analysis_run(self) -> None:
        """Removes the selected protein pair from the list of protein pairs to analyze."""
        self._view.ui.list_pred_analysis_multi_overview.takeItem(
            self._view.ui.list_pred_analysis_multi_overview.currentRow(),
        )
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.lbl_pred_analysis_multi_images,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
            self._view.ui.btn_pred_analysis_multi_remove.show()
            self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_multi_start.show()
            self._view.ui.lbl_pred_analysis_multi_images.show()
            self._view.ui.cb_pred_analysis_multi_images.show()
        # if self._view.ui.list_pred_analysis_multi_overview.count() == 0:
        #
        #     self._view.ui.btn_pred_analysis_multi_back_pred_setup.show()
        #     self._view.ui.btn_pred_analysis_multi_remove.hide()

    def check_multi_pred_analysis_if_same_no_of_chains_selected(self) -> None:
        """Checks if the same number of chains were selected."""
        self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(False)
        if self.no_of_selected_chains == len(self._view.ui.list_pred_analysis_multi_model_chains.selectedItems()):
            self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(True)

        prot_1_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self._view.ui.list_pred_analysis_multi_overview.count()):
            if analysis_name == self._view.ui.list_pred_analysis_multi_overview.item(tmp_row).text():
                self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_next_4)
                return

    def check_multi_pred_analysis_if_prot_structs_are_filled(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self._view.ui.box_pred_analysis_multi_prot_struct_1.itemText(
            self._view.ui.box_pred_analysis_multi_prot_struct_1.currentIndex(),
        )
        prot_2 = self._view.ui.box_pred_analysis_multi_prot_struct_2.itemText(
            self._view.ui.box_pred_analysis_multi_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self._view.ui.btn_pred_analysis_multi_next_2.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_multi_next_2.setEnabled(False)

    def count_multi_pred_analysis_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains in protein 1."""
        self.no_of_selected_chains = len(self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self._view.ui.btn_pred_analysis_multi_next_3.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_multi_next_3.setEnabled(False)

    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Distance analysis">
    def structure_analysis_add(self) -> None:
        """Shows the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_batch_overview,
            self._view.ui.list_analysis_batch_overview,
            self._view.ui.lbl_analysis_batch_prot_struct_1,
            self._view.ui.box_analysis_batch_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs,
            self._view.ui.lbl_analysis_batch_prot_struct_2,
            self._view.ui.box_analysis_batch_prot_struct_2,
            self._view.ui.btn_analysis_batch_back,
            self._view.ui.btn_analysis_batch_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_batch_remove,
            self._view.ui.btn_analysis_batch_add,
            self._view.ui.lbl_analysis_batch_ref_chains,
            self._view.ui.list_analysis_batch_ref_chains,
            self._view.ui.btn_analysis_batch_back_2,
            self._view.ui.btn_analysis_batch_next_2,
            self._view.ui.lbl_analysis_batch_model_chains,
            self._view.ui.list_analysis_batch_model_chains,
            self._view.ui.btn_analysis_batch_back_3,
            self._view.ui.btn_analysis_batch_next_3,
            self._view.ui.lbl_analysis_batch_images,
            self._view.ui.cb_analysis_batch_images,
            self._view.ui.btn_analysis_batch_start,
        ]

        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_analysis_batch_prot_struct_1.clear()
        self._view.ui.lbl_analysis_batch_prot_struct_2.clear()
        self._view.ui.lbl_analysis_batch_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_analysis_batch_prot_struct_2.setText("Protein structure 2")
        self.fill_protein_boxes_batch()
        if self._view.ui.list_analysis_batch_overview.count() > 0:
            try:
                self._view.ui.list_analysis_batch_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def structure_analysis_back(self) -> None:
        """Hides the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_batch_overview,
            self._view.ui.list_analysis_batch_overview,
            self._view.ui.btn_analysis_batch_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_batch_remove,
            self._view.ui.lbl_analysis_batch_prot_struct_1,
            self._view.ui.lbl_analysis_batch_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs,
            self._view.ui.lbl_analysis_batch_ref_chains,
            self._view.ui.list_analysis_batch_ref_chains,
            self._view.ui.btn_analysis_batch_back_2,
            self._view.ui.btn_analysis_batch_next_2,
            self._view.ui.box_analysis_batch_prot_struct_1,
            self._view.ui.box_analysis_batch_prot_struct_2,
            self._view.ui.btn_analysis_batch_back,
            self._view.ui.btn_analysis_batch_next,
            self._view.ui.lbl_analysis_batch_model_chains,
            self._view.ui.list_analysis_batch_model_chains,
            self._view.ui.btn_analysis_batch_back_3,
            self._view.ui.btn_analysis_batch_next_3,
            self._view.ui.lbl_analysis_batch_images,
            self._view.ui.cb_analysis_batch_images,
            self._view.ui.btn_analysis_batch_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_analysis_batch_overview.count() > 0:
            self._view.ui.btn_analysis_batch_remove.show()
            self._view.ui.btn_analysis_batch_remove.setEnabled(False)
            self._view.ui.btn_analysis_batch_start.show()
            self._view.ui.lbl_analysis_batch_images.show()
            self._view.ui.cb_analysis_batch_images.show()

    def structure_analysis_next_2(self) -> None:
        """Shows the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_batch_overview,
            self._view.ui.list_analysis_batch_overview,
            self._view.ui.lbl_analysis_batch_prot_struct_1,
            self._view.ui.lbl_analysis_batch_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs,
            self._view.ui.lbl_analysis_batch_ref_chains,
            self._view.ui.list_analysis_batch_ref_chains,
            self._view.ui.lbl_analysis_batch_model_chains,
            self._view.ui.list_analysis_batch_model_chains,
            self._view.ui.btn_analysis_batch_back_3,
            self._view.ui.btn_analysis_batch_next_3,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_batch_remove,
            self._view.ui.btn_analysis_batch_add,
            self._view.ui.box_analysis_batch_prot_struct_1,
            self._view.ui.box_analysis_batch_prot_struct_2,
            self._view.ui.btn_analysis_batch_back,
            self._view.ui.btn_analysis_batch_next,
            self._view.ui.btn_analysis_batch_back_2,
            self._view.ui.btn_analysis_batch_next_2,
            self._view.ui.lbl_analysis_batch_images,
            self._view.ui.cb_analysis_batch_images,
            self._view.ui.btn_analysis_batch_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_analysis_batch_model_chains.clear()
        self._view.ui.list_analysis_batch_ref_chains.setEnabled(False)
        self._view.ui.btn_analysis_batch_next_3.setEnabled(False)

        tmp_protein = self._current_project.search_protein(self._view.ui.box_analysis_batch_prot_struct_2.currentText())
        for tmp_chain in tmp_protein.chains:
            if tmp_chain.chain_type == "protein_chain":
                self._view.ui.list_analysis_batch_model_chains.addItem(tmp_chain.chain_letter)
        if len(self._view.ui.list_analysis_batch_ref_chains.selectedItems()) == 1:
            self._view.ui.lbl_analysis_batch_model_chains.setText(
                f"Select 1 chain in protein structure {self._view.ui.lbl_analysis_batch_prot_struct_2.text()}.",
            )
        else:
            self._view.ui.lbl_analysis_batch_model_chains.setText(
                f"Select {len(self._view.ui.list_analysis_batch_ref_chains.selectedItems())} chains in "
                f"protein structure {self._view.ui.lbl_analysis_batch_prot_struct_2.text()}.",
            )

    def structure_analysis_back_2(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_batch_overview,
            self._view.ui.list_analysis_batch_overview,
            self._view.ui.lbl_analysis_batch_prot_struct_1,
            self._view.ui.box_analysis_batch_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs,
            self._view.ui.lbl_analysis_batch_prot_struct_2,
            self._view.ui.box_analysis_batch_prot_struct_2,
            self._view.ui.btn_analysis_batch_back,
            self._view.ui.btn_analysis_batch_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_batch_remove,
            self._view.ui.btn_analysis_batch_add,
            self._view.ui.lbl_analysis_batch_ref_chains,
            self._view.ui.list_analysis_batch_ref_chains,
            self._view.ui.btn_analysis_batch_back_2,
            self._view.ui.btn_analysis_batch_next_2,
            self._view.ui.lbl_analysis_batch_model_chains,
            self._view.ui.list_analysis_batch_model_chains,
            self._view.ui.btn_analysis_batch_back_3,
            self._view.ui.btn_analysis_batch_next_3,
            self._view.ui.lbl_analysis_batch_images,
            self._view.ui.cb_analysis_batch_images,
            self._view.ui.btn_analysis_batch_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_analysis_batch_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_analysis_batch_prot_struct_2.setText("Protein structure 2")

    def structure_analysis_next_3(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self._view.ui.btn_analysis_batch_remove,
            self._view.ui.btn_analysis_batch_add,
            self._view.ui.lbl_analysis_batch_overview,
            self._view.ui.list_analysis_batch_overview,
            self._view.ui.lbl_analysis_batch_images,
            self._view.ui.cb_analysis_batch_images,
            self._view.ui.btn_analysis_batch_start,
        ]
        gui_elements_to_hide = [
            self._view.ui.box_analysis_batch_prot_struct_1,
            self._view.ui.box_analysis_batch_prot_struct_2,
            self._view.ui.lbl_analysis_batch_prot_struct_1,
            self._view.ui.lbl_analysis_batch_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs,
            self._view.ui.lbl_analysis_batch_ref_chains,
            self._view.ui.list_analysis_batch_ref_chains,
            self._view.ui.lbl_analysis_batch_model_chains,
            self._view.ui.list_analysis_batch_model_chains,
            self._view.ui.btn_analysis_batch_back,
            self._view.ui.btn_analysis_batch_next,
            self._view.ui.btn_analysis_batch_back_2,
            self._view.ui.btn_analysis_batch_next_2,
            self._view.ui.btn_analysis_batch_back_3,
            self._view.ui.btn_analysis_batch_next_3,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self._view.ui.lbl_analysis_batch_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_analysis_batch_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_analysis_batch_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_analysis_batch_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self._view.ui.list_analysis_batch_overview.addItem(item)
        self._view.ui.btn_analysis_batch_remove.setEnabled(False)

    def structure_analysis_back_3(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_batch_overview,
            self._view.ui.list_analysis_batch_overview,
            self._view.ui.lbl_analysis_batch_prot_struct_1,
            self._view.ui.lbl_analysis_batch_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs,
            self._view.ui.lbl_analysis_batch_ref_chains,
            self._view.ui.list_analysis_batch_ref_chains,
            self._view.ui.btn_analysis_batch_back_2,
            self._view.ui.btn_analysis_batch_next_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_batch_remove,
            self._view.ui.btn_analysis_batch_add,
            self._view.ui.box_analysis_batch_prot_struct_1,
            self._view.ui.box_analysis_batch_prot_struct_2,
            self._view.ui.btn_analysis_batch_back,
            self._view.ui.btn_analysis_batch_next,
            self._view.ui.btn_analysis_batch_back_3,
            self._view.ui.btn_analysis_batch_next_3,
            self._view.ui.lbl_analysis_batch_images,
            self._view.ui.cb_analysis_batch_images,
            self._view.ui.btn_analysis_batch_start,
            self._view.ui.lbl_analysis_batch_model_chains,
            self._view.ui.list_analysis_batch_model_chains,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_analysis_batch_ref_chains.setEnabled(True)

        # tmp_protein = self._current_project.search_protein(self._view.ui.box_analysis_batch_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self._view.ui.list_analysis_batch_ref_chains.addItem(tmp_chain.chain_letter)

    def structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_analysis_batch_remove.setEnabled(True)

    def fill_protein_boxes_batch(self) -> None:
        """Fills the combo boxes with the protein names."""
        proteins = []
        for tmp_protein in self._current_project.proteins:
            proteins.append(tmp_protein.get_molecule_object())
        proteins.insert(0, "")
        self._view.ui.box_analysis_batch_prot_struct_1.clear()
        self._view.ui.box_analysis_batch_prot_struct_2.clear()
        gui_utils.fill_combo_box(self._view.ui.box_analysis_batch_prot_struct_1, proteins)
        gui_utils.fill_combo_box(self._view.ui.box_analysis_batch_prot_struct_2, proteins)

    def remove_analysis_run(self) -> None:
        """Removes the selected protein pair from the list of protein pairs to analyze."""
        self._view.ui.list_analysis_batch_overview.takeItem(self._view.ui.list_analysis_batch_overview.currentRow())
        if self._view.ui.list_analysis_batch_overview.count() == 0:
            gui_elements_to_show = [
                self._view.ui.lbl_analysis_batch_overview,
                self._view.ui.list_analysis_batch_overview,
                self._view.ui.btn_analysis_batch_add,
            ]

            gui_elements_to_hide = [
                self._view.ui.btn_analysis_batch_remove,
                self._view.ui.lbl_analysis_batch_prot_struct_1,
                self._view.ui.lbl_analysis_batch_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs,
                self._view.ui.lbl_analysis_batch_ref_chains,
                self._view.ui.list_analysis_batch_ref_chains,
                self._view.ui.btn_analysis_batch_back_2,
                self._view.ui.btn_analysis_batch_next_2,
                self._view.ui.box_analysis_batch_prot_struct_1,
                self._view.ui.box_analysis_batch_prot_struct_2,
                self._view.ui.btn_analysis_batch_back,
                self._view.ui.btn_analysis_batch_next,
                self._view.ui.lbl_analysis_batch_model_chains,
                self._view.ui.list_analysis_batch_model_chains,
                self._view.ui.btn_analysis_batch_back_3,
                self._view.ui.btn_analysis_batch_next_3,
                self._view.ui.lbl_analysis_batch_images,
                self._view.ui.cb_analysis_batch_images,
                self._view.ui.btn_analysis_batch_start,
            ]

            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_analysis_batch_remove.hide()
        else:
            if self._view.ui.list_analysis_batch_overview.count() > 0:
                try:
                    self._view.ui.list_analysis_batch_overview.currentItem().setSelected(False)
                except AttributeError:
                    constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")
        self._view.ui.btn_analysis_batch_remove.setEnabled(False)

    def check_if_same_no_of_chains_selected_batch(self) -> None:
        """Checks if the same number of proteins were selected."""
        self._view.ui.btn_analysis_batch_next_3.setEnabled(False)
        if self.no_of_selected_chains == len(self._view.ui.list_analysis_batch_model_chains.selectedItems()):
            self._view.ui.btn_analysis_batch_next_3.setEnabled(True)

        prot_1_name = self._view.ui.lbl_analysis_batch_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_analysis_batch_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_analysis_batch_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_analysis_batch_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self._view.ui.list_analysis_batch_overview.count()):
            if analysis_name == self._view.ui.list_analysis_batch_overview.item(tmp_row).text():
                self._view.ui.btn_analysis_batch_next_3.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_analysis_batch_next_3)
                return

    def check_if_prot_structs_are_filled_batch(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self._view.ui.box_analysis_batch_prot_struct_1.itemText(
            self._view.ui.box_analysis_batch_prot_struct_1.currentIndex(),
        )
        prot_2 = self._view.ui.box_analysis_batch_prot_struct_2.itemText(
            self._view.ui.box_analysis_batch_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self._view.ui.btn_analysis_batch_next.setEnabled(True)
        else:
            self._view.ui.btn_analysis_batch_next.setEnabled(False)

    def count_batch_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains of protein 1."""
        self.no_of_selected_chains = len(self._view.ui.list_analysis_batch_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self._view.ui.btn_analysis_batch_next_2.setEnabled(True)
        else:
            self._view.ui.btn_analysis_batch_next_2.setEnabled(False)

    # </editor-fold>

    # <editor-fold desc="Analysis images">
    def display_image_analysis_page(self) -> None:
        """Displays the analysis image work area."""
        # get all protein pairs without images
        self._view.ui.list_analysis_images_struct_analysis.clear()
        self._view.ui.list_analysis_images_creation_struct_analysis.clear()
        for tmp_protein_pair in self._current_project.protein_pairs:
            if len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) == 0:
                self._view.ui.list_analysis_images_struct_analysis.addItem(tmp_protein_pair.name)
        self.last_sidebar_button = styles.color_sidebar_buttons(
            self.last_sidebar_button,
            self._view.ui.btn_image_analysis_page,
        )
        tools.switch_page(self._view.ui.stackedWidget, self._view.ui.lbl_page_title, 23, "Analysis Images")
        self._view.ui.btn_add_analysis_images_struct_analysis.setEnabled(False)
        self._view.ui.btn_remove_analysis_images_creation_struct_analysis.setEnabled(False)
        self._view.ui.btn_start_automatic_image_creation.setEnabled(False)

    def analysis_images_enable_add(self) -> None:
        """Enables the add button."""
        self._view.ui.btn_add_analysis_images_struct_analysis.setEnabled(True)

    def analysis_images_enable_remove(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_remove_analysis_images_creation_struct_analysis.setEnabled(True)

    def add_protein_pair_to_image_creation_queue(self) -> None:
        """Adds a protein pair from the list of protein pairs to make images of."""
        protein_pair_to_add = self._view.ui.list_analysis_images_struct_analysis.currentItem().text()
        self._view.ui.list_analysis_images_creation_struct_analysis.addItem(protein_pair_to_add)
        self._view.ui.list_analysis_images_struct_analysis.takeItem(
            self._view.ui.list_analysis_images_struct_analysis.currentRow(),
        )
        self._view.ui.btn_add_analysis_images_struct_analysis.setEnabled(False)
        self.analysis_images_check_if_creation_can_start()

    def remove_protein_pair_from_image_creation_queue(self) -> None:
        """Removes a protein pair from the list of protein pairs to make images of."""
        protein_pair_to_remove = self._view.ui.list_analysis_images_creation_struct_analysis.currentItem()
        self._view.ui.list_analysis_images_creation_struct_analysis.takeItem(
            self._view.ui.list_analysis_images_creation_struct_analysis.currentRow(),
        )
        self._view.ui.list_analysis_images_struct_analysis.addItem(protein_pair_to_remove)
        self._view.ui.btn_remove_analysis_images_creation_struct_analysis.setEnabled(False)
        self.analysis_images_check_if_creation_can_start()

    def analysis_images_check_if_creation_can_start(self) -> None:
        """Checks if the list of protein pairs which get images are empty."""
        if self._view.ui.list_analysis_images_creation_struct_analysis.count() > 0:
            self._view.ui.btn_start_automatic_image_creation.setEnabled(True)
        else:
            self._view.ui.btn_start_automatic_image_creation.setEnabled(False)

    # </editor-fold>

    # </editor-fold>

    def update_status(self, message: str) -> None:
        """Updates the status bar of the main view with a custom message."""
        self._view.status_bar.showMessage(message)
