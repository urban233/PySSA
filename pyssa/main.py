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
"""Module for the main window of the pyssa plugin"""
import logging
import os
import shutil
import subprocess
import sys
import webbrowser
import pathlib
import PyQt5.QtCore
import numpy as np
import pymol
import csv
import copy
from pyssa.internal.data_structures.data_classes import current_session
from pyssa.util import protein_pair_util
from pyssa.gui.ui.dialogs import dialog_settings_global
from pyssa.gui.ui.dialogs import dialog_startup
from pyssa.util import constants, input_validator, gui_page_management, tools, global_variables, gui_utils
from pyssa.gui.ui.styles import styles
from PyQt5.QtGui import QIcon

from pymol import cmd
# TODO: fix import statements so that they do not import a class!
from urllib.request import urlopen
from urllib.error import URLError
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from pyssa.gui.ui.forms.auto_generated.auto_main_window import Ui_MainWindow
from pyssa.internal.thread import workers
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import filesystem_io
from pyssa.gui.ui.dialogs import dialog_distance_plot
from pyssa.gui.ui.dialogs import dialog_distance_histogram
from pyssa.gui.ui.dialogs import dialog_about
from pyssa.gui.ui.dialogs import dialog_add_models
from pyssa.gui.ui.dialogs import dialog_display_docs
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures import project_watcher
from pyssa.internal.data_structures import settings
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.internal.data_structures.data_classes import stage
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.io_pyssa import path_util
from pyssa.util import pyssa_keys
from pyssa.util import prediction_util
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair

# setup logger
# logging.basicConfig(level=logging.DEBUG)
# global variables
global_var_project_dict = {0: project.Project("", pathlib.Path(""))}
global_var_list_widget_row = 0


class MainWindow(QMainWindow):
    """This class contains all information about the MainWindow in the
    application

    Attributes:
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
    def __init__(self):
        # ----- Initialize the ui building process
        super().__init__()
        # build ui object
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.lbl_page_title.setText("Home")
        self.setMinimumWidth(580)
        self.setMinimumHeight(200)

        # <editor-fold desc="Setup app settings">
        self.app_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        if not os.path.exists(constants.SETTINGS_FULL_FILEPATH):
            dialog = dialog_startup.DialogStartup()
            dialog.exec_()
            # checks if the cancel button was pressed
            if dialog_startup.global_var_terminate_app == 1:
                os.remove(constants.SETTINGS_FULL_FILEPATH)
                sys.exit()
            self.app_settings.app_launch = 1
            self.app_settings.workspace_path = pathlib.Path(dialog_startup.global_var_startup_workspace)
            self.app_settings.serialize_settings()
        try:
            self.app_settings = self.app_settings.deserialize_settings()
        except ValueError:
            print("The settings file is corrupted. Please restore the settings!")
            gui_utils.error_dialog_settings("The settings file is corrupted. Please restore the settings!", "",
                                            self.app_settings)
        # </editor-fold>

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
        self.block_box_expert_install = basic_boxes.no_buttons("Local Colabfold installation", "An installation process is currently running.", QMessageBox.Information)
        # type, name, session
        self.current_session = current_session.CurrentSession("", "", "")

        # </editor-fold>

        # sets up the status bar
        self._setup_statusbar()
        tools.create_directory(constants.SETTINGS_DIR, "scratch")
        self._setup_default_configuration()

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

        self.worker_analysis = workers.AnalysisWorkerPool(
            self.ui.list_analysis_batch_overview, self.ui.cb_analysis_images,
            self.status_bar, self.app_project, self.app_settings, self._init_batch_analysis_page)
        self.worker_analysis.signals.finished.connect(self.post_analysis_process)
        self.worker_analysis.setAutoDelete(True)

        self.worker_image_creation = workers.BatchImageWorkerPool(
            self.ui.list_analysis_images_struct_analysis, self.ui.list_analysis_images_creation_struct_analysis,
            self.status_bar, self.app_project)
        self.worker_image_creation.signals.finished.connect(self.post_image_creation_process)
        self.worker_image_creation.setAutoDelete(True)

        self.block_box_analysis = basic_boxes.no_buttons("Analysis", "A analysis is currently running, please wait.", QMessageBox.Information)
        self.block_box_images = basic_boxes.no_buttons("Analysis Images", "Images getting created, please wait.", QMessageBox.Information)
        # configure gui element properties
        self.ui.txt_results_aligned_residues.setAlignment(QtCore.Qt.AlignRight)
        self.ui.table_pred_mono_prot_to_predict.setSizeAdjustPolicy(PyQt5.QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.ui.table_pred_mono_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self.ui.table_pred_multi_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)

        self.pymol_session_specs = {
            pyssa_keys.SESSION_SPEC_PROTEIN: [0, ""],
            pyssa_keys.SESSION_SPEC_COLOR: [0, ""],
            pyssa_keys.SESSION_SPEC_REPRESENTATION: [0, ""],
            pyssa_keys.SESSION_SPEC_BG_COLOR: [0, ""],
        }

        # helper attributes
        self.used_proteins: list[str] = ["", "", "", ""]

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
        self.last_sidebar_button = PyQt5.QtWidgets.QPushButton()

        # connections
        self._connect_all_gui_elements()
        # create tooltips
        self._create_all_tooltips()
        self._project_watcher.show_valid_options(self.ui)
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

        # fixme: should the pdf documentation be accessible through the pyssa gui?
        self.ui.action_help_docs_pdf.setVisible(False)

        # sets additional parameters
        self.ui.lbl_logo.setPixmap(PyQt5.QtGui.QPixmap(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowTitle(f"PySSA {constants.VERSION_NUMBER}")
        constants.PYSSA_LOGGER.info("PySSA started.")

    # <editor-fold desc="GUI page management functions">
    def _create_local_pred_monomer_management(self):
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
                }
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
                }
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
                }
            ),
            # prediction stage (with advanced configurations)
            stage.Stage(
                {
                    "label_advanced_config": self.ui.lbl_pred_mono_advanced_config,
                    "button_advanced_config": self.ui.btn_pred_mono_advanced_config,
                },
                {
                    "predict_button": self.ui.btn_pred_mono_predict,
                }
            ),
        ]
        self.local_pred_monomer_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_local_pred_multimer_management(self):
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
                }
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
                }
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
                }
            ),
            # prediction stage (with advanced configurations)
            stage.Stage(
                {
                    "label_advanced_config": self.ui.lbl_pred_multi_advanced_config,
                    "button_advanced_config": self.ui.btn_pred_multi_advanced_config,
                },
                {
                    "predict_button": self.ui.btn_pred_multi_predict,
                }
            ),
        ]
        self.local_pred_multimer_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_single_analysis_management(self):
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
                }
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
                }
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
                }
            ),
        ]
        self.single_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_batch_analysis_management(self):
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
                }
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
                }
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
                }
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
                }
            ),
        ]
        self.batch_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_results_management(self):
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
                }
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
                }
            ),
        ]
        self.results_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_monomer_prediction_analysis_management(self):
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
                }
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
                }
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
                }
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
                }
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
                }
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
                }
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
                }
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
                    "back_button":self.ui.btn_pred_analysis_mono_back_pred_setup,
                    "start_button": self.ui.btn_pred_analysis_mono_start,
                }
            ),
        ]
        self.monomer_prediction_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    def _create_multimer_prediction_analysis_management(self):
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
                }
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
                }
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
                }
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
                }
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
                }
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
                }
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
                }
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
                }
            ),
        ]
        self.multimer_prediction_analysis_management = gui_page_management.GuiPageManagement(tmp_stages)

    # </editor-fold>

    def _setup_statusbar(self):
        """This function sets up the status bar and fills it with the current workspace

        """
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage(str(self.workspace_path))

    def _setup_default_configuration(self):
        self.ui.lbl_current_project_name.setText("")
        # menu
        self.ui.action_install_from_file.setVisible(False)
        self.ui.action_add_multiple_models.setVisible(False)
        self.ui.action_file_save_as.setVisible(False)
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
        self.ui.btn_prediction_only_start.setEnabled(False)
        self.ui.btn_prediction_only_next.setEnabled(False)
        self.ui.lbl_prediction_only_status_protein_name.setText("")
        # sequence vs .pdb page
        self.ui.btn_s_v_p_start.setEnabled(False)
        self.ui.list_s_v_p_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        # single analysis page
        self.ui.lbl_analysis_model_chains.hide()
        self.ui.list_analysis_model_chains.hide()
        self.ui.btn_analysis_back.hide()
        self.ui.btn_analysis_start.hide()
        # batch analysis page

        # results page

        # image page

    def _connect_all_gui_elements(self):
        """This function connects all gui elements with their corresponding slots

        """
        # <editor-fold desc="Menu">
        self.ui.action_file_quit.triggered.connect(self.quit_app)
        self.ui.action_file_restore_settings.triggered.connect(self.restore_settings)
        self.ui.action_settings_edit_all.triggered.connect(self.open_settings_global)
        self.ui.action_add_multiple_models.triggered.connect(self.open_add_models)
        self.ui.action_install_from_file.triggered.connect(self.install_local_colabfold_from_file)
        self.ui.action_help_docs.triggered.connect(self.open_documentation)
        self.ui.action_help_docs_pdf.triggered.connect(self.open_documentation_pdf)
        self.ui.action_help_about.triggered.connect(self.open_about)

        # </editor-fold>

        # <editor-fold desc="Side Menu">
        self.ui.btn_new_page.clicked.connect(self.display_new_page)
        self.ui.btn_open_page.clicked.connect(self.display_open_page)
        self.ui.btn_delete_page.clicked.connect(self.display_delete_page)
        self.ui.btn_save_project.clicked.connect(self.save_project)
        self.ui.btn_edit_page.clicked.connect(self.display_edit_page)
        self.ui.btn_view_page.clicked.connect(self.display_view_page)
        self.ui.btn_use_page.clicked.connect(self.display_use_page)
        self.ui.btn_import_project.clicked.connect(self.import_project)
        self.ui.btn_export_project.clicked.connect(self.export_current_project)
        self.ui.btn_close_project.clicked.connect(self.close_project)
        self.ui.btn_pred_local_monomer_page.clicked.connect(self.display_local_pred_mono)
        self.ui.btn_pred_local_multimer_page.clicked.connect(self.display_local_pred_multi)
        self.ui.btn_prediction_abort.clicked.connect(self.abort_prediction)
        self.ui.btn_pred_analysis_monomer_page.clicked.connect(self.display_monomer_pred_analysis)
        self.ui.btn_pred_analysis_multimer_page.clicked.connect(self.display_multimer_pred_analysis)
        self.ui.btn_single_analysis_page.clicked.connect(self.display_single_analysis_page)
        self.ui.btn_batch_analysis_page.clicked.connect(self.display_job_analysis_page)
        self.ui.btn_image_analysis_page.clicked.connect(self.display_image_analysis_page)
        self.ui.btn_results_page.clicked.connect(self.display_results_page)
        self.ui.btn_analysis_abort.clicked.connect(self.abort_analysis)
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
        self.ui.btn_open_open_project.clicked.connect(self.open_project)
        self.ui.list_open_projects.doubleClicked.connect(self.open_project)
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

        # </editor-fold>

        # <editor-fold desc="View project page">
        self.ui.btn_view_project_show.clicked.connect(self.view_sequence)
        self.ui.btn_view_project_show_structure.clicked.connect(self.view_structure)
        self.ui.list_view_project_proteins.doubleClicked.connect(self.view_sequence)

        # </editor-fold>

        # <editor-fold desc="Use project page">
        self.ui.txt_use_project_name.textChanged.connect(self.validate_use_project_name)
        self.ui.btn_use_next.clicked.connect(self.show_protein_selection_for_use)
        self.ui.txt_use_search.textChanged.connect(self.validate_use_search)
        self.ui.btn_use_add_available_protein_structures.clicked.connect(self.add_protein_structure_to_new_project)
        self.ui.list_use_available_protein_structures.doubleClicked.connect(self.add_protein_structure_to_new_project)
        self.ui.btn_use_remove_selected_protein_structures.clicked.connect(self.remove_protein_structure_to_new_project)
        self.ui.list_use_selected_protein_structures.doubleClicked.connect(self.remove_protein_structure_to_new_project)
        self.ui.btn_use_back.clicked.connect(self.hide_protein_selection_for_use)
        self.ui.btn_use_create_new_project.clicked.connect(self.create_use_project)

        # </editor-fold>

        # <editor-fold desc="Monomer local prediction page">
        self.ui.btn_pred_mono_seq_to_predict.clicked.connect(self.local_pred_mono_show_protein_name)
        self.ui.btn_pred_mono_seq_to_predict_remove.clicked.connect(self.local_pred_mono_remove_protein_to_predict)
        self.ui.btn_pred_mono_next.clicked.connect(self.local_pred_mono_show_protein_sequence)
        self.ui.btn_pred_mono_back.clicked.connect(self.local_pred_mono_show_protein_overview)
        self.ui.btn_pred_mono_add_protein.clicked.connect(self.local_pred_mono_add_protein_to_predict)
        self.ui.btn_pred_mono_back_2.clicked.connect(self.local_pred_mono_show_protein_name)
        self.ui.txt_pred_mono_prot_name.textChanged.connect(self.local_pred_mono_validate_protein_name)
        self.ui.txt_pred_mono_seq_name.textChanged.connect(self.local_pred_mono_validate_protein_sequence)
        self.ui.btn_pred_mono_advanced_config.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_mono_predict.clicked.connect(self.predict_local_monomer)

        # </editor-fold>

        # <editor-fold desc="Multimer prediction page">
        self.ui.btn_pred_multi_prot_to_predict_add.clicked.connect(self.local_pred_multi_show_protein_name)
        self.ui.btn_pred_multi_prot_to_predict_remove.clicked.connect(self.local_pred_multi_remove_protein_to_predict)
        self.ui.btn_pred_multi_next.clicked.connect(self.local_pred_multi_show_protein_sequence)
        self.ui.btn_pred_multi_back.clicked.connect(self.local_pred_multi_show_protein_overview)
        self.ui.btn_pred_multi_prot_to_predict_add_2.clicked.connect(self.local_pred_multi_add_protein_to_predict)
        self.ui.btn_pred_multi_back_2.clicked.connect(self.local_pred_multi_show_protein_name)
        self.ui.txt_pred_multi_prot_name.textChanged.connect(self.local_pred_multi_validate_protein_name)
        self.ui.txt_pred_multi_prot_seq.textChanged.connect(self.local_pred_multi_validate_protein_sequence)
        self.ui.btn_pred_multi_prot_seq_add.clicked.connect(self.local_pred_multi_add_sequence_to_list)
        self.ui.btn_pred_multi_prot_seq_overview_remove.clicked.connect(self.local_pred_multi_remove_sequence_to_list)
        self.ui.btn_pred_multi_advanced_config.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_multi_predict.clicked.connect(self.predict_local_multimer)
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
        self.ui.btn_pred_analysis_mono_seq_to_predict.clicked.connect(self.mono_pred_analysis_show_protein_name)
        self.ui.btn_pred_analysis_mono_seq_to_predict_remove.clicked.connect(self.mono_pred_analysis_remove_protein_to_predict)
        self.ui.btn_pred_analysis_mono_next.clicked.connect(self.mono_pred_analysis_show_protein_sequence)
        self.ui.btn_pred_analysis_mono_back.clicked.connect(self.mono_pred_analysis_show_protein_overview)
        self.ui.btn_pred_analysis_mono_add_protein.clicked.connect(self.mono_pred_analysis_add_protein_to_predict)
        self.ui.btn_pred_analysis_mono_back_2.clicked.connect(self.mono_pred_analysis_show_protein_name)
        self.ui.txt_pred_analysis_mono_prot_name.textChanged.connect(self.mono_pred_analysis_validate_protein_name)
        self.ui.txt_pred_analysis_mono_seq_name.textChanged.connect(self.mono_pred_analysis_validate_protein_sequence)
        self.ui.btn_pred_mono_advanced_config_2.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_analysis_mono_go_analysis_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self.ui.btn_pred_analysis_mono_add.clicked.connect(self.show_mono_pred_analysis_stage_1)
        self.ui.btn_pred_analysis_mono_remove.clicked.connect(self.remove_mono_pred_analysis_analysis_run)
        self.ui.btn_pred_analysis_mono_back_3.clicked.connect(self.show_mono_pred_analysis_stage_0)
        self.ui.btn_pred_analysis_mono_next_2.clicked.connect(self.show_mono_pred_analysis_stage_2)
        self.ui.btn_pred_analysis_mono_back_4.clicked.connect(self.show_mono_pred_analysis_stage_1)
        self.ui.btn_pred_analysis_mono_next_3.clicked.connect(self.show_mono_pred_analysis_stage_3)
        self.ui.btn_pred_analysis_mono_back_5.clicked.connect(self.show_mono_pred_analysis_stage_2)
        self.ui.btn_pred_analysis_mono_next_4.clicked.connect(self.show_mono_pred_analysis_stage_0)
        self.ui.box_pred_analysis_mono_prot_struct_1.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled)
        self.ui.box_pred_analysis_mono_prot_struct_2.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled)
        self.ui.list_pred_analysis_mono_ref_chains.itemSelectionChanged.connect(
            self.count_mono_pred_analysis_selected_chains_for_prot_struct_1)
        self.ui.list_pred_analysis_mono_model_chains.itemSelectionChanged.connect(
            self.check_mono_pred_analysis_if_same_no_of_chains_selected)
        self.ui.btn_pred_analysis_mono_back_pred_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self.ui.btn_pred_analysis_mono_start.clicked.connect(self.start_process_batch)

        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Monomer Prediction + Analysis page">
        # <editor-fold desc="Prediction section">
        self.ui.btn_pred_analysis_multi_prot_to_predict_add.clicked.connect(self.multi_pred_analysis_show_protein_name)
        self.ui.btn_pred_analysis_multi_prot_to_predict_remove.clicked.connect(self.multi_pred_analysis_remove_protein_to_predict)
        self.ui.btn_pred_analysis_multi_next.clicked.connect(self.multi_pred_analysis_show_protein_sequence)
        self.ui.btn_pred_analysis_multi_back.clicked.connect(self.multi_pred_analysis_show_protein_overview)
        self.ui.btn_pred_analysis_multi_prot_seq_add.clicked.connect(self.multi_pred_analysis_add_sequence_to_list)
        self.ui.btn_pred_analysis_multi_prot_seq_overview_remove.clicked.connect(self.multi_pred_analysis_remove_sequence_to_list)
        self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.clicked.connect(self.multi_pred_analysis_add_protein_to_predict)
        self.ui.btn_pred_analysis_multi_back_2.clicked.connect(self.multi_pred_analysis_show_protein_name)

        self.ui.txt_pred_analysis_multi_prot_name.textChanged.connect(self.multi_pred_analysis_validate_protein_name)
        self.ui.txt_pred_analysis_multi_prot_seq.textChanged.connect(self.multi_pred_analysis_validate_protein_sequence)
        self.ui.btn_pred_analysis_multi_advanced_config.clicked.connect(self.show_prediction_configuration)
        self.ui.btn_pred_analysis_multi_go_analysis_setup.clicked.connect(self.switch_multimer_pred_analysis_tab)

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self.ui.btn_pred_analysis_multi_add.clicked.connect(self.show_multi_pred_analysis_stage_1)
        self.ui.btn_pred_analysis_multi_remove.clicked.connect(self.remove_multi_pred_analysis_analysis_run)
        self.ui.btn_pred_analysis_multi_back_3.clicked.connect(self.show_multi_pred_analysis_stage_0)
        self.ui.btn_pred_analysis_multi_next_2.clicked.connect(self.show_multi_pred_analysis_stage_2)
        self.ui.btn_pred_analysis_multi_back_4.clicked.connect(self.show_multi_pred_analysis_stage_1)
        self.ui.btn_pred_analysis_multi_next_3.clicked.connect(self.show_multi_pred_analysis_stage_3)
        self.ui.btn_pred_analysis_multi_back_5.clicked.connect(self.show_multi_pred_analysis_stage_2)
        self.ui.btn_pred_analysis_multi_next_4.clicked.connect(self.show_multi_pred_analysis_stage_0)
        self.ui.box_pred_analysis_multi_prot_struct_1.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled)
        self.ui.box_pred_analysis_multi_prot_struct_2.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled)
        self.ui.list_pred_analysis_multi_ref_chains.itemSelectionChanged.connect(
            self.count_multi_pred_analysis_selected_chains_for_prot_struct_1)
        self.ui.list_pred_analysis_multi_model_chains.itemSelectionChanged.connect(
            self.check_multi_pred_analysis_if_same_no_of_chains_selected)
        self.ui.btn_pred_analysis_multi_back_pred_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self.ui.btn_pred_analysis_multi_start.clicked.connect(self.start_process_batch)

        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Batch analysis page">
        self.ui.btn_analysis_batch_add.clicked.connect(self.show_batch_analysis_stage_1)
        self.ui.btn_analysis_batch_remove.clicked.connect(self.remove_analysis_run)
        self.ui.btn_analysis_batch_back.clicked.connect(self.show_batch_analysis_stage_0)
        self.ui.btn_analysis_batch_next.clicked.connect(self.show_batch_analysis_stage_2)
        self.ui.btn_analysis_batch_back_2.clicked.connect(self.show_batch_analysis_stage_1)
        self.ui.btn_analysis_batch_next_2.clicked.connect(self.show_batch_analysis_stage_3)
        self.ui.btn_analysis_batch_back_3.clicked.connect(self.show_batch_analysis_stage_2)
        self.ui.btn_analysis_batch_next_3.clicked.connect(self.show_batch_analysis_stage_0)
        self.ui.box_analysis_batch_prot_struct_1.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch)
        self.ui.box_analysis_batch_prot_struct_2.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch)
        self.ui.list_analysis_batch_ref_chains.itemSelectionChanged.connect(
            self.count_batch_selected_chains_for_prot_struct_1)
        self.ui.list_analysis_batch_model_chains.itemSelectionChanged.connect(
            self.check_if_same_no_of_chains_selected_batch)
        self.ui.btn_analysis_batch_start.clicked.connect(self.start_process_batch)

        # </editor-fold>

        # <editor-fold desc="Analysis images page">
        self.ui.btn_add_analysis_images_struct_analysis.clicked.connect(self.add_protein_pair_to_image_creation_queue)
        self.ui.list_analysis_images_struct_analysis.doubleClicked.connect(self.add_protein_pair_to_image_creation_queue)
        self.ui.btn_remove_analysis_images_creation_struct_analysis.clicked.connect(self.remove_protein_pair_from_image_creation_queue)
        self.ui.list_analysis_images_creation_struct_analysis.doubleClicked.connect(self.remove_protein_pair_from_image_creation_queue)
        self.ui.btn_start_automatic_image_creation.clicked.connect(self.start_automatic_image_creation)

        # </editor-fold>

        # <editor-fold desc="Results page">
        self.ui.cb_results_analysis_options.currentIndexChanged.connect(self.load_results)
        self.ui.btn_color_rmsd.clicked.connect(self.color_protein_pair_by_rmsd)
        self.ui.btn_view_struct_alignment.clicked.connect(self.display_structure_alignment)
        self.ui.btn_view_distance_plot.clicked.connect(self.display_distance_plot)
        self.ui.btn_view_distance_histogram.clicked.connect(self.display_distance_histogram)
        self.ui.btn_view_interesting_region.clicked.connect(self.display_interesting_region)
        self.ui.btn_view_distance_table.clicked.connect(self.display_distance_table)

        # </editor-fold>

        # <editor-fold desc="Manage">
        self.ui.box_manage_choose_color.activated.connect(self.choose_manage_color_selected_protein)
        self.ui.box_manage_choose_representation.activated.connect(self.choose_manage_representation)
        self.ui.box_manage_choose_bg_color.activated.connect(self.choose_manage_bg_color)

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

        # </editor-fold>

        # <editor-fold desc="Hotspots page">
        self.ui.list_hotspots_choose_protein.currentItemChanged.connect(self.open_protein)
        self.ui.btn_hotspots_resi_show.clicked.connect(self.show_resi_sticks)
        self.ui.btn_hotspots_resi_hide.clicked.connect(self.hide_resi_sticks)
        self.ui.btn_hotspots_resi_zoom.clicked.connect(self.zoom_resi_position)

        # </editor-fold>

    def _create_all_tooltips(self):
        # menu

        # side menu

        # new project page

        # open project page

        # delete project page

        # edit project page

        # view project page

        # use project page

        # new sequence page

        # sequence vs .pdb page
        self.ui.btn_prediction_load_reference.setToolTip("Open reference pdb file")
        self.ui.btn_prediction_start.setToolTip("Predict with Colab Notebook")
        self.ui.txt_prediction_load_reference.setToolTip("Reference file path")
        # self.ui.txt_prediction_chain_ref.setToolTip("Enter chain(s) of reference")
        self.ui.txt_prediction_chain_model.setToolTip("Enter chain(s) of model")
        # single analysis page
        self.ui.btn_analysis_start.setToolTip("Start analysis process")
        self.status_bar.setToolTip("Status information: Current process")
        # batch analysis page
        self.status_bar.setToolTip("Status information: Current process")
        # results page

        # image page
        self.ui.btn_save_scene.setToolTip("Create new PyMOL scene")
        self.ui.btn_update_scene.setToolTip("Overwrite current scene")
        self.ui.btn_save_image.setToolTip("Save current viewpoint as png file")
        self.ui.cb_ray_tracing.setToolTip("Enable ray-tracing")
        self.ui.cb_transparent_bg.setToolTip("Enable transparent background")
        self.ui.box_representation.setToolTip("Choose a representation")
        self.ui.box_bg_color.setToolTip("Choose a background color")
        self.ui.box_renderer.setToolTip("Choose a ray-tracing renderer")
        self.ui.box_ray_trace_mode.setToolTip("Choose a ray-trace mode")

    # <editor-fold desc="Page init functions">
    def _init_fill_combo_boxes(self):
        """This function fills all combo boxes of the plugin

        """
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

    def _init_new_page(self):
        """This function clears all text fields and hides everything which is needed

        """
        self.ui.txt_new_project_name.clear()
        self.ui.txt_new_choose_reference.clear()
        self.ui.lbl_new_status_project_name.setText("")
        self.ui.lbl_new_status_choose_reference.setText("")
        self.ui.cb_new_add_reference.setCheckState(0)
        self.ui.btn_new_create_project.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_new_create_project)

    def _init_use_page(self):
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

    def _init_edit_page(self):
        self.ui.list_edit_project_proteins.clear()
        gui_elements_to_hide = [
            self.ui.lbl_edit_clean_new_prot,
            self.ui.btn_edit_clean_new_prot,
            self.ui.lbl_edit_clean_update_prot,
            self.ui.btn_edit_clean_update_prot,
            self.ui.label_12,
            self.ui.btn_edit_project_delete,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_edit_project_proteins)
        # self.project_scanner.scan_project_for_valid_proteins(self.ui.list_edit_project_proteins)

    def _init_local_pred_mono_page(self):
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

    def _init_local_pred_multi_page(self):
        # clears everything
        self.ui.txt_pred_multi_prot_name.clear()
        self.ui.txt_pred_multi_prot_seq.clear()
        self.ui.list_pred_multi_prot_seq_overview.clear()
        # sets up defaults: Prediction
        self.ui.btn_pred_multi_next.setEnabled(False)
        self.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        self.ui.lbl_pred_multi_prot_name_status.setText("")
        self.ui.lbl_pred_multi_prot_seq_status.setText("")
        # self.local_pred_multimer_management.show_stage_x(0)
        # self.local_pred_multimer_management.disable_all_next_buttons()
        # self.local_pred_multimer_management.clear_all_text_boxes()
        # self.local_pred_multimer_management.set_empty_string_in_label()

    def _init_mono_pred_analysis_page(self):
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

    def _init_multi_pred_analysis_page(self):
        # <editor-fold desc="Prediction section">
        # clears everything
        self.ui.txt_pred_analysis_multi_prot_name.clear()
        self.ui.txt_pred_analysis_multi_prot_seq.clear()
        self.ui.list_pred_analysis_multi_prot_seq_overview.clear()
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

    def _init_sequence_vs_pdb_page(self):
        """This function clears all text fields and hides everything which is needed

        """
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

    def _init_single_analysis_page(self):
        """This function clears all text fields and hides everything which is needed

        """
        self.single_analysis_management.show_stage_x(0)
        self.single_analysis_management.disable_all_next_buttons()
        self.single_analysis_management.set_empty_string_in_label()
        self.ui.lbl_analysis_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_analysis_prot_struct_2.setText("Protein structure 2")
        self.ui.box_analysis_prot_struct_1.setCurrentIndex(0)
        self.ui.box_analysis_prot_struct_2.setCurrentIndex(0)

    def _init_results_page(self):
        """This function clears all text fields and hides everything which is needed

        """
        # stage 1
        self.ui.list_results_interest_regions.clear()
        self.ui.txt_results_rmsd.clear()
        self.ui.txt_results_aligned_residues.clear()

    def _init_analysis_image_page(self):
        self.ui.list_analysis_images_struct_analysis.setEnabled(True)
        self.ui.list_analysis_images_creation_struct_analysis.setEnabled(True)
        self.display_image_analysis_page()

    def _init_image_page(self):
        """This function clears all text fields and hides everything which is needed

        """
        # stage 1
        self.ui.box_representation.setCurrentIndex(0)
        self.ui.box_bg_color.setCurrentIndex(0)
        self.ui.box_renderer.setCurrentIndex(0)
        self.ui.box_ray_trace_mode.setCurrentIndex(0)
        self.ui.box_ray_texture.setCurrentIndex(0)
        self.ui.cb_ray_tracing.setChecked(False)
        self.ui.cb_transparent_bg.setChecked(False)

    def _init_batch_analysis_page(self):
        # sets up defaults: Batch
        self.batch_analysis_management.show_stage_x(0)
        self.ui.list_analysis_batch_overview.clear()
        self.ui.btn_analysis_batch_remove.hide()

    def _init_all_pages(self):
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
    def display_home_page(self):
        """This function displays the homepage of the plugin

        """
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 0, "Home")

    def display_prediction_and_analysis_page(self):
        """This function displays the prediction + analysis work area

        """
        self.ui.list_widget_projects.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        workspace_projects: list[str] = os.listdir(self.workspace_path)
        workspace_projects.sort()
        for project in workspace_projects:
            self.ui.list_widget_projects.addItem(project)

        # regular opening of work area
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 2, "Sequence vs .pdb")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_prediction_analysis_page)

    def display_sequence_vs_pdb_page(self):
        """This function displays the sequence vs .pdb page

        """
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 10, "Sequence vs .pdb")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_prediction_analysis_page)

    def display_single_analysis_page(self):
        """This function displays the single analysis work area

        """
        self.fill_protein_structure_boxes()
        self.ui.list_analysis_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.list_analysis_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        # regular work area opening
        self._init_single_analysis_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 3, "Single Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_single_analysis_page)

    def display_job_analysis_page(self):
        """This function displays the job analysis work area

        """
        self.ui.list_analysis_batch_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.list_analysis_batch_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        # regular work area opening
        self._init_batch_analysis_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 4, "Structure Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_batch_analysis_page)

    def display_results_page(self):
        """This function displays the results work area

        """
        results = []
        results.insert(0, "")
        current_results_index = 0
        i = 0
        for tmp_protein_pair in self.app_project.protein_pairs:
            results.append(tmp_protein_pair.name)
            if tmp_protein_pair.name == self.results_name:
                current_results_index = i
            i += 1
        self.ui.cb_results_analysis_options.clear()
        gui_utils.fill_combo_box(self.ui.cb_results_analysis_options, results)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 5, "Results")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_results_page)
        self.ui.cb_results_analysis_options.setCurrentIndex(current_results_index + 1)

    def display_image_analysis_page(self):
        """This function displays the analysis image work area

        """
        # get all protein pairs without images
        self.ui.list_analysis_images_struct_analysis.clear()
        self.ui.list_analysis_images_creation_struct_analysis.clear()
        for tmp_protein_pair in self.app_project.protein_pairs:
            if len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) == 0:
                self.ui.list_analysis_images_struct_analysis.addItem(tmp_protein_pair.name)
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_image_analysis_page)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 23, "Analysis Images")

    def display_image_page(self):
        """This function displays the image work area

        """
        self._init_image_page()
        if self.ui.box_renderer.currentText() == "":
            self.ui.cb_ray_tracing.hide()
            self.ui.label_26.hide()
        else:
            self.ui.cb_ray_tracing.show()
            self.ui.label_26.show()
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_image_page)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 6, "Image")

    def display_new_page(self):
        """This function displays the new project work area

        """
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

    def display_open_page(self):
        """This function displays the open project work area

        """
        self.ui.txt_open_search.clear()
        self.ui.txt_open_selected_project.clear()
        if safeguard.Safeguard.check_filepath(self.workspace_path):
            self.ui.list_open_projects.clear()
            # pre-process
            self.status_bar.showMessage(self.workspace.text())
            try:
                tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_open_projects)
            except PermissionError:
                gui_utils.error_dialog_settings("The settings file is corrupted. Please restore the settings!", "",
                                                self.app_settings)
                self.display_home_page()
                return
            tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 8, "Open existing project")
            self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_open_page)
        else:
            gui_utils.error_dialog_settings("The settings file is corrupted. Please restore the settings!", "",
                                            self.app_settings)
            self.display_home_page()

    def display_delete_page(self):
        """This function displays the "delete" project work area

        """
        self.ui.txt_delete_search.clear()
        self.ui.txt_delete_selected_projects.clear()
        self.ui.list_delete_projects.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_delete_projects)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 9, "Delete existing project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_delete_page)

    def display_edit_page(self):
        """This function displays the edit project page

        """
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        self._init_edit_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 13, "Edit proteins of current project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_edit_page)

    def display_view_page(self):
        """This function displays the edit project page

        """
        self.ui.list_view_project_proteins.clear()
        self.ui.txtedit_view_sequence.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        # list all proteins from pdb directory
        gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_view_project_proteins)
        # self.project_scanner.scan_project_for_valid_proteins(list_view_project_proteins=self.ui.list_view_project_proteins)

        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 11, "View proteins of current project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_view_page)

    def display_local_pred_mono(self):
        self._init_local_pred_mono_page()
        self.local_pred_monomer_management.show_stage_x(0)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 19, "Local Monomer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_pred_local_monomer_page)

    def display_local_pred_multi(self):
        self.local_pred_multimer_management.show_stage_x(0)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 20, "Local Multimer Prediction")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_pred_local_multimer_page)

    def display_use_page(self):
        self._init_use_page()
        self.ui.list_use_available_protein_structures.clear()
        self.ui.list_use_selected_protein_structures.clear()
        valid_projects = tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_use_existing_projects)
        # filesystem operations
        gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_use_selected_protein_structures)
        # self.project_scanner.scan_project_for_valid_proteins(self.ui.list_use_selected_protein_structures)
        protein_dict, protein_names = tools.scan_workspace_for_non_duplicate_proteins(self.workspace_path)
        global_variables.global_var_workspace_proteins = protein_dict
        # this for-loop is necessary for eliminating all proteins which are in the current project from the ones which
        # are available
        for i in range(self.ui.list_use_selected_protein_structures.count()):
            self.ui.list_use_selected_protein_structures.setCurrentRow(i)
            tmp_prot_name = self.ui.list_use_selected_protein_structures.currentItem().text()
            if tmp_prot_name in protein_names:
                protein_names.remove(tmp_prot_name)

        self.ui.list_use_available_protein_structures.addItems(protein_names)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 14, "Use existing project")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_use_page)

    def display_hotspots_page(self):
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
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_hotspots_page)

    def display_monomer_pred_analysis(self):
        self._init_mono_pred_analysis_page()
        self.monomer_prediction_analysis_management.show_stage_x(0)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 21, "Monomer Prediction + Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_pred_analysis_monomer_page)

    def display_multimer_pred_analysis(self):
        self._init_multi_pred_analysis_page()
        self.multimer_prediction_analysis_management.show_stage_x(0)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 22, "Multimer Prediction + Analysis")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_pred_analysis_multimer_page)

    def display_manage_pymol_session(self):
        pymol_objs = cmd.get_object_list()
        pymol_objs.insert(0, "")
        for tmp_object in pymol_objs:
            self.ui.box_manage_choose_protein.addItem(tmp_object)
        self.ui.box_manage_choose_protein.setCurrentIndex(self.pymol_session_specs[pyssa_keys.SESSION_SPEC_PROTEIN][0])
        self.ui.box_manage_choose_color.setCurrentIndex(self.pymol_session_specs[pyssa_keys.SESSION_SPEC_COLOR][0])
        self.ui.box_manage_choose_representation.setCurrentIndex(self.pymol_session_specs[pyssa_keys.SESSION_SPEC_REPRESENTATION][0])
        self.ui.box_manage_choose_bg_color.setCurrentIndex(self.pymol_session_specs[pyssa_keys.SESSION_SPEC_BG_COLOR][0])
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 24, "Manage PyMOL session")
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_manage_session)

    # </editor-fold>

    def restore_settings(self):
        """This function deletes the old settings.xml and creates a new one,
        with the default values.

        """
        out = gui_utils.warning_dialog_restore_settings("Are you sure you want to restore all settings?")
        if out:
            tools.restore_default_settings(self.app_settings)
            self.status_bar.showMessage("Settings were successfully restored.")
            logging.info("Settings were successfully restored.")
        else:
            self.status_bar.showMessage("Settings were not modified.")
            logging.info("Settings were not modified.")

    def quit_app(self):
        """This function closes the entire plugin.

        """
        self.close()

    def open_settings_global(self):
        """This function open the dialog for the global settings.

        """
        # TODO: Wsl2 und Localcolabfold btn -> Is it on machine? -> Don't show btn
        # wsl2 btn
        # subprocess.run("(gcm wsl).Version"): # gives version of wsl
        # PS C:\Users\hannah> wsl hostname -i 127.0.1.1


        # local colabfold # TODO: make it possible that it grap on btn
        # "/home/$USER/.pyssa/colabfold_batch/bin/colabfold_batch"
        # if os.path.exists("/home/$USER/.pyssa/colabfold_batch/bin/colabfold_batch"):
        #     self.ui.btn_install_local_prediction.hide()
        #
        # else:
        #     self.ui.btn_install_local_prediction.show()
        dialog = dialog_settings_global.DialogSettingsGlobal()
        dialog.exec_()
        self.app_settings = self.app_settings.deserialize_settings()
        self.workspace_path = self.app_settings.workspace_path
        self._setup_statusbar()

    @staticmethod
    def open_documentation():
        """This function opens the official plugin documentation as HTML page.

        """
        dialog = dialog_display_docs.WebViewDialog(str(constants.DOCS_HTML))
        dialog.exec_()

    @staticmethod
    def open_documentation_pdf():
        """This function opens the official plugin documentation as PDF.

        """
        # webbrowser.open_new(
        #     f"file://{os.getcwd()}/docs/pyssa/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")
        # opens the documentation of the os
        if sys.platform.startswith("darwin"):
            # macOS path
            webbrowser.open_new(f"file://{constants.path_list[1]}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")
        elif sys.platform.startswith("linux"):
            # Linux path
            webbrowser.open_new(f"file://{constants.path_list[0]}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")
        elif sys.platform.startswith("win32"):
            # Windows path
            webbrowser.open_new(f"file://{constants.path_list[2]}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")

    @staticmethod
    def open_about():
        """This function opens the about dialog.

        """
        dialog = dialog_about.DialogAbout()
        dialog.exec_()

    def open_add_models(self):
        """This function opens the add models dialog.

        """
        dialog = dialog_add_models.DialogAddModels()
        dialog.exec_()

        if len(dialog_add_models.global_var_pdb_files) > 0:
            for pdb_path in dialog_add_models.global_var_pdb_files:
                pdb_path_info = QtCore.QFileInfo(pdb_path)
                pdb_name = pdb_path_info.baseName()
                # save project folder in current workspace
                # mkdir project
                folder_paths = [
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/pdb"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/alignment_files"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/distance_csv"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/images"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/images/interesting_regions"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/plots"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/plots/distance_histogram"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/plots/distance_plot"),
                    pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/sessions/"),
                ]
                for path in folder_paths:
                    os.mkdir(path)
                reference_name = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")
                shutil.copy(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{reference_name[0]}",
                            f"{folder_paths[1]}/{reference_name[0]}")
                shutil.copy(pdb_path, f"{folder_paths[1]}/{pdb_name}.pdb")

    def post_install_local_colabfold_from_file(self):
        self.block_box_expert_install.destroy(True)
        self.app_settings.local_colabfold = 1
        self.app_settings.serialize_settings()
        basic_boxes.ok("Local Colabfold installation", "Installation is finished!", QMessageBox.Information)

    def install_local_colabfold_from_file(self):
        """This function installs the local colabfold from a .tar file

        """
        file_dialog = QFileDialog()
        download_path = PyQt5.QtCore.QStandardPaths.standardLocations(PyQt5.QtCore.QStandardPaths.DownloadLocation)[0]
        file_dialog.setDirectory(download_path)
        file_path, _ = file_dialog.getOpenFileName(self, "Select the UbuntuColabfold.tar file", "", "Colabfold WSL tar (UbuntuColabfold.tar)")
        if file_path:
            install_worker = workers.ColabfoldInstallerWorkerPool(True)
            install_worker.signals.finished.connect(self.post_install_local_colabfold_from_file)
            if basic_boxes.yes_or_no("Local Colabfold installation", "Are you sure that you want to install Local Colabfold?", QMessageBox.Question) is True:
                install_worker.install = True
                install_worker.local_install = True, file_path
                self.threadpool.start(install_worker)
                self.block_box_expert_install.exec_()
            else:
                # logical message: the user does NOT want to install local colabfold
                basic_boxes.ok("Local Colabfold installation", "Installation process aborted.", QMessageBox.Information)
                return

    # <editor-fold desc="New project page functions">
    def show_add_reference(self):
        """This function shows the reference input section

        """
        # checkbox is checked
        test = self.ui.cb_new_add_reference.checkState()
        if self.ui.cb_new_add_reference.checkState() == 2:
            self.ui.txt_new_choose_reference.clear()
            self.ui.txt_new_choose_reference.setStyleSheet("background-color: white")
            self.ui.lbl_new_choose_reference.show()
            self.ui.txt_new_choose_reference.show()
            self.ui.btn_new_choose_reference.show()
            self.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_new_create_project)
            # check internet connectivity
            timeout: float = 3
            try:
                urlopen('https://www.google.com', timeout=timeout)
            except URLError:
                self.ui.txt_new_choose_reference.setEnabled(False)
                self.ui.lbl_new_status_choose_reference.setText("You cannot enter a PDB ID (no internet).")
                return
            self.ui.txt_new_choose_reference.setEnabled(True)

        else:
            self.ui.lbl_new_choose_reference.hide()
            self.ui.txt_new_choose_reference.hide()
            self.ui.btn_new_choose_reference.hide()
            self.ui.lbl_new_status_choose_reference.setText("")
            self.ui.btn_new_create_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_new_create_project)

    def load_reference_in_project(self):
        """This function loads a reference in a new project

        """
        try:
            # open file dialog
            file_name = QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                 QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
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

    def validate_reference_in_project(self):
        """This function checks if the entered reference protein is
        valid or not.

        """
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

    def validate_project_name(self):
        """This function validates the input of the project name in real-time

        """
        input_validator.InputValidator.validate_project_name(
            self.ui.list_new_projects, self.ui.txt_new_project_name,
            self.ui.lbl_new_status_project_name, self.ui.btn_new_create_project,
            self.ui.cb_new_add_reference
        )

    def create_new_project(self):
        """This function creates a new project based on the plugin New ... page

        """
        if self.app_settings.wsl_install == 0:
            basic_boxes.ok("Create new project", "Please install local colabfold to create a project!", QMessageBox.Warning)
            return
        elif self.app_settings.local_colabfold == 0:
            basic_boxes.ok("Create new project", "Please install local colabfold to create a project!",
                           QMessageBox.Warning)
            return
        self._project_watcher.on_home_page = False
        self.ui.lbl_current_project_name.setText(self.ui.txt_new_project_name.text())
        self.status_bar.showMessage(f"Current project path: {self.workspace_path}/{self.ui.txt_new_project_name.text()}")
        # save project folder in current workspace
        self.app_project = project.Project(self.ui.txt_new_project_name.text(), self.workspace_path)
        # self.app_project.create_project_tree()
        # save reference .pdb
        if self.ui.cb_new_add_reference.checkState() == 2 and self.ui.btn_new_create_project.isEnabled() is True:
            if len(self.ui.txt_new_choose_reference.text()) == 4:
                # PDB ID as input
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                pdb_id = self.ui.txt_new_choose_reference.text().upper()
                try:
                    cmd.fetch(pdb_id, type="pdb", path=self.scratch_path)
                except pymol.CmdException:
                    tools.clean_scratch_folder()
                    # TODO: add message that fetching the reference failed
                    return

                tmp_ref_protein = protein.Protein(molecule_object=pdb_id,
                                                  pdb_filepath=path_util.FilePath(pathlib.Path(f"{self.scratch_path}/{pdb_id}.pdb")))
            else:
                # local pdb file as input
                pdb_filepath = path_util.FilePath(pathlib.Path(self.ui.txt_new_choose_reference.text()))
                tmp_ref_protein = protein.Protein(molecule_object=pdb_filepath.get_filename(),
                                                  pdb_filepath=pdb_filepath)
            self.app_project.add_existing_protein(tmp_ref_protein)
        self.ui.cb_new_add_reference.setCheckState(0)

        self.app_project.serialize_project(pathlib.Path(f"{self.workspace_path}/{self.app_project.get_project_name()}.xml"))
        constants.PYSSA_LOGGER.info(f"Created the project {self.app_project.get_project_name()}.")
        protein_names = []
        for tmp_protein in self.app_project.proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        constants.PYSSA_LOGGER.debug(f"These are the proteins {protein_names}.")
        # shows options which can be done with the data in the project folder
        self._project_watcher.current_project = self.app_project
        self.project_scanner.project = self.app_project
        constants.PYSSA_LOGGER.info(f"{self._project_watcher.current_project.get_project_name()} is the current project.")
        self._project_watcher.show_valid_options(self.ui)
        self.display_view_page()

    # </editor-fold>

    # <editor-fold desc="Open project page functions">
    def validate_open_search(self):
        """This function validates the input of the project name in real-time

        """
        if self.ui.list_open_projects.currentItem() is not None:
            self.ui.list_open_projects.currentItem().setSelected(False)
        # set color for lineEdit
        input_validator.InputValidator.validate_search_input(
            self.ui.list_open_projects, self.ui.txt_open_search,
            self.ui.lbl_open_status_search, self.ui.txt_open_selected_project
        )

    def select_project_from_open_list(self):
        try:
            self.ui.txt_open_selected_project.setText(self.ui.list_open_projects.currentItem().text())
        except AttributeError:
            self.ui.txt_open_selected_project.setText("")

    def activate_open_button(self):
        """This function is used to activate the open button

        """
        if self.ui.txt_open_selected_project.text() == "":
            self.ui.btn_open_open_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_open_open_project)
        else:
            self.ui.btn_open_open_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_open_open_project)

    def open_project(self):
        """This function opens an existing project

        """
        # show project management options in side menu
        tmp_project_path = pathlib.Path(f"{self.workspace_path}/{self.ui.list_open_projects.currentItem().text()}")
        self.app_project = project.Project.deserialize_project(tmp_project_path, self.app_settings)
        constants.PYSSA_LOGGER.info(f"Opening the project {self.app_project.get_project_name()}.")
        self._project_watcher.current_project = self.app_project
        self.project_scanner.project = self.app_project
        constants.PYSSA_LOGGER.info(f"{self._project_watcher.current_project.get_project_name()} is the current project.")
        self.ui.lbl_current_project_name.setText(self.app_project.get_project_name())
        self._project_watcher.on_home_page = False
        self._project_watcher.show_valid_options(self.ui)
        cmd.reinitialize()
        self.display_view_page()

    # </editor-fold>

    # <editor-fold desc="Delete project page functions">
    def select_project_from_delete_list(self):
        """This function selects a project from the project list on the delete page

        """
        try:
            self.ui.txt_delete_selected_projects.setText(self.ui.list_delete_projects.currentItem().text())
        except AttributeError:
            self.ui.txt_delete_selected_projects.setText("")

    def activate_delete_button(self):
        """This function is used to activate the open button

        """
        if self.ui.txt_delete_selected_projects.text() == "":
            self.ui.btn_delete_delete_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_delete_delete_project)
        else:
            self.ui.btn_delete_delete_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_delete_delete_project)

    def validate_delete_search(self):
        """This function validates the input of the project name in real-time

        """

        if self.ui.list_delete_projects.currentItem() is not None:
            self.ui.list_delete_projects.currentItem().setSelected(False)
        # set color for lineEdit
        input_validator.InputValidator.validate_search_input(
            self.ui.list_delete_projects, self.ui.txt_delete_search,
            self.ui.lbl_delete_status_search, self.ui.txt_delete_selected_projects
        )

    def delete_project(self):
        """This function deletes an existing project

        """
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
        else:
            constants.PYSSA_LOGGER.info("No project has been deleted. No changes were made.")
            return

    # </editor-fold>

    # <editor-fold desc="Save project functions">
    def save_project(self):
        """This function saves the "project" which is currently only the pymol session

        """
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_save_project)
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        self.app_project.serialize_project(self.app_project.get_project_xml_path())

    # </editor-fold>

    # <editor-fold desc="Edit project page functions">
    def check_for_cleaning(self):
        try:
            protein_name = self.ui.list_edit_project_proteins.currentItem().text()
        except AttributeError:
            return
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        cmd.reinitialize()
        tmp_protein = self.app_project.search_protein(
            protein_name.replace(".pdb", ""))
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

    def clean_protein_new(self):
        tmp_protein = self.app_project.search_protein(self.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""))
        clean_tmp_protein = tmp_protein.clean_protein(new_protein=True)
        constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", clean_tmp_protein.get_molecule_object())
        self.app_project.add_existing_protein(clean_tmp_protein)
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        self._init_edit_page()

    def clean_protein_update(self):
        tmp_protein = self.app_project.search_protein(self.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""))
        if basic_boxes.yes_or_no("Clean protein",
                                 "Are you sure you want to clean this protein?\n"
                                 "This will remove all organic and solvent components!",
                                 QMessageBox.Information):
            tmp_protein.clean_protein()
            constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", tmp_protein.get_molecule_object())
            self.app_project.serialize_project(self.app_project.get_project_xml_path())
            self._init_edit_page()
        else:
            constants.PYSSA_LOGGER.info("No protein has been cleaned.")
            return

    def delete_protein(self):
        protein_name = self.ui.list_edit_project_proteins.currentItem().text()
        response = gui_utils.warning_message_protein_gets_deleted()
        if response:
            try:
                self.app_project.delete_specific_protein(protein_name)
            except ValueError:
                constants.PYSSA_LOGGER.error("The protein %s could not be deleted, because it is not in the project.", protein_name)
            self.ui.list_edit_project_proteins.clear()
            gui_utils.fill_list_view_with_protein_names(self.app_project, self.ui.list_edit_project_proteins)
            self._project_watcher.show_valid_options(self.ui)
        else:
            constants.PYSSA_LOGGER.info("No protein was deleted.")

    # </editor-fold>

    # <editor-fold desc="View project page functions">
    def view_sequence(self):
        tmp_protein_basename = self.ui.list_view_project_proteins.currentItem().text()
        tmp_protein_sequences = self.app_project.search_protein(tmp_protein_basename).get_protein_sequences()
        self.ui.txtedit_view_sequence.clear()
        for tmp_sequence in tmp_protein_sequences:
            self.ui.txtedit_view_sequence.append("".join(tmp_sequence.sequence))
        # fixme: experimental sequence viewer gui
        # dialog = dialog_sequence_viewer.SequenceViewer(tmp_protein_sequences, tmp_protein_filename)
        # dialog.exec_()

    def view_structure(self):
        protein_name = self.ui.list_view_project_proteins.currentItem().text()
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        cmd.reinitialize()
        try:
            self.app_project.search_protein(protein_name).load_protein_pymol_session()
            constants.PYSSA_LOGGER.info("Loaded PyMOL session of protein %s", protein_name)
        except pymol.CmdException:
            constants.PYSSA_LOGGER.error("Error while loading protein in PyMOL!")
            return
        self.current_session = current_session.CurrentSession("protein", protein_name, self.app_project.search_protein(protein_name).load_protein_pymol_session())

    # </editor-fold>

    # <editor-fold desc="Use project page functions">
    def validate_use_project_name(self):
        """This function validates the input of the project name in real-time

        """
        input_validator.InputValidator.validate_project_name(
            self.ui.list_use_existing_projects, self.ui.txt_use_project_name,
            self.ui.lbl_use_status_project_name, self.ui.btn_use_next
        )

    def validate_use_search(self):
        """This function validates the input of the project name in real-time

        """
        message = "Protein structure does not exists."
        input_validator.InputValidator.validate_search_input(
            self.ui.list_use_available_protein_structures, self.ui.txt_use_search,
            self.ui.lbl_use_status_search, status_message=message
        )

    def add_protein_structure_to_new_project(self):
        prot_to_add = self.ui.list_use_available_protein_structures.currentItem().text()
        self.ui.list_use_selected_protein_structures.addItem(prot_to_add)
        self.ui.list_use_available_protein_structures.takeItem(self.ui.list_use_available_protein_structures.currentRow())

    def remove_protein_structure_to_new_project(self):
        prot_to_remove = self.ui.list_use_selected_protein_structures.currentItem()
        self.ui.list_use_selected_protein_structures.takeItem(self.ui.list_use_selected_protein_structures.currentRow())
        self.ui.list_use_available_protein_structures.addItem(prot_to_remove)

    def show_protein_selection_for_use(self):
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

    def hide_protein_selection_for_use(self):
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

    def create_use_project(self):
        self.ui.lbl_current_project_name.setText(self.ui.txt_use_project_name.text())
        self.status_bar.showMessage(
            f"Current project path: {self.workspace_path}/{self.ui.txt_use_project_name.text()}")
        # save project folder in current workspace
        # existing_project = self.app_project
        new_project = project.Project(self.ui.txt_use_project_name.text(), self.workspace_path)
        # new_project.create_project_tree()
        self.app_project = new_project
        # copy proteins in new project
        proteins_to_copy = []
        for i in range(self.ui.list_use_selected_protein_structures.count()):
            self.ui.list_use_selected_protein_structures.setCurrentRow(i)
            proteins_to_copy.append(self.ui.list_use_selected_protein_structures.currentItem().text())
        for tmp_protein in proteins_to_copy:
            # protein_path = global_variables.global_var_workspace_proteins[tmp_protein] # TODO: global_var needs to be removed
            # shutil.copy(protein_path, f"{self.app_project.folder_paths["proteins"]}/{tmp_protein}")
            protein_infos = (tools.scan_workspace_for_non_duplicate_proteins(self.workspace_path))[0]
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


            #
            # new_protein = protein.Protein(molecule_object=tmp_protein,
            #                               pdb_filepath=path_util.FilePath(workspace_scanner.scan_workspace_for_non_duplicate_proteins().get(tmp_protein)))
            #new_protein.serialize_protein()
            self.app_project.add_existing_protein(tmp_protein_obj)
        self.app_project.serialize_project(pathlib.Path(f"{self.workspace_path}/{self.app_project.get_project_name()}.xml"))
        # shows options which can be done with the data in the project folder
        self._project_watcher.current_project = self.app_project
        self._project_watcher.on_home_page = False
        self._project_watcher.show_valid_options(self.ui)
        self.project_scanner.project = self.app_project
        self._init_use_page()
        constants.PYSSA_LOGGER.info(f"The project {self.app_project.get_project_name()} was successfully created through a use.")
        self.display_view_page()

    # </editor-fold>

    # <editor-fold desc="Import, Export functions">
    def import_project(self):
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_import_project)
        file_dialog = QFileDialog()
        desktop_path = PyQt5.QtCore.QStandardPaths.standardLocations(PyQt5.QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getOpenFileName(self, "Select a project file to import", "", "XML Files (*.xml)")
        if file_path:
            tmp_project = project.Project("", self.workspace_path)
            tmp_project = tmp_project.deserialize_project(file_path, self.app_settings)
            tmp_project.set_workspace_path(self.workspace_path)
            if len(tmp_project.proteins) <= 1:
                if self.app_settings.wsl_install == 0:
                    basic_boxes.ok("Create new project", "Please install local colabfold to import this project!",
                                   QMessageBox.Warning)
                    return
                elif self.app_settings.local_colabfold == 0:
                    basic_boxes.ok("Create new project", "Please install local colabfold to import this project!",
                                   QMessageBox.Warning)
                    return
            new_filepath = pathlib.Path(f"{self.workspace_path}/{tmp_project.get_project_name()}.xml")
            tmp_project.serialize_project(new_filepath)
            self.app_project = self.app_project.deserialize_project(new_filepath, self.app_settings)
            constants.PYSSA_LOGGER.info(f"Opening the project {self.app_project.get_project_name()}.")
            self._project_watcher.current_project = self.app_project
            self.project_scanner.project = self.app_project
            constants.PYSSA_LOGGER.info(
                f"{self._project_watcher.current_project.get_project_name()} is the current project.")
            self.ui.lbl_current_project_name.setText(self.app_project.get_project_name())
            self._project_watcher.on_home_page = False
            self._project_watcher.show_valid_options(self.ui)
            self.display_view_page()

    def export_current_project(self):
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button, self.ui.btn_export_project)
        file_dialog = QFileDialog()
        desktop_path = PyQt5.QtCore.QStandardPaths.standardLocations(PyQt5.QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(self, "Save current project", "", "XML Files (*.xml)")
        if file_path:
            self.app_project.serialize_project(file_path)

    # </editor-fold>

    # <editor-fold desc="Close project functions">
    def close_project(self):
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        cmd.reinitialize()
        self._project_watcher.on_home_page = True
        self._project_watcher.current_project = project.Project("", pathlib.Path(""))
        self._project_watcher.show_valid_options(self.ui)
        self.ui.lbl_current_project_name.setText("")
        self._init_all_pages()
        self.results_name = ""
        constants.PYSSA_LOGGER.info(f"The project {self.app_project.get_project_name()} was closed")
        self.display_home_page()

    # </editor-fold>

    # <editor-fold desc="Monomer Local Prediction functions">
    def local_pred_mono_validate_protein_name(self):
        """This function validates the input of the project name in real-time

        """
        # TODO: does not work as expected
        if not safeguard.Safeguard.check_if_value_is_in_table_v_header(self.ui.txt_pred_mono_prot_name.text(),
                                                                       self.ui.table_pred_mono_prot_to_predict):
            self.ui.lbl_pred_mono_prot_name_status.setText("Protein name already used.")
        tools.validate_protein_name(self.ui.txt_pred_mono_prot_name,
                                    self.ui.lbl_pred_mono_prot_name_status,
                                    self.ui.btn_pred_mono_next)

    def show_prediction_configuration(self):
        config = dialog_advanced_prediction_configurations.DialogAdvancedPredictionConfigurations(self.prediction_configuration)
        config.exec_()
        self.prediction_configuration.amber_force_field = config.prediction_config.amber_force_field
        self.prediction_configuration.templates = config.prediction_config.templates

    def setup_defaults_monomer_prediction(self):
        # clears everything
        self.ui.txt_pred_mono_prot_name.clear()
        self.ui.txt_pred_mono_seq_name.clear()
        # sets up defaults: Prediction
        self.ui.btn_pred_mono_next.setEnabled(False)
        self.ui.btn_pred_mono_add_protein.setEnabled(False)
        self.ui.lbl_pred_mono_prot_name_status.setText("")
        self.ui.lbl_pred_mono_seq_name_status.setText("")

    def local_pred_mono_validate_protein_sequence(self):
        """This function validates the input of the protein sequence in real-time

        """
        tools.validate_protein_sequence(self.ui.txt_pred_mono_seq_name,
                                        self.ui.lbl_pred_mono_seq_name_status,
                                        self.ui.btn_pred_mono_add_protein)

    def local_pred_mono_show_protein_overview(self):
        if self.ui.table_pred_mono_prot_to_predict.rowCount() == 0:
            self.local_pred_monomer_management.show_stage_x(0)
        else:
            gui_elements_to_show = [
                self.ui.btn_pred_mono_seq_to_predict,
                self.ui.btn_pred_mono_seq_to_predict_remove,
            ]
            self.local_pred_monomer_management.show_gui_elements_stage_x(
                [0, 3], [1, 2], show_specific_elements=gui_elements_to_show
            )

    def local_pred_mono_show_protein_name(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_mono_seq_to_predict,
            self.ui.btn_pred_mono_seq_to_predict_remove,
        ]
        self.local_pred_monomer_management.show_gui_elements_stage_x(
            [0, 1], [2, 3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_mono_prot_name,
                                  self.ui.lbl_pred_mono_prot_name)

    def local_pred_mono_show_protein_sequence(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_mono_seq_to_predict,
            self.ui.btn_pred_mono_seq_to_predict_remove,
            self.ui.btn_pred_mono_back,
            self.ui.btn_pred_mono_next,
        ]
        self.local_pred_monomer_management.show_gui_elements_stage_x(
            [0, 1, 2], [3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.disable_text_box(self.ui.txt_pred_mono_prot_name,
                                   self.ui.lbl_pred_mono_prot_name)

    def local_pred_mono_check_if_table_is_empty(self):
        if self.ui.table_pred_mono_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_mono_predict)
            self.ui.btn_pred_mono_predict.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_mono_predict)
            self.ui.btn_pred_mono_predict.setEnabled(True)

    def local_pred_mono_add_protein_to_predict(self):
        self.ui.table_pred_mono_prot_to_predict.setRowCount(self.ui.table_pred_mono_prot_to_predict.rowCount()+1)
        self.ui.table_pred_mono_prot_to_predict.insertRow(self.ui.table_pred_mono_prot_to_predict.rowCount()+1)
        self.ui.table_pred_mono_prot_to_predict.setItem(self.ui.table_pred_mono_prot_to_predict.rowCount() - 1, 0,
                                                        QTableWidgetItem("A"))
        self.ui.table_pred_mono_prot_to_predict.setItem(self.ui.table_pred_mono_prot_to_predict.rowCount()-1, 1,
                                                        QTableWidgetItem(self.ui.txt_pred_mono_seq_name.toPlainText()))
        self.ui.table_pred_mono_prot_to_predict.setVerticalHeaderItem(self.ui.table_pred_mono_prot_to_predict.rowCount()-1,
                                                                      QTableWidgetItem(self.ui.txt_pred_mono_prot_name.text()))
        self.ui.table_pred_mono_prot_to_predict.resizeColumnsToContents()
        self.local_pred_mono_check_if_table_is_empty()
        self.local_pred_mono_show_protein_overview()
        self.setup_defaults_monomer_prediction()

    def local_pred_mono_remove_protein_to_predict(self):
        self.ui.table_pred_mono_prot_to_predict.removeRow(self.ui.table_pred_mono_prot_to_predict.currentRow())
        self.local_pred_mono_check_if_table_is_empty()
        self.local_pred_mono_show_protein_overview()

    # def local_pred_mono_validate_protein_name(self):
    #     """This function validates the input of the project name in real-time
    #
    #     """
    #     tools.validate_protein_name(self.ui.txt_local_pred_mono_protein_name,
    #                                 self.ui.lbl_local_pred_mono_status_protein_name,
    #                                 self.ui.btn_local_pred_mono_next)
    #
    # def local_pred_mono_show_prediction_configuration(self):
    #     config = dialog_advanced_prediction_configurations.DialogAdvancedPredictionConfigurations(self.prediction_configuration)
    #     config.exec_()
    #     self.prediction_configuration.amber_force_field = config.prediction_config.amber_force_field
    #     self.prediction_configuration.templates = config.prediction_config.templates
    #
    # def local_pred_mono_validate_protein_sequence(self):
    #     """This function validates the input of the protein sequence in real-time
    #
    #     """
    #     tools.validate_protein_sequence(self.ui.txt_local_pred_mono_prot_seq,
    #                                     self.ui.lbl_local_pred_mono_status_prot_seq,
    #                                     self.ui.btn_local_pred_mono_next_2)
    #
    # def local_pred_mono_show_protein_sequence(self):
    #     gui_elements_hide = [
    #         self.ui.btn_local_pred_mono_next,
    #     ]
    #     gui_elements_show = [
    #         self.ui.lbl_local_pred_mono_prot_seq,
    #         self.ui.lbl_local_pred_mono_status_prot_seq,
    #         self.ui.txt_local_pred_mono_prot_seq,
    #         self.ui.btn_local_pred_mono_back,
    #         self.ui.btn_local_pred_mono_next_2,
    #     ]
    #     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
    #     gui_utils.disable_text_box(self.ui.txt_local_pred_mono_protein_name,
    #                                self.ui.lbl_local_pred_mono_protein_name)
    #
    # def local_pred_mono_hide_protein_sequence(self):
    #     gui_elements_hide = [
    #         self.ui.lbl_local_pred_mono_prot_seq,
    #         self.ui.lbl_local_pred_mono_status_prot_seq,
    #         self.ui.txt_local_pred_mono_prot_seq,
    #         self.ui.btn_local_pred_mono_back,
    #         self.ui.btn_local_pred_mono_next_2,
    #     ]
    #     gui_elements_show = [
    #         self.ui.btn_local_pred_mono_next,
    #     ]
    #     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
    #     gui_utils.enable_text_box(self.ui.txt_local_pred_mono_protein_name,
    #                               self.ui.lbl_local_pred_mono_protein_name)
    #
    # def local_pred_mono_show_advanced_config(self):
    #     gui_elements_show = [
    #         self.ui.btn_local_pred_mono_back_2,
    #         self.ui.btn_local_pred_mono_predict,
    #         self.ui.lbl_local_pred_mono_advanced_config,
    #         self.ui.btn_local_pred_mono_advanced_config,
    #     ]
    #     gui_elements_hide = [
    #         self.ui.btn_local_pred_mono_back,
    #         self.ui.btn_local_pred_mono_next_2,
    #     ]
    #     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
    #     gui_utils.disable_text_box(self.ui.txt_local_pred_mono_prot_seq,
    #                                self.ui.lbl_local_pred_mono_prot_seq)
    #     self.ui.btn_local_pred_mono_predict.setEnabled(True)
    #     styles.color_button_ready(self.ui.btn_local_pred_mono_predict)
    #
    # def local_pred_mono_hide_advanced_config(self):
    #     gui_elements_hide = [
    #         self.ui.btn_local_pred_mono_back_2,
    #         self.ui.btn_local_pred_mono_predict,
    #         self.ui.lbl_local_pred_mono_advanced_config,
    #         self.ui.btn_local_pred_mono_advanced_config,
    #     ]
    #     gui_elements_show = [
    #         self.ui.btn_local_pred_mono_back,
    #         self.ui.btn_local_pred_mono_next_2,
    #     ]
    #     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
    #     gui_utils.enable_text_box(self.ui.txt_local_pred_mono_prot_seq,
    #                               self.ui.lbl_local_pred_mono_prot_seq)
    #     self.ui.btn_local_pred_mono_predict.setEnabled(False)
    #     styles.color_button_not_ready(self.ui.btn_local_pred_mono_predict)

    def post_prediction_process(self):
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
        basic_boxes.ok("Structure prediction", "All structure predictions are done. Go to View to check the new proteins.",
                       QMessageBox.Information)
        constants.PYSSA_LOGGER.info("All structure predictions are done.")
        self._project_watcher.show_valid_options(self.ui)
        self._init_local_pred_mono_page()

    def predict_local_monomer(self):
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
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        worker = workers.PredictionWorkerPool(self.ui.table_pred_mono_prot_to_predict,
                                              self.prediction_configuration, self.app_project)
        worker.signals.finished.connect(self.post_prediction_process)
        constants.PYSSA_LOGGER.info("Thread started for prediction process.")
        self.threadpool.start(worker)
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
        # self._thread_controller.thread_worker_pairs.get(constants.PREDICTION_TASK).setup_and_run_thread(display_msg_box2)
        # gui_elements_to_show = [
        #     self.ui.btn_prediction_abort,
        # ]
        # gui_elements_to_hide = [
        #     self.ui.btn_use_page,
        #     self.ui.btn_close_project,
        # ]
        # gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        # #self._project_watcher.show_valid_options(self.ui)
        self.display_view_page()

    def abort_prediction(self):
        constants.PYSSA_LOGGER.info("Structure prediction process was aborted manually.")
        subprocess.run(["wsl", "--shutdown"])
        constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
        filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
        constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
        basic_boxes.ok("Abort prediction", "The structure prediction was aborted.", QMessageBox.Information)
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_prediction_abort)
        self._project_watcher.show_valid_options(self.ui)

    # </editor-fold>

    # <editor-fold desc="Multimer Local Prediction functions">
    def local_pred_multi_validate_protein_name(self):
        """This function validates the input of the project name in real-time

        """
        tools.validate_protein_name(self.ui.txt_pred_multi_prot_name,
                                    self.ui.lbl_pred_multi_prot_name_status,
                                    self.ui.btn_pred_multi_next)

    def local_pred_multi_validate_protein_sequence(self):
        """This function validates the input of the protein sequence in real-time

        """
        tools.validate_protein_sequence(self.ui.txt_pred_multi_prot_seq,
                                        self.ui.lbl_pred_multi_prot_seq_status,
                                        self.ui.btn_pred_multi_prot_seq_add)

    def local_pred_multi_show_protein_overview(self):
        if self.ui.table_pred_multi_prot_to_predict.rowCount() == 0:
            self.local_pred_multimer_management.show_stage_x(0)
        else:
            gui_elements_to_show = [
                self.ui.btn_pred_multi_prot_to_predict_add,
                self.ui.btn_pred_multi_prot_to_predict_remove,
            ]
            self.local_pred_multimer_management.show_gui_elements_stage_x(
                [0, 3], [1, 2], show_specific_elements=gui_elements_to_show
            )

    def local_pred_multi_show_protein_name(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.btn_pred_multi_prot_to_predict_remove,
        ]
        self.local_pred_multimer_management.show_gui_elements_stage_x(
            [0, 1], [2, 3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_multi_prot_name,
                                  self.ui.lbl_pred_multi_prot_name)

    def local_pred_multi_show_protein_sequence(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.btn_pred_multi_back,
            self.ui.btn_pred_multi_next,
        ]
        self.local_pred_multimer_management.show_gui_elements_stage_x(
            [0, 1, 2], [3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.disable_text_box(self.ui.txt_pred_multi_prot_name,
                                   self.ui.lbl_pred_multi_prot_name)

    def local_pred_multi_add_sequence_to_list(self):
        self.ui.list_pred_multi_prot_seq_overview.addItem(QListWidgetItem(self.ui.txt_pred_multi_prot_seq.toPlainText()))
        self.local_pred_multi_check_if_list_is_empty()

    def local_pred_multi_remove_sequence_to_list(self):
        self.ui.list_pred_multi_prot_seq_overview.takeItem(self.ui.list_pred_multi_prot_seq_overview.currentRow())
        self.local_pred_multi_check_if_list_is_empty()

    def local_pred_multi_check_if_list_is_empty(self):
        if self.ui.list_pred_multi_prot_seq_overview.count() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_multi_prot_to_predict_add_2)
            self.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_multi_prot_to_predict_add_2)
            self.ui.btn_pred_multi_prot_to_predict_add_2.setEnabled(True)

    def local_pred_multi_check_if_table_is_empty(self):
        if self.ui.table_pred_multi_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_multi_predict)
            self.ui.btn_pred_multi_predict.setEnabled(False)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_multi_predict)
            self.ui.btn_pred_multi_predict.setEnabled(True)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(True)

    def local_pred_multi_add_protein_to_predict(self):
        for i in range(self.ui.list_pred_multi_prot_seq_overview.count()):
            self.ui.table_pred_multi_prot_to_predict.setRowCount(
                self.ui.table_pred_multi_prot_to_predict.rowCount() + 1)
            self.ui.table_pred_multi_prot_to_predict.insertRow(self.ui.table_pred_multi_prot_to_predict.rowCount() + 1)
            tmp_chain_seq = (constants.chain_dict.get(i), self.ui.list_pred_multi_prot_seq_overview.item(i).text())
            self.ui.table_pred_multi_prot_to_predict.setItem(self.ui.table_pred_multi_prot_to_predict.rowCount() - 1, 0,
                                                             QTableWidgetItem(tmp_chain_seq[0]))
            self.ui.table_pred_multi_prot_to_predict.setItem(self.ui.table_pred_multi_prot_to_predict.rowCount() - 1, 1,
                                                             QTableWidgetItem(tmp_chain_seq[1]))
            name_item = QTableWidgetItem(self.ui.txt_pred_multi_prot_name.text())
            self.ui.table_pred_multi_prot_to_predict.setVerticalHeaderItem(self.ui.table_pred_multi_prot_to_predict.rowCount()-1, name_item)
        self.ui.table_pred_multi_prot_to_predict.resizeColumnsToContents()
        self.local_pred_multi_check_if_table_is_empty()
        self.local_pred_multi_show_protein_overview()
        self._init_local_pred_multi_page()

    def local_pred_multi_remove_protein_to_predict(self):
        self.ui.table_pred_multi_prot_to_predict.removeRow(self.ui.table_pred_multi_prot_to_predict.currentRow())
        prot_name = self.ui.table_pred_multi_prot_to_predict.verticalHeaderItem(self.ui.table_pred_multi_prot_to_predict.currentRow()).text()
        for i in range(self.ui.table_pred_multi_prot_to_predict.rowCount()):
            if self.ui.table_pred_multi_prot_to_predict.verticalHeaderItem(i).text() == prot_name:
                self.ui.table_pred_multi_prot_to_predict.setItem(i, 0, QTableWidgetItem(constants.chain_dict.get(i)))
        self.local_pred_multi_check_if_table_is_empty()
        self.local_pred_multi_show_protein_overview()

    def predict_local_multimer(self):
        constants.PYSSA_LOGGER.info("Begin multimer prediction process.")
        worker = workers.PredictionWorkerPool(self.ui.table_pred_multi_prot_to_predict,
                                              self.prediction_configuration, self.app_project)
        worker.signals.finished.connect(self.post_prediction_process)
        constants.PYSSA_LOGGER.info("Thread started for prediction process.")
        self.threadpool.start(worker)
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
        self.display_view_page()

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
    # <editor-fold desc="Prediction section">
    def mono_pred_analysis_validate_protein_name(self):
        """This function validates the input of the project name in real-time

        """
        # TODO: does not work as expected
        if not safeguard.Safeguard.check_if_value_is_in_table_v_header(self.ui.txt_pred_analysis_mono_prot_name.text(),
                                                                       self.ui.table_pred_analysis_mono_prot_to_predict):
            self.ui.lbl_pred_analysis_mono_prot_name_status.setText("Protein name already used.")
        tools.validate_protein_name(self.ui.txt_pred_analysis_mono_prot_name,
                                    self.ui.lbl_pred_analysis_mono_prot_name_status,
                                    self.ui.btn_pred_analysis_mono_next)

    def setup_defaults_monomer_prediction_analysis(self):
        # clears everything
        self.ui.txt_pred_analysis_mono_prot_name.clear()
        self.ui.txt_pred_analysis_mono_seq_name.clear()
        # sets up defaults: Prediction
        self.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self.ui.btn_pred_analysis_mono_add_protein.setEnabled(False)
        self.ui.lbl_pred_analysis_mono_prot_name_status.setText("")
        self.ui.lbl_pred_analysis_mono_seq_name_status.setText("")

    def mono_pred_analysis_validate_protein_sequence(self):
        """This function validates the input of the protein sequence in real-time

        """
        tools.validate_protein_sequence(self.ui.txt_pred_analysis_mono_seq_name,
                                        self.ui.lbl_pred_analysis_mono_seq_name_status,
                                        self.ui.btn_pred_analysis_mono_add_protein)

    def mono_pred_analysis_show_protein_overview(self):
        if self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() == 0:
            self.monomer_prediction_analysis_management.show_stage_x(0)
        else:
            gui_elements_to_show = [
                self.ui.btn_pred_analysis_mono_seq_to_predict,
                self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            ]
            self.monomer_prediction_analysis_management.show_gui_elements_stage_x(
                [0, 3], [1, 2], show_specific_elements=gui_elements_to_show
            )

    def mono_pred_analysis_show_protein_name(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_seq_to_predict,
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
        ]
        self.monomer_prediction_analysis_management.show_gui_elements_stage_x(
            [0, 1], [2, 3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_mono_prot_name,
                                  self.ui.lbl_pred_analysis_mono_prot_name)

    def mono_pred_analysis_show_protein_sequence(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_analysis_mono_seq_to_predict,
            self.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self.ui.btn_pred_analysis_mono_back,
            self.ui.btn_pred_analysis_mono_next,
        ]
        self.monomer_prediction_analysis_management.show_gui_elements_stage_x(
            [0, 1, 2], [3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_mono_prot_name,
                                   self.ui.lbl_pred_analysis_mono_prot_name)

    def mono_pred_analysis_check_if_table_is_empty(self):
        if self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_analysis_mono_go_analysis_setup)
            self.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_analysis_mono_go_analysis_setup)
            self.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(True)

    def mono_pred_analysis_add_protein_to_predict(self):
        self.ui.table_pred_analysis_mono_prot_to_predict.setRowCount(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()+1)
        self.ui.table_pred_analysis_mono_prot_to_predict.insertRow(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()+1)
        self.ui.table_pred_analysis_mono_prot_to_predict.setItem(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1, 0,
                                                        QTableWidgetItem("A"))
        self.ui.table_pred_analysis_mono_prot_to_predict.setItem(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()-1, 1,
                                                        QTableWidgetItem(self.ui.txt_pred_analysis_mono_seq_name.toPlainText()))
        self.ui.table_pred_analysis_mono_prot_to_predict.setVerticalHeaderItem(self.ui.table_pred_analysis_mono_prot_to_predict.rowCount()-1,
                                                                      QTableWidgetItem(self.ui.txt_pred_analysis_mono_prot_name.text()))
        self.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        self.mono_pred_analysis_check_if_table_is_empty()
        self.mono_pred_analysis_show_protein_overview()
        self.setup_defaults_monomer_prediction()

    def mono_pred_analysis_remove_protein_to_predict(self):
        self.ui.table_pred_analysis_mono_prot_to_predict.removeRow(self.ui.table_pred_analysis_mono_prot_to_predict.currentRow())
        self.mono_pred_analysis_check_if_table_is_empty()
        self.mono_pred_analysis_show_protein_overview()

    # </editor-fold>

    def switch_monomer_pred_analysis_tab(self):
        if self.ui.tabWidget.currentIndex() == 0:
            self.ui.tabWidget.setCurrentIndex(1)
            self.show_mono_pred_analysis_stage_0()
            self.mono_pred_analysis_show_protein_overview()
        else:
            self.ui.list_pred_analysis_mono_overview.clear()
            self.ui.tabWidget.setCurrentIndex(0)

    # <editor-fold desc="Analysis section">
    def show_mono_pred_analysis_stage_0(self):
        gui_page_management.show_analysis_page_stage_0(self.monomer_prediction_analysis_management,
                                                       self.ui.list_pred_analysis_mono_ref_chains,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_2,
                                                       self.ui.list_pred_analysis_mono_model_chains,
                                                       self.ui.list_pred_analysis_mono_overview,
                                                       self.ui.btn_pred_analysis_mono_remove,
                                                       self.ui.btn_pred_analysis_mono_add,
                                                       4)

    def show_mono_pred_analysis_stage_1(self):
        gui_page_management.show_analysis_page_stage_1(self.monomer_prediction_analysis_management,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_2,
                                                       self.ui.btn_pred_analysis_mono_remove,
                                                       self.ui.btn_pred_analysis_mono_add,
                                                       4,
                                                       self.fill_mono_pred_analysis_protein_boxes)

    def show_mono_pred_analysis_stage_2(self):
        gui_page_management.show_analysis_page_stage_2(self.app_project,
                                                       self.monomer_prediction_analysis_management,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_2,
                                                       self.ui.box_pred_analysis_mono_prot_struct_1,
                                                       self.ui.box_pred_analysis_mono_prot_struct_2,
                                                       self.ui.lbl_pred_analysis_mono_ref_chains,
                                                       self.ui.list_pred_analysis_mono_ref_chains,
                                                       self.ui.btn_pred_analysis_mono_next_2,
                                                       self.ui.btn_pred_analysis_mono_next_3,
                                                       self.ui.btn_pred_analysis_mono_back_3,
                                                       4,
                                                       table_prot_to_predict=self.ui.table_pred_analysis_mono_prot_to_predict,
                                                       state=constants.PREDICTION_ANALYSIS)

    def show_mono_pred_analysis_stage_3(self):
        gui_page_management.show_analysis_page_stage_3(self.app_project,
                                                       self.monomer_prediction_analysis_management,
                                                       self.ui.list_pred_analysis_mono_overview,
                                                       self.ui.btn_pred_analysis_mono_add,
                                                       self.ui.btn_pred_analysis_mono_remove,
                                                       self.ui.btn_pred_analysis_mono_next_2,
                                                       self.ui.btn_pred_analysis_mono_back_3,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_mono_prot_struct_2,
                                                       self.ui.box_pred_analysis_mono_prot_struct_1,
                                                       self.ui.box_pred_analysis_mono_prot_struct_2,
                                                       self.ui.btn_pred_analysis_mono_next_3,
                                                       self.ui.btn_pred_analysis_mono_back_4,
                                                       self.ui.lbl_pred_analysis_mono_model_chains,
                                                       self.ui.list_pred_analysis_mono_model_chains,
                                                       self.ui.btn_pred_analysis_mono_next_4,
                                                       self.no_of_selected_chains,
                                                       4,
                                                       table_prot_to_predict=self.ui.table_pred_analysis_mono_prot_to_predict,
                                                       state=constants.PREDICTION_ANALYSIS)

    def fill_mono_pred_analysis_protein_boxes(self):
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

    def remove_mono_pred_analysis_analysis_run(self):
        self.ui.list_pred_analysis_mono_overview.takeItem(self.ui.list_pred_analysis_mono_overview.currentRow())
        if self.ui.list_pred_analysis_mono_overview.count() == 0:
            self.monomer_prediction_analysis_management.show_stage_x(4)
            self.ui.btn_pred_analysis_mono_remove.hide()

    def check_mono_pred_analysis_if_same_no_of_chains_selected(self):
        self.ui.btn_pred_analysis_mono_next_4.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_pred_analysis_mono_next_4)
        if self.no_of_selected_chains == len(self.ui.list_pred_analysis_mono_model_chains.selectedItems()):
            styles.color_button_ready(self.ui.btn_pred_analysis_mono_next_4)
            self.ui.btn_pred_analysis_mono_next_4.setEnabled(True)

    def check_mono_pred_analysis_if_prot_structs_are_filled(self):
        prot_1 = self.ui.box_pred_analysis_mono_prot_struct_1.itemText(self.ui.box_pred_analysis_mono_prot_struct_1.currentIndex())
        prot_2 = self.ui.box_pred_analysis_mono_prot_struct_2.itemText(self.ui.box_pred_analysis_mono_prot_struct_2.currentIndex())
        if prot_1 != "" and prot_2 != "":
            self.ui.btn_pred_analysis_mono_next_2.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_mono_next_2.setEnabled(False)

    def count_mono_pred_analysis_selected_chains_for_prot_struct_1(self):
        self.no_of_selected_chains = len(self.ui.list_pred_analysis_mono_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self.ui.btn_pred_analysis_mono_next_3.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_mono_next_3.setEnabled(False)

    # </editor-fold>

    def post_prediction_analysis_process(self):
        # TODO: Wedding
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
        basic_boxes.ok("Structure prediction", "All structure predictions are done. Go to View to check the new proteins.",
                       QMessageBox.Information)
        constants.PYSSA_LOGGER.info("All structure predictions are done.")
        self._project_watcher.show_valid_options(self.ui)

    def start_monomer_prediction_analysis(self):
        # TODO: Wedding
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        worker = workers.PredictionWorkerPool(self.ui.table_pred_analysis_mono_prot_to_predict,
                                              self.prediction_configuration, self.app_project)
        worker.signals.finished.connect(self.post_prediction_process)
        constants.PYSSA_LOGGER.info("Thread started for prediction process.")
        self.threadpool.start(worker)
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
        self.display_view_page()
    # </editor-fold>

    # <editor-fold desc="Multimer Prediction + Analysis functions">
    # <editor-fold desc="Prediction section">
    def multi_pred_analysis_validate_protein_name(self):
        """This function validates the input of the project name in real-time

        """
        tools.validate_protein_name(self.ui.txt_pred_analysis_multi_prot_name,
                                    self.ui.lbl_pred_analysis_multi_prot_name_status,
                                    self.ui.btn_pred_analysis_multi_next)

    def multi_pred_analysis_validate_protein_sequence(self):
        """This function validates the input of the protein sequence in real-time

        """
        tools.validate_protein_sequence(self.ui.txt_pred_analysis_multi_prot_seq,
                                        self.ui.lbl_pred_analysis_multi_prot_seq_status,
                                        self.ui.btn_pred_analysis_multi_prot_seq_add)

    def multi_pred_analysis_show_protein_overview(self):
        if self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 0:
            self.multimer_prediction_analysis_management.show_stage_x(0)
        else:
            gui_elements_to_show = [
                self.ui.btn_pred_multi_prot_to_predict_add,
                self.ui.btn_pred_multi_prot_to_predict_remove,
            ]
            self.multimer_prediction_analysis_management.show_gui_elements_stage_x(
                [0, 3], [1, 2], show_specific_elements=gui_elements_to_show
            )

    def multi_pred_analysis_show_protein_name(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.btn_pred_multi_prot_to_predict_remove,
        ]
        self.multimer_prediction_analysis_management.show_gui_elements_stage_x(
            [0, 1], [2, 3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.enable_text_box(self.ui.txt_pred_analysis_multi_prot_name,
                                  self.ui.lbl_pred_analysis_multi_prot_name)

    def multi_pred_analysis_show_protein_sequence(self):
        gui_elements_to_hide = [
            self.ui.btn_pred_multi_prot_to_predict_add,
            self.ui.btn_pred_multi_prot_to_predict_remove,
            self.ui.btn_pred_analysis_multi_back,
            self.ui.btn_pred_analysis_multi_next,
        ]
        self.multimer_prediction_analysis_management.show_gui_elements_stage_x(
            [0, 1, 2], [3], hide_specific_elements=gui_elements_to_hide)
        gui_utils.disable_text_box(self.ui.txt_pred_analysis_multi_prot_name,
                                   self.ui.lbl_pred_analysis_multi_prot_name)

    def multi_pred_analysis_add_sequence_to_list(self):
        self.ui.list_pred_analysis_multi_prot_seq_overview.addItem(
            QListWidgetItem(self.ui.txt_pred_analysis_multi_prot_seq.toPlainText()))
        self.multi_pred_analysis_check_if_list_is_empty()

    def multi_pred_analysis_remove_sequence_to_list(self):
        self.ui.list_pred_analysis_multi_prot_seq_overview.takeItem(self.ui.list_pred_analysis_multi_prot_seq_overview.currentRow())
        self.multi_pred_analysis_check_if_list_is_empty()

    def multi_pred_analysis_check_if_list_is_empty(self):
        if self.ui.list_pred_analysis_multi_prot_seq_overview.count() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_prot_to_predict_add_2)
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_analysis_multi_prot_to_predict_add_2)
            self.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(True)

    def multi_pred_analysis_check_if_table_is_empty(self):
        if self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self.ui.btn_pred_multi_predict)
            self.ui.btn_pred_multi_predict.setEnabled(False)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(False)
        else:
            styles.color_button_ready(self.ui.btn_pred_multi_predict)
            self.ui.btn_pred_multi_predict.setEnabled(True)
            self.ui.btn_pred_multi_prot_to_predict_remove.setEnabled(True)

    def multi_pred_analysis_add_protein_to_predict(self):
        for i in range(self.ui.list_pred_analysis_multi_prot_seq_overview.count()):
            self.ui.table_pred_analysis_multi_prot_to_predict.setRowCount(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1)
            self.ui.table_pred_analysis_multi_prot_to_predict.insertRow(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1)
            tmp_chain_seq = (constants.chain_dict.get(i), self.ui.list_pred_analysis_multi_prot_seq_overview.item(i).text())
            self.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1, 0,
                QTableWidgetItem(tmp_chain_seq[0]))
            self.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1, 1,
                QTableWidgetItem(tmp_chain_seq[1]))
            name_item = QTableWidgetItem(self.ui.txt_pred_analysis_multi_prot_name.text())
            self.ui.table_pred_analysis_multi_prot_to_predict.setVerticalHeaderItem(
                self.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1, name_item)
        self.ui.table_pred_analysis_multi_prot_to_predict.resizeColumnsToContents()
        self.multi_pred_analysis_check_if_table_is_empty()
        self.multi_pred_analysis_show_protein_overview()
        self._init_multi_pred_analysis_page()

    def multi_pred_analysis_remove_protein_to_predict(self):
        self.ui.table_pred_analysis_multi_prot_to_predict.removeRow(self.ui.table_pred_analysis_multi_prot_to_predict.currentRow())
        prot_name = self.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(
            self.ui.table_pred_analysis_multi_prot_to_predict.currentRow()).text()
        for i in range(self.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            if self.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text() == prot_name:
                self.ui.table_pred_analysis_multi_prot_to_predict.setItem(i, 0,
                                                                 QTableWidgetItem(constants.chain_dict.get(i)))
        self.multi_pred_analysis_check_if_table_is_empty()
        self.multi_pred_analysis_show_protein_overview()

    # </editor-fold>

    def switch_multimer_pred_analysis_tab(self):
        if self.ui.tabWidget_2.currentIndex() == 0:
            self.ui.tabWidget_2.setCurrentIndex(1)
            self.show_multi_pred_analysis_stage_0()
            self.multi_pred_analysis_show_protein_overview()
        else:
            self.ui.list_pred_analysis_multi_overview.clear()
            self.ui.tabWidget_2.setCurrentIndex(0)

    # <editor-fold desc="Analysis section">
    def show_multi_pred_analysis_stage_0(self):
        gui_page_management.show_analysis_page_stage_0(self.multimer_prediction_analysis_management,
                                                       self.ui.list_pred_analysis_multi_ref_chains,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_2,
                                                       self.ui.list_pred_analysis_multi_model_chains,
                                                       self.ui.list_pred_analysis_multi_overview,
                                                       self.ui.btn_pred_analysis_multi_remove,
                                                       self.ui.btn_pred_analysis_multi_add,
                                                       4)

    def show_multi_pred_analysis_stage_1(self):
        gui_page_management.show_analysis_page_stage_1(self.multimer_prediction_analysis_management,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_2,
                                                       self.ui.btn_pred_analysis_multi_remove,
                                                       self.ui.btn_pred_analysis_multi_add,
                                                       4,
                                                       self.fill_multi_pred_analysis_protein_boxes)

    def show_multi_pred_analysis_stage_2(self):
        gui_page_management.show_analysis_page_stage_2(self.app_project,
                                                       self.multimer_prediction_analysis_management,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_2,
                                                       self.ui.box_pred_analysis_multi_prot_struct_1,
                                                       self.ui.box_pred_analysis_multi_prot_struct_2,
                                                       self.ui.lbl_pred_analysis_multi_ref_chains,
                                                       self.ui.list_pred_analysis_multi_ref_chains,
                                                       self.ui.btn_pred_analysis_multi_next_2,
                                                       self.ui.btn_pred_analysis_multi_next_3,
                                                       self.ui.btn_pred_analysis_multi_back_3,
                                                       4,
                                                       table_prot_to_predict=self.ui.table_pred_analysis_multi_prot_to_predict,
                                                       state=constants.PREDICTION_ANALYSIS)

    def show_multi_pred_analysis_stage_3(self):
        gui_page_management.show_analysis_page_stage_3(self.app_project,
                                                       self.multimer_prediction_analysis_management,
                                                       self.ui.list_pred_analysis_multi_overview,
                                                       self.ui.btn_pred_analysis_multi_add,
                                                       self.ui.btn_pred_analysis_multi_remove,
                                                       self.ui.btn_pred_analysis_multi_next_2,
                                                       self.ui.btn_pred_analysis_multi_back_3,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_1,
                                                       self.ui.lbl_pred_analysis_multi_prot_struct_2,
                                                       self.ui.box_pred_analysis_multi_prot_struct_1,
                                                       self.ui.box_pred_analysis_multi_prot_struct_2,
                                                       self.ui.btn_pred_analysis_multi_next_3,
                                                       self.ui.btn_pred_analysis_multi_back_4,
                                                       self.ui.lbl_pred_analysis_multi_model_chains,
                                                       self.ui.list_pred_analysis_multi_model_chains,
                                                       self.ui.btn_pred_analysis_multi_next_4,
                                                       self.no_of_selected_chains,
                                                       4,
                                                       table_prot_to_predict=self.ui.table_pred_analysis_multi_prot_to_predict,
                                                       state=constants.PREDICTION_ANALYSIS)

    def fill_multi_pred_analysis_protein_boxes(self):
        protein_names = []
        prediction_runs = prediction_util.get_prediction_name_and_seq_from_table(self.ui.table_pred_analysis_multi_prot_to_predict)
        for tmp_prediction_run in prediction_runs:
            protein_names.append(tmp_prediction_run.name)
        for tmp_protein in self.app_project.proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        self.ui.box_pred_analysis_multi_prot_struct_1.clear()
        self.ui.box_pred_analysis_multi_prot_struct_2.clear()
        gui_utils.fill_combo_box(self.ui.box_pred_analysis_multi_prot_struct_1, protein_names)
        gui_utils.fill_combo_box(self.ui.box_pred_analysis_multi_prot_struct_2, protein_names)

    def remove_multi_pred_analysis_analysis_run(self):
        self.ui.list_pred_analysis_multi_overview.takeItem(self.ui.list_pred_analysis_multi_overview.currentRow())
        if self.ui.list_pred_analysis_multi_overview.count() == 0:
            self.multimer_prediction_analysis_management.show_stage_x(4)
            self.ui.btn_pred_analysis_multi_remove.hide()

    def check_multi_pred_analysis_if_same_no_of_chains_selected(self):
        self.ui.btn_pred_analysis_multi_next_4.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_pred_analysis_multi_next_4)
        if self.no_of_selected_chains == len(self.ui.list_pred_analysis_multi_model_chains.selectedItems()):
            styles.color_button_ready(self.ui.btn_pred_analysis_multi_next_4)
            self.ui.btn_pred_analysis_multi_next_4.setEnabled(True)

    def check_multi_pred_analysis_if_prot_structs_are_filled(self):
        prot_1 = self.ui.box_pred_analysis_multi_prot_struct_1.itemText(self.ui.box_pred_analysis_multi_prot_struct_1.currentIndex())
        prot_2 = self.ui.box_pred_analysis_multi_prot_struct_2.itemText(self.ui.box_pred_analysis_multi_prot_struct_2.currentIndex())
        if prot_1 != "" and prot_2 != "":
            self.ui.btn_pred_analysis_multi_next_2.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_multi_next_2.setEnabled(False)

    def count_multi_pred_analysis_selected_chains_for_prot_struct_1(self):
        self.no_of_selected_chains = len(self.ui.list_pred_analysis_multi_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self.ui.btn_pred_analysis_multi_next_3.setEnabled(True)
        else:
            self.ui.btn_pred_analysis_multi_next_3.setEnabled(False)

    # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Structure Analysis functions">
    def show_batch_analysis_stage_0(self):
        gui_page_management.show_analysis_page_stage_0(self.batch_analysis_management,
                                                       self.ui.list_analysis_batch_ref_chains,
                                                       self.ui.lbl_analysis_batch_prot_struct_1,
                                                       self.ui.lbl_analysis_batch_prot_struct_2,
                                                       self.ui.list_analysis_batch_model_chains,
                                                       self.ui.list_analysis_batch_overview,
                                                       self.ui.btn_analysis_batch_remove,
                                                       self.ui.btn_analysis_batch_add,
                                                       0)

    def show_batch_analysis_stage_1(self):
        gui_page_management.show_analysis_page_stage_1(self.batch_analysis_management,
                                                       self.ui.lbl_analysis_batch_prot_struct_1,
                                                       self.ui.lbl_analysis_batch_prot_struct_2,
                                                       self.ui.btn_analysis_batch_remove,
                                                       self.ui.btn_analysis_batch_add,
                                                       0,
                                                       self.fill_protein_boxes_batch)

    def show_batch_analysis_stage_2(self):
        gui_page_management.show_analysis_page_stage_2(self.app_project,
                                                       self.batch_analysis_management,
                                                       self.ui.lbl_analysis_batch_prot_struct_1,
                                                       self.ui.lbl_analysis_batch_prot_struct_2,
                                                       self.ui.box_analysis_batch_prot_struct_1,
                                                       self.ui.box_analysis_batch_prot_struct_2,
                                                       self.ui.lbl_analysis_batch_ref_chains,
                                                       self.ui.list_analysis_batch_ref_chains,
                                                       self.ui.btn_analysis_batch_next,
                                                       self.ui.btn_analysis_batch_next_2,
                                                       self.ui.btn_analysis_batch_back,
                                                       0)

    def show_batch_analysis_stage_3(self):
        gui_page_management.show_analysis_page_stage_3(self.app_project,
                                                       self.batch_analysis_management,
                                                       self.ui.list_analysis_batch_overview,
                                                       self.ui.btn_analysis_batch_add,
                                                       self.ui.btn_analysis_batch_remove,
                                                       self.ui.btn_analysis_batch_next,
                                                       self.ui.btn_analysis_batch_back,
                                                       self.ui.lbl_analysis_batch_prot_struct_1,
                                                       self.ui.lbl_analysis_batch_prot_struct_2,
                                                       self.ui.box_analysis_batch_prot_struct_1,
                                                       self.ui.box_analysis_batch_prot_struct_2,
                                                       self.ui.btn_analysis_batch_next_2,
                                                       self.ui.btn_analysis_batch_back_2,
                                                       self.ui.lbl_analysis_batch_model_chains,
                                                       self.ui.list_analysis_batch_model_chains,
                                                       self.ui.btn_analysis_batch_next_3,
                                                       self.no_of_selected_chains,
                                                       0)

    def fill_protein_boxes_batch(self):
        proteins = []
        for tmp_protein in self.app_project.proteins:
            proteins.append(tmp_protein.get_molecule_object())
        proteins.insert(0, "")
        self.ui.box_analysis_batch_prot_struct_1.clear()
        self.ui.box_analysis_batch_prot_struct_2.clear()
        gui_utils.fill_combo_box(self.ui.box_analysis_batch_prot_struct_1, proteins)
        gui_utils.fill_combo_box(self.ui.box_analysis_batch_prot_struct_2, proteins)

    def remove_analysis_run(self):
        self.ui.list_analysis_batch_overview.takeItem(self.ui.list_analysis_batch_overview.currentRow())
        if self.ui.list_analysis_batch_overview.count() == 0:
            self.batch_analysis_management.show_stage_x(0)
            self.ui.btn_analysis_batch_remove.hide()

    def check_if_same_no_of_chains_selected_batch(self):
        self.ui.btn_analysis_batch_next_3.setEnabled(False)
        styles.color_button_not_ready(self.ui.btn_analysis_batch_next_3)
        if self.no_of_selected_chains == len(self.ui.list_analysis_batch_model_chains.selectedItems()):
            styles.color_button_ready(self.ui.btn_analysis_batch_next_3)
            self.ui.btn_analysis_batch_next_3.setEnabled(True)

    def check_if_prot_structs_are_filled_batch(self):
        prot_1 = self.ui.box_analysis_batch_prot_struct_1.itemText(self.ui.box_analysis_batch_prot_struct_1.currentIndex())
        prot_2 = self.ui.box_analysis_batch_prot_struct_2.itemText(self.ui.box_analysis_batch_prot_struct_2.currentIndex())
        if prot_1 != "" and prot_2 != "":
            self.ui.btn_analysis_batch_next.setEnabled(True)
        else:
            self.ui.btn_analysis_batch_next.setEnabled(False)

    def count_batch_selected_chains_for_prot_struct_1(self):
        self.no_of_selected_chains = len(self.ui.list_analysis_batch_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self.ui.btn_analysis_batch_next_2.setEnabled(True)
        else:
            self.ui.btn_analysis_batch_next_2.setEnabled(False)

    def post_analysis_process(self):
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
        self.block_box_analysis.destroy(True)
        basic_boxes.ok("Structure analysis", "All structure analysis' are done. Go to results to check the new results.",
                       QMessageBox.Information)
        constants.PYSSA_LOGGER.info("All structure analysis' are done.")
        self._project_watcher.show_valid_options(self.ui)
        self._init_batch_analysis_page()

    def start_process_batch(self):
        """This function contains the main analysis algorithm for the
        Protein structure comparison.

        """
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
        #         basic_boxes.ok("Single Analysis", f"The structure analysis: {analysis_data[3]} already exists!", QMessageBox.Critical)
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
        self.worker_analysis = workers.AnalysisWorkerPool(
            self.ui.list_analysis_batch_overview, self.ui.cb_analysis_images,
            self.status_bar, self.app_project, self.app_settings, self._init_batch_analysis_page)
        constants.PYSSA_LOGGER.info("Thread started for analysis process.")
        self.threadpool.start(self.worker_analysis)
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

    def abort_analysis(self):
        # TODO: abort analysis does not work!
        self.worker_analysis.__del__()
        basic_boxes.ok("Abort structure analysis", "The structure analysis was aborted.", QMessageBox.Information)
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_analysis_abort)
        self.last_sidebar_button = styles.color_sidebar_buttons(self.last_sidebar_button,
                                                                self.ui.btn_analysis_abort)

    # </editor-fold>

    # <editor-fold desc="Analysis Images">
    def add_protein_pair_to_image_creation_queue(self):
        protein_pair_to_add = self.ui.list_analysis_images_struct_analysis.currentItem().text()
        self.ui.list_analysis_images_creation_struct_analysis.addItem(protein_pair_to_add)
        self.ui.list_analysis_images_struct_analysis.takeItem(self.ui.list_analysis_images_struct_analysis.currentRow())

    def remove_protein_pair_from_image_creation_queue(self):
        protein_pair_to_remove = self.ui.list_analysis_images_creation_struct_analysis.currentItem()
        self.ui.list_analysis_images_creation_struct_analysis.takeItem(self.ui.list_analysis_images_creation_struct_analysis.currentRow())
        self.ui.list_analysis_images_struct_analysis.addItem(protein_pair_to_remove)

    def post_image_creation_process(self):
        self.app_project.serialize_project(self.app_project.get_project_xml_path())
        constants.PYSSA_LOGGER.info("Project has been saved to XML file.")
        self.block_box_images.destroy(True)
        basic_boxes.ok("Analysis Images",
                       "All images of all analysis' have been created. Go to results to check the new results.",
                       QMessageBox.Information)
        constants.PYSSA_LOGGER.info("All images of all analysis' have been created.")
        self._init_analysis_image_page()
        self.display_view_page()
        self._project_watcher.show_valid_options(self.ui)

    def start_automatic_image_creation(self):
        constants.PYSSA_LOGGER.info("Begin image creation process.")
        self.worker_image_creation = workers.BatchImageWorkerPool(
            self.ui.list_analysis_images_struct_analysis, self.ui.list_analysis_images_creation_struct_analysis,
            self.status_bar, self.app_project)
        constants.PYSSA_LOGGER.info("Thread started for image creation process.")
        self.threadpool.start(self.worker_image_creation)
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
    def show_analysis_results_options(self):
        self.results_management.show_stage_x(0)

    def show_results_interactions(self, gui_elements_to_show=None, gui_elements_to_hide=None):
        if gui_elements_to_hide is not None:
            self.results_management.show_gui_elements_stage_x([0, 1], [],
                                                              show_specific_elements=[self.ui.lbl_results_analysis_options,
                                                                                      self.ui.cb_results_analysis_options],
                                                              hide_specific_elements=gui_elements_to_hide)
        else:
            self.results_management.show_gui_elements_stage_x([0, 1], [],
                                                              show_specific_elements=[
                                                                  self.ui.lbl_results_analysis_options,
                                                                  self.ui.cb_results_analysis_options])

    def post_load_results(self):
        print("Results were loaded.")

    def load_results(self):
        shutil.rmtree(constants.CACHE_DIR)
        os.mkdir(constants.CACHE_DIR)
        os.mkdir(constants.CACHE_IMAGES)
        self.results_name = self.ui.cb_results_analysis_options.currentText()
        if self.results_name == "":
            self.show_analysis_results_options()
            return
        # TODO: check if the function below works correctly
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        # session_msg = basic_boxes.yes_or_no("PyMOL Session", "Do you want to save your current pymol session?",
        #                                     QMessageBox.Information)
        # if session_msg is True:
        #     self.current_session.session = pymol_io.convert_pymol_session_to_base64_string(self.current_session.name)
        #     self.app_project.save_pymol_session(self.current_session)

        # <editor-fold desc="Worker variant">
        # worker = workers.ResultsWorkerPool(
        #    self.app_project, self.ui, self.show_analysis_results_options, self.show_results_interactions
        # )
        # worker.signals.finished.connect(self.post_load_results)
        # self.threadpool.start(worker)
        # print("Worker started.")
        #
        # </editor-fold>

        # <editor-fold desc="Main Thread variant">
        gui_elements_to_hide = []
        # TODO: implement image check for xml format
        # if not os.path.exists(pathlib.Path(f"{current_results_path}/images")):
        #     # no images where made
        #     gui_elements_to_hide.append(self.ui.lbl_results_structure_alignment)
        #     gui_elements_to_hide.append(self.ui.btn_view_struct_alignment)
        #     gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
        #     gui_elements_to_hide.append(self.ui.list_results_interest_regions)
        #     gui_elements_to_hide.append(self.ui.btn_view_interesting_region)
        # elif os.path.exists(pathlib.Path(f"{current_results_path}/images")):
        #     if not os.path.exists(pathlib.Path(f"{current_results_path}/images/structure_alignment.png")):
        #         gui_elements_to_hide.append(self.ui.lbl_results_structure_alignment)
        #         gui_elements_to_hide.append(self.ui.btn_view_struct_alignment)
        #     elif not os.path.exists(pathlib.Path(f"{current_results_path}/images/interesting_")):
        #         gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
        #         gui_elements_to_hide.append(self.ui.list_results_interest_regions)
        #         gui_elements_to_hide.append(self.ui.btn_view_interesting_region)

        tmp_protein_pair = self.app_project.search_protein_pair(self.results_name)
        distance_data: dict[str, np.ndarray] = tmp_protein_pair.distance_analysis.analysis_results.distance_data
        distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])

        # check if histogram can be created
        distance_list.sort()
        x, y = np.histogram(distance_list, bins=np.arange(0, distance_list[len(distance_list) - 1], 0.25))
        if x.size != y.size:
            x = np.resize(x, (1, y.size))
        # this conversion is needed for the pyqtgraph library!
        x = x.tolist()
        try:
            x = x[0]
        except IndexError:
            # histogram could not be created
            gui_elements_to_hide.append(self.ui.lbl_results_distance_histogram)
            gui_elements_to_hide.append(self.ui.btn_view_distance_histogram)

        filesystem_io.XmlDeserializer(self.app_project.get_project_xml_path()).deserialize_analysis_images(
            tmp_protein_pair.name, tmp_protein_pair.distance_analysis.analysis_results)
        if len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) != 0 and len(tmp_protein_pair.distance_analysis.analysis_results.interesting_regions_images) != 0:
            # if both image types were made during analysis
            tmp_protein_pair.distance_analysis.analysis_results.create_image_png_files_from_base64()
            self.ui.list_results_interest_regions.clear()
            for tmp_filename in os.listdir(constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                self.ui.list_results_interest_regions.addItem(tmp_filename)
        elif len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) != 0 and len(tmp_protein_pair.distance_analysis.analysis_results.interesting_regions_images) == 0:
            # only struct align image were made
            tmp_protein_pair.distance_analysis.analysis_results.create_image_png_files_from_base64()
            self.ui.lbl_results_structure_alignment.show()
            self.ui.btn_view_struct_alignment.show()
            gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
            gui_elements_to_hide.append(self.ui.list_results_interest_regions)
            gui_elements_to_hide.append(self.ui.btn_view_interesting_region)
        else:
            # no images were made
            gui_elements_to_hide.append(self.ui.lbl_results_structure_alignment)
            gui_elements_to_hide.append(self.ui.btn_view_struct_alignment)
            gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
            gui_elements_to_hide.append(self.ui.list_results_interest_regions)
            gui_elements_to_hide.append(self.ui.btn_view_interesting_region)
        if gui_elements_to_hide:
            self.show_results_interactions(gui_elements_to_hide=gui_elements_to_hide)
        else:
            self.show_results_interactions()
        self.ui.list_results_interest_regions.sortItems()
        self.ui.txt_results_rmsd.setText(str(tmp_protein_pair.distance_analysis.analysis_results.rmsd))
        self.ui.txt_results_aligned_residues.setText(str(tmp_protein_pair.distance_analysis.analysis_results.aligned_aa))
        cmd.reinitialize()
        tmp_protein_pair.load_pymol_session()
        self.current_session = current_session.CurrentSession("protein_pair", tmp_protein_pair.name, tmp_protein_pair.pymol_session)

        # </editor-fold>

    def color_protein_pair_by_rmsd(self):
        """This function colors the residues of the reference and the model protein in 5 colors
                depending on their distance to the reference

        Args:
            alignment_filename:
                filename of the alignment_file
        Returns:

        """
        tmp_protein_pair: 'protein_pair.ProteinPair' = self.app_project.search_protein_pair(self.results_name)

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
                atom_info = protein_pair_util.get_chain_and_position(tmp_protein_pair.distance_analysis.analysis_results.distance_data, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{tmp_protein_pair.protein_1.get_molecule_object()}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{tmp_protein_pair.protein_2.get_molecule_object()}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_1, atom1)
                cmd.color(color_1, atom2)
                # cmd.do(f"color {color_1}, {atom1}")
                # cmd.do(f"color {color_1}, {atom2}")
                i += 1

            elif distance_value <= cutoff_2:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{tmp_protein_pair.protein_1.get_molecule_object()}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{tmp_protein_pair.protein_2.get_molecule_object()}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_2, atom1)
                cmd.color(color_2, atom2)
                i += 1

            elif distance_value <= cutoff_3:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{tmp_protein_pair.protein_1.get_molecule_object()}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}/CA"
                atom2: str = f"/{tmp_protein_pair.protein_2.get_molecule_object()}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}/CA"
                # coloring
                cmd.color(color_3, atom1)
                cmd.color(color_3, atom2)
                i += 1

            elif distance_value <= cutoff_4:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{tmp_protein_pair.protein_1.get_molecule_object()}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{tmp_protein_pair.protein_2.get_molecule_object()}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_4, atom1)
                cmd.color(color_4, atom2)
                i += 1

            elif distance_value <= cutoff_5:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{tmp_protein_pair.protein_1.get_molecule_object()}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{tmp_protein_pair.protein_2.get_molecule_object()}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_5, atom1)
                cmd.color(color_5, atom2)
                i += 1

            elif distance_value > cutoff_5:
                atom_info = protein_pair_util.get_chain_and_position(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data, i)
                # create two atoms for the get_distance command
                atom1: str = f"/{tmp_protein_pair.protein_1.get_molecule_object()}//" \
                             f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                atom2: str = f"/{tmp_protein_pair.protein_2.get_molecule_object()}//" \
                             f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                # coloring
                cmd.color(color_6, f"({atom1})")
                cmd.color(color_6, f"({atom2})")
                i += 1

        # hide unnecessary representations
        cmd.hide("cartoon", tmp_protein_pair.protein_1.get_molecule_object())
        cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}_CA")
        cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}_CA")

    def change_interesting_regions(self):
        """This function is used to switch between projects within a job.

        """
        global global_var_list_widget_row
        global global_var_project_dict

        if gui_utils.warning_switch_pymol_session("") is True:
            try:
                file_path = global_var_project_dict[global_var_list_widget_row].get_protein_pairs_path()
                cmd.save(f"{file_path}/sessions/session_file_model_s.pse")
                tools.quick_log_and_display("info", "Saving the pymol session was successful.",
                                            self.status_bar, "Saving was successful.")
            except KeyError:
                tools.quick_log_and_display("error", "No project has been opened.", self.status_bar,
                                            "No project has been opened.")
            except pymol.CmdException:
                tools.quick_log_and_display("error", "Unexpected Error from PyMOL while saving the "
                                                     "current pymol session", self.status_bar,
                                            "Unexpected Error from PyMOL")

        current_row = self.ui.project_list.currentRow()
        global_var_list_widget_row = current_row
        results_path = global_var_project_dict[current_row].get_protein_pairs_path()
        dir_content = os.listdir(f"{results_path}/images/interesting_regions")
        self.ui.cb_interesting_regions.clear()
        for tmp_file in dir_content:
            self.ui.cb_interesting_regions.addItem(tmp_file)
        cmd.load(global_var_project_dict[current_row].get_session_file())

    def display_structure_alignment(self):
        """This function opens a window which displays the image of the structure alignment.

        """
        png_dialog = QtWidgets.QDialog(self)
        label = QtWidgets.QLabel(self)
        file_path = pathlib.Path(
            f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/{self.results_name}")
        self.ui.cb_results_analysis_options.currentText()
        pixmap = QtGui.QPixmap(f"{constants.CACHE_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{self.ui.cb_results_analysis_options.currentText()}")
        # TO-DO: Create setting for min. image size
        pixmap = pixmap.scaled(450, 450, transformMode=PyQt5.QtCore.Qt.SmoothTransformation)
        label.setPixmap(pixmap)
        label.setScaledContents(True)
        png_dialog_layout = QHBoxLayout()
        png_dialog_layout.addWidget(label)
        png_dialog.setLayout(png_dialog_layout)
        png_dialog.setWindowTitle("Image of: structure alignment")
        png_dialog.show()

    def display_distance_plot(self):
        """This function opens a window which displays the distance plot.

        """
        protein_pair_of_analysis = self.app_project.search_protein_pair(self.ui.cb_results_analysis_options.currentText())
        dialog = dialog_distance_plot.DialogDistancePlot(protein_pair_of_analysis)
        dialog.exec_()

    def display_distance_histogram(self):
        """This function opens a window which displays the distance histogram.

        """
        # item = self.ui.project_list.selectedItems()
        # if item is None:
        #     raise ValueError
        protein_pair_of_analysis = self.app_project.search_protein_pair(
            self.ui.cb_results_analysis_options.currentText())
        dialog = dialog_distance_histogram.DialogDistanceHistogram(protein_pair_of_analysis)
        dialog.exec_()

    def display_interesting_region(self):
        """This function displays an image of an interesting region.

        """
        png_dialog = QtWidgets.QDialog(self)
        label = QtWidgets.QLabel(self)
        global global_var_project_dict
        file_path = constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR
        #file_name = self.ui.cb_interesting_regions.currentText()
        file_name = self.ui.list_results_interest_regions.currentItem().text()
        pixmap = QtGui.QPixmap(f"{constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{file_name}")
        # TO-DO: Create setting for min. image size
        pixmap = pixmap.scaled(450, 450, transformMode=PyQt5.QtCore.Qt.SmoothTransformation)
        label.setPixmap(pixmap)
        label.setScaledContents(True)
        png_dialog_layout = QHBoxLayout()
        png_dialog_layout.addWidget(label)
        png_dialog.setLayout(png_dialog_layout)
        png_dialog.setWindowTitle(f"Image of: {file_name}")
        png_dialog.show()

    def display_distance_table(self):
        """This function displays the distances in a table.

        """
        csv_model = QtGui.QStandardItemModel()
        csv_model.setColumnCount(7)
        labels = ["Residue pair no.", "Reference Chain", "Reference Position", "Reference Residue",
                  "Model Chain", "Model Position", "Model Residue", "Distance",
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
            [distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX], distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN],
             distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION],
             distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI],
             distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN],
             distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION],
             distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI], distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]])
        distance_data_array_transpose = distance_data_array.transpose()
        with open(csv_filepath, mode='w', newline='') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerows(distance_data_array_transpose)

        file.close()

        with open(csv_filepath, 'r', encoding="utf-8") as csv_file:
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
        table_dialog_layout = QHBoxLayout()
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
        table_dialog.setWindowTitle("Distances of Structure Alignment")
        table_dialog.show()

    # </editor-fold>

    # <editor-fold desc="Manage page functions">
    def choose_manage_color_selected_protein(self):
        input = self.ui.box_manage_choose_protein.currentText()
        tmp_protein = self.app_project.search_protein(input)
        tmp_protein.color_protein_in_pymol(self.ui.box_manage_choose_color.currentText(), f"/{tmp_protein.get_molecule_object()}")

    def choose_manage_representation(self):
        """This function sets the representation.

        """
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

    def choose_manage_bg_color(self):
        """This function sets the background color

        """
        if self.ui.box_manage_choose_bg_color.currentIndex() == 0:
            print("Please select a background color.")
            self.status_bar.showMessage("Please select a background color.")
        elif self.ui.box_manage_choose_bg_color.currentIndex() == 1:
            cmd.bg_color("black")
        elif self.ui.box_manage_choose_bg_color.currentIndex() == 2:
            cmd.bg_color("white")
        else:
            print("Missing implementation!")

    # </editor-fold>

    # <editor-fold desc="Image page functions">
    def show_representation(self):
        """This function sets the representation.

        """
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

    def choose_bg_color(self):
        """This function sets the background color

        """
        if self.ui.box_bg_color.currentIndex() == 0:
            print("Please select a background color.")
            self.status_bar.showMessage("Please select a background color.")
        elif self.ui.box_bg_color.currentIndex() == 1:
            cmd.bg_color("black")
        elif self.ui.box_bg_color.currentIndex() == 2:
            cmd.bg_color("white")
        else:
            print("Missing implementation!")

    def choose_renderer(self):
        """This function sets the renderer.

        """
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

    def choose_ray_trace_mode(self):
        """This function sets the ray-trace mode.

        """
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

    def choose_ray_texture(self):
        """This function sets the ray texture.

        """
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

    def decide_transparent_bg(self):
        """This function sets the transparent background.

        """
        if self.ui.cb_transparent_bg.isChecked():
            cmd.set("ray_opaque_background", "off")
        else:
            cmd.set("ray_opaque_background", "on")

    @staticmethod
    def update_scene():
        """This function updates the current selected PyMOL scene.

        """
        cmd.scene(key="auto_generated", action="update")

    def save_scene(self):
        """This function saves the current view as a new PyMOL scene.

        """
        # returns tuple with (name, bool)
        scene_name = QtWidgets.QInputDialog.getText(self, "Save Scene", "Enter scene name:")
        if scene_name[1]:
            cmd.scene(key=scene_name[0], action="append")

    def preview_image(self):
        """This function previews the image

        """
        if self.ui.cb_ray_tracing.isChecked():
            self.status_bar.showMessage("Preview ray-traced image ...")
            cmd.ray(2400, 2400, renderer=int(self.renderer))
            self.status_bar.showMessage("Finished preview of ray-traced image.")
        else:
            self.status_bar.showMessage("Preview draw image ...")
            cmd.draw(2400, 2400)
            self.status_bar.showMessage("Finished preview of drawn image.")

    def save_image(self):
        """This function saves the image as a png file.

        """
        if self.ui.cb_ray_tracing.isChecked():
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image",
                                                             filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display("info", "No file has been selected.",
                                                self.status_bar, "No file has been selected.")
                    return
                self.status_bar.showMessage("Creating ray-traced image ...")
                cmd.ray(2400, 2400, renderer=int(self.renderer))
                cmd.png(full_file_name[0], dpi=300)
                self.status_bar.showMessage("Finished image creation.")
            except FileExistsError:
                tools.quick_log_and_display("error", "File exists already.",
                                            self.status_bar, "File exists already.")
            except pymol.CmdException:
                tools.quick_log_and_display("error", "Unexpected Error from PyMOL while saving the "
                                                     "an image", self.status_bar,
                                            "Unexpected Error from PyMOL")
        else:
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image",
                                                             filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display("info", "No file has been selected.",
                                                self.status_bar, "No file has been selected.")
                    return
                self.status_bar.showMessage("Creating draw image ...")
                cmd.draw(2400, 2400)
                cmd.png(full_file_name[0], dpi=300)
                self.status_bar.showMessage("Finished image creation.")
            except FileExistsError:
                tools.quick_log_and_display("error", "File exists already.",
                                            self.status_bar, "File exists already.")
            except pymol.CmdException:
                tools.quick_log_and_display("error", "Unexpected Error from PyMOL while saving the "
                                                     "an image", self.status_bar,
                                            "Unexpected Error from PyMOL")

    # </editor-fold>

    # <editor-fold desc="Hotspots page functions">
    def open_protein(self):
        tools.ask_to_save_pymol_session(self.app_project, self.current_session)
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find("_vs_") == -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.load_protein_pymol_session()
            tmp_protein.pymol_selection.set_selections_without_chains_ca()
            tmp_model = cmd.get_model(tmp_protein.pymol_selection.selection_string)
            first_amoino_acid_no = tmp_model.atom[0].resi
            self.ui.sp_hotspots_resi_no.setMinimum(int(first_amoino_acid_no))
            tmp_sequence = cmd.get_fastastr('all')
            self.ui.sp_hotspots_resi_no.setMaximum(len(tmp_sequence))
            gui_elements_to_show = [
                self.ui.lbl_hotspots_resi_no,
                self.ui.sp_hotspots_resi_no,
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
            tmp_protein_pair = self.app_project.get_specific_protein_pair(input)
            prot1_no = int(tools.check_first_seq_no(tmp_protein_pair.protein_1))
            prot2_no = int(tools.check_first_seq_no(tmp_protein_pair.protein_2))
            if prot1_no > prot2_no:
                self.ui.sp_hotspots_resi_no.setMinimum(prot2_no)
                self.used_proteins[0] = tmp_protein_pair.protein_1.get_molecule_object()
                self.used_proteins[1] = tmp_protein_pair.protein_2.get_molecule_object()
            else:
                self.ui.sp_hotspots_resi_no.setMinimum(prot1_no)
                self.used_proteins[0] = tmp_protein_pair.protein_1.get_molecule_object()
                self.used_proteins[1] = tmp_protein_pair.protein_2.get_molecule_object()

            prot1_seq_len = int(tools.check_seq_length(tmp_protein_pair.protein_1))
            prot2_seq_len = int(tools.check_seq_length(tmp_protein_pair.protein_2))
            if prot1_seq_len > prot2_seq_len:
                self.ui.sp_hotspots_resi_no.setMaximum(prot2_seq_len)
                self.used_proteins[3] = tmp_protein_pair.protein_2.get_molecule_object()
                self.used_proteins[2] = tmp_protein_pair.protein_1.get_molecule_object()
            else:
                self.ui.sp_hotspots_resi_no.setMaximum(prot1_seq_len)
                self.used_proteins[2] = tmp_protein_pair.protein_1.get_molecule_object()
                self.used_proteins[3] = tmp_protein_pair.protein_2.get_molecule_object()

            gui_elements_to_show = [
                self.ui.lbl_hotspots_resi_no,
                self.ui.sp_hotspots_resi_no,
                self.ui.lbl_hotspots_resi_show,
                self.ui.btn_hotspots_resi_show,
                self.ui.lbl_hotspots_resi_hide,
                self.ui.btn_hotspots_resi_hide,
                self.ui.lbl_hotspots_resi_zoom,
                self.ui.btn_hotspots_resi_zoom,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            tmp_protein_pair.load_pymol_session()

    def show_resi_sticks(self):
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find("_vs_") == -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein.get_molecule_object()}///{self.ui.sp_hotspots_resi_no.text()}/")
            tmp_protein.show_resi_as_balls_and_sticks()
        else:
            tmp_protein_pair = self.app_project.search_protein_pair(input)
            resi_no = int(self.ui.sp_hotspots_resi_no.text())
            if resi_no < int(self.used_proteins[0]):
                # use prot 2
                tmp_protein_pair.protein_2.pymol_selection.set_custom_selection(
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}///{self.ui.sp_hotspots_resi_no.text()}/")
                tmp_protein_pair.protein_2.show_resi_as_balls_and_sticks()
            elif resi_no < int(self.used_proteins[1]):
                # use prot 1
                tmp_protein_pair.protein_1.pymol_selection.set_custom_selection(
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}///{self.ui.sp_hotspots_resi_no.text()}/")
                tmp_protein_pair.protein_1.show_resi_as_balls_and_sticks()
            elif resi_no > int(self.used_proteins[2]):
                # use prot 2
                tmp_protein_pair.protein_2.pymol_selection.set_custom_selection(
                    f"/{tmp_protein_pair.protein_2.get_molecule_object()}///{self.ui.sp_hotspots_resi_no.text()}/")
                tmp_protein_pair.protein_2.show_resi_as_balls_and_sticks()
            elif resi_no > int(self.used_proteins[3]):
                # use prot 1
                tmp_protein_pair.protein_1.pymol_selection.set_custom_selection(
                    f"/{tmp_protein_pair.protein_1.get_molecule_object()}///{self.ui.sp_hotspots_resi_no.text()}/")
                tmp_protein_pair.protein_1.show_resi_as_balls_and_sticks()

    def hide_resi_sticks(self):
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find("_vs_") == -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein.get_molecule_object()}///{self.ui.sp_hotspots_resi_no.text()}/")
            tmp_protein.hide_resi_as_balls_and_sticks()

    def zoom_resi_position(self):
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find("_vs_") == -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein.get_molecule_object()}///{self.ui.sp_hotspots_resi_no.text()}/")
            tmp_protein.zoom_resi_protein_position()

    # </editor-fold>


if __name__ == '__main__':
    app = QApplication(sys.argv)
    styles.set_stylesheet(app)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())
