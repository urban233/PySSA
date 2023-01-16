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
import json
import logging
import os
import shutil
import sys
import webbrowser
import pathlib
import subprocess
import PyQt5.QtCore
import numpy as np
import pymol
import pyqtgraph as pg

from pyssa.gui.ui.dialogs import dialog_sequence_viewer
from pyssa.gui.ui.dialogs import dialog_settings_global
from pyssa.gui.utilities import global_variables
from pyssa.gui.ui.dialogs import dialog_startup
from pyssa.gui.utilities import constants

from PyQt5.QtGui import QIcon
from pymol import Qt
from pymol import cmd
# TODO: fix import statements so that they do not import a class!
from urllib.request import urlopen
from urllib.error import URLError
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5 import QtCore
from pyssa.gui.ui.forms.auto_generated.auto_main_window import Ui_MainWindow
from pyssa.gui.data_structures import project
from pyssa.gui.data_structures import structure_analysis
from pyssa.gui.data_structures.data_classes import protein_analysis_info
from pyssa.gui.data_structures import project_watcher
from pyssa.gui.data_structures import data_transformer
from pyssa.gui.data_structures import safeguard
from pyssa.gui.data_structures.data_classes import prediction_configuration
from pyssa.gui.ui.dialogs import dialog_distance_plot
from pyssa.gui.ui.dialogs import dialog_distance_histogram
from pyssa.gui.ui.dialogs import dialog_about
from pyssa.gui.ui.dialogs import dialog_add_models
from pyssa.gui.ui.dialogs import dialog_add_model
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.dialogs import dialog_add_sequence_monomer
from pyssa.gui.utilities import gui_utils
from pyssa.gui.utilities import tools
from pyssa.gui.utilities import styles
from pyssa.gui.utilities import gui_page_management
from pyssa.gui.utilities import input_validator
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.utilities.data_classes import stage
from pyssa.pymol_protein_tools import protein
from pyssa.gui.data_structures import settings

# setup logger
logging.basicConfig(level=logging.DEBUG)
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
        self.setMinimumWidth(550)
        self.setMinimumHeight(200)

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

        # ----- All class attributes are listed here
        self.app_project = project.Project("", pathlib.Path(""))
        self._project_watcher = project_watcher.ProjectWatcher(self.app_project, no_of_pdb_files=None)
        self.scratch_path = constants.SCRATCH_DIR
        self.workspace_path = self.app_settings.workspace_path
        self.workspace = Qt.QtWidgets.QLabel(f"Current Workspace: {self.workspace_path}")
        self.status_bar = Qt.QtWidgets.QStatusBar()
        self.prediction_configuration = prediction_configuration.PredictionConfiguration(True, "pdb70")

        self.results_name = ""
        self.no_of_selected_chains = 0
        self.plot_dialog = Qt.QtWidgets.QDialog(self)
        self.view_box = None

        # -- Gui page management vars
        self.local_pred_monomer_management: gui_page_management.GuiPageManagement
        self.local_pred_multimer_management: gui_page_management.GuiPageManagement
        self.single_analysis_management: gui_page_management.GuiPageManagement
        self.batch_analysis_management: gui_page_management.GuiPageManagement
        self.results_management: gui_page_management.GuiPageManagement

        # sets up the status bar
        self._setup_statusbar()
        tools.create_directory(constants.SETTINGS_DIR, "scratch")
        self._setup_default_configuration()

        # management
        self._create_local_pred_monomer_management()
        self._create_local_pred_multimer_management()
        self._create_single_analysis_management()
        self._create_batch_analysis_management()
        self._create_results_management()

        # configure gui element properties
        self.ui.txt_results_aligned_residues.setAlignment(QtCore.Qt.AlignRight)
        self.ui.table_pred_mono_prot_to_predict.setSizeAdjustPolicy(PyQt5.QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.ui.table_pred_mono_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)
        self.ui.table_pred_multi_prot_to_predict.horizontalHeader().setDefaultAlignment(QtCore.Qt.AlignLeft)

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
        # connections
        self._connect_all_gui_elements()
        # create tooltips
        self._create_all_tooltips()
        self._project_watcher.show_valid_options(self.ui)
        # setting additional parameters
        self.ui.lbl_logo.setPixmap(PyQt5.QtGui.QPixmap(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowTitle(f"PySSA {constants.VERSION_NUMBER}")

    # ----- Functions for GuiPageManagement obj creation
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

    # def _create_local_pred_monomer_management(self):
    #     # gui element management
    #     tmp_stages = [
    #         # protein name stage
    #         stage.Stage(
    #             {
    #                 "label_protein_name": self.ui.lbl_local_pred_mono_protein_name,
    #                 "text_field_protein_name": self.ui.txt_local_pred_mono_protein_name,
    #                 "label_protein_name_status": self.ui.lbl_local_pred_mono_status_protein_name,
    #              },
    #             {
    #                 "next_button": self.ui.btn_local_pred_mono_next,
    #             }
    #         ),
    #         # protein sequence stage
    #         stage.Stage(
    #             {
    #                 "label_protein_sequence": self.ui.lbl_local_pred_mono_prot_seq,
    #                 "text_field_protein_sequence": self.ui.txt_local_pred_mono_prot_seq,
    #                 "label_protein_sequence_status": self.ui.lbl_local_pred_mono_status_prot_seq,
    #             },
    #             {
    #                 "back_button": self.ui.btn_local_pred_mono_back,
    #                 "next_button": self.ui.btn_local_pred_mono_next_2,
    #             }
    #         ),
    #         # prediction stage (with advanced configurations)
    #         stage.Stage(
    #             {
    #                 "label_advanced_config": self.ui.lbl_local_pred_mono_advanced_config,
    #                 "button_advanced_config": self.ui.btn_local_pred_mono_advanced_config,
    #             },
    #             {
    #                 "back_button": self.ui.btn_local_pred_mono_back_2,
    #                 "predict_button": self.ui.btn_local_pred_mono_predict,
    #             }
    #         ),
    #     ]
    #     self.local_pred_monomer_management = gui_page_management.GuiPageManagement(tmp_stages)

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

    # def _create_local_pred_multimer_management(self):
    #     # gui element management
    #     tmp_stages = [
    #         # prediction mode: stage 0
    #         stage.Stage(
    #             {
    #                 "label_prediction_mode": self.ui.lbl_local_pred_prediction_mode,
    #              },
    #             {
    #                 "single_button": self.ui.btn_local_pred_multi_single,
    #                 "batch_button": self.ui.btn_local_pred_multi_batch,
    #             }
    #         ),
    #         # protein name: stage 1
    #         stage.Stage(
    #             {
    #                 "label_protein_name": self.ui.lbl_local_pred_multi_protein_name,
    #                 "text_field_protein_name": self.ui.txt_local_pred_multi_protein_name,
    #                 "label_protein_name_status": self.ui.lbl_local_pred_multi_status_protein_name,
    #              },
    #             {
    #                 "back_button": self.ui.btn_local_pred_multi_back_prediction_mode,
    #                 "next_button": self.ui.btn_local_pred_multi_next,
    #             }
    #         ),
    #         # single protein sequence: stage 2
    #         stage.Stage(
    #             {
    #                 "label_protein_sequence": self.ui.lbl_local_pred_multi_prot_seq_single,
    #                 "text_field_protein_sequence": self.ui.txt_local_pred_multi_prot_seq,
    #                 "label_protein_sequence_status": self.ui.lbl_local_pred_multi_status_prot_seq,
    #                 "button_add_sequence": self.ui.btn_local_pred_multi_add_seq_single,
    #                 "label_protein_sequence_overview": self.ui.lbl_local_pred_multi_prot_overview,
    #                 "table_protein_sequence_overview": self.ui.table_local_pred_multi_prot_overview,
    #             },
    #             {
    #                 "back_button": self.ui.btn_local_pred_multi_back,
    #                 "next_button": self.ui.btn_local_pred_multi_next_2,
    #             }
    #         ),
    #         # batch protein sequence: stage 3
    #         stage.Stage(
    #             {
    #                 "label_protein_sequence_overview": self.ui.lbl_local_pred_multi_prot_overview,
    #                 "table_protein_sequence_overview": self.ui.table_local_pred_multi_prot_overview,
    #                 "label_protein_sequence_batch": self.ui.lbl_local_pred_multi_prot_seq_batch,
    #                 "button_add_sequence_batch": self.ui.btn_local_pred_multi_add_seq_batch
    #             },
    #             {
    #                 "back_button": self.ui.btn_local_pred_multi_back,
    #                 "next_button": self.ui.btn_local_pred_multi_next_2,
    #             }
    #         ),
    #         # single prediction: stage 4 (with advanced configurations)
    #         stage.Stage(
    #             {
    #                 "label_advanced_config": self.ui.lbl_local_pred_multi_advanced_config,
    #                 "button_advanced_config": self.ui.btn_local_pred_multi_advanced_config,
    #             },
    #             {
    #                 "back_button": self.ui.btn_local_pred_multi_back_2,
    #                 "predict_button": self.ui.btn_local_pred_multi_predict,
    #             }
    #         ),
    #         # batch prediction: stage 5 (with advanced configurations)
    #         stage.Stage(
    #             {
    #                 "label_advanced_config_batch": self.ui.lbl_local_pred_multi_advanced_config,
    #                 "button_advanced_config_batch": self.ui.btn_local_pred_multi_advanced_config,
    #             },
    #             {
    #                 "back_button": self.ui.btn_local_pred_multi_back_3,
    #                 "predict_button": self.ui.btn_local_pred_multi_predict,
    #             }
    #         ),
    #     ]
    #     self.local_pred_multimer_management = gui_page_management.GuiPageManagement(tmp_stages)

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

    def _setup_statusbar(self):
        """This function sets up the status bar and fills it with the current workspace

        """
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage(str(self.workspace_path))

    def _setup_default_configuration(self):
        self.ui.lbl_current_project_name.setText("")
        # menu
        self.ui.action_add_model.setVisible(False)
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
        # menu
        self.ui.action_file_quit.triggered.connect(self.quit_app)
        self.ui.action_file_restore_settings.triggered.connect(self.restore_settings)
        self.ui.action_settings_edit_all.triggered.connect(self.open_settings_global)
        self.ui.action_add_multiple_models.triggered.connect(self.open_add_models)
        self.ui.action_add_model.triggered.connect(self.open_add_model)
        self.ui.action_help_docs.triggered.connect(self.open_documentation)
        self.ui.action_help_docs_pdf.triggered.connect(self.open_documentation_pdf)
        self.ui.action_help_about.triggered.connect(self.open_about)
        # side menu
        self.ui.btn_new_page.clicked.connect(self.display_new_page)
        self.ui.btn_open_page.clicked.connect(self.display_open_page)
        self.ui.btn_delete_page.clicked.connect(self.display_delete_page)
        self.ui.btn_save_project.clicked.connect(self.save_project)
        self.ui.btn_edit_page.clicked.connect(self.display_view_page)
        self.ui.btn_view_page.clicked.connect(self.display_view_page)
        self.ui.btn_use_page.clicked.connect(self.display_use_page)
        self.ui.btn_close_project.clicked.connect(self.close_project)
        self.ui.btn_pred_local_monomer_page.clicked.connect(self.display_local_pred_mono)
        self.ui.btn_pred_local_multimer_page.clicked.connect(self.display_local_pred_multi)
        #self.ui.btn_prediction_page.clicked.connect(self.display_sequence_vs_pdb_page)
        self.ui.btn_single_analysis_page.clicked.connect(self.display_single_analysis_page)
        self.ui.btn_batch_analysis_page.clicked.connect(self.display_job_analysis_page)
        self.ui.btn_results_page.clicked.connect(self.display_results_page)
        self.ui.btn_image_page.clicked.connect(self.display_image_page)
        self.ui.btn_hotspots_page.clicked.connect(self.display_hotspots_page)
        # new project page
        self.ui.btn_new_choose_reference.clicked.connect(self.load_reference_in_project)
        self.ui.txt_new_project_name.textChanged.connect(self.validate_project_name)
        self.ui.txt_new_choose_reference.textChanged.connect(self.validate_reference_in_project)
        self.ui.cb_new_add_reference.stateChanged.connect(self.show_add_reference)
        self.ui.btn_new_create_project.clicked.connect(self.create_new_project)
        # open project page
        self.ui.btn_open_open_project.clicked.connect(self.open_project)
        self.ui.list_open_projects.doubleClicked.connect(self.open_project)
        self.ui.txt_open_search.textChanged.connect(self.validate_open_search)
        self.ui.txt_open_selected_project.textChanged.connect(self.activate_open_button)
        self.ui.list_open_projects.currentItemChanged.connect(self.select_project_from_open_list)
        # delete project page
        self.ui.btn_delete_delete_project.clicked.connect(self.delete_project)
        self.ui.txt_delete_search.textChanged.connect(self.validate_delete_search)
        self.ui.txt_delete_selected_projects.textChanged.connect(self.activate_delete_button)
        self.ui.list_delete_projects.currentItemChanged.connect(self.select_project_from_delete_list)
        # edit project page
        self.ui.btn_edit_page.clicked.connect(self.display_edit_page)
        self.ui.list_edit_project_proteins.currentItemChanged.connect(self.check_for_cleaning)
        self.ui.btn_edit_clean_new_prot.clicked.connect(self.clean_protein_new)
        self.ui.btn_edit_clean_update_prot.clicked.connect(self.clean_protein_update)
        self.ui.btn_edit_project_delete.clicked.connect(self.delete_protein)
        # view project page
        self.ui.btn_view_project_show.clicked.connect(self.view_sequence)
        self.ui.btn_view_project_show_structure.clicked.connect(self.view_structure)
        self.ui.list_view_project_proteins.doubleClicked.connect(self.view_sequence)
        # use project page
        self.ui.txt_use_project_name.textChanged.connect(self.validate_use_project_name)
        self.ui.btn_use_next.clicked.connect(self.show_protein_selection_for_use)
        self.ui.txt_use_search.textChanged.connect(self.validate_use_search)
        self.ui.btn_use_add_available_protein_structures.clicked.connect(self.add_protein_structure_to_new_project)
        self.ui.list_use_available_protein_structures.doubleClicked.connect(self.add_protein_structure_to_new_project)
        self.ui.btn_use_remove_selected_protein_structures.clicked.connect(self.remove_protein_structure_to_new_project)
        self.ui.list_use_selected_protein_structures.doubleClicked.connect(self.remove_protein_structure_to_new_project)
        self.ui.btn_use_back.clicked.connect(self.hide_protein_selection_for_use)
        self.ui.btn_use_create_new_project.clicked.connect(self.create_use_project)
        # sequence vs .pdb page
        #self.ui.btn_s_v_p_start.clicked.connect(self.predict)
        # monomer local prediction page
        self.ui.btn_pred_mono_seq_to_predict.clicked.connect(self.local_pred_mono_show_protein_name)
        self.ui.btn_pred_mono_seq_to_predict_remove.clicked.connect(self.local_pred_mono_remove_protein_to_predict)
        self.ui.btn_pred_mono_next.clicked.connect(self.local_pred_mono_show_protein_sequence)
        self.ui.btn_pred_mono_back.clicked.connect(self.local_pred_mono_show_protein_overview)
        self.ui.btn_pred_mono_add_protein.clicked.connect(self.local_pred_mono_add_protein_to_predict)
        self.ui.btn_pred_mono_back_2.clicked.connect(self.local_pred_mono_show_protein_name)
        self.ui.txt_pred_mono_prot_name.textChanged.connect(self.local_pred_mono_validate_protein_name)
        self.ui.txt_pred_mono_seq_name.textChanged.connect(self.local_pred_mono_validate_protein_sequence)
        self.ui.btn_pred_mono_advanced_config.clicked.connect(self.show_prediction_configuration)
        # TODO: predict_monomer connection is missing!

        # multimer local prediction page
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
        # self.ui.btn_pred_mono_advanced_config.clicked.connect(self.local_pred_mono_show_prediction_configuration)


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
        self.ui.btn_analysis_next.clicked.connect(self.show_single_analysis_stage_1)
        self.ui.btn_analysis_next_2.clicked.connect(self.show_single_analysis_stage_2)
        self.ui.btn_analysis_back.clicked.connect(self.show_single_analysis_stage_0)
        self.ui.btn_analysis_back_2.clicked.connect(self.show_single_analysis_stage_1)
        self.ui.btn_analysis_start.clicked.connect(self.start_process)
        self.ui.box_analysis_prot_struct_1.currentIndexChanged.connect(self.check_if_prot_structs_are_filled)
        self.ui.box_analysis_prot_struct_2.currentIndexChanged.connect(self.check_if_prot_structs_are_filled)
        self.ui.list_analysis_ref_chains.itemSelectionChanged.connect(self.count_selected_chains_for_prot_struct_1)
        self.ui.list_analysis_model_chains.itemSelectionChanged.connect(self.check_if_same_no_of_chains_selected)
        # batch analysis page
        self.ui.btn_analysis_batch_add.clicked.connect(self.show_batch_analysis_stage_1)
        self.ui.btn_analysis_batch_remove.clicked.connect(self.remove_analysis_run)
        self.ui.btn_analysis_batch_back.clicked.connect(self.show_batch_analysis_stage_0)
        self.ui.btn_analysis_batch_next.clicked.connect(self.show_batch_analysis_stage_2)
        self.ui.btn_analysis_batch_back_2.clicked.connect(self.show_batch_analysis_stage_1)
        self.ui.btn_analysis_batch_next_2.clicked.connect(self.show_batch_analysis_stage_3)
        self.ui.btn_analysis_batch_back_3.clicked.connect(self.show_batch_analysis_stage_2)
        self.ui.btn_analysis_batch_next_3.clicked.connect(self.show_batch_analysis_stage_0)
        self.ui.box_analysis_batch_prot_struct_1.currentIndexChanged.connect(self.check_if_prot_structs_are_filled_batch)
        self.ui.box_analysis_batch_prot_struct_2.currentIndexChanged.connect(self.check_if_prot_structs_are_filled_batch)
        self.ui.list_analysis_batch_ref_chains.itemSelectionChanged.connect(self.count_batch_selected_chains_for_prot_struct_1)
        self.ui.list_analysis_batch_model_chains.itemSelectionChanged.connect(self.check_if_same_no_of_chains_selected_batch)
        self.ui.btn_analysis_batch_start.clicked.connect(self.start_process_batch)
        # results page
        self.ui.cb_results_analysis_options.currentIndexChanged.connect(self.load_results)
        self.ui.btn_view_struct_alignment.clicked.connect(self.display_structure_alignment)
        self.ui.btn_view_distance_plot.clicked.connect(self.display_distance_plot)
        self.ui.btn_view_distance_histogram.clicked.connect(self.display_distance_histogram)
        # self.ui.btn_view_interesting_region.clicked.connect(self.tmp_dialog_change)
        self.ui.btn_view_distance_table.clicked.connect(self.display_distance_table)
        # image page
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
        # hotspots page
        self.ui.list_hotspots_choose_protein.currentItemChanged.connect(self.open_protein)
        self.ui.btn_hotspots_resi_show.clicked.connect(self.show_resi_sticks)
        self.ui.btn_hotspots_resi_hide.clicked.connect(self.hide_resi_sticks)
        self.ui.btn_hotspots_resi_zoom.clicked.connect(self.zoom_resi_position)

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

    # def _init_hide_ui_elements(self):
    #     """This function hides all UI elements which need to be hidden during the
    #     plugin startup
    #
    #     TODO: this function should be obsolete with the introduction of the project_watcher?!
    #     """
    #     gui_elements_to_hide = [
    #         self.ui.lbl_new_choose_reference,
    #         self.ui.txt_new_choose_reference,
    #         self.ui.btn_new_choose_reference,
    #         self.ui.lbl_prediction_load_reference,
    #         self.ui.txt_prediction_load_reference,
    #         self.ui.btn_prediction_load_reference,
    #         self.ui.btn_prediction_next_2,
    #         self.ui.btn_prediction_back_2,
    #         # self.ui.lbl_prediction_ref_chains.hide()
    #         self.ui.lbl_prediction_model_chains,
    #         self.ui.txt_prediction_chain_model,
    #         self.ui.btn_prediction_back_3,
    #         # self.ui.btn_prediction_start.hide()
    #         self.ui.btn_hotspots_page,
    #         # sidebar elements
    #         # self.ui.lbl_prediction.hide()
    #         self.ui.lbl_analysis,
    #         self.ui.lbl_handle_pymol_session,
    #         self.ui.btn_save_project,
    #         # self.ui.btn_prediction_only_page.hide()
    #         # self.ui.btn_prediction_page.hide()
    #         self.ui.btn_single_analysis_page,
    #         self.ui.btn_batch_analysis_page,
    #         self.ui.btn_results_page,
    #         self.ui.btn_image_page,
    #     ]
    #     gui_utils.hide_gui_elements(gui_elements_to_hide)

    def _init_fill_combo_boxes(self):
        """This function fills all combo boxes of the plugin

        """
        item_list_representation = [
            "",
            "cartoon",
            "ribbon",
        ]
        gui_utils.fill_combo_box(self.ui.box_representation, item_list_representation)
        # combo box BgColor
        item_list_bg_color = [
            "",
            "black",
            "white",
        ]
        gui_utils.fill_combo_box(self.ui.box_bg_color, item_list_bg_color)
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

    def _init_new_page(self):
        """This function clears all text fields and hides everything which is needed

        """
        self.ui.txt_new_project_name.clear()
        self.ui.txt_new_choose_reference.clear()
        self.ui.lbl_new_status_project_name.setText("")
        self.ui.lbl_new_status_choose_reference.setText("")

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
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        tools.scan_project_for_valid_proteins(f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}",
                                              self.ui.list_edit_project_proteins)

    def _init_local_pred_mono_page(self):
        # clears everything
        self.ui.txt_pred_mono_prot_name.clear()
        self.ui.txt_pred_mono_seq_name.clear()
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

    def _init_all_pages(self):
        self._init_local_pred_mono_page()
        self._init_local_pred_multi_page()
        self._init_sequence_vs_pdb_page()
        self._init_single_analysis_page()
        self._init_batch_analysis_page()
        self._init_results_page()
        self._init_image_page()

    # Slots
    # def handle_side_menu(self):
    #     """This function is used to hide and show the side menu
    #
    #     """
    #     width = self.ui.side_menu_container.width()
    #     if width == 0:
    #         # runs if sidebar will be opened (current status: closed)
    #         new_width = 170
    #         self.setMinimumWidth(650)
    #         self.setMaximumWidth(12000)
    #         self.ui.btn_side_menu.setText("<-")
    #     else:
    #         # runs if sidebar will be closed (current status: open)
    #         new_width = 0
    #         # self.setFixedWidth(650 - 170)
    #         self.setMinimumWidth(650-170)
    #         self.setMaximumWidth(480)
    #         # self.ui.main_body.setMinimumWidth(650-170)
    #         self.ui.btn_side_menu.setText("->")
    #         self.ui.btn_side_menu.setToolTip("Expand")
    #     self.ui.side_menu_container.setFixedWidth(new_width)
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

    def display_sequence_vs_pdb_page(self):
        """This function displays the sequence vs .pdb page

        """
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 10, "Sequence vs .pdb")

    def display_single_analysis_page(self):
        """This function displays the single analysis work area

        """
        self.fill_protein_structure_boxes()
        self.ui.list_analysis_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.list_analysis_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        # regular work area opening
        self._init_single_analysis_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 3, "Single Analysis")

    def display_job_analysis_page(self):
        """This function displays the job analysis work area

        """
        self.ui.list_analysis_batch_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.list_analysis_batch_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        # regular work area opening
        self._init_batch_analysis_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 4, "Structure Analysis")

    def display_results_page(self):
        """This function displays the results work area

        """
        results = os.listdir(self.app_project.get_results_path())
        results.insert(0, "")
        self.ui.cb_results_analysis_options.clear()
        gui_utils.fill_combo_box(self.ui.cb_results_analysis_options, results)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 5, "Results")
        self.show_analysis_results_options()

    def display_image_page(self):
        """This function displays the image work area

        """
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

    def display_open_page(self):
        """This function displays the open project work area

        """
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
        else:
            gui_utils.error_dialog_settings("The settings file is corrupted. Please restore the settings!", "",
                                            self.app_settings)
            self.display_home_page()

    def display_delete_page(self):
        """This function displays the "delete" project work area

        """
        self.ui.list_delete_projects.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_delete_projects)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 9, "Delete existing project")

    def display_edit_page(self):
        """This function displays the edit project page

        """
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        self._init_edit_page()
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 13, "Edit proteins of current project")

    def display_view_page(self):
        """This function displays the edit project page

        """
        self.ui.list_view_project_proteins.clear()
        self.ui.txtedit_view_sequence.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 11, "View proteins of current project")

        # list all proteins from pdb directory
        tools.scan_project_for_valid_proteins(f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}",
                                              self.ui.list_view_project_proteins)

    def display_local_pred_mono(self):
        self.local_pred_monomer_management.show_stage_x(0)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 19, "Local Monomer Prediction")

    def display_local_pred_multi(self):
        self.local_pred_multimer_management.show_stage_x(0)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 20, "Local Multimer Prediction")

    def display_use_page(self):
        self._init_use_page()
        self.ui.list_use_available_protein_structures.clear()
        self.ui.list_use_selected_protein_structures.clear()
        valid_projects = tools.scan_workspace_for_valid_projects(self.workspace_path, self.ui.list_use_existing_projects)
        # filesystem operations
        tools.scan_project_for_valid_proteins(pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}"),
                                              self.ui.list_use_selected_protein_structures)
        protein_dict, protein_names = tools.scan_workspace_for_non_duplicate_proteins(valid_projects,
                                                                                      self.ui.lbl_current_project_name.text(),
                                                                                      self.workspace_path,
                                                                                      self.ui.list_use_available_protein_structures)
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

    def display_hotspots_page(self):
        self.ui.list_hotspots_choose_protein.clear()
        tools.scan_project_for_valid_proteins(self.app_project.project_path, self.ui.list_hotspots_choose_protein)
        tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 18, "Hotspots")

    # def __check_start_possibility(self):
    #     """This function is used to determine if the Start button can be
    #     enabled for the single analysis.
    #
    #     """
    #     i = 0
    #     j = 0
    #     if len(str(self.ui.txt_analysis_project_name.text())) > 0:
    #         i += 1
    #     if len(str(self.ui.txt_analysis_load_reference.text())) > 0:
    #         i += 1
    #     if len(str(self.ui.txt_analysis_load_model.toPlainText())) > 0:
    #         i += 1
    #     if self.ui.cb_analysis_chain_info.isChecked() and len(self.ui.txt_analysis_chain_ref.text()) > 0 and len(
    #             self.ui.txt_analysis_chain_model.text()) > 0:
    #         i += 1
    #     if self.ui.cb_analysis_chain_info.isChecked():
    #         j = 4
    #     else:
    #         j = 3
    #     if i == j:
    #         self.ui.btn_analysis_start.setEnabled(True)
    #         with open('styles/styles_start_button_ready.css', 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_analysis_start.setStyleSheet(button_style)
    #     else:
    #         self.ui.btn_analysis_start.setEnabled(False)
    #         with open('styles/styles_start_button_not_ready.css', 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_analysis_start.setStyleSheet(button_style)
    #
    # def __check_start_possibility_batch(self):
    #     """This function is used to determine if the Start button can be
    #     enabled for the single analysis.
    #
    #     """
    #     i = 0
    #     j = 0
    #     if len(str(self.ui.txt_batch_job_name.text())) > 0:
    #         i += 1
    #     if len(str(self.ui.txt_batch_load_reference.text())) > 0:
    #         i += 1
    #     if len(str(self.ui.txt_batch_load_model.toPlainText())) > 0:
    #         i += 1
    #     if self.ui.cb_batch_chain_info.isChecked() and len(self.ui.txt_batch_chain_ref.text()) > 0 and len(
    #             self.ui.txt_batch_chain_model.text()) > 0:
    #         i += 1
    #     if self.ui.cb_batch_chain_info.isChecked():
    #         j = 4
    #     else:
    #         j = 3
    #     if i == j:
    #         self.ui.btn_batch_start.setEnabled(True)
    #         with open('styles/styles_start_button_ready.css', 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_batch_start.setStyleSheet(button_style)
    #     else:
    #         self.ui.btn_batch_start.setEnabled(False)
    #         with open('styles/styles_start_button_not_ready.css', 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_batch_start.setStyleSheet(button_style)
    #
    # def __check_start_possibility_prediction(self):
    #     """This function is used to determine if the Start button can be
    #     enabled for the prediction + analysis.
    #
    #     """
    #     # creates path for specific stylesheets
    #     if sys.platform.startswith("darwin"):
    #         # macOS path
    #         styles_path_btn_ready = f"{constants.path_list[1]}/styles/styles_start_button_ready.css"
    #         styles_path_btn_not_ready = f"{constants.path_list[1]}/styles/styles_start_button_not_ready.css"
    #     elif sys.platform.startswith("linux"):
    #         # Linux path
    #         styles_path_btn_ready = f"{constants.path_list[0]}/styles/styles_start_button_ready.css"
    #         styles_path_btn_not_ready = f"{constants.path_list[0]}/styles/styles_start_button_not_ready.css"
    #     elif sys.platform.startswith("win32"):
    #         # Windows path
    #         styles_path_btn_ready = f"{constants.path_list[2]}/styles/styles_start_button_ready.css"
    #         styles_path_btn_not_ready = f"{constants.path_list[2]}/styles/styles_start_button_not_ready.css"
    #
    #     i = 0
    #     j = 0
    #     if len(str(self.ui.txt_prediction_project_name.text())) > 0:
    #         i += 1
    #     if len(str(self.ui.txt_prediction_load_reference.text())) > 0:
    #         i += 1
    #     if self.ui.cb_prediction_chain_info.isChecked() and len(self.ui.txt_prediction_chain_ref.text()) > 0 and len(
    #             self.ui.txt_prediction_chain_model.text()) > 0:
    #         i += 1
    #     if self.ui.cb_prediction_chain_info.isChecked():
    #         j = 3
    #     else:
    #         j = 2
    #     if i == j:
    #         self.ui.btn_prediction_start.setEnabled(True)
    #         with open(styles_path_btn_ready, 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_prediction_start.setStyleSheet(button_style)
    #     else:
    #         self.ui.btn_prediction_start.setEnabled(False)
    #         with open(styles_path_btn_not_ready, 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_prediction_start.setStyleSheet(button_style)
    #
    # def __check_start_possibility_prediction_only(self):
    #     """This function is used to determine if the Start button can be
    #     enabled for the single analysis.
    #
    #     """
    #     if len(str(self.ui.txt_prediction_only_notebook_url.text())) > 0:
    #         self.ui.btn_prediction_only_start.setEnabled(True)
    #         with open('styles/styles_start_button_ready.css', 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_prediction_only_start.setStyleSheet(button_style)
    #     else:
    #         self.ui.btn_prediction_only_start.setEnabled(False)
    #         with open('styles/styles_start_button_not_ready.css', 'r') as style_sheet_file:
    #             button_style = style_sheet_file.read()
    #
    #             # Set the stylesheet of the application
    #             self.ui.btn_prediction_only_start.setStyleSheet(button_style)

    # @SLOT
    # Menu
    # def open(self):
    #     """This function opens a project.xml file and fills the input boxes with the right values
    #
    #     """
    #     # open file dialog
    #     try:
    #         global global_var_project_dict
    #         file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self,
    #                                                              "Open project file",
    #                                                              Qt.QtCore.QDir.homePath(),
    #                                                              "Plugin Project File (*.xml)")
    #         if file_name == ("", ""):
    #             tools.quick_log_and_display("info", "No file has been selected.",
    #                                         self.status_bar, "No file has been selected.")
    #             return
    #         # clear current project list and interesting regions
    #         self.ui.project_list.clear()
    #         self.ui.cb_interesting_regions.clear()
    #
    #         xml_file = minidom.parse(file_name[0])
    #         element = xml_file.getElementsByTagName('project_name')
    #         job_file_info = Qt.QtCore.QFileInfo(file_name[0])
    #         job_name = job_file_info.baseName()
    #
    #         if len(element) == 0:
    #             # job
    #             i = 0
    #             loop_element = ["", ""]
    #             while loop_element != []:
    #                 loop_element = xml_file.getElementsByTagName(f'project_{i}')
    #                 if loop_element != []:
    #                     path = loop_element[0].getAttribute('value')
    #                     path_split = path.split("/")
    #                     project_name = path_split[(len(path_split) - 1)]
    #                     global_var_project_dict[i] = project.Project(project_name,
    #                                                                              f"{self.workspace_path}/{job_name}")
    #                     index = len("project_of_")
    #                     model_name = project_name[index:]
    #                     # TO-DO: insert modelname
    #                     global_var_project_dict[i].set_pdb_model(model_name)
    #                     self.ui.project_list.addItem(global_var_project_dict[i].get_project_name())
    #                 i += 1
    #
    #         else:
    #             # project
    #             file_name_file_info = Qt.QtCore.QFileInfo(file_name[0])
    #             global_var_project_dict[0] = project.Project(file_name_file_info.baseName(),
    #                                                                                    self.workspace_path)
    #             global_var_project_dict[0].set_pdb_model(xml_file.getElementsByTagName('pdb_model')[0].getAttribute('value'))
    #             # add filename to project list (results tab)
    #             self.ui.project_list.addItem(global_var_project_dict[0].get_project_name())
    #
    #         # fill combo box of interesting regions
    #         results_path = global_var_project_dict[0].get_results_path()
    #         dir_content = os.listdir(f"{results_path}/images/interesting_regions")
    #         for tmp_file in dir_content:
    #             self.ui.cb_interesting_regions.addItem(tmp_file)
    #         self.ui.project_list.setCurrentRow(0)
    #
    #     except FileNotFoundError:
    #         print("File could not be opened.")
    #     except ValueError:
    #         print("No file has been selected.")
    #
    # def save_as(self):
    #     """This function saves the current pymol session.
    #
    #     """
    #     try:
    #         file_path = Qt.QtWidgets.QFileDialog.getSaveFileName(self,
    #                                                              "Save PyMOL session",
    #                                                              f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}",
    #                                                              "PyMOL session file (.pse)")
    #         if file_path == ("", ""):
    #             tools.quick_log_and_display("info", "No file has been created.", self.status_bar,
    #                                         "No file has been created.")
    #         cmd.save(f"{file_path[0]}.pse")
    #         tools.quick_log_and_display("info", "Saving the pymol session was successful.",
    #                                     self.status_bar, "Saving was successful.")
    #     except FileExistsError:
    #         tools.quick_log_and_display("warning", "File already exists!", self.status_bar,
    #                                     "File already exists!")
    #     except pymol.CmdException:
    #         tools.quick_log_and_display("error", "Unexpected Error from PyMOL while saving the "
    #                                              "current pymol session", self.status_bar,
    #                                     "Unexpected Error from PyMOL")
    #
    # def save(self):
    #     """This function saves the current pymol session.
    #
    #     """
    #     global global_var_project_dict
    #     try:
    #         file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
    #         cmd.save(f"{file_path}/sessions/session_file_model_s.pse")
    #         tools.quick_log_and_display("info", "Saving the pymol session was successful.",
    #                                     self.status_bar, "Saving was successful.")
    #     except KeyError:
    #         tools.quick_log_and_display("error", "No project has been opened.", self.status_bar,
    #                                     "No project has been opened.")
    #     except pymol.CmdException:
    #         tools.quick_log_and_display("error", "Unexpected Error from PyMOL while saving the "
    #                                              "current pymol session", self.status_bar,
    #                                     "Unexpected Error from PyMOL")
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

    # def prepare_model_pdbs(self):
    #     """This function extracts and moves the .pdb files directly from the
    #     download directory to another directory
    #
    #     """
    #     try:
    #         # path to download directory
    #         tmp_path = f"{os.path.expanduser('~')}/Downloads"
    #         os.chdir(tmp_path)
    #         # path to store all pdb files (should be defined in the Settings tab!)
    #         target_dir = f"{os.path.expanduser('~')}/Documents/data"
    #
    #         # tmp_list contains a list with all zips starting with prediction
    #         tmp_list = tools.filter_prediction_zips(tmp_path)
    #
    #         # create temporary dir_settings to handle unzipping
    #         tmp_dir = f"{tmp_path}/tmp"
    #         if not os.path.exists(tmp_dir):
    #             os.mkdir(tmp_dir)
    #
    #         for file in tmp_list:
    #             tools.extract_and_move_model_pdb(tmp_path, tmp_dir, file, target_dir)
    #
    #         shutil.rmtree(tmp_dir)
    #         self.status_bar.showMessage("Preparing the .pdb files was successful.")
    #     except Exception:
    #         self.status_bar.showMessage("Preparing the .pdb files failed!")

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
        #webbrowser.open_new(f"file://{os.getcwd()}/docs/pymol_plugin/build/html/index.html")
        # opens the documentation of the os
        if sys.platform.startswith("darwin"):
            # macOS path
            webbrowser.open_new(f"file://{constants.path_list[1]}/docs/pymol_plugin/build/html/index.html")
        elif sys.platform.startswith("linux"):
            # Linux path
            webbrowser.open_new(f"file://{constants.path_list[0]}/docs/pymol_plugin/build/html/index.html")
        elif sys.platform.startswith("win32"):
            # Windows path
            webbrowser.open_new(f"file://{constants.path_list[2]}/docs/pymol_plugin/build/html/index.html")

    @staticmethod
    def open_documentation_pdf():
        """This function opens the official plugin documentation as PDF.

        """
        # webbrowser.open_new(
        #     f"file://{os.getcwd()}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")
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
                pdb_path_info = Qt.QtCore.QFileInfo(pdb_path)
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

    def open_add_model(self):
        """This function opens the add model dialog.

        """
        dialog = dialog_add_model.DialogAddModel()
        dialog.exec_()
        if dialog_add_model.global_var_add_model[1] is True:
            file_info = Qt.QtCore.QFileInfo(dialog_add_model.global_var_add_model[0])
            shutil.copy(dialog_add_model.global_var_add_model[0],
                        f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{file_info.baseName()}")
            print("Model was copied to project.")
            self.ui.lbl_analysis.show()
            self.ui.btn_single_analysis_page.show()
            self.ui.action_add_model.setVisible(False)
            self.ui.action_add_multiple_models.setVisible(False)
        else:
            print("No model was added.")

    # ----- Functions for New project page
    def show_add_reference(self):
        """This function shows the reference input section

        """
        # checkbox is checked
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
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self.ui.txt_new_choose_reference.setText(str(file_name[0]))
            self.ui.txt_new_choose_reference.setEnabled(False)
            self.ui.txt_new_choose_reference.setStyleSheet("background-color: #33C065")
            self.ui.btn_new_create_project.setEnabled(True)
            styles.color_button_ready(self.ui.btn_new_create_project)
        except ValueError:
            print("No file has been selected.")

    def validate_reference_in_project(self):
        """This function checks if the entered reference protein is
        valid or not.

        """
        if len(self.ui.txt_new_choose_reference.text()) == 0:
            self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
            self.ui.lbl_new_status_choose_reference.setText("")
            self.ui.btn_new_create_project.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_new_create_project)
        elif len(self.ui.txt_new_choose_reference.text()) < 4:
            self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
            styles.color_button_not_ready(self.ui.btn_new_create_project)
            self.ui.btn_new_create_project.setEnabled(False)
        # checks if a pdb id was entered
        elif len(self.ui.txt_new_choose_reference.text()) == 4:
            pdb_id = self.ui.txt_new_choose_reference.text().upper()
            try:
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                cmd.fetch(pdb_id, type="pdb", path=self.scratch_path)
                os.remove(f"{self.scratch_path}/{pdb_id}.pdb")
                cmd.reinitialize()
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #33C065")
                self.ui.btn_new_create_project.setEnabled(True)
                styles.color_button_ready(self.ui.btn_new_create_project)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
                styles.color_button_not_ready(self.ui.btn_new_create_project)
                return
            except FileNotFoundError:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
                self.ui.lbl_new_status_choose_reference.setText("Invalid PDB ID.")
                self.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_new_create_project)
                return
        else:
            if self.ui.txt_new_choose_reference.text().find("/") == -1:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_new_create_project.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_new_create_project)

            elif self.ui.txt_new_choose_reference.text().find("\\") == -1:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
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
        self.ui.lbl_current_project_name.setText(self.ui.txt_new_project_name.text())
        self.status_bar.showMessage(f"Current project path: {self.workspace_path}/{self.ui.txt_new_project_name.text()}")
        # save project folder in current workspace
        self.app_project = project.Project(self.ui.txt_new_project_name.text(), self.workspace_path)
        self.app_project.create_project_tree()
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

                shutil.move(f"{self.scratch_path}/{pdb_id}.pdb", self.app_project.get_pdb_path())
                tmp_ref_protein = protein.Protein(pdb_id,
                                                  filepath=pathlib.Path(self.app_project.get_pdb_path()),
                                                  export_data_dir=self.scratch_path)
                tmp_ref_protein.clean_pdb_file()
            else:
                # local pdb file as input
                shutil.copy(self.ui.txt_new_choose_reference.text(), self.app_project.get_pdb_path())
                protein_file_info = QtCore.QFileInfo(self.ui.txt_new_choose_reference.text())
                pdb_id = protein_file_info.baseName()
                tmp_ref_protein = protein.Protein(pdb_id,
                                                  filepath=pathlib.Path(self.app_project.get_pdb_path()),
                                                  export_data_dir=self.scratch_path)
                cmd.load(self.ui.txt_new_choose_reference.text(), object=pdb_id)
            tmp_ref_protein.set_chains()
            tmp_ref_protein.set_sequence()
            chains = cmd.get_chains(pdb_id)
            cmd.reinitialize()
            self.app_project.add_existing_protein(tmp_ref_protein)
            tmp_ref_protein.serialize_protein(self.app_project.get_objects_proteins_path(), pdb_id)
            for chain in chains:
                self.ui.list_s_v_p_ref_chains.addItem(chain)
        self.app_project.serialize_project(self.app_project.project_path, "project")
        # shows options which can be done with the data in the project folder
        self._project_watcher.current_project = self.app_project
        self._project_watcher.show_valid_options(self.ui)
        self.display_view_page()

    # ----- Functions for Open project page
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
        # self._init_project_management()
        tmp_project_path = pathlib.Path(f"{self.workspace_path}/{self.ui.list_open_projects.currentItem().text()}")
        self.app_project = project.Project.deserialize_project(tmp_project_path)
        self._project_watcher.current_project = self.app_project
        # if self.app_project.get_number_of_proteins() > 0:
        #     tmp_proteins = os.listdir(self.app_project.get_pdb_path())
        #     self.app_project.proteins.clear()
        #     self.app_project.protein_pairs.clear()
        #     for single_protein in tmp_proteins:
        #         tmp_protein_name = single_protein.replace(".pdb", "")
        #         json_path = pathlib.Path(f"{self.app_project.get_objects_proteins_path()}/{tmp_protein_name}.json")
        #         self.app_project.add_existing_protein(protein.Protein.deserialize_protein(json_path))
        self.ui.lbl_current_project_name.setText(self.app_project.get_project_name())
        self._project_watcher.show_valid_options(self.ui)
        self.display_view_page()
        # check if directory in empty
        # TO-DO: code below does not work with new project class!
        # results_paths = ["results/alignment_files"]
        # dir_content = os.listdir(f"{self.workspace_path}/{self.ui.list_open_projects.currentItem().text()}/{results_paths[0]}")
        # self.ui.lbl_handle_pymol_session.show()
        # self.ui.btn_image_page.show()
        # if not dir_content:
        #     # no analysis was done
        #     tmp_content = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")
        #     if len(tmp_content) == 1:
        #         # if no model is in the pdb dir
        #         self.display_image_page()
        #         self.ui.btn_save_project.show()
        #         self.ui.action_file_save_as.setVisible(True)
        #         self.ui.action_add_model.setVisible(True)
        #         self.ui.action_add_multiple_models.setVisible(True)
        #     elif len(tmp_content) == 2:
        #         # if a model is in the pdb dir
        #         self.display_single_analysis_page()
        #         self.ui.lbl_analysis.show()
        #         self.ui.btn_single_analysis_page.show()
        #         self.ui.btn_save_project.show()
        #         self.ui.action_file_save_as.setVisible(True)
        #         self.ui.action_add_model.setVisible(False)
        #         self.ui.action_add_multiple_models.setVisible(False)
        #     elif len(tmp_content) == 0:
        #         # check if data is corrupted
        #         response = tools.check_results_for_integrity(self.workspace_path,
        #                                                      self.ui.lbl_current_project_name.text())
        #         if response[0]:
        #             self.ui.lbl_handle_pymol_session.hide()
        #             self.ui.btn_image_page.hide()
        #             self.display_prediction_only_page()
        #
        #     else:
        #         gui_utils.error_project_data_is_invalid(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")
        #     return
        # else:
        #     # check if data is corrupted
        #     response = tools.check_results_for_integrity(self.workspace_path, self.ui.lbl_current_project_name.text())
        #     if response[0]:
        #         # an analysis happened
        #         # show results gui elements
        #         self.ui.lbl_analysis.show()
        #         self.ui.btn_results_page.show()
        #         self.display_results_page()
        #         self.ui.btn_save_project.show()
        #         self.ui.action_file_save_as.setVisible(True)
        #     else:
        #         gui_utils.error_project_data_is_invalid(response[1])
        #         self._init_hide_ui_elements()
        #         self.ui.lbl_current_project_name.clear()
        #         self.ui.btn_save_project.hide()
        #     return

    # ----- Functions for Delete project page
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

        if response is True:
            shutil.rmtree(pathlib.Path(f"{self.workspace_path}/{self.ui.txt_delete_selected_projects.text()}"))
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
        else:
            return

    # ----- Functions for Save project
    def save_project(self):
        """This function saves the "project" which is currently only the pymol session

        """
        pass
        # if not os.path.exists(pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/wsl/")):
        #     os.mkdir(pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/wsl/"))
        # if not os.path.exists(constants.WSL_STORAGE_PATH):
        #     os.mkdir(constants.WSL_STORAGE_PATH)
        # print(subprocess.run(["wsl", "--import", constants.WSL_DISTRO_NAME, str(constants.WSL_STORAGE_PATH),
        #                       str(constants.WSL_DISTRO_IMPORT_PATH)]))
        # print(subprocess.run(["C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe", pathlib.Path(
        #     f"{os.path.expanduser('~')}/github_repos/tmpPySSA/pyssa/scripts/install_wsl_env.ps1")]))

        # sequence = "MAQAKHKQRKRLKSSCKRHPLYVDFSDVGWNDWIVAPPGYHAFYCHGECPFPLADHLNSTNHAIVQTLVNSVNSKIPKACCVPTELSAISMLYLDENEKVVLKNYQDMVVEGCGCRMAQAKHKQRKRLKSSCKRHPLYVDFSDVGWNDWIVAPPGYHAFYCHGECPFPLADHLNSTNHAIVQTLVNSVNSKIPKACCVPTELSAISMLYLDENEKVVLKNYQDMVVEGCGCR"
        # print(subprocess.run(["wsl", "curl", "-X", "POST", "--data", f"{sequence}", "https://api.esmatlas.com/foldSequence/v1/pdb/", "-OutFile", "protein.pdb"]))

        # new_project = project.Project(self.ui.txt_use_project_name.text(), pathlib.Path(self.workspace_path))
        # new_project.serialize_project("C:\\Users\\martin\\github_repos\\tmpPySSA", "testproject")
        #
        # opened_project = project.Project("", pathlib.Path("")).deserialize_project("C:\\Users\\martin\\github_repos\\tmpPySSA\\testproject.json")
        # opened_project.create_folder_paths()
        # # old saving process
        # session_file = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/sessions/")
        # try:
        #     cmd.save(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/sessions/{session_file[0]}")
        #     self.status_bar.showMessage("Saved project successfully.")
        # except pymol.CmdException:
        #     self.status_bar.showMessage("PyMOL internal error while saving.")
        # except IndexError:
        #     cmd.save(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/sessions/auto_generated_session_file.pse")
        #     self.status_bar.showMessage("Saved project successfully.")

    # ----- Functions for Edit project page
    def check_for_cleaning(self):
        try:
            self.ui.list_edit_project_proteins.currentItem().text()
        except AttributeError:
            return
        cmd.reinitialize()
        tmp_protein = self.app_project.search_protein(
            self.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""))
        tmp_protein.load_protein()
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

    def clean_protein_new(self):
        tmp_protein = self.app_project.search_protein(self.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""))
        clean_tmp_protein = tmp_protein.clean_protein(new_protein=True)
        self.app_project.add_existing_protein(clean_tmp_protein)
        self.app_project.serialize_project(self.app_project.project_path, "project")
        self._init_edit_page()

    def clean_protein_update(self):
        tmp_protein = self.app_project.search_protein(self.ui.list_edit_project_proteins.currentItem().text().replace(".pdb", ""))
        if basic_boxes.yes_or_no("Clean protein",
                                 "Are you sure you want to clean this protein?\n"
                                 "This will remove all organic and solvent components!",
                                 QMessageBox.Information):
            tmp_protein.clean_protein()
            self._init_edit_page()
        else:
            return

    def delete_protein(self):
        protein_name = self.ui.list_edit_project_proteins.currentItem().text()
        project_path = f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}"
        response = gui_utils.warning_message_protein_gets_deleted()
        if response:
            tools.remove_pdb_file(f"{project_path}/pdb/{protein_name}")
            self.ui.list_edit_project_proteins.clear()
            tools.scan_project_for_valid_proteins(
                pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}"),
                self.ui.list_edit_project_proteins)
            self._project_watcher.show_valid_options(self.ui)
        else:
            print("Nothing happened.")

    # ----- Functions for View project page
    def view_sequence(self):
        project_path = f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}"
        tmp_protein_filename = self.ui.list_view_project_proteins.currentItem().text()
        sequence = tools.get_sequence_from_pdb_file(f"{project_path}/pdb/{tmp_protein_filename}")
        self.ui.txtedit_view_sequence.clear()
        self.ui.txtedit_view_sequence.append(sequence)
        # fixme: experimental sequence viewer gui
        dialog = dialog_sequence_viewer.SequenceViewer(sequence, tmp_protein_filename)
        dialog.exec_()

    def view_structure(self):
        protein_name = self.ui.list_view_project_proteins.currentItem().text()
        protein_path = pathlib.Path(f"{self.app_project.get_pdb_path()}/{protein_name}")
        # TODO: ask if the session should be saved
        cmd.reinitialize()
        try:
            cmd.load(protein_path)
        except pymol.CmdException:
            print("Error while loading protein in PyMOL!")

    # ----- Functions for Use project page
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
        existing_project = self.app_project
        new_project = project.Project(self.ui.txt_use_project_name.text(), self.workspace_path)
        new_project.create_project_tree()
        self.app_project = new_project
        # copy proteins in new project
        prots_to_copy = []
        for i in range(self.ui.list_use_selected_protein_structures.count()):
            self.ui.list_use_selected_protein_structures.setCurrentRow(i)
            prots_to_copy.append(self.ui.list_use_selected_protein_structures.currentItem().text())
        for tmp_protein in prots_to_copy:
            protein_path = global_variables.global_var_workspace_proteins[tmp_protein]
            shutil.copy(protein_path, f"{self.app_project.get_pdb_path()}/{tmp_protein}")
            new_protein = protein.Protein(tmp_protein, filepath=pathlib.Path(self.app_project.get_pdb_path()))
            new_protein.serialize_protein(self.app_project.get_objects_proteins_path(), tmp_protein)
            self.app_project.add_existing_protein(new_protein)
        self.app_project.serialize_project(self.app_project.project_path, "project")
        # shows options which can be done with the data in the project folder
        self._project_watcher.current_project = self.app_project
        self._project_watcher.show_valid_options(self.ui)
        self._init_use_page()
        self.display_view_page()

    # ----- Functions for Close project
    def close_project(self):
        self._project_watcher.current_project = project.Project("", pathlib.Path(""))
        self._project_watcher.show_valid_options(self.ui)
        self.ui.lbl_current_project_name.setText("")
        self._init_all_pages()
        self.display_home_page()

    # Prediction + Analysis
    def show_prediction_load_reference(self):
        """Shows the text field and tool button for the load reference functionality

        """
        gui_elements_show = [
            self.ui.lbl_prediction_load_reference,
            self.ui.txt_prediction_load_reference,
            self.ui.btn_prediction_load_reference,
            self.ui.btn_prediction_next_2,
            self.ui.btn_prediction_back_2,
        ]
        gui_elements_hide = [
            self.ui.btn_prediction_next_1
        ]
        gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
        self.ui.txt_prediction_project_name.setEnabled(False)
        self.ui.list_widget_projects.setEnabled(False)
        self.ui.btn_prediction_next_2.setEnabled(False)
        self.ui.txt_prediction_load_reference.setEnabled(True)
        self.ui.txt_prediction_project_name.setStyleSheet("background-color: white")

    def hide_prediction_load_reference(self):
        """Hides the text field and tool button for the load reference functionality

        """
        gui_elements_show = [
            self.ui.btn_prediction_next_1,
        ]
        gui_elements_hide = [
            self.ui.lbl_prediction_load_reference,
            self.ui.txt_prediction_load_reference,
            self.ui.btn_prediction_load_reference,
            self.ui.btn_prediction_next_2,
            self.ui.btn_prediction_back_2,
        ]
        gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
        self.ui.txt_prediction_project_name.setEnabled(True)
        self.ui.list_widget_projects.setEnabled(True)
        self.ui.txt_prediction_load_reference.clear()
        # color green
        self.ui.txt_prediction_project_name.setStyleSheet("background-color: #33C065")

    def load_reference_for_prediction(self):
        """This function loads a reference for the analysis part

        """
        try:
            # open file dialog
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self.ui.txt_prediction_load_reference.setText(str(file_name[0]))
            self.ui.txt_prediction_load_reference.setEnabled(False)
            self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #33C065")
            self.ui.btn_prediction_next_2.setEnabled(True)
            styles.color_button_ready(self.ui.btn_prediction_next_2)
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the reference failed!")
        except ValueError:
            print("No file has been selected.")

    def validate_reference_for_prediction(self):
        """This function checks if the entered reference protein is
        valid or not.

        """
        if len(self.ui.txt_prediction_load_reference.text()) == 0:
            self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
            self.ui.btn_prediction_next_2.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_prediction_next_2)
        elif len(self.ui.txt_prediction_load_reference.text()) < 4:
            self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
            self.ui.btn_prediction_next_2.setEnabled(False)
            styles.color_button_not_ready(self.ui.btn_prediction_next_2)
        # checks if a pdb id was entered
        elif len(self.ui.txt_prediction_load_reference.text()) == 4:
            pdb_id = self.ui.txt_prediction_load_reference.text().upper()
            try:
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                cmd.fetch(pdb_id, type="pdb", path=self.scratch_path)
                os.remove(f"{self.scratch_path}/{pdb_id}.pdb")
                cmd.reinitialize()
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #33C065")
                self.ui.btn_prediction_next_2.setEnabled(True)
                styles.color_button_ready(self.ui.btn_prediction_next_2)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_prediction_next_2)
                return
            except FileNotFoundError:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_prediction_next_2)
                self.ui.lbl_prediction_status_load_reference.setText("Invalid PDB ID.")
                return
        else:
            if self.ui.txt_prediction_load_reference.text().find("/") == -1:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_prediction_next_2)

            elif self.ui.txt_prediction_load_reference.text().find("\\") == -1:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles.color_button_not_ready(self.ui.btn_prediction_next_2)

    def enable_chain_information_input_for_prediction(self):
        """This function enables the text boxes to enter the chains for the
        reference and the model

        """
        try:
            self.ui.txt_prediction_chain_ref.setEnabled(self.ui.cb_prediction_chain_info.checkState())
            self.ui.txt_prediction_chain_model.setEnabled(self.ui.cb_prediction_chain_info.checkState())
            self.status_bar.showMessage("Enter the chain information.")
            self.__check_start_possibility_prediction()
        except Exception:
            self.status_bar.showMessage("Unexpected Error.")

    # def show_prediction_chain_info_reference(self):
    #     """This function shows all gui elements for the chain information
    #     selection of the reference protein
    #
    #     """
    #     self.ui.list_widget_ref_chains.setStyleSheet("border-style: solid;"
    #                                                  "border-width: 2px;"
    #                                                  "border-radius: 8px;"
    #                                                  "border-color: #DCDBE3;")
    #     # show and hide relevant gui elements
    #     self.ui.btn_prediction_back_2.hide()
    #     self.ui.btn_prediction_next_2.hide()
    #     self.ui.lbl_prediction_ref_chains.show()
    #     self.ui.list_widget_ref_chains.show()
    #     self.ui.btn_prediction_back_3.show()
    #     self.ui.btn_prediction_start.show()
    #     # disable important gui elements
    #     self.ui.txt_prediction_load_reference.setEnabled(False)
    #     self.ui.btn_prediction_load_reference.setEnabled(False)
    #     # colors white
    #     self.ui.txt_prediction_load_reference.setStyleSheet("background-color: white")
    #
    #     # fill chains list widget
    #     if len(self.ui.txt_prediction_load_reference.text()) == 4:
    #         pdb_id = self.ui.txt_prediction_load_reference.text().upper()
    #         tmp_ref_protein = core.Protein(pdb_id, export_data_dir=self.scratch_path)
    #         tmp_ref_protein.clean_pdb_file()
    #         chains = cmd.get_chains(pdb_id)
    #     else:
    #         cmd.load(self.ui.txt_prediction_load_reference.text(), object="reference_protein")
    #         chains = cmd.get_chains("reference_protein")
    #     for chain in chains:
    #         self.ui.list_widget_ref_chains.addItem(chain)
    #     styles.color_button_ready(self.ui.btn_prediction_start)
    #     cmd.reinitialize()

    # def hide_prediction_chain_info_reference(self):
    #     """Hides the gui elements for the reference chains
    #
    #     """
    #     self.ui.list_widget_ref_chains.clear()
    #     # show and hide relevant gui elements
    #     self.ui.btn_prediction_back_2.show()
    #     self.ui.btn_prediction_next_2.show()
    #     self.ui.lbl_prediction_ref_chains.hide()
    #     self.ui.list_widget_ref_chains.hide()
    #     self.ui.btn_prediction_back_3.hide()
    #     self.ui.btn_prediction_start.hide()
    #     # disable important gui elements
    #     self.ui.txt_prediction_load_reference.setEnabled(True)
    #     self.ui.btn_prediction_load_reference.setEnabled(True)
    #     # colors green
    #     self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #33C065")

    # def check_prediction_if_txt_prediction_chain_ref_is_filled(self):
    #     """This function checks if any chains are in the text field for the
    #     reference.
    #
    #     """
    #     self.__check_start_possibility_prediction()
    #
    # def check_prediction_if_txt_prediction_chain_model_is_filled(self):
    #     """This function checks if any chains are in the text field for the
    #     model.
    #
    #     """
    #     self.__check_start_possibility_prediction()
    #
    # def check_prediction_if_txt_prediction_project_name_is_filled(self):
    #     """This function checks if the project name text field is filled.
    #
    #     """
    #     self.__check_start_possibility_prediction()
    #
    # def check_prediction_if_txt_prediction_load_reference_is_filled(self):
    #     """This function checks if a reference pdb file is selected.
    #
    #     """
    #     self.__check_start_possibility_prediction()

    def enable_predict_button(self):
        self.ui.btn_prediction_only_start.setEnabled(True)
        styles.color_button_ready(self.ui.btn_prediction_only_start)

    # def predict(self):
    #     """This function opens a webbrowser with a colab notebook, to run the
    #     prediction. In addition, it runs the entire analysis after the
    #     prediction.
    #
    #     """
    #     ref_chain_items = self.ui.list_widget_ref_chains.selectedItems()
    #     ref_chains = []
    #     for chain in ref_chain_items:
    #         ref_chains.append(chain.text())
    #     # global global_var_abort_prediction
    #     # global_var_abort_prediction = False
    #     # check if a prediction is already finished
    #     if os.path.isfile(f"{global_variables.global_var_settings_obj.get_prediction_path()}/prediction.zip"):
    #         self.status_bar.showMessage(
    #             f"Warning! | Current Workspace: {self.workspace_path}")
    #         check = gui_utils.warning_message_prediction_exists(
    #             f"The prediction is here: {global_variables.global_var_settings_obj.get_prediction_path()}/prediction.zip ",
    #             f"{global_variables.global_var_settings_obj.get_prediction_path()}/prediction.zip")
    #         if not check:
    #             return
    #     # creates project without xml creation and model adding these come after the prediction
    #     project_obj = project.Project(self.ui.txt_prediction_project_name.text(),
    #                                                         self.workspace_path)
    #     project_obj.create_project_tree()
    #     project_obj.set_pdb_file(self.ui.txt_prediction_load_reference.text())
    #     project_obj.set_pdb_id(self.ui.txt_prediction_load_reference.text())
    #     project_obj.set_ref_chains(ref_chains)
    #     project_obj.set_model_chains((self.ui.txt_prediction_chain_model.text()))
    #     # gets reference filename and filepath
    #     if len(self.ui.txt_prediction_load_reference.text()) == 4:
    #         tmp_protein = core.Protein(self.ui.txt_prediction_load_reference.text(),
    #                                    export_data_dir=project_obj.get_pdb_path())
    #         tmp_protein.clean_pdb_file()
    #         REFERENCE_OBJ_NAME = self.ui.txt_prediction_load_reference.text()
    #         REFERENCE_DIR = project_obj.get_pdb_path()
    #     else:
    #         ref_file_info = Qt.QtCore.QFileInfo(self.ui.txt_prediction_load_reference.text())
    #         REFERENCE_OBJ_NAME = ref_file_info.baseName()
    #         REFERENCE_DIR = ref_file_info.canonicalPath()
    #     # starting the default web browser to display the colab notebook
    #     self.status_bar.showMessage("Opening Google colab notebook ...")
    #     if self.ui.action_settings_model_w_off_colab_notebook.isChecked():
    #         webbrowser.open_new(constants.OFFICIAL_NOTEBOOK_URL)
    #     else:
    #         webbrowser.open_new(constants.NOTEBOOK_URL)
    #
    #     # # waiting for the colab notebook to finish
    #     # archive = "prediction.zip"
    #     # source_path = global_variables.global_var_settings_obj.get_prediction_path()
    #     # FILE_NAME = f"{source_path}/{archive}"
    #     # # flag = False
    #     # # while flag == False:
    #     # #     print("AlphaFold is still running ...")
    #     # #     time.sleep(5)
    #     # #     # time.sleep(120)
    #     # # TO-DO: loop doesn't work correctly
    #     # while os.path.isfile(FILE_NAME) is False:
    #     #     print("AlphaFold is still running ...")
    #     #     # time.sleep(5)
    #     #     time.sleep(20)
    #     #     # time.sleep(120)
    #     #     # global global_var_abort_prediction
    #     #     # if global_var_abort_prediction:
    #     #     #     return
    #
    #     # alphafold specific process
    #     archive = "prediction.zip"
    #     source_path = global_variables.global_var_settings_obj.get_prediction_path()
    #     filename = f"{source_path}/{archive}"
    #     while os.path.isfile(filename) is False:
    #         print("AlphaFold is still running ...")
    #         time.sleep(20)
    #     # move prediction.zip in scratch folder
    #     shutil.copy(filename, f"{self.scratch_path}/{archive}")
    #     shutil.unpack_archive(f"{self.scratch_path}/{archive}", self.scratch_path, "zip")
    #     shutil.copy(f"{self.scratch_path}/prediction/selected_prediction.pdb",
    #                 f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/selected_prediction.pdb")
    #     shutil.rmtree(f"{self.scratch_path}/prediction")
    #     os.remove(f"{self.scratch_path}/{archive}")
    #
    #     # ----------------------------------------------------------------- #
    #     # start of the analysis algorithm
    #     self.status_bar.showMessage("Protein structure analysis started ...")
    #
    #     # # extracts and moves the prediction.pdb to the workspace/pdb folder
    #     # tools.extract_and_move_model_pdb(
    #     #     str(source_path), f"{str(source_path)}/tmp", archive, project.get_pdb_path())
    #
    #     # gets model filename and filepath
    #     PREDICTION_NAME = tools.get_prediction_file_name(project_obj.get_pdb_path())
    #     full_model_file_path = f"{project_obj.get_pdb_path()}/{PREDICTION_NAME[0]}"
    #     model_file_info = Qt.QtCore.QFileInfo(full_model_file_path)
    #     MODEL_OBJ_NAME = model_file_info.baseName()
    #     MODEL_DIR = model_file_info.canonicalPath()
    #
    #     # set model in project object
    #     project_obj.set_pdb_model(full_model_file_path)
    #     project_obj.create_xml_file()
    #
    #     # create the Protein object for the reference
    #     reference_protein: list[core.Protein] = [core.Protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]
    #
    #     # create model Protein object
    #     model_proteins: list[core.Protein] = [core.Protein(MODEL_OBJ_NAME, MODEL_DIR)]
    #     # sets the filepath of the model in the project xml file
    #     export_dir = project_obj.get_results_path()
    #     structure_analysis_obj = structure_analysis.StructureAnalysis(
    #         reference_protein, model_proteins,
    #         project_obj.get_ref_chains().split(","), project_obj.get_model_chains().split(","),
    #         export_dir, cycles=global_variables.global_var_settings_obj.get_cycles(),
    #         cutoff=global_variables.global_var_settings_obj.get_cutoff(),
    #     )
    #     structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.ref_chains,
    #                                                      structure_analysis_obj.reference_protein)
    #     structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.model_chains,
    #                                                      structure_analysis_obj.model_proteins)
    #     structure_analysis_obj.do_analysis_in_pymol(structure_analysis_obj.create_protein_pairs(),
    #                                             self.status_bar, "2")

    # ----- Functions for Monomer Local Prediction
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
        self._init_local_pred_mono_page()

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

    def predict_local_monomer(self):
        # creating tmp directories in scratch folder to organize prediction inputs and outputs
        # TODO: is there a more elegant way to do it?
        if not os.path.exists(pathlib.Path(f"{self.scratch_path}/local_predictions")):
            os.mkdir(pathlib.Path(f"{self.scratch_path}/local_predictions"))
        if not os.path.exists(constants.PREDICTION_FASTA_DIR):
            os.mkdir(constants.PREDICTION_FASTA_DIR)
        if not os.path.exists(constants.PREDICTION_PDB_DIR):
            os.mkdir(constants.PREDICTION_PDB_DIR)
        # create fasta file and move to scratch fasta dir
        fasta_file = open(f"{constants.PREDICTION_FASTA_DIR}/{self.ui.txt_local_pred_mono_protein_name.text()}.fasta", "w")
        fasta_file.write(f">{self.ui.txt_local_pred_mono_protein_name.text()}\n")
        fasta_file.write(self.ui.txt_local_pred_mono_prot_seq.toPlainText())
        fasta_file.close()
        user_name = os.getlogin()
        fasta_path = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/fasta"
        pdb_path = f"/mnt/c/Users/{user_name}/.pyssa/scratch/local_predictions/pdb"

        try:
            subprocess.run(["C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe", pathlib.Path(
                f"{os.path.expanduser('~')}/github_repos/tmpPySSA/pyssa/scripts/convert_dos_to_unix.ps1")])
            subprocess.run(["wsl", f"/mnt/c/Users/{user_name}/github_repos/tmpPySSA/pyssa/scripts/colabfold_predict.sh",
                            fasta_path, pdb_path])
            subprocess.run(["wsl", "--shutdown"])
        except OSError:
            shutil.rmtree(pathlib.Path(f"{self.scratch_path}/local_predictions"))
            return
        # shutil.unpack_archive(f"{self.scratch_path}/{archive}",
        #                       f"{self.scratch_path}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}", "zip")
        # prediction_results: list[str] = os.listdir(
        #     pathlib.Path(f"{constants.SCRATCH_DIR}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}"))
        # for filename in prediction_results:
        #     check = filename.find("_relaxed_rank_1")
        #     if check != -1:
        #         src = pathlib.Path(f"{constants.SCRATCH_DIR}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}/{filename}")
        #         dest = pathlib.Path(
        #             f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}")
        #         shutil.copy(src, dest)
        #         os.rename(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}",
        #                   f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{self.ui.txt_prediction_only_protein_name.text()}.pdb")
        #         break
        # shutil.rmtree(f"{self.scratch_path}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}")
        # os.remove(f"{self.scratch_path}/{archive}")
        # try:
        #     cmd.load(
        #         f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{self.ui.txt_prediction_only_protein_name.text()}.pdb")
        # except pymol.CmdException:
        #     print("Loading the model failed.")
        #     return
        # self._project_watcher.show_valid_options(self.ui)
        # self.display_view_page()

    # ----- Functions for Multimer Local Prediction
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

    # ----- Functions for Single Analysis
    def show_single_analysis_stage_0(self):
        self.single_analysis_management.show_stage_x(0)
        self.ui.lbl_analysis_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_analysis_prot_struct_2.setText("Protein structure 2")

    def show_single_analysis_stage_1(self):
        self.single_analysis_management.show_stage_x(1)
        self.ui.lbl_analysis_prot_struct_1.setText(self.ui.box_analysis_prot_struct_1.currentText())
        self.ui.lbl_analysis_prot_struct_2.setText(self.ui.box_analysis_prot_struct_2.currentText())
        tools.add_chains_from_pdb_file_to_list(f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}",
                                               self.ui.box_analysis_prot_struct_1.currentText(),
                                               self.ui.list_analysis_ref_chains)
        self.ui.btn_analysis_next_2.setEnabled(True)
        self.ui.lbl_analysis_ref_chains.setText(
            f"Select chains in protein structure {self.ui.lbl_analysis_prot_struct_1.text()}.")

    def show_single_analysis_stage_2(self):
        self.single_analysis_management.show_stage_x(2)
        tools.add_chains_from_pdb_file_to_list(f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}",
                                               self.ui.box_analysis_prot_struct_2.currentText(),
                                               self.ui.list_analysis_model_chains)
        gui_elements_to_show = [
            self.ui.lbl_analysis_ref_chains,
            self.ui.list_analysis_ref_chains,
            self.ui.lbl_analysis_model_chains,
            self.ui.list_analysis_model_chains,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        self.ui.btn_analysis_start.setEnabled(False)
        if self.no_of_selected_chains == 1:
            self.ui.lbl_analysis_model_chains.setText(
                f"Please select {self.no_of_selected_chains} chains in protein structure {self.ui.lbl_analysis_prot_struct_2.text()}.")
            self.ui.list_analysis_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.SingleSelection)
        elif self.no_of_selected_chains > 1:
            self.ui.lbl_analysis_model_chains.setText(
                f"Please select {self.no_of_selected_chains} chains in protein structure {self.ui.lbl_analysis_prot_struct_2.text()}.")
            self.ui.list_analysis_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        else:
            gui_elements_to_hide = [
                self.ui.lbl_analysis_ref_chains,
                self.ui.list_analysis_ref_chains,
                self.ui.lbl_analysis_model_chains,
                self.ui.list_analysis_model_chains,
            ]
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self.ui.btn_analysis_start.setEnabled(True)

    def fill_protein_structure_boxes(self):
        proteins = tools.scan_project_for_valid_proteins(f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}")
        proteins.insert(0, "")
        self.ui.box_analysis_prot_struct_1.clear()
        self.ui.box_analysis_prot_struct_2.clear()
        gui_utils.fill_combo_box(self.ui.box_analysis_prot_struct_1, proteins)
        gui_utils.fill_combo_box(self.ui.box_analysis_prot_struct_2, proteins)

    def count_batch_selected_chains_for_prot_struct_1(self):
        self.no_of_selected_chains = len(self.ui.list_analysis_batch_ref_chains.selectedItems())

    def check_if_prot_structs_are_filled(self):
        prot_1 = self.ui.box_analysis_prot_struct_1.itemText(self.ui.box_analysis_prot_struct_1.currentIndex())
        prot_2 = self.ui.box_analysis_prot_struct_2.itemText(self.ui.box_analysis_prot_struct_2.currentIndex())
        if prot_1 != "" and prot_2 != "":
            self.ui.btn_analysis_next.setEnabled(True)
        else:
            self.ui.btn_analysis_next.setEnabled(False)

    def count_selected_chains_for_prot_struct_1(self):
        self.no_of_selected_chains = len(self.ui.list_analysis_ref_chains.selectedItems())

    def check_if_same_no_of_chains_selected(self):
        self.ui.btn_analysis_start.setEnabled(False)
        if self.no_of_selected_chains == len(self.ui.list_analysis_model_chains.selectedItems()):
            self.ui.btn_analysis_start.setEnabled(True)

    def start_process(self):
        """This function contains the main analysis algorithm for the
        Protein structure comparison.

        """
        self.ui.btn_analysis_start.setEnabled(False)
        self.status_bar.showMessage("Protein structure analysis started ...")
        cmd.reinitialize()
        data_transformer_analysis = data_transformer.DataTransformer(self.ui)
        transformed_analysis_data = data_transformer_analysis.transform_to_analysis(self.app_project)

        if not os.path.exists(transformed_analysis_data[2]):
            os.mkdir(transformed_analysis_data[2])
        else:
            basic_boxes.ok("Single Analysis", "A structure analysis already exists!", QMessageBox.Critical)
            self._init_single_analysis_page()
            return

        structure_analysis_obj = structure_analysis.StructureAnalysis(
            reference_protein=[transformed_analysis_data[0]], model_proteins=[transformed_analysis_data[1]],
            ref_chains=transformed_analysis_data[0].chains, model_chains=transformed_analysis_data[1].chains,
            export_dir=transformed_analysis_data[2], cycles=self.app_settings.get_cycles(),
            cutoff=self.app_settings.get_cutoff(),
        )
        if self.ui.cb_analysis_images.isChecked():
            structure_analysis_obj.response_create_images = True
        structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.ref_chains,
                                                             structure_analysis_obj.reference_protein)
        structure_analysis_obj.create_selection_for_proteins(structure_analysis_obj.model_chains,
                                                             structure_analysis_obj.model_proteins)
        protein_pairs = structure_analysis_obj.create_protein_pairs()
        structure_analysis_obj.do_analysis_in_pymol(protein_pairs, self.status_bar)
        protein_pairs[0].name = transformed_analysis_data[3]
        protein_pairs[0].cutoff = self.app_settings.cutoff
        self.app_project.add_protein_pair(protein_pairs[0])
        protein_pairs[0].serialize_protein_pair(self.app_project.get_objects_protein_pairs_path())
        self.app_project.serialize_project(self.app_project.project_path, "project")
        self.app_project = project.Project.deserialize_project(self.app_project.project_path)
        self._project_watcher.show_valid_options(self.ui)
        self._init_single_analysis_page()

    # ----- Functions for Batch
    def show_batch_analysis_stage_0(self):
        if self.ui.lbl_analysis_batch_prot_struct_1.text() != "Protein structure 1":
            prot_1_name = self.ui.lbl_analysis_batch_prot_struct_1.text().replace(".pdb", "")
            prot_1_chains = []
            for chain in self.ui.list_analysis_batch_ref_chains.selectedItems():
                prot_1_chains.append(chain.text())
            prot_1_chains = ','.join([str(elem) for elem in prot_1_chains])
            prot_2_name = self.ui.lbl_analysis_batch_prot_struct_2.text().replace(".pdb", "")
            prot_2_chains = []
            for chain in self.ui.list_analysis_batch_model_chains.selectedItems():
                prot_2_chains.append(chain.text())
            prot_2_chains = ','.join([str(elem) for elem in prot_2_chains])
            analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
            item = QListWidgetItem(analysis_name)
            self.ui.list_analysis_batch_overview.addItem(item)
        if self.ui.list_analysis_batch_overview.count() == 0:
            self.batch_analysis_management.show_stage_x(0)
        else:
            gui_elements_to_show = [
                self.ui.btn_analysis_batch_add,
                self.ui.btn_analysis_batch_remove,
            ]
            self.batch_analysis_management.show_gui_elements_stage_x(
                [0, 4], [1, 2, 3], show_specific_elements=gui_elements_to_show
            )

    def show_batch_analysis_stage_1(self):
        gui_elements_to_hide = [
            self.ui.btn_analysis_batch_add,
            self.ui.btn_analysis_batch_remove,
        ]
        self.batch_analysis_management.show_gui_elements_stage_x(
            [0, 1], [2, 3, 4], hide_specific_elements=gui_elements_to_hide)
        self.fill_protein_boxes_batch()
        self.ui.lbl_analysis_batch_prot_struct_1.setText("Protein structure 1")
        self.ui.lbl_analysis_batch_prot_struct_2.setText("Protein structure 2")

    def show_batch_analysis_stage_2(self):
        self.batch_analysis_management.show_gui_elements_stage_x(
            [0, 1, 2], [3, 4], hide_specific_elements=[self.ui.box_analysis_batch_prot_struct_1,
                                                       self.ui.box_analysis_batch_prot_struct_2,
                                                       self.ui.btn_analysis_batch_next,
                                                       self.ui.btn_analysis_batch_back]
        )
        self.ui.lbl_analysis_batch_prot_struct_1.setText(self.ui.box_analysis_batch_prot_struct_1.currentText())
        self.ui.lbl_analysis_batch_prot_struct_2.setText(self.ui.box_analysis_batch_prot_struct_2.currentText())
        tools.add_chains_from_pdb_file_to_list(f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}",
                                               self.ui.box_analysis_batch_prot_struct_1.currentText(),
                                               self.ui.list_analysis_batch_ref_chains)
        self.ui.lbl_analysis_batch_ref_chains.setText(
            f"Select chains in protein structure {self.ui.lbl_analysis_batch_prot_struct_1.text()}.")

    def show_batch_analysis_stage_3(self):
        self.ui.btn_analysis_batch_next_3.setEnabled(False)
        tools.add_chains_from_pdb_file_to_list(f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}",
                                               self.ui.box_analysis_batch_prot_struct_2.currentText(),
                                               self.ui.list_analysis_batch_model_chains)
        gui_elements_to_hide = [
            self.ui.box_analysis_batch_prot_struct_1,
            self.ui.box_analysis_batch_prot_struct_2,
            self.ui.btn_analysis_batch_next,
            self.ui.btn_analysis_batch_back,
            self.ui.btn_analysis_batch_next_2,
            self.ui.btn_analysis_batch_back_2,
        ]
        if self.no_of_selected_chains == 1:
            # only one chain was selected
            self.ui.lbl_analysis_batch_model_chains.setText(
                f"Please select {self.no_of_selected_chains} chains in protein structure {self.ui.lbl_analysis_batch_prot_struct_2.text()}.")
            self.ui.list_analysis_batch_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.SingleSelection)
            self.batch_analysis_management.show_gui_elements_stage_x(
                [0, 1, 2, 3], [4], hide_specific_elements=gui_elements_to_hide
            )
        elif self.no_of_selected_chains > 1:
            # multiple chains were selected
            self.ui.lbl_analysis_batch_model_chains.setText(
                f"Please select {self.no_of_selected_chains} chains in protein structure {self.ui.lbl_analysis_batch_prot_struct_2.text()}.")
            self.ui.list_analysis_batch_model_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
            self.batch_analysis_management.show_gui_elements_stage_x(
                [0, 1, 2, 3], [4], hide_specific_elements=gui_elements_to_hide
            )
        else:
            # no chains were selected
            gui_elements_to_show = [
                self.ui.btn_analysis_batch_add,
                self.ui.btn_analysis_batch_remove,
            ]
            self.batch_analysis_management.show_gui_elements_stage_x(
                [0, 4], [1, 2, 3], show_specific_elements=gui_elements_to_show
            )
            prot_1_name = self.ui.lbl_analysis_batch_prot_struct_1.text().replace(".pdb", "")
            prot_2_name = self.ui.lbl_analysis_batch_prot_struct_2.text().replace(".pdb", "")
            analysis_name = f"{prot_1_name}_vs_{prot_2_name}"
            item = QListWidgetItem(analysis_name)
            self.ui.list_analysis_batch_overview.addItem(item)

    def fill_protein_boxes_batch(self):
        proteins = tools.scan_project_for_valid_proteins(
            f"{self.workspace_path}\\{self.ui.lbl_current_project_name.text()}")
        proteins.insert(0, "")
        self.ui.box_analysis_batch_prot_struct_1.clear()
        self.ui.box_analysis_batch_prot_struct_2.clear()
        gui_utils.fill_combo_box(self.ui.box_analysis_batch_prot_struct_1, proteins)
        gui_utils.fill_combo_box(self.ui.box_analysis_batch_prot_struct_2, proteins)

    def remove_analysis_run(self):
        self.ui.list_analysis_batch_overview.takeItem(self.ui.list_analysis_batch_overview.currentRow())
        if self.ui.list_analysis_batch_overview.count() == 0:
            self.batch_analysis_management.show_stage_x(0)

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

    def start_process_batch(self):
        """This function contains the main analysis algorithm for the
        Protein structure comparison.

        """
        batch_analysis = []
        for row_no in range(self.ui.list_analysis_batch_overview.count()):
            tmp_batch_analysis = self.ui.list_analysis_batch_overview.item(row_no).text()
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

            tmp_prot_1 = protein_analysis_info.ProteinAnalysisInfo(prot_1_name, prot_1_chains, tmp_batch_analysis)
            tmp_prot_2 = protein_analysis_info.ProteinAnalysisInfo(prot_2_name, prot_2_chains, tmp_batch_analysis)
            batch_analysis.append((tmp_prot_1, tmp_prot_2))

        transformer = data_transformer.DataTransformer(self.ui)
        # contains analysis-ready data format: list(tuple(prot_1, prot_2, export_dir), ...)
        batch_data = transformer.transform_data_for_analysis(self.app_project, batch_analysis)

        for analysis_data in batch_data:
            if not os.path.exists(analysis_data[2]):
                os.mkdir(analysis_data[2])
            else:
                basic_boxes.ok("Single Analysis", f"The structure analysis: {analysis_data[3]} already exists!", QMessageBox.Critical)
                self._init_batch_analysis_page()
                return

            cmd.reinitialize()
            structure_analysis_obj = structure_analysis.StructureAnalysis(
                reference_protein=[analysis_data[0]], model_proteins=[analysis_data[1]],
                ref_chains=analysis_data[0].chains, model_chains=analysis_data[1].chains,
                export_dir=analysis_data[2], cycles=self.app_settings.get_cycles(),
                cutoff=self.app_settings.get_cutoff(),
            )
            if self.ui.cb_analysis_images.isChecked():
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
            self._project_watcher.show_valid_options(self.ui)
            self._init_batch_analysis_page()

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

    # ----- Functions for Results
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

    def load_results(self):
        self.results_name = self.ui.cb_results_analysis_options.currentText()
        if self.results_name == "":
            self.show_analysis_results_options()
            return
        current_results_path = pathlib.Path(f"{self.app_project.get_results_path()}/{self.results_name}")
        gui_elements_to_hide = []
        if not os.path.exists(pathlib.Path(f"{current_results_path}/images")):
            # no images where made
            gui_elements_to_hide.append(self.ui.lbl_results_structure_alignment)
            gui_elements_to_hide.append(self.ui.btn_view_struct_alignment)
            gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
            gui_elements_to_hide.append(self.ui.list_results_interest_regions)
            gui_elements_to_hide.append(self.ui.btn_view_interesting_region)
        elif os.path.exists(pathlib.Path(f"{current_results_path}/images")):
            if not os.path.exists(pathlib.Path(f"{current_results_path}/images/structure_alignment.png")):
                gui_elements_to_hide.append(self.ui.lbl_results_structure_alignment)
                gui_elements_to_hide.append(self.ui.btn_view_struct_alignment)
            elif not os.path.exists(pathlib.Path(f"{current_results_path}/images/interesting_")):
                gui_elements_to_hide.append(self.ui.lbl_results_interest_regions)
                gui_elements_to_hide.append(self.ui.list_results_interest_regions)
                gui_elements_to_hide.append(self.ui.btn_view_interesting_region)

        # check if histogram can be created
        # read csv file
        file_path = pathlib.Path(
            f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/{self.results_name}")
        path = f"{file_path}/distance_csv/distances.csv"
        distance_list = []
        with open(path, 'r', encoding="utf-8") as csv_file:
            i = 0
            for line in csv_file:
                cleaned_line = line.replace("\n", "")
                if cleaned_line.split(",")[8] != 'distance':
                    distance_list.append(float(cleaned_line.split(",")[8]))
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

        if gui_elements_to_hide:
            self.show_results_interactions(gui_elements_to_hide=gui_elements_to_hide)
        else:
            self.show_results_interactions()

        try:
            rmsd_file = open(pathlib.Path(f"{current_results_path}/rmsd.json"), "r", encoding="utf-8")
        except FileNotFoundError:
            print(f"There is no valid protein pair json file under: {pathlib.Path(f'{current_results_path}')}")
            return
        rmsd_dict = json.load(rmsd_file)
        self.ui.txt_results_rmsd.setText(str(rmsd_dict.get("rmsd")))
        self.ui.txt_results_aligned_residues.setText(str(rmsd_dict.get("aligned_residues")))
        session_filepath = pathlib.Path(f"{current_results_path}/sessions/session_file_model_s.pse")
        cmd.reinitialize()
        cmd.load(str(session_filepath))

    def change_interesting_regions(self):
        """This function is used to switch between projects within a job.

        """
        global global_var_list_widget_row
        global global_var_project_dict

        if gui_utils.warning_switch_pymol_session("") is True:
            try:
                file_path = global_var_project_dict[global_var_list_widget_row].get_results_path()
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
        results_path = global_var_project_dict[current_row].get_results_path()
        dir_content = os.listdir(f"{results_path}/images/interesting_regions")
        self.ui.cb_interesting_regions.clear()
        for tmp_file in dir_content:
            self.ui.cb_interesting_regions.addItem(tmp_file)
        cmd.load(global_var_project_dict[current_row].get_session_file())

    def display_structure_alignment(self):
        """This function opens a window which displays the image of the structure alignment.

        """
        png_dialog = Qt.QtWidgets.QDialog(self)
        label = Qt.QtWidgets.QLabel(self)
        file_path = pathlib.Path(
            f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/{self.results_name}")
        pixmap = Qt.QtGui.QPixmap(f"{file_path}/images/structure_alignment.png")
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
        protein_pair_of_analysis = self.app_project.get_specific_protein_pair(self.ui.cb_results_analysis_options.currentText())
        dialog = dialog_distance_plot.DialogDistancePlot(protein_pair_of_analysis)
        dialog.exec_()

        # item = self.ui.project_list.selectedItems()
        # if item is None:
        #     raise ValueError

        # # plot_dialog = Qt.QtWidgets.QDialog(self)
        # plot_dialog_layout = QHBoxLayout()
        # graph_widget = pg.PlotWidget()
        #
        # # read csv file
        # file_path = f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results"
        # model_name = self.ui.lbl_current_project_name.text()
        #
        # path = f"{file_path}/distance_csv/distances.csv"
        # distance_list = []
        # cutoff_line = []
        # with open(path, 'r', encoding="utf-8") as csv_file:
        #     for line in csv_file:
        #         cleaned_line = line.replace("\n", "")
        #         if cleaned_line.split(",")[8] != 'distance':
        #             distance_list.append(float(cleaned_line.split(",")[8]))
        #             cutoff_line.append(global_variables.global_var_settings_obj.get_cutoff())
        # # creates actual distance plot line
        # graph_widget.plotItem.plot(distance_list, pen=pg.mkPen(color="#4B91F7", width=6),
        #                            symbol="o", symbolSize=10, symbolBrush=('b'))
        # self.view_box = graph_widget.plotItem.getViewBox()
        # self.view_box.setRange(xRange=[23, 40])
        # # creates cutoff line
        # graph_widget.plotItem.plot(cutoff_line, pen=pg.mkPen(color="#f83021", width=6))
        # # styling the plot
        # graph_widget.setBackground('w')
        # graph_widget.setTitle(f"Distance Plot of {model_name}", size="23pt")
        # styles = {'font-size': '14px'}
        # ax_label_y = "Distance in Å"
        # graph_widget.setLabel('left', ax_label_y, **styles)
        # graph_widget.setLabel('bottom', "Residue pair no.", **styles)
        # graph_widget.plotItem.showGrid(x=True, y=True)
        # plot_dialog_layout.addWidget(graph_widget)
        # self.plot_dialog.setLayout(plot_dialog_layout)
        # self.plot_dialog.setWindowTitle("Distance Plot")
        # self.plot_dialog.show()

    def display_distance_histogram(self):
        """This function opens a window which displays the distance histogram.

        """
        # item = self.ui.project_list.selectedItems()
        # if item is None:
        #     raise ValueError
        protein_pair_of_analysis = self.app_project.get_specific_protein_pair(
            self.ui.cb_results_analysis_options.currentText())
        dialog = dialog_distance_histogram.DialogDistanceHistogram(protein_pair_of_analysis)
        dialog.exec_()

        # plot_dialog = Qt.QtWidgets.QDialog(self)
        # plot_dialog_layout = QHBoxLayout()
        # graph_widget = pg.PlotWidget()
        # view_box = graph_widget.plotItem.getViewBox()
        # # read csv file
        # path = pathlib.Path(f"{protein_pair_of_analysis.results_dir}/distance_csv/distances.csv")
        # distance_list = []
        # with open(path, 'r', encoding="utf-8") as csv_file:
        #     i = 0
        #     for line in csv_file:
        #         cleaned_line = line.replace("\n", "")
        #         if cleaned_line.split(",")[8] != 'distance':
        #             distance_list.append(float(cleaned_line.split(",")[8]))
        # distance_list.sort()
        # length = len(distance_list)
        # max_distance = distance_list[length-1]
        # x, y = np.histogram(distance_list, bins=np.arange(0, max_distance, 0.25))
        #
        # if x.size != y.size:
        #     x = np.resize(x, (1, y.size))
        # # this conversion is needed for the pyqtgraph library!
        # x = x.tolist()
        # try:
        #     x = x[0]
        # except IndexError:
        #     # Error got raised where the distances where all 0
        #     tools.quick_log_and_display("error", "The histogram could not be created.",
        #                                 self.status_bar, "The histogram could not be created. "
        #                                                  " Check the distance table!")
        #     return
        #
        # y = y.tolist()
        # color = Qt.QtGui.QColor.fromRgb(255, 128, 128)
        # # creates bar chart item
        # graph_bar_item = pg.BarGraphItem(x0=0, y=y, height=0.2, width=x,
        #                                  pen=pg.mkPen(color="#4B91F7"), brush=pg.mkBrush(color="#4B91F7"))
        # # creates y-labels for bar chart
        # y_labels = []
        # for i in range(len(y)):
        #     try:
        #         label = f"{y[i]} - {y[i+1]}"
        #     except IndexError:
        #         # detects if a last label is necessary
        #         label = f"{y[i]} - {y[i]+0.25}"
        #     y_labels.append(label)
        # y_values = y
        # ticks = []
        # for i, item in enumerate(y_labels):
        #     ticks.append((y_values[i], item))
        # ticks = [ticks]
        #
        # # styling the plot
        # graph_widget.setBackground('w')
        # graph_widget.setTitle(f"Distance Histogram of {protein_pair_of_analysis.name}", size="23pt")
        # styles = {'font-size': '14px'}
        # ax_label_x = "Distance in Å"
        # graph_widget.setLabel('left', ax_label_x, **styles)
        # graph_widget.setLabel('bottom', "Frequency of α-C atoms distance", **styles)
        # graph_widget.addItem(graph_bar_item)
        # bar_ax = graph_widget.getAxis('left')
        # bar_ax.setTicks(ticks)
        # view_box.invertY(True)
        # text = pg.TextItem(text='', border='w', fill=(0, 0, 255, 100), anchor=(0, 0))
        # graph_widget.addItem(text)
        # plot_dialog_layout.addWidget(graph_widget)
        # plot_dialog.setLayout(plot_dialog_layout)
        # plot_dialog.setWindowTitle("Distance Histogram")
        # plot_dialog.show()

    def display_interesting_region(self):
        """This function displays an image of an interesting region.

        """
        png_dialog = Qt.QtWidgets.QDialog(self)
        label = Qt.QtWidgets.QLabel(self)
        global global_var_project_dict
        file_path = pathlib.Path(
            f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/{self.results_name}")
        #file_name = self.ui.cb_interesting_regions.currentText()
        file_name = ""
        pixmap = Qt.QtGui.QPixmap(f"{file_path}/images/interesting_regions/{file_name}")
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
        csv_model = Qt.QtGui.QStandardItemModel()
        csv_model.setColumnCount(7)
        labels = ["Residue pair no.", "Reference Chain", "Reference Position", "Reference Residue",
                  "Model Chain", "Model Position", "Model Residue", "Distance",
                  ]
        csv_model.setHorizontalHeaderLabels(labels)
        table_dialog = Qt.QtWidgets.QDialog(self)
        table_view = Qt.QtWidgets.QTableView()
        table_view.setModel(csv_model)

        file_path = pathlib.Path(
            f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/{self.results_name}")
        path = f"{file_path}/distance_csv/distances.csv"

        with open(path, 'r', encoding="utf-8") as csv_file:
            i = 0
            for line in csv_file:
                tmp_list = line.split(",")
                tmp_list.pop(0)
                standard_item_list = []
                if i == 0:
                    for tmp in tmp_list:
                        tmp_item = Qt.QtGui.QStandardItem(tmp)
                        standard_item_list.append(tmp_item)
                else:
                    pair_no_item = Qt.QtGui.QStandardItem()
                    pair_no_item.setData(int(tmp_list[0]), role=QtCore.Qt.DisplayRole)
                    ref_chain_item = Qt.QtGui.QStandardItem()
                    ref_chain_item.setData(str(tmp_list[1]), role=QtCore.Qt.DisplayRole)
                    ref_pos_item = Qt.QtGui.QStandardItem()
                    ref_pos_item.setData(int(tmp_list[2]), role=QtCore.Qt.DisplayRole)
                    ref_resi_item = Qt.QtGui.QStandardItem()
                    ref_resi_item.setData(str(tmp_list[3]), role=QtCore.Qt.DisplayRole)
                    model_chain_item = Qt.QtGui.QStandardItem()
                    model_chain_item.setData(str(tmp_list[4]), role=QtCore.Qt.DisplayRole)
                    model_pos_item = Qt.QtGui.QStandardItem()
                    model_pos_item.setData(int(tmp_list[5]), role=QtCore.Qt.DisplayRole)
                    model_resi_item = Qt.QtGui.QStandardItem()
                    model_resi_item.setData(str(tmp_list[6]), role=QtCore.Qt.DisplayRole)
                    distance_item = Qt.QtGui.QStandardItem()
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

    # ----- Functions for Image page
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
            print("Please select a renderer.")
            self.status_bar.showMessage("Please select a renderer.")
        elif self.ui.box_renderer.currentIndex() == 1:
            self.renderer = "-1"
        elif self.ui.box_renderer.currentIndex() == 2:
            self.renderer = "0"
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
        scene_name = Qt.QtWidgets.QInputDialog.getText(self, "Save Scene", "Enter scene name:")
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
            save_dialog = Qt.QtWidgets.QFileDialog()
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
            save_dialog = Qt.QtWidgets.QFileDialog()
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

    # ----- Functions for Hotspots page
    def open_protein(self):
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find(".pdb") != -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.load_protein()
            tmp_sequence = cmd.get_fastastr('all')
            self.ui.sp_hotspots_resi_no.setMaximum(len(tmp_sequence))
        else:
            # protein pair is selected
            tmp_protein_pair = self.app_project.get_specific_protein_pair(input)
            tmp_protein_pair.load_protein_pair()

    def show_resi_sticks(self):
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find(".pdb") != -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.selection = f"/{tmp_protein.molecule_object}///{self.ui.sp_hotspots_resi_no.text()}/"
            tmp_protein.show_resi_as_balls_and_sticks()
        # protein pair is selected

    def hide_resi_sticks(self):
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find(".pdb") != -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.selection = f"/{tmp_protein.molecule_object}///{self.ui.sp_hotspots_resi_no.text()}/"
            tmp_protein.hide_resi_as_balls_and_sticks()

    def zoom_resi_position(self):
        input = self.ui.list_hotspots_choose_protein.currentItem().text()
        if input.find(".pdb") != -1:
            # one protein is selected
            tmp_protein = self.app_project.search_protein(input.replace(".pdb", ""))
            tmp_protein.selection = f"/{tmp_protein.molecule_object}///{self.ui.sp_hotspots_resi_no.text()}/"
            tmp_protein.zoom_resi_protein_position()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    styles.set_stylesheet(app)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())
