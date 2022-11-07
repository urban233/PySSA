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
import sys
import time
import webbrowser
from pathlib import Path
from xml.dom import minidom

import PyQt5.QtCore
import numpy as np
import pymol
import pyqtgraph as pg
from PyQt5 import QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtWidgets import QHBoxLayout
from pymol import Qt
from pymol import cmd
from urllib.request import urlopen
from urllib.error import URLError

import utils.project_utils
from alternatives.auto_main_window_alt import Ui_MainWindow
from dialogs import dialog_about
from dialogs import dialog_add_models
from dialogs import dialog_add_model
from dialogs import dialog_startup
from dialogs import dialog_distance_plot
from pymolproteintools import core
from utils import gui_utils
from utils import job_utils
from utils import project_constants
from utils import structure_analysis_utils
from utils import tools
from utils import global_utils
from utils import styles_utils
from utils import project_utils

# setup logger
logging.basicConfig(level=logging.DEBUG)
# global variables
global_var_project_dict = {0: utils.project_utils.Project("", "")}
global_var_list_widget_row = 0


class MainWindow(QMainWindow):
    """This class contains all information about the new design MainWindow in the
    application

    """

    def __init__(self):
        super().__init__()
        # build ui object
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.lbl_page_title.setText("Home")
        self.setMinimumWidth(550)

        # checks if the plugin launched for the first time
        if global_utils.global_var_settings_obj.get_app_launch() == 0:
            dialog = dialog_startup.DialogStartup()
            dialog.exec_()
            # checks if the cancel button was pressed
            if dialog_startup.global_var_terminate_app == 1:
                global_utils.global_var_settings_obj.delete_settings_xml()
                sys.exit()
            global_utils.global_var_settings_obj.set_workspace_path(dialog_startup.global_var_startup_workspace)
            # sets the launch_app var to 1 to indicate that the plugin has successfully overcome the first launch
            global_utils.global_var_settings_obj.set_app_launch(1)
            global_utils.global_var_settings_obj.save_settings_to_xml()

        self.plot_dialog = Qt.QtWidgets.QDialog(self)
        self.view_box = None
        # sets up the status bar
        self._setup_statusbar()
        # # sets up settings.xml
        # if not os.path.exists(self.SETTINGS):
        #     settings = utils.settings_utils.SettingsXml(self.SETTINGS)
        #     settings.create_settings_xml_file()
        # settings = utils.settings_utils.SettingsXml(self.SETTINGS)
        # self.tmp_settings = settings.load_xml_in_memory()
        #
        # # safeguard workspace path
        # sg_1 = tools.safeguard_filepath_xml(self.tmp_settings, 'workspacePath', 'value')
        # # safeguard pdb path
        # sg_2 = tools.safeguard_filepath_xml(self.tmp_settings, 'pdbPath', 'value')
        # # safeguard zip path
        # sg_3 = tools.safeguard_filepath_xml(self.tmp_settings, 'zipPath', 'value')
        # # safeguard cycles value
        # sg_4 = tools.safeguard_numerical_value_xml(self.tmp_settings, 'cyclesValue', 'value', 'int')
        # # safeguard cutoff value
        # sg_5 = tools.safeguard_numerical_value_xml(self.tmp_settings, 'cutoffValue', 'value', 'float')
        #
        # if sg_1 is False or sg_2 is False or sg_3 is False or sg_4 is False or sg_5 is False:
        #     self.status_bar.showMessage("The settings.xml is corrupted! Please fix this issue first!")
        #     gui_utils.error_dialog_settings("The settings.xml is corrupted! Please fix this issue first!",
        #                                     "Check the log for more info.")
        # else:
        #     self.workspace_path = utils.settings_utils.SettingsXml.get_path(self.tmp_settings,
        #                                                                     "workspacePath",
        #                                                                     "value")

        # startup_dialog = dialog_startup.DialogStartup()
        # startup_dialog.exec_()

        self.scratch_path = f"{project_constants.SETTINGS_DIR}/scratch"
        tools.create_directory(project_constants.SETTINGS_DIR, "scratch")

        self.ui.action_add_model.setVisible(False)
        self.ui.action_add_multiple_models.setVisible(False)
        self.ui.action_file_save_as.setVisible(False)

        self.ui.btn_s_v_p_start.setEnabled(False)
        self.ui.btn_prediction_only_start.setEnabled(False)
        self.ui.lbl_current_project_name.setText("")
        self.ui.btn_new_create_project.setEnabled(False)
        self.ui.cb_new_add_reference.setCheckable(False)
        self.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
        self.ui.txt_open_search.textChanged.connect(self.validate_open_search)
        self.ui.lbl_open_status_search.setText("")
        self.ui.txt_delete_search.textChanged.connect(self.validate_delete_search)
        self.ui.lbl_delete_status_search.setText("")
        self.ui.btn_open_open_project.setEnabled(False)
        self.ui.txt_open_selected_project.textChanged.connect(self.activate_open_button)
        self.ui.list_open_projects.currentItemChanged.connect(self.select_project_from_open_list)
        self.ui.btn_delete_delete_project.setEnabled(False)
        self.ui.txt_delete_selected_projects.textChanged.connect(self.activate_delete_button)
        self.ui.list_delete_projects.currentItemChanged.connect(self.select_project_from_delete_list)
        self.ui.list_s_v_p_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)

        self.ui.list_new_seq_notebooks.addItem(project_constants.OFFICIAL_NOTEBOOK_NAME)
        self.ui.list_new_seq_notebooks.currentItemChanged.connect(self.enable_predict_button)

        self.ui.btn_new_create_project.clicked.connect(self.create_new_project)
        self.ui.btn_delete_delete_project.clicked.connect(self.delete_project)
        self.ui.btn_open_open_project.clicked.connect(self.open_project)
        self.ui.list_open_projects.doubleClicked.connect(self.open_project)
        self.ui.btn_save_project.clicked.connect(self.save_project)
        self.ui.action_add_multiple_models.triggered.connect(self.open_add_models)
        self.ui.action_add_model.triggered.connect(self.open_add_model)
        self.ui.btn_analysis_next.clicked.connect(self.analysis_next_step)
        self.ui.btn_analysis_back.clicked.connect(self.analysis_back_step)

        self.ui.lbl_analysis_model_chains.hide()
        self.ui.list_analysis_model_chains.hide()
        self.ui.btn_analysis_back.hide()
        self.ui.btn_analysis_start.hide()

        # setup defaults for pages
        self._init_hide_ui_elements()
        self._init_fill_combo_boxes()
        self._init_new_page()
        self._init_new_sequence_page()
        self._init_sequence_vs_pdb_page()
        self._init_single_analysis_page()
        self._init_batch_page()

        # connections
        self._connect_sidebar_buttons()
        self._connect_menu_entries()
        self._connect_new_buttons()
        self._connect_new_sequence_buttons()
        self._connect_sequence_vs_pdb_buttons()
        self._connect_single_analysis_buttons()
        self._connect_batch_buttons()
        self._connect_results_buttons()
        self._connect_image_buttons()

        self._connect_new_text_input()
        self._connect_batch_text_fields()

        self._connect_image_combo_boxes()

        self._connect_new_checkbox()
        self._connect_batch_checkbox()
        self._connect_image_checkbox()

        # create tooltips
        self._create_tooltips_sequence_vs_pdb()
        self._create_tooltips_single_analysis()
        self._create_tooltips_batch()
        self._create_tooltips_image()

        # setting additional parameters
        self.setWindowTitle("PySSA v0.9.2")

    def _setup_statusbar(self):
        """This function sets up the status bar and fills it with the current workspace

        """
        self.status_bar = Qt.QtWidgets.QStatusBar()
        self.setStatusBar(self.status_bar)
        self.workspace_path = global_utils.global_var_settings_obj.get_workspace_path()
        self.workspace = Qt.QtWidgets.QLabel(f"Current Workspace: {self.workspace_path}")
        self.status_bar.showMessage(self.workspace.text())

    def _init_hide_ui_elements(self):
        """This function hides all UI elements which need to be hidden during the
        plugin startup

        """
        self.ui.lbl_new_choose_reference.hide()
        self.ui.txt_new_choose_reference.hide()
        self.ui.btn_new_choose_reference.hide()
        self.ui.lbl_prediction_load_reference.hide()
        self.ui.txt_prediction_load_reference.hide()
        self.ui.btn_prediction_load_reference.hide()
        self.ui.btn_prediction_next_2.hide()
        self.ui.btn_prediction_back_2.hide()
        # self.ui.lbl_prediction_ref_chains.hide()
        self.ui.lbl_prediction_model_chains.hide()
        self.ui.txt_prediction_chain_model.hide()
        self.ui.btn_prediction_back_3.hide()
        # self.ui.btn_prediction_start.hide()
        self.ui.btn_hotspots_page.hide()
        # sidebar elements
        self.ui.lbl_prediction.hide()
        self.ui.lbl_analysis.hide()
        self.ui.lbl_handle_pymol_session.hide()
        self.ui.btn_save_project.hide()
        self.ui.btn_prediction_only_page.hide()
        self.ui.btn_prediction_page.hide()
        self.ui.btn_single_analysis_page.hide()
        self.ui.btn_job_analysis_page.hide()
        self.ui.btn_results_page.hide()
        self.ui.btn_image_page.hide()

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
        self.ui.lbl_new_status_project_name.setText("")
        self.ui.lbl_new_status_choose_reference.setText("")

    def _init_new_sequence_page(self):
        # sets up defaults: Prediction
        self.ui.btn_prediction_only_start.setEnabled(False)

    def _init_sequence_vs_pdb_page(self):
        # sets up defaults: Prediction + Analysis
        self.ui.label_18.hide()
        self.ui.txt_prediction_project_name.hide()
        self.ui.lbl_prediction_status_project_name.hide()
        self.ui.list_widget_projects.hide()
        self.ui.btn_prediction_next_1.hide()

        self.ui.lbl_prediction_ref_chains.show()
        self.ui.list_widget_ref_chains.show()
        self.ui.list_widget_ref_chains.setStyleSheet("border-style: solid;"
                                                     "border-width: 2px;"
                                                     "border-radius: 8px;"
                                                     "border-color: #DCDBE3;")
        self.ui.btn_prediction_start.show()
        self.ui.btn_prediction_start.setEnabled(False)
        self.ui.btn_prediction_next_1.setEnabled(False)
        self.ui.list_widget_ref_chains.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ui.lbl_prediction_status_project_name.setText("")
        self.ui.lbl_prediction_status_load_reference.setText("")

    def _init_single_analysis_page(self):
        # sets up defaults: Single Analysis
        self.ui.btn_analysis_start.setEnabled(False)

    def _init_batch_page(self):
        # sets up defaults: Batch
        self.ui.txt_batch_chain_ref.setEnabled(False)
        self.ui.txt_batch_chain_model.setEnabled(False)
        self.ui.btn_batch_start.setEnabled(False)
        self.ui.progress_bar_batch.setProperty("value", 0)

    def _connect_sidebar_buttons(self):
        """This function connects all buttons on the sidebar
        to their slots

        """
        self.ui.btn_new_page.clicked.connect(self.display_new_page)
        self.ui.btn_open_page.clicked.connect(self.display_open_page)
        self.ui.btn_delete_page.clicked.connect(self.display_delete_page)
        self.ui.btn_prediction_only_page.clicked.connect(self.display_prediction_only_page)
        self.ui.btn_prediction_page.clicked.connect(self.display_sequence_vs_pdb_page)
        self.ui.btn_single_analysis_page.clicked.connect(self.display_single_analysis_page)
        self.ui.btn_job_analysis_page.clicked.connect(self.display_job_analysis_page)
        self.ui.btn_results_page.clicked.connect(self.display_results_page)
        self.ui.btn_image_page.clicked.connect(self.display_image_page)

    def _connect_menu_entries(self):
        """This function connects all menu entries with their
        slots

        """
        self.ui.action_file_open.triggered.connect(self.open)
        self.ui.action_file_save_as.triggered.connect(self.save_as)
        self.ui.action_file_save.triggered.connect(self.save)
        self.ui.action_file_quit.triggered.connect(self.quit_app)
        self.ui.action_file_restore_settings.triggered.connect(self.restore_settings)
        self.ui.action_settings_edit_all.triggered.connect(self.open_settings_global)
        self.ui.action_help_docs.triggered.connect(self.open_documentation)
        self.ui.action_help_docs_pdf.triggered.connect(self.open_documentation_pdf)
        self.ui.action_help_about.triggered.connect(self.open_about)

    def _connect_new_buttons(self):
        self.ui.btn_new_choose_reference.clicked.connect(self.load_reference_in_project)

    def _connect_new_sequence_buttons(self):
        """This function connects all buttons on the new
        sequence page with their slots

        """
        self.ui.btn_prediction_only_start.clicked.connect(self.predict_only)

    def _connect_sequence_vs_pdb_buttons(self):
        """This function connects all buttons on the sequence
         vs pdb page with their slots

        """
        # self.ui.txt_prediction_project_name.textChanged.connect(self.validate_project_name)
        # self.ui.txt_prediction_project_name.returnPressed.connect(self.show_prediction_load_reference)
        # self.ui.btn_prediction_next_1.clicked.connect(self.show_prediction_load_reference)
        # self.ui.btn_prediction_back_2.clicked.connect(self.hide_prediction_load_reference)
        # self.ui.btn_prediction_load_reference.clicked.connect(self.load_reference_for_prediction)
        # self.ui.txt_prediction_load_reference.textChanged.connect(self.validate_reference_for_prediction)
        # self.ui.txt_prediction_load_reference.returnPressed.connect(self.show_prediction_chain_info_reference)
        # self.ui.btn_prediction_next_2.clicked.connect(self.show_prediction_chain_info_reference)
        # self.ui.btn_prediction_back_3.clicked.connect(self.hide_prediction_chain_info_reference)
        self.ui.btn_s_v_p_start.clicked.connect(self.predict)

    def _connect_single_analysis_buttons(self):
        self.ui.btn_analysis_start.clicked.connect(self.start_process)

    def _connect_batch_buttons(self):
        self.ui.btn_batch_load_reference.clicked.connect(self.load_reference_for_batch)
        self.ui.btn_batch_load_model.clicked.connect(self.load_model_for_batch)
        self.ui.btn_batch_start.clicked.connect(self.start_process_batch)

    def _connect_results_buttons(self):
        self.ui.btn_view_struct_alignment.clicked.connect(self.display_structure_alignment)
        self.ui.btn_view_distance_plot.clicked.connect(self.display_distance_plot)
        self.ui.btn_view_distance_histogram.clicked.connect(self.display_distance_histogram)
        self.ui.btn_view_distance_table.clicked.connect(self.display_distance_table)
        #self.ui.btn_view_interesting_region.clicked.connect(self.tmp_dialog_change)

    def _connect_image_buttons(self):
        self.ui.btn_update_scene.clicked.connect(self.update_scene)
        self.ui.btn_save_scene.clicked.connect(self.save_scene)
        self.ui.btn_save_image.clicked.connect(self.save_image)
        self.ui.btn_preview_image.clicked.connect(self.preview_image)

    def _connect_new_text_input(self):
        self.ui.txt_new_project_name.textChanged.connect(self.validate_project_name)
        self.ui.txt_new_choose_reference.textChanged.connect(self.validate_reference_in_project)

    def _connect_batch_text_fields(self):
        self.ui.txt_batch_job_name.textChanged.connect(
            self.check_batch_if_txt_batch_job_name_is_filled)
        self.ui.txt_batch_chain_ref.textChanged.connect(
            self.check_batch_if_txt_batch_chain_ref_is_filled)
        self.ui.txt_batch_chain_model.textChanged.connect(
            self.check_batch_if_txt_batch_chain_model_is_filled)

    def _connect_image_combo_boxes(self):
        self.ui.box_representation.activated.connect(self.show_representation)
        self.ui.box_bg_color.activated.connect(self.choose_bg_color)
        self.ui.box_renderer.activated.connect(self.choose_renderer)
        self.ui.box_ray_trace_mode.activated.connect(self.choose_ray_trace_mode)
        self.ui.box_ray_texture.activated.connect(self.choose_ray_texture)

    def _connect_new_checkbox(self):
        self.ui.cb_new_add_reference.stateChanged.connect(self.show_add_reference)

    def _connect_batch_checkbox(self):
        self.ui.cb_batch_chain_info.stateChanged.connect(
            self.enable_chain_information_input_for_batch)

    def _connect_image_checkbox(self):
        self.ui.cb_transparent_bg.stateChanged.connect(self.decide_transparent_bg)

    def _create_tooltips_sequence_vs_pdb(self):
        # for buttons
        self.ui.btn_prediction_load_reference.setToolTip("Open reference pdb file")
        self.ui.btn_prediction_start.setToolTip("Predict with Colab Notebook")

        # for text fields
        self.ui.txt_prediction_load_reference.setToolTip("Reference file path")
        # self.ui.txt_prediction_chain_ref.setToolTip("Enter chain(s) of reference")
        self.ui.txt_prediction_chain_model.setToolTip("Enter chain(s) of model")

    def _create_tooltips_single_analysis(self):
        # for buttons
        self.ui.btn_analysis_start.setToolTip("Start analysis process")

        # for statusbar
        self.status_bar.setToolTip("Status information: Current process")

    def _create_tooltips_batch(self):
        # for buttons
        self.ui.btn_batch_load_reference.setToolTip("Open reference pdb file")
        self.ui.btn_batch_load_model.setToolTip("Open model pdb files")
        self.ui.btn_batch_start.setToolTip("Start batch process")

        # for text fields
        self.ui.txt_batch_load_reference.setToolTip("Reference file path")
        self.ui.txt_batch_load_model.setToolTip("Model file paths")
        self.ui.txt_batch_chain_ref.setToolTip("Enter chain(s) of reference")
        self.ui.txt_batch_chain_model.setToolTip("Enter chain(s) of models")

        # for checkbox
        self.ui.cb_batch_chain_info.setToolTip("Enable input of chains")

        # for statusbar
        self.status_bar.setToolTip("Status information: Current process")

    def _create_tooltips_image(self):
        # for buttons
        self.ui.btn_save_scene.setToolTip("Create new PyMOL scene")
        self.ui.btn_update_scene.setToolTip("Overwrite current scene")
        self.ui.btn_save_image.setToolTip("Save current viewpoint as png file")

        # for checkboxes
        self.ui.cb_ray_tracing.setToolTip("Enable ray-tracing")
        self.ui.cb_transparent_bg.setToolTip("Enable transparent background")

        # for combo-boxes
        self.ui.box_representation.setToolTip("Choose a representation")
        self.ui.box_bg_color.setToolTip("Choose a background color")
        self.ui.box_renderer.setToolTip("Choose a ray-tracing renderer")
        self.ui.box_ray_trace_mode.setToolTip("Choose a ray-trace mode")

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

    def display_prediction_only_page(self):
        """This function displays the prediction only work area

        """
        self.ui.stackedWidget.setCurrentIndex(1)
        self.ui.lbl_page_title.setText("New Sequence")
        item = self.ui.list_new_seq_notebooks.findItems("AlphaFold",
                                                        Qt.QtCore.Qt.MatchContains |
                                                        Qt.QtCore.Qt.MatchExactly
                                                        )
        self.ui.list_new_seq_notebooks.setCurrentItem(item[0])

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
        self.ui.stackedWidget.setCurrentIndex(2)
        self.ui.lbl_page_title.setText("Sequence vs .pdb")

    def display_sequence_vs_pdb_page(self):
        # regular opening of work area
        self.ui.stackedWidget.setCurrentIndex(10)
        self.ui.lbl_page_title.setText("Sequence vs .pdb")

    def display_single_analysis_page(self):
        """This function displays the single analysis work area

        """
        # pre-process
        # fill chains list widget
        self.ui.list_analysis_ref_chains.clear()
        self.ui.list_analysis_model_chains.clear()

        pdb_files = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")

        # TODO: should there be sub-dirs to determine if it is the ref or the model?
        cmd.load(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{pdb_files[0]}", object="reference_protein")
        cmd.load(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{pdb_files[1]}", object="model_protein")
        ref_chains = cmd.get_chains("reference_protein")
        model_chains = cmd.get_chains("model_protein")
        for chain in ref_chains:
            self.ui.list_analysis_ref_chains.addItem(chain)
        for chain in model_chains:
            self.ui.list_analysis_model_chains.addItem(chain)
        styles_utils.color_button_ready(self.ui.btn_analysis_start)
        cmd.reinitialize()
        # regular work area opening
        self.ui.stackedWidget.setCurrentIndex(3)
        self.ui.lbl_page_title.setText("Single Analysis")

    def display_job_analysis_page(self):
        """This function displays the job analysis work area

        """
        self.ui.stackedWidget.setCurrentIndex(4)
        self.ui.lbl_page_title.setText("Job Analysis")
        global global_var_work_area_history
        global_var_work_area_history.append(4)

    def display_results_page(self):
        """This function displays the results work area

        """
        self.ui.stackedWidget.setCurrentIndex(5)
        self.ui.lbl_page_title.setText("Results")

    def display_image_page(self):
        """This function displays the image work area

        """
        self.ui.stackedWidget.setCurrentIndex(6)
        self.ui.lbl_page_title.setText("Image")

    def display_new_page(self):
        """This function displays the new project work area

        """
        self.ui.list_new_projects.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        tools.scan_workspace_for_vaild_projects(self.workspace_path, self.ui.list_new_projects)
        # workspace_projects: list[str] = os.listdir(self.workspace_path)
        # valid_directories = []
        # # iterates over possible project directories
        # for directory in workspace_projects:
        #     try:
        #         directory_content = os.listdir(f"{self.workspace_path}/{directory}")
        #         # iterates over the content in a single project directory
        #         for content in directory_content:
        #             if content == "project.xml":
        #                 valid_directories.append(directory)
        #     except NotADirectoryError:
        #         print(f"This: {directory} is not a directory.")
        #
        # valid_directories.sort()
        # for project in valid_directories:
        #     self.ui.list_new_projects.addItem(project)

        self.ui.stackedWidget.setCurrentIndex(7)
        self.ui.lbl_page_title.setText("Create new project")

    def display_open_page(self):
        """This function displays the open project work area

        """
        self.ui.list_open_projects.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        tools.scan_workspace_for_vaild_projects(self.workspace_path, self.ui.list_open_projects)
        # workspace_projects: list[str] = os.listdir(self.workspace_path)
        # workspace_projects.sort()
        # for project in workspace_projects:
        #     self.ui.list_open_projects.addItem(project)

        self.ui.stackedWidget.setCurrentIndex(8)
        self.ui.lbl_page_title.setText("Open existing project")

    def display_delete_page(self):
        """This function displays the "delete" project work area

        """
        self.ui.list_delete_projects.clear()
        # pre-process
        self.status_bar.showMessage(self.workspace.text())
        tools.scan_workspace_for_vaild_projects(self.workspace_path, self.ui.list_delete_projects)
        # workspace_projects: list[str] = os.listdir(self.workspace_path)
        # workspace_projects.sort()
        # for project in workspace_projects:
        #     self.ui.list_delete_projects.addItem(project)

        self.ui.stackedWidget.setCurrentIndex(9)
        self.ui.lbl_page_title.setText("Delete existing project")

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
    #         styles_path_btn_ready = f"{project_constants.path_list[1]}/styles/styles_start_button_ready.css"
    #         styles_path_btn_not_ready = f"{project_constants.path_list[1]}/styles/styles_start_button_not_ready.css"
    #     elif sys.platform.startswith("linux"):
    #         # Linux path
    #         styles_path_btn_ready = f"{project_constants.path_list[0]}/styles/styles_start_button_ready.css"
    #         styles_path_btn_not_ready = f"{project_constants.path_list[0]}/styles/styles_start_button_not_ready.css"
    #     elif sys.platform.startswith("win32"):
    #         # Windows path
    #         styles_path_btn_ready = f"{project_constants.path_list[2]}/styles/styles_start_button_ready.css"
    #         styles_path_btn_not_ready = f"{project_constants.path_list[2]}/styles/styles_start_button_not_ready.css"
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
    def open(self):
        """This function opens a project.xml file and fills the input boxes with the right values

        """
        # open file dialog
        try:
            global global_var_project_dict
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self,
                                                                 "Open project file",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "Plugin Project File (*.xml)")
            if file_name == ("", ""):
                tools.quick_log_and_display("info", "No file has been selected.",
                                            self.status_bar, "No file has been selected.")
                return
            # clear current project list and interesting regions
            self.ui.project_list.clear()
            self.ui.cb_interesting_regions.clear()

            xml_file = minidom.parse(file_name[0])
            element = xml_file.getElementsByTagName('project_name')
            job_file_info = Qt.QtCore.QFileInfo(file_name[0])
            job_name = job_file_info.baseName()

            if len(element) == 0:
                # job
                i = 0
                loop_element = ["", ""]
                while loop_element != []:
                    loop_element = xml_file.getElementsByTagName(f'project_{i}')
                    if loop_element != []:
                        path = loop_element[0].getAttribute('value')
                        path_split = path.split("/")
                        project_name = path_split[(len(path_split) - 1)]
                        global_var_project_dict[i] = utils.project_utils.Project(project_name,
                                                                                 f"{self.workspace_path}/{job_name}")
                        index = len("project_of_")
                        model_name = project_name[index:]
                        # TODO: insert modelname
                        global_var_project_dict[i].set_pdb_model(model_name)
                        self.ui.project_list.addItem(global_var_project_dict[i].get_project_name())
                    i += 1

            else:
                # project
                file_name_file_info = Qt.QtCore.QFileInfo(file_name[0])
                global_var_project_dict[0] = utils.project_utils.Project(file_name_file_info.baseName(),
                                                                         self.workspace_path)
                global_var_project_dict[0].set_pdb_model(xml_file.getElementsByTagName('pdb_model')[0].getAttribute('value'))
                # add filename to project list (results tab)
                self.ui.project_list.addItem(global_var_project_dict[0].get_project_name())

            # fill combo box of interesting regions
            results_path = global_var_project_dict[0].get_results_path()
            dir_content = os.listdir(f"{results_path}/images/interesting_regions")
            for tmp_file in dir_content:
                self.ui.cb_interesting_regions.addItem(tmp_file)
            self.ui.project_list.setCurrentRow(0)

        except FileNotFoundError:
            print("File could not be opened.")
        except ValueError:
            print("No file has been selected.")

    def save_as(self):
        """This function saves the current pymol session.

        """
        try:
            file_path = Qt.QtWidgets.QFileDialog.getSaveFileName(self,
                                                                 "Save PyMOL session",
                                                                 f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}",
                                                                 "PyMOL session file (.pse)")
            if file_path == ("", ""):
                tools.quick_log_and_display("info", "No file has been created.", self.status_bar,
                                            "No file has been created.")
            cmd.save(f"{file_path[0]}.pse")
            tools.quick_log_and_display("info", "Saving the pymol session was successful.",
                                        self.status_bar, "Saving was successful.")
        except FileExistsError:
            tools.quick_log_and_display("warning", "File already exists!", self.status_bar,
                                        "File already exists!")
        except pymol.CmdException:
            tools.quick_log_and_display("error", "Unexpected Error from PyMOL while saving the "
                                                 "current pymol session", self.status_bar,
                                        "Unexpected Error from PyMOL")

    def save(self):
        """This function saves the current pymol session.

        """
        global global_var_project_dict
        try:
            file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
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

    def restore_settings(self):
        """This function deletes the old settings.xml and creates a new one,
        with the default values.

        """
        out = gui_utils.warning_dialog_restore_settings("Are you sure you want to restore all settings?", "")
        if out:
            tools.restore_default_settings()
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
        tools.open_global_settings()
        self._setup_statusbar()

    @staticmethod
    def open_documentation():
        """This function opens the official plugin documentation as HTML page.

        """
        #webbrowser.open_new(f"file://{os.getcwd()}/docs/pymol_plugin/build/html/index.html")
        # opens the documentation of the os
        if sys.platform.startswith("darwin"):
            # macOS path
            webbrowser.open_new(f"file://{project_constants.path_list[1]}/docs/pymol_plugin/build/html/index.html")
        elif sys.platform.startswith("linux"):
            # Linux path
            webbrowser.open_new(f"file://{project_constants.path_list[0]}/docs/pymol_plugin/build/html/index.html")
        elif sys.platform.startswith("win32"):
            # Windows path
            webbrowser.open_new(f"file://{project_constants.path_list[2]}/docs/pymol_plugin/build/html/index.html")

    @staticmethod
    def open_documentation_pdf():
        """This function opens the official plugin documentation as PDF.

        """
        # webbrowser.open_new(
        #     f"file://{os.getcwd()}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")
        # opens the documentation of the os
        if sys.platform.startswith("darwin"):
            # macOS path
            webbrowser.open_new(f"file://{project_constants.path_list[1]}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")
        elif sys.platform.startswith("linux"):
            # Linux path
            webbrowser.open_new(f"file://{project_constants.path_list[0]}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")
        elif sys.platform.startswith("win32"):
            # Windows path
            webbrowser.open_new(f"file://{project_constants.path_list[2]}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")

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
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/pdb"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/alignment_files"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/distance_csv"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/images"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/images/interesting_regions"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/plots"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/plots/distance_histogram"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/plots/distance_plot"),
                    Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}_{pdb_name}/results/sessions/"),
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

    # new project
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
            styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
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
            styles_utils.color_button_ready(self.ui.btn_new_create_project)

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
            styles_utils.color_button_ready(self.ui.btn_new_create_project)
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
            styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
        elif len(self.ui.txt_new_choose_reference.text()) < 4:
            self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
            styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
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
                styles_utils.color_button_ready(self.ui.btn_new_create_project)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
                styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
                return
            except FileNotFoundError:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
                self.ui.lbl_new_status_choose_reference.setText("Invalid PDB ID.")
                self.ui.btn_new_create_project.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
                return
        else:
            if self.ui.txt_new_choose_reference.text().find("/") == -1:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_new_create_project.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_new_create_project)

            elif self.ui.txt_new_choose_reference.text().find("\\") == -1:
                self.ui.txt_new_choose_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_new_create_project.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_new_create_project)

    def validate_project_name(self):
        """This function validates the input of the project name in real-time

        """

        if self.ui.list_new_projects.currentItem() is not None:
            self.ui.list_new_projects.currentItem().setSelected(False)
        # set color for lineEdit
        self.ui.txt_new_project_name.setStyleSheet("background-color: #FC5457")
        if len(self.ui.txt_new_project_name.text()) == 0:
            self.ui.lbl_new_status_project_name.setText("")
            self.ui.cb_new_add_reference.setCheckable(False)
            self.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
            self.ui.btn_new_create_project.setEnabled(False)
            styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
            return
        elif len(self.ui.txt_new_project_name.text()) > 20:
            self.ui.lbl_new_status_project_name.setText("Project name is too long (max. 20 characters).")
            self.ui.btn_new_create_project.setEnabled(False)
            styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
            return
        else:
            regex = Qt.QtCore.QRegularExpression()
            # TODO: has no dash in regex!
            regex.setPattern("\\w{20}")
            validator = QtGui.QRegularExpressionValidator(regex)
            for i in range(len(self.ui.txt_new_project_name.text())):
                result = validator.validate(self.ui.txt_new_project_name.text(), i)
                if result[0] > 0:
                    self.ui.txt_new_project_name.setStyleSheet("background-color: #33C065")
                    self.ui.lbl_new_status_project_name.setText("")
                    self.ui.cb_new_add_reference.setCheckable(True)
                    self.ui.cb_new_add_reference.setStyleSheet("color: black;")
                    self.ui.btn_new_create_project.setEnabled(True)
                    styles_utils.color_button_ready(self.ui.btn_new_create_project)
                else:
                    self.ui.txt_new_project_name.setStyleSheet("background-color: #FC5457")
                    self.ui.lbl_new_status_project_name.setText("Invalid character.")
                    self.ui.cb_new_add_reference.setCheckable(False)
                    self.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
                    self.ui.btn_new_create_project.setEnabled(False)
                    styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
                    return
            item = self.ui.list_new_projects.findItems(self.ui.txt_new_project_name.text(),
                                                          Qt.QtCore.Qt.MatchContains |
                                                          Qt.QtCore.Qt.MatchExactly
                                                          )
            if len(item) != 0:
                self.ui.list_new_projects.setCurrentItem(item[0])
                self.ui.txt_new_project_name.setStyleSheet("background-color: #FC5457")
                self.ui.lbl_new_status_project_name.setText("Project name already exists.")
                self.ui.cb_new_add_reference.setCheckable(False)
                self.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
                self.ui.btn_new_create_project.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
            # else:
            #     self.ui.list_widget_projects.currentItem().setSelected(False)
            #     self.ui.txt_prediction_project_name.setStyleSheet("background-color: green")
            #     self.ui.btn_prediction_next_1.setEnabled(True)
            #     styles_utils.color_button_ready(self.ui.btn_prediction_next_1)
            print("Check successful.")

    def create_new_project(self):
        """This function creates a new project based on the plugin New ... page

        """
        self._init_hide_ui_elements()
        self.ui.lbl_current_project_name.setText(self.ui.txt_new_project_name.text())
        self.status_bar.showMessage(f"Current project path: {self.workspace_path}/{self.ui.txt_new_project_name.text()}")
        # save project folder in current workspace
        new_project = project_utils.Project(self.ui.txt_new_project_name.text(), self.workspace_path)
        new_project.create_project_tree()
        new_project.save_project_to_xml()

        # # mkdir project
        # folder_paths = [
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/pdb"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/alignment_files"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/distance_csv"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/images"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/images/interesting_regions"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/plots"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/plots/distance_histogram"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/plots/distance_plot"),
        #     Path(f"{self.workspace_path}/{self.ui.txt_new_project_name.text()}/results/sessions/"),
        # ]
        # for path in folder_paths:
        #     os.mkdir(path)

        # save reference .pdb
        if self.ui.cb_new_add_reference.checkState() == 2 and self.ui.btn_new_create_project.isEnabled() == True:
            if len(self.ui.txt_new_choose_reference.text()) == 4:
                # PDB ID as input
                # the pdb file gets saved in a scratch directory where it gets deleted immediately
                pdb_id = self.ui.txt_new_choose_reference.text().upper()
                cmd.fetch(pdb_id, type="pdb", path=self.scratch_path)
                shutil.move(f"{self.scratch_path}/{pdb_id}.pdb", new_project.get_pdb_path())
                cmd.reinitialize()
            else:
                # local pdb file as input
                shutil.copy(self.ui.txt_new_choose_reference.text(), new_project.get_pdb_path())

        # shows options which can be done with the data in the project folder
        self.ui.lbl_prediction.show()
        if self.ui.txt_new_choose_reference.text() == "":
            # no reference pdb given
            self.ui.btn_prediction_only_page.show()
            self.ui.btn_prediction_page.hide()
            # change to page with new sequence
            self.ui.btn_prediction_only_page.click()
            # clean up new project page
            self.ui.txt_new_project_name.clear()
            self.ui.cb_new_add_reference.setCheckState(0)
        else:
            # reference pdb given
            self.ui.action_add_model.setVisible(True)
            self.ui.action_add_multiple_models.setVisible(True)
            self.ui.btn_prediction_page.show()
            self.ui.btn_prediction_only_page.hide()
            # change to page with sequence vs pdb
            self.display_sequence_vs_pdb_page()
            # fill chains list widget
            if len(self.ui.txt_new_choose_reference.text()) == 4:
                pdb_id = self.ui.txt_new_choose_reference.text().upper()
                tmp_ref_protein = core.Protein(pdb_id, export_data_dir=self.scratch_path)
                tmp_ref_protein.clean_pdb_file()
                chains = cmd.get_chains(pdb_id)
            else:
                cmd.load(self.ui.txt_new_choose_reference.text(), object="reference_protein")
                chains = cmd.get_chains("reference_protein")
            for chain in chains:
                self.ui.list_s_v_p_ref_chains.addItem(chain)
            styles_utils.color_button_ready(self.ui.btn_s_v_p_start)
            cmd.reinitialize()
            self.ui.btn_s_v_p_start.setEnabled(True)
            # clean up new project page
            self.ui.txt_new_project_name.clear()
            self.ui.txt_new_choose_reference.clear()
            self.ui.cb_new_add_reference.setCheckState(0)
            self.ui.btn_new_create_project.isEnabled()

    # open
    def validate_open_search(self):
        """This function validates the input of the project name in real-time

        """

        if self.ui.list_open_projects.currentItem() is not None:
            self.ui.list_open_projects.currentItem().setSelected(False)
        # set color for lineEdit
        self.ui.txt_open_search.setStyleSheet("background-color: white")
        if len(self.ui.txt_open_search.text()) == 0:
            self.ui.lbl_open_status_search.setText("")
            self.ui.txt_open_selected_project.setText("")
            return
        else:
            # regex = Qt.QtCore.QRegularExpression()
            # # TODO: has no dash in regex!
            # regex.setPattern("\\w{20}")
            # validator = QtGui.QRegularExpressionValidator(regex)
            # for i in range(len(self.ui.txt_open_search.text())):
            #     result = validator.validate(self.ui.txt_open_search.text(), i)
            #     if result[0] > 0:
            #         self.ui.txt_open_search.setStyleSheet("background-color: #33C065")
            #         self.ui.lbl_new_status_project_name.setText("")
            #         self.ui.cb_new_add_reference.setCheckable(True)
            #         self.ui.cb_new_add_reference.setStyleSheet("color: black;")
            #         self.ui.btn_new_create_project.setEnabled(True)
            #         styles_utils.color_button_ready(self.ui.btn_new_create_project)
            #     else:
            #         self.ui.txt_open_search.setStyleSheet("background-color: #FC5457")
            #         self.ui.lbl_new_status_project_name.setText("Invalid character.")
            #         self.ui.cb_new_add_reference.setCheckable(False)
            #         self.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
            #         self.ui.btn_new_create_project.setEnabled(False)
            #         styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
            #         return
            item = self.ui.list_open_projects.findItems(self.ui.txt_open_search.text(),
                                                          Qt.QtCore.Qt.MatchContains |
                                                          Qt.QtCore.Qt.MatchExactly
                                                          )
            if len(item) != 0:
                self.ui.list_open_projects.setCurrentItem(item[0])
                self.ui.txt_open_selected_project.setText(self.ui.list_open_projects.currentItem().text())
                self.ui.lbl_open_status_search.setText("")
            else:
                self.ui.txt_open_selected_project.setText("")
                self.ui.txt_open_search.setStyleSheet("background-color: #FC5457")
                self.ui.lbl_open_status_search.setText("Project name does not exists.")

            # else:
            #     self.ui.list_widget_projects.currentItem().setSelected(False)
            #     self.ui.txt_prediction_project_name.setStyleSheet("background-color: green")
            #     self.ui.btn_prediction_next_1.setEnabled(True)
            #     styles_utils.color_button_ready(self.ui.btn_prediction_next_1)
            print("Check successful.")

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
            styles_utils.color_button_not_ready(self.ui.btn_open_open_project)
        else:
            self.ui.btn_open_open_project.setEnabled(True)
            styles_utils.color_button_ready(self.ui.btn_open_open_project)

    def open_project(self):
        """This function opens an existing project

        """
        self._init_hide_ui_elements()
        self.ui.lbl_current_project_name.setText(self.ui.list_open_projects.currentItem().text())
        # check if directory in empty
        results_paths = ["results/alignment_files"]
        dir_content = os.listdir(f"{self.workspace_path}/{self.ui.list_open_projects.currentItem().text()}/{results_paths[0]}")
        self.ui.lbl_handle_pymol_session.show()
        self.ui.btn_image_page.show()
        if not dir_content:
            # no analysis was done
            tmp_content = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")
            if len(tmp_content) == 1:
                # if no model is in the pdb dir
                self.display_image_page()
                self.ui.btn_save_project.show()
                self.ui.action_file_save_as.setVisible(True)
                self.ui.action_add_model.setVisible(True)
                self.ui.action_add_multiple_models.setVisible(True)
            elif len(tmp_content) == 2:
                # if a model is in the pdb dir
                self.display_single_analysis_page()
                self.ui.lbl_analysis.show()
                self.ui.btn_single_analysis_page.show()
                self.ui.btn_save_project.show()
                self.ui.action_file_save_as.setVisible(True)
                self.ui.action_add_model.setVisible(False)
                self.ui.action_add_multiple_models.setVisible(False)
            elif len(tmp_content) == 0:
                # check if data is corrupted
                response = tools.check_results_for_integrity(self.workspace_path,
                                                             self.ui.lbl_current_project_name.text())
                if response[0]:
                    self.ui.lbl_handle_pymol_session.hide()
                    self.ui.btn_image_page.hide()
                    self.ui.lbl_prediction.show()
                    self.ui.btn_prediction_only_page.show()
                    self.display_prediction_only_page()

            else:
                gui_utils.error_project_data_is_invalid(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")
            return
        else:
            # check if data is corrupted
            response = tools.check_results_for_integrity(self.workspace_path, self.ui.lbl_current_project_name.text())
            if response[0]:
                # an analysis happened
                # show results gui elements
                self.ui.lbl_analysis.show()
                self.ui.btn_results_page.show()
                self.display_results_page()
                self.ui.btn_save_project.show()
                self.ui.action_file_save_as.setVisible(True)
            else:
                gui_utils.error_project_data_is_invalid(response[1])
                self._init_hide_ui_elements()
                self.ui.lbl_current_project_name.clear()
                self.ui.btn_save_project.hide()
            return

    # delete
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
            styles_utils.color_button_not_ready(self.ui.btn_delete_delete_project)
        else:
            self.ui.btn_delete_delete_project.setEnabled(True)
            styles_utils.color_button_ready(self.ui.btn_delete_delete_project)

    def validate_delete_search(self):
        """This function validates the input of the project name in real-time

        """

        if self.ui.list_delete_projects.currentItem() is not None:
            self.ui.list_delete_projects.currentItem().setSelected(False)
        # set color for lineEdit
        self.ui.txt_delete_search.setStyleSheet("background-color: white")
        if len(self.ui.txt_delete_search.text()) == 0:
            self.ui.lbl_delete_status_search.setText("")
            self.ui.txt_delete_selected_projects.setText("")
            return
        else:
            # regex = Qt.QtCore.QRegularExpression()
            # # TODO: has no dash in regex!
            # regex.setPattern("\\w{20}")
            # validator = QtGui.QRegularExpressionValidator(regex)
            # for i in range(len(self.ui.txt_delete_search.text())):
            #     result = validator.validate(self.ui.txt_delete_search.text(), i)
            #     if result[0] > 0:
            #         self.ui.txt_delete_search.setStyleSheet("background-color: #33C065")
            #         self.ui.lbl_new_status_project_name.setText("")
            #         self.ui.cb_new_add_reference.setCheckable(True)
            #         self.ui.cb_new_add_reference.setStyleSheet("color: black;")
            #         self.ui.btn_new_create_project.setEnabled(True)
            #         styles_utils.color_button_ready(self.ui.btn_new_create_project)
            #     else:
            #         self.ui.txt_delete_search.setStyleSheet("background-color: #FC5457")
            #         self.ui.lbl_new_status_project_name.setText("Invalid character.")
            #         self.ui.cb_new_add_reference.setCheckable(False)
            #         self.ui.cb_new_add_reference.setStyleSheet("color: #E1E1E1;")
            #         self.ui.btn_new_create_project.setEnabled(False)
            #         styles_utils.color_button_not_ready(self.ui.btn_new_create_project)
            #         return
            item = self.ui.list_delete_projects.findItems(self.ui.txt_delete_search.text(),
                                                          Qt.QtCore.Qt.MatchContains |
                                                          Qt.QtCore.Qt.MatchExactly
                                                          )
            if len(item) != 0:
                self.ui.list_delete_projects.setCurrentItem(item[0])
                self.ui.txt_delete_selected_projects.setText(self.ui.list_delete_projects.currentItem().text())
                self.ui.lbl_delete_status_search.setText("")
            else:
                self.ui.txt_delete_selected_projects.setText("")
                self.ui.txt_delete_search.setStyleSheet("background-color: #FC5457")
                self.ui.lbl_delete_status_search.setText("Project name does not exists.")

            # else:
            #     self.ui.list_widget_projects.currentItem().setSelected(False)
            #     self.ui.txt_prediction_project_name.setStyleSheet("background-color: green")
            #     self.ui.btn_prediction_next_1.setEnabled(True)
            #     styles_utils.color_button_ready(self.ui.btn_prediction_next_1)
            print("Check successful.")

    def delete_project(self):
        """This function deletes an existing project

        """
        # popup message which warns the user that the selected project gets deleted
        response: bool = gui_utils.warning_message_project_gets_deleted()

        if response is True:
            shutil.rmtree(f"{self.workspace_path}/{self.ui.txt_delete_selected_projects.text()}")
            if self.ui.txt_delete_selected_projects.text() == self.ui.lbl_current_project_name.text():
                self.ui.lbl_current_project_name.clear()
                self._init_hide_ui_elements()
            self.ui.txt_delete_selected_projects.clear()
            # update list
            self.ui.list_delete_projects.clear()
            # pre-process
            self.status_bar.showMessage(self.workspace.text())
            self.ui.list_delete_projects.clear()
            self.status_bar.showMessage(self.workspace.text())
            tools.scan_workspace_for_vaild_projects(self.workspace_path, self.ui.list_delete_projects)
        else:
            return

    # save
    def save_project(self):
        """This function saves the "project" which is currently only the pymol session

        """
        session_file = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/sessions/")
        try:
            cmd.save(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/sessions/{session_file[0]}")
            self.status_bar.showMessage("Saved project successfully.")
        except pymol.CmdException:
            self.status_bar.showMessage("PyMOL internal error while saving.")
        except IndexError:
            cmd.save(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results/sessions/auto_generated_session_file.pse")
            self.status_bar.showMessage("Saved project successfully.")

    # Prediction
    def check_prediction_only_if_txt_notebook_url_is_filled(self):
        """This function checks if a reference pdb file is selected.

        """
        self.__check_start_possibility_prediction_only()

    def predict_only(self):
        """This function is used to predict with any google colab notebook.

        """
        notebook_name = self.ui.list_new_seq_notebooks.currentItem().text()
        if notebook_name == project_constants.OFFICIAL_NOTEBOOK_NAME:
            webbrowser.open_new(project_constants.OFFICIAL_NOTEBOOK_URL)
            self.ui.btn_prediction_only_start.setEnabled(False)
            styles_utils.color_button_not_ready(self.ui.btn_prediction_only_start)
            # alphafold specific process
            archive = "prediction.zip"
            source_path = global_utils.global_var_settings_obj.get_prediction_path()
            filename = f"{source_path}/{archive}"
            while os.path.isfile(filename) is False:
                print("Prediction is still running ...")
                time.sleep(20)
            # move prediction.zip in scratch folder
            shutil.copy(filename, f"{self.scratch_path}/{archive}")
            shutil.unpack_archive(f"{self.scratch_path}/{archive}", self.scratch_path, "zip")
            shutil.copy(f"{self.scratch_path}/prediction/selected_prediction.pdb",
                        f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/selected_prediction.pdb")
            shutil.rmtree(f"{self.scratch_path}/prediction")
            os.remove(f"{self.scratch_path}/{archive}")
            try:
                cmd.load(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/selected_prediction.pdb")
            except pymol.CmdException:
                print("Loading the model failed.")
                return
        else:
            return
        self.ui.lbl_prediction.hide()
        self.ui.btn_prediction_only_page.hide()
        self.ui.btn_save_project.show()
        self.ui.lbl_handle_pymol_session.show()
        self.ui.btn_image_page.show()
        self.display_image_page()

    # Prediction + Analysis
    # def validate_project_name(self):
    #     """This function validates the input of the project name in real-time
    #
    #     """
    #     if self.ui.list_widget_projects.currentItem() is not None:
    #         self.ui.list_widget_projects.currentItem().setSelected(False)
    #     # set color for lineEdit
    #     self.ui.txt_prediction_project_name.setStyleSheet("background-color: #FC5457")
    #     if len(self.ui.txt_prediction_project_name.text()) == 0:
    #         self.ui.btn_prediction_next_1.setEnabled(False)
    #         styles_utils.color_button_not_ready(self.ui.btn_prediction_next_1)
    #         return
    #     elif len(self.ui.txt_prediction_project_name.text()) > 20:
    #         self.ui.btn_prediction_next_1.setEnabled(False)
    #         styles_utils.color_button_not_ready(self.ui.btn_prediction_next_1)
    #         self.ui.lbl_prediction_status_project_name.setText("Project name is too long (max. 20 characters).")
    #         return
    #     else:
    #         regex = Qt.QtCore.QRegularExpression()
    #         # TODO: has no dash in regex!
    #         regex.setPattern("\\w{20}")
    #         validator = QtGui.QRegularExpressionValidator(regex)
    #         for i in range(len(self.ui.txt_prediction_project_name.text())):
    #             result = validator.validate(self.ui.txt_prediction_project_name.text(), i)
    #             if result[0] > 0:
    #                 self.ui.txt_prediction_project_name.setStyleSheet("background-color: #33C065")
    #                 self.ui.btn_prediction_next_1.setEnabled(True)
    #                 styles_utils.color_button_ready(self.ui.btn_prediction_next_1)
    #             else:
    #                 self.ui.txt_prediction_project_name.setStyleSheet("background-color: #FC5457")
    #                 self.ui.btn_prediction_next_1.setEnabled(False)
    #                 styles_utils.color_button_not_ready(self.ui.btn_prediction_next_1)
    #                 self.ui.lbl_prediction_status_project_name.setText("Invalid character.")
    #                 return
    #         item = self.ui.list_widget_projects.findItems(self.ui.txt_prediction_project_name.text(),
    #                                                       Qt.QtCore.Qt.MatchContains |
    #                                                       Qt.QtCore.Qt.MatchExactly
    #                                                       )
    #         if len(item) != 0:
    #             self.ui.list_widget_projects.setCurrentItem(item[0])
    #             self.ui.txt_prediction_project_name.setStyleSheet("background-color: #FC5457")
    #             self.ui.btn_prediction_next_1.setEnabled(False)
    #             styles_utils.color_button_not_ready(self.ui.btn_prediction_next_1)
    #             self.ui.lbl_prediction_status_project_name.setText("Project name already exists.")
    #         # else:
    #         #     self.ui.list_widget_projects.currentItem().setSelected(False)
    #         #     self.ui.txt_prediction_project_name.setStyleSheet("background-color: green")
    #         #     self.ui.btn_prediction_next_1.setEnabled(True)
    #         #     styles_utils.color_button_ready(self.ui.btn_prediction_next_1)
    #         print("Check successful.")

    def show_prediction_load_reference(self):
        """Shows the text field and tool button for the load reference functionality

        """
        self.ui.lbl_prediction_load_reference.show()
        self.ui.txt_prediction_load_reference.show()
        self.ui.btn_prediction_load_reference.show()
        self.ui.txt_prediction_project_name.setEnabled(False)
        self.ui.list_widget_projects.setEnabled(False)
        self.ui.btn_prediction_next_1.hide()
        self.ui.btn_prediction_next_2.show()
        self.ui.btn_prediction_next_2.setEnabled(False)
        self.ui.btn_prediction_back_2.show()
        self.ui.txt_prediction_load_reference.setEnabled(True)
        self.ui.txt_prediction_project_name.setStyleSheet("background-color: white")

    def hide_prediction_load_reference(self):
        """Hides the text field and tool button for the load reference functionality

        """
        self.ui.lbl_prediction_load_reference.hide()
        self.ui.txt_prediction_load_reference.hide()
        self.ui.btn_prediction_load_reference.hide()
        self.ui.txt_prediction_project_name.setEnabled(True)
        self.ui.list_widget_projects.setEnabled(True)
        self.ui.btn_prediction_next_1.show()
        self.ui.btn_prediction_next_2.hide()
        self.ui.btn_prediction_back_2.hide()
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
            styles_utils.color_button_ready(self.ui.btn_prediction_next_2)
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
            styles_utils.color_button_not_ready(self.ui.btn_prediction_next_2)
        elif len(self.ui.txt_prediction_load_reference.text()) < 4:
            self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
            self.ui.btn_prediction_next_2.setEnabled(False)
            styles_utils.color_button_not_ready(self.ui.btn_prediction_next_2)
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
                styles_utils.color_button_ready(self.ui.btn_prediction_next_2)
            # if the id does not exist an exception gets raised
            except pymol.CmdException:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_prediction_next_2)
                return
            except FileNotFoundError:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_prediction_next_2)
                self.ui.lbl_prediction_status_load_reference.setText("Invalid PDB ID.")
                return
        else:
            if self.ui.txt_prediction_load_reference.text().find("/") == -1:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_prediction_next_2)

            elif self.ui.txt_prediction_load_reference.text().find("\\") == -1:
                self.ui.txt_prediction_load_reference.setStyleSheet("background-color: #FC5457")
                self.ui.btn_prediction_next_2.setEnabled(False)
                styles_utils.color_button_not_ready(self.ui.btn_prediction_next_2)

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
    #     styles_utils.color_button_ready(self.ui.btn_prediction_start)
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
        styles_utils.color_button_ready(self.ui.btn_prediction_only_start)

    def predict(self):
        """This function opens a webbrowser with a colab notebook, to run the
        prediction. In addition, it runs the entire analysis after the
        prediction.

        """
        ref_chain_items = self.ui.list_widget_ref_chains.selectedItems()
        ref_chains = []
        for chain in ref_chain_items:
            ref_chains.append(chain.text())
        # global global_var_abort_prediction
        # global_var_abort_prediction = False
        # check if a prediction is already finished
        if os.path.isfile(f"{global_utils.global_var_settings_obj.get_prediction_path()}/prediction.zip"):
            self.status_bar.showMessage(
                f"Warning! | Current Workspace: {self.workspace_path}")
            check = gui_utils.warning_message_prediction_exists(
                f"The prediction is here: {global_utils.global_var_settings_obj.get_prediction_path()}/prediction.zip ",
                f"{global_utils.global_var_settings_obj.get_prediction_path()}/prediction.zip")
            if not check:
                return
        # creates project without xml creation and model adding these come after the prediction
        project = utils.project_utils.Project(self.ui.txt_prediction_project_name.text(),
                                              self.workspace_path)
        project.create_project_tree()
        project.set_pdb_file(self.ui.txt_prediction_load_reference.text())
        project.set_pdb_id(self.ui.txt_prediction_load_reference.text())
        # TODO: check if the new list type conflicts with the analysis
        project.set_ref_chains(ref_chains)
        project.set_model_chains((self.ui.txt_prediction_chain_model.text()))
        # gets reference filename and filepath
        if len(self.ui.txt_prediction_load_reference.text()) == 4:
            tmp_protein = core.Protein(self.ui.txt_prediction_load_reference.text(),
                                       export_data_dir=project.get_pdb_path())
            tmp_protein.clean_pdb_file()
            REFERENCE_OBJ_NAME = self.ui.txt_prediction_load_reference.text()
            REFERENCE_DIR = project.get_pdb_path()
        else:
            ref_file_info = Qt.QtCore.QFileInfo(self.ui.txt_prediction_load_reference.text())
            REFERENCE_OBJ_NAME = ref_file_info.baseName()
            REFERENCE_DIR = ref_file_info.canonicalPath()
        # starting the default web browser to display the colab notebook
        self.status_bar.showMessage("Opening Google colab notebook ...")
        if self.ui.action_settings_model_w_off_colab_notebook.isChecked():
            webbrowser.open_new(project_constants.OFFICIAL_NOTEBOOK_URL)
        else:
            webbrowser.open_new(project_constants.NOTEBOOK_URL)

        # # waiting for the colab notebook to finish
        # archive = "prediction.zip"
        # source_path = global_utils.global_var_settings_obj.get_prediction_path()
        # FILE_NAME = f"{source_path}/{archive}"
        # # flag = False
        # # while flag == False:
        # #     print("AlphaFold is still running ...")
        # #     time.sleep(5)
        # #     # time.sleep(120)
        # # TODO: loop doesn't work correctly
        # while os.path.isfile(FILE_NAME) is False:
        #     print("AlphaFold is still running ...")
        #     # time.sleep(5)
        #     time.sleep(20)
        #     # time.sleep(120)
        #     # global global_var_abort_prediction
        #     # if global_var_abort_prediction:
        #     #     return

        # alphafold specific process
        archive = "prediction.zip"
        source_path = global_utils.global_var_settings_obj.get_prediction_path()
        filename = f"{source_path}/{archive}"
        while os.path.isfile(filename) is False:
            print("AlphaFold is still running ...")
            time.sleep(20)
        # move prediction.zip in scratch folder
        shutil.copy(filename, f"{self.scratch_path}/{archive}")
        shutil.unpack_archive(f"{self.scratch_path}/{archive}", self.scratch_path, "zip")
        shutil.copy(f"{self.scratch_path}/prediction/selected_prediction.pdb",
                    f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/selected_prediction.pdb")
        shutil.rmtree(f"{self.scratch_path}/prediction")
        os.remove(f"{self.scratch_path}/{archive}")

        # ----------------------------------------------------------------- #
        # start of the analysis algorithm
        self.status_bar.showMessage("Protein structure analysis started ...")

        # # extracts and moves the prediction.pdb to the workspace/pdb folder
        # tools.extract_and_move_model_pdb(
        #     str(source_path), f"{str(source_path)}/tmp", archive, project.get_pdb_path())

        # gets model filename and filepath
        PREDICTION_NAME = tools.get_prediction_file_name(project.get_pdb_path())
        full_model_file_path = f"{project.get_pdb_path()}/{PREDICTION_NAME[0]}"
        model_file_info = Qt.QtCore.QFileInfo(full_model_file_path)
        MODEL_OBJ_NAME = model_file_info.baseName()
        MODEL_DIR = model_file_info.canonicalPath()

        # set model in project object
        project.set_pdb_model(full_model_file_path)
        project.create_xml_file()

        # create the Protein object for the reference
        reference_protein: list[core.Protein] = [core.Protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]

        # create model Protein object
        model_proteins: list[core.Protein] = [core.Protein(MODEL_OBJ_NAME, MODEL_DIR)]
        # sets the filepath of the model in the project xml file
        export_dir = project.get_results_path()
        structure_analysis = structure_analysis_utils.StructureAnalysis(
            reference_protein, model_proteins,
            project.get_ref_chains().split(","), project.get_model_chains().split(","),
            export_dir, cycles=global_utils.global_var_settings_obj.get_cycles(),
            cutoff=global_utils.global_var_settings_obj.get_cutoff(),
        )
        structure_analysis.create_selection_for_proteins(structure_analysis.ref_chains,
                                                         structure_analysis.reference_protein)
        structure_analysis.create_selection_for_proteins(structure_analysis.model_chains,
                                                         structure_analysis.model_proteins)
        structure_analysis.do_analysis_in_pymol(structure_analysis.create_protein_pairs(),
                                                self.status_bar, "2")

    # Single Analysis
    def analysis_next_step(self):
        self.ui.lbl_analysis_model_chains.show()
        self.ui.list_analysis_model_chains.show()
        self.ui.btn_analysis_back.show()
        self.ui.btn_analysis_start.show()
        self.ui.btn_analysis_next.hide()
        self.ui.lbl_analysis_ref_chains.setStyleSheet("color: #CACACA")
        self.ui.list_analysis_ref_chains.setEnabled(False)

    def analysis_back_step(self):
        self.ui.lbl_analysis_model_chains.hide()
        self.ui.list_analysis_model_chains.hide()
        self.ui.btn_analysis_back.hide()
        self.ui.btn_analysis_start.hide()
        self.ui.btn_analysis_next.show()
        self.ui.lbl_analysis_ref_chains.setStyleSheet("color: black")
        self.ui.list_analysis_ref_chains.setEnabled(True)

    def load_reference_for_analysis(self):
        """This function opens a file dialog to choose a .pdb file as
        reference and displays the path in a text box

        """
        try:
            # open file dialog
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
            # display path in text box
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self.ui.txt_analysis_load_reference.setText(str(file_name[0]))
            self.status_bar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_prediction()
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the reference failed!")
        except ValueError:
            print("No file has been selected.")
            self.__check_start_possibility()

    def load_model_for_analysis(self):
        """This function opens a file dialog to choose a .pdb file as
        model and displays the path in a text box

        """
        try:
            # open file dialog
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Model",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self.ui.txt_analysis_load_model.setText(str(file_name[0]))
            self.status_bar.showMessage("Loading the model was successful.")
            self.__check_start_possibility()
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the model failed!")
        except ValueError:
            print("No file has been selected.")

    def enable_chain_information_input_for_analysis(self):
        """This function enables the text boxes to enter the chains for the
        reference and the model

        """
        try:
            self.ui.txt_analysis_chain_ref.setEnabled(self.ui.cb_analysis_chain_info.checkState())
            self.ui.txt_analysis_chain_model.setEnabled(self.ui.cb_analysis_chain_info.checkState())
            self.status_bar.showMessage("Enter the chain information.")
            self.__check_start_possibility()
        except Exception:
            self.status_bar.showMessage("Unexpected Error.")

    def check_analysis_if_txt_analysis_project_name_is_filled(self):
        """This function checks if a project name is entered into the text field

        """
        self.__check_start_possibility()

    def check_analysis_if_txt_analysis_chain_ref_is_filled(self):
        """This function checks if any chains are in the text field for the
        reference.

        """
        self.__check_start_possibility()

    def check_analysis_if_txt_analysis_chain_model_is_filled(self):
        """This function checks if any chains are in the text field for the
        model.

        """
        self.__check_start_possibility()

    def start_process(self):
        """This function contains the main analysis algorithm for the
        Protein structure comparison.

        """
        self.ui.btn_analysis_start.setEnabled(False)
        self.status_bar.showMessage("Protein structure analysis started ...")
        cmd.reinitialize()

        model_full_filepath = self.ui.txt_analysis_load_model.toPlainText()
        model_file_info = Qt.QtCore.QFileInfo(model_full_filepath)
        MODEL_OBJ_NAME = model_file_info.baseName()
        MODEL_DIR = model_file_info.canonicalPath()

        project = utils.project_utils.Project(self.ui.txt_analysis_project_name.text(),
                                              self.workspace_path)
        if project.create_project_tree() is False:
            self.ui.btn_analysis_start.setEnabled(True)
            return
        project.set_pdb_file(self.ui.txt_analysis_load_reference.text())
        project.set_pdb_id(self.ui.txt_analysis_load_reference.text())
        project.set_pdb_model(self.ui.txt_analysis_load_model.toPlainText())
        project.set_ref_chains(self.ui.txt_analysis_chain_ref.text())
        project.set_model_chains((self.ui.txt_analysis_chain_model.text()))
        project.create_xml_file()

        # gets reference filename and filepath
        if len(self.ui.txt_analysis_load_reference.text()) == 4:
            tmp_protein = core.Protein(self.ui.txt_analysis_load_reference.text(),
                                       export_data_dir=project.get_pdb_path())
            tmp_protein.clean_pdb_file()
            REFERENCE_OBJ_NAME = self.ui.txt_analysis_load_reference.text()
            REFERENCE_DIR = project.get_pdb_path()
        else:
            ref_file_info = Qt.QtCore.QFileInfo(self.ui.txt_analysis_load_reference.text())
            REFERENCE_OBJ_NAME = ref_file_info.baseName()
            REFERENCE_DIR = ref_file_info.canonicalPath()

        reference_protein: list[core.Protein] = [core.Protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]
        model_proteins: list[core.Protein] = [core.Protein(MODEL_OBJ_NAME, MODEL_DIR)]
        export_dir = project.get_results_path()
        structure_analysis = structure_analysis_utils.StructureAnalysis(
            reference_protein, model_proteins,
            project.get_ref_chains().split(","), project.get_model_chains().split(","),
            export_dir, cycles=global_utils.global_var_settings_obj.get_cycles(),
            cutoff=global_utils.global_var_settings_obj.get_cutoff(),
        )
        structure_analysis.create_selection_for_proteins(structure_analysis.ref_chains,
                                                         structure_analysis.reference_protein)
        structure_analysis.create_selection_for_proteins(structure_analysis.model_chains,
                                                         structure_analysis.model_proteins)
        protein_pairs = structure_analysis.create_protein_pairs()
        structure_analysis.do_analysis_in_pymol(protein_pairs,
                                                self.status_bar, self.ui.progress_bar_analysis)
        del structure_analysis

    # Batch
    def load_reference_for_batch(self):
        """This function opens a file dialog to choose a .pdb file as
        reference and displays the path in a text box

        """
        try:
            # open file dialog
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
            # display path in text box
            if file_name == ("", ""):
                raise ValueError
            # display path in text box
            self.ui.txt_batch_load_reference.setText(str(file_name[0]))
            self.status_bar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_batch()
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the reference failed!")
        except ValueError:
            print("No file has been selected.")
            self.__check_start_possibility()

    def load_model_for_batch(self):
        """This function loads multiple files as models.

        """
        try:
            # open file dialog
            file_names = Qt.QtWidgets.QFileDialog.getOpenFileNames(self, "Open Models",
                                                                   Qt.QtCore.QDir.homePath(),
                                                                   "PDB Files (*.pdb)")
            if file_names == ([""], ""):
                raise ValueError
            # display path in text box
            for file in file_names[0]:
                self.ui.txt_batch_load_model.append(str(file))
            self.status_bar.showMessage("Loading the models was successful.")
            self.__check_start_possibility_batch()
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the models failed!")
        except ValueError:
            print("No file has been selected.")
            self.__check_start_possibility()

    def enable_chain_information_input_for_batch(self):
        """This function enables the text boxes to enter the chains for the
        reference and the model

        """
        try:
            self.ui.txt_batch_chain_ref.setEnabled(
                self.ui.cb_batch_chain_info.checkState())
            self.ui.txt_batch_chain_model.setEnabled(
                self.ui.cb_batch_chain_info.checkState())
            self.status_bar.showMessage("Enter the chain information.")
            self.__check_start_possibility_batch()
        except Exception:
            self.status_bar.showMessage("Unexpected Error.")

    def check_batch_if_txt_batch_job_name_is_filled(self):
        """This function checks if a job name was entered in the text field.

        """
        self.__check_start_possibility_batch()

    def check_batch_if_txt_batch_chain_ref_is_filled(self):
        """This function checks if any chains are in the text field for the
        reference.

        """
        self.__check_start_possibility_batch()

    def check_batch_if_txt_batch_chain_model_is_filled(self):
        """This function checks if any chains are in the text field for the
        model.

        """
        self.__check_start_possibility_batch()

    def start_process_batch(self):
        """This function contains the main analysis algorithm for the
        Protein structure comparison.

        """
        self.status_bar.showMessage("Checking user input ...")
        job = job_utils.Job(self.ui.txt_batch_job_name.text(), self.workspace_path)
        raw_models_input = self.ui.txt_batch_load_model.toPlainText()
        models = raw_models_input.split("\n")

        # runs analysis with project creation
        for model in models:
            cmd.reinitialize()
            model_file_info = Qt.QtCore.QFileInfo(model)
            MODEL_OBJ_NAME = model_file_info.baseName()
            MODEL_DIR = model_file_info.canonicalPath()

            project = utils.project_utils.Project(f"project_of_{MODEL_OBJ_NAME}",
                                                  f"{self.workspace_path}/{job.get_job_name()}")
            project.create_project_tree()
            project.set_job_name(job.get_job_name())
            project.set_pdb_file(self.ui.txt_batch_load_reference.text())
            project.set_pdb_id(self.ui.txt_batch_load_reference.text())
            project.set_pdb_model(MODEL_OBJ_NAME)
            project.set_ref_chains(self.ui.txt_batch_chain_ref.text())
            project.set_model_chains((self.ui.txt_batch_chain_model.text()))
            project.create_xml_file()
            job.add_project_to_job(project)

            # gets reference filename and filepath
            if len(self.ui.txt_batch_load_reference.text()) == 4:
                tmp_protein = core.Protein(self.ui.txt_batch_load_reference.text(),
                                           export_data_dir=project.get_pdb_path())
                tmp_protein.clean_pdb_file()
                REFERENCE_OBJ_NAME = self.ui.txt_batch_load_reference.text()
                REFERENCE_DIR = project.get_pdb_path()
            else:
                ref_file_info = Qt.QtCore.QFileInfo(self.ui.txt_batch_load_reference.text())
                REFERENCE_OBJ_NAME = ref_file_info.baseName()
                REFERENCE_DIR = ref_file_info.canonicalPath()
            reference_protein: list[core.Protein] = [core.Protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]
            model_proteins: list[core.Protein] = [core.Protein(MODEL_OBJ_NAME, MODEL_DIR)]
            export_dir = project.get_results_path()
            structure_analysis = structure_analysis_utils.StructureAnalysis(
                reference_protein, model_proteins,
                project.get_ref_chains().split(","), project.get_model_chains().split(","),
                export_dir, cycles=global_utils.global_var_settings_obj.get_cycles(),
                cutoff=global_utils.global_var_settings_obj.get_cutoff(),
            )
            structure_analysis.create_selection_for_proteins(structure_analysis.ref_chains,
                                                             structure_analysis.reference_protein)
            structure_analysis.create_selection_for_proteins(structure_analysis.model_chains,
                                                             structure_analysis.model_proteins)
            structure_analysis.do_analysis_in_pymol(structure_analysis.create_protein_pairs(),
                                                    self.status_bar, self.ui.progress_bar_batch)
        job.create_xml_file()

    # Results
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
        global global_var_project_dict
        try:
            file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
        except KeyError:
            tools.quick_log_and_display("error", "No project has been opened.", self.status_bar,
                                        "No project has been opened.")
            return
        pixmap = Qt.QtGui.QPixmap(f"{file_path}/images/structure_alignment.png")
        # TODO: Create setting for min. image size
        pixmap = pixmap.scaled(450, 450, transformMode=PyQt5.QtCore.Qt.SmoothTransformation)
        label.setPixmap(pixmap)
        label.setScaledContents(True)
        png_dialog_layout = QHBoxLayout()
        png_dialog_layout.addWidget(label)
        png_dialog.setLayout(png_dialog_layout)
        png_dialog.setWindowTitle("Image of: structure alignment")
        png_dialog.show()

    def update_distance_plot(self):
        """This function updates the distance plot

        """
        from_aa = int(self.ui.sp_distance_plot_from.text())
        to_aa = int(self.ui.sp_distance_plot_to.text())
        from_range = float(self.ui.dsp_distance_plot_from_range.text())
        to_range = float(self.ui.dsp_distance_plot_to_range.text())
        self.view_box.setRange(xRange=[from_aa, to_aa], yRange=[from_range, to_range])
        # if self.ui.cb_sync_with_pymol.isChecked():
        #     print("Sync is active.")
        #     # TODO: handle objects better
        #     pdb_list = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")
        #     pdb_name = pdb_list[0].replace(".pdb", "")
        #     zoom_selection = f"/{pdb_name}///{from_aa}-{to_aa}/CA"
        #     cmd.select("zoom_sele", zoom_selection)
        #     cmd.zoom("zoom_sele")

    def display_distance_plot(self):
        """This function opens a window which displays the distance plot.

        """
        file_path = f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results"
        model_name = self.ui.lbl_current_project_name.text()
        global_utils.global_var_tmp_project_info.clear()
        global_utils.global_var_tmp_project_info.append(file_path)
        global_utils.global_var_tmp_project_info.append(model_name)

        dialog = dialog_distance_plot.DialogDistancePlot()
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
        #             cutoff_line.append(global_utils.global_var_settings_obj.get_cutoff())
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

        plot_dialog = Qt.QtWidgets.QDialog(self)
        plot_dialog_layout = QHBoxLayout()
        graph_widget = pg.PlotWidget()

        # read csv file
        file_path = f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/results"
        model_name = self.ui.lbl_current_project_name.text()
        path = f"{file_path}/distance_csv/distances.csv"
        distance_list = []
        with open(path, 'r', encoding="utf-8") as csv_file:
            i = 0
            for line in csv_file:
                cleaned_line = line.replace("\n", "")
                if cleaned_line.split(",")[8] != 'distance':
                    distance_list.append(float(cleaned_line.split(",")[8]))
        distance_list.sort()
        length = len(distance_list)
        max_distance = distance_list[length-1]
        x, y = np.histogram(distance_list, bins=np.arange(0, max_distance, 0.25))

        if x.size != y.size:
            x = np.resize(x, (1, y.size))
        # this conversion is needed for the pyqtgraph library!
        x = x.tolist()
        try:
            x = x[0]
        except IndexError:
            # Error got raised where the distances where all 0
            tools.quick_log_and_display("error", "The histogram could not be created.",
                                        self.status_bar, "The histogram could not be created. "
                                                         " Check the distance table!")
            return

        y = y.tolist()
        color = Qt.QtGui.QColor.fromRgb(255, 128, 128)
        # creates bar chart item
        graph_bar_item = pg.BarGraphItem(x0=0, y=y, height=0.2, width=x,
                                         pen=pg.mkPen(color="#4B91F7"), brush=pg.mkBrush(color="#4B91F7"))
        # creates y-labels for bar chart
        y_labels = []
        # TODO: last histogram bar has no label
        for i in range((len(y)-1)):
            label = f"{y[i]} - {y[i+1]}"
            y_labels.append(label)
        y_values = y
        ticks = []
        for i, item in enumerate(y_labels):
            ticks.append((y_values[i], item))
        ticks = [ticks]

        # styling the plot
        graph_widget.setBackground('w')
        graph_widget.setTitle(f"Distance Histogram of {model_name}", size="23pt")
        styles = {'font-size': '14px'}
        ax_label_x = "Distance in angstrom"
        graph_widget.setLabel('left', ax_label_x, **styles)
        graph_widget.setLabel('bottom', "Frequency of α-C atoms distance", **styles)
        graph_widget.addItem(graph_bar_item)
        bar_ax = graph_widget.getAxis('left')
        bar_ax.setTicks(ticks)
        plot_dialog_layout.addWidget(graph_widget)
        plot_dialog.setLayout(plot_dialog_layout)
        plot_dialog.setWindowTitle("Distance Histogram")
        plot_dialog.show()

    def display_interesting_region(self):
        """This function displays an image of an interesting region.

        """
        png_dialog = Qt.QtWidgets.QDialog(self)
        label = Qt.QtWidgets.QLabel(self)
        global global_var_project_dict
        try:
            file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
        except KeyError:
            tools.quick_log_and_display("error", "No project has been opened.", self.status_bar,
                                        "No project has been opened.")
            return
        file_name = self.ui.cb_interesting_regions.currentText()
        pixmap = Qt.QtGui.QPixmap(f"{file_path}/images/interesting_regions/{file_name}")
        # TODO: Create setting for min. image size
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
                  "Model Chain", "Model Pos", "Model Residue", "Distance",
                  ]
        csv_model.setHorizontalHeaderLabels(labels)
        table_dialog = Qt.QtWidgets.QDialog(self)
        table_view = Qt.QtWidgets.QTableView()
        table_view.setModel(csv_model)

        global global_var_project_dict
        try:
            file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
        except KeyError:
            tools.quick_log_and_display("error", "No project has been opened.", self.status_bar,
                                        "No project has been opened.")
            return
        path = f"{file_path}/distance_csv/distances.csv"
        with open(path, 'r', encoding="utf-8") as csv_file:
            i = 0
            for line in csv_file:
                tmp_list = line.split(",")
                tmp_list.pop(0)
                standard_item_list = []
                for tmp in tmp_list:
                    tmp_item = Qt.QtGui.QStandardItem(tmp)
                    standard_item_list.append(tmp_item)
                csv_model.insertRow(i, standard_item_list)
                i += 1
            csv_file.close()
        csv_model.removeRow(0)
        table_view.setAlternatingRowColors(True)
        table_view.resizeColumnsToContents()
        table_view.verticalHeader().setVisible(False)

        table_dialog_layout = QHBoxLayout()
        table_dialog_layout.addWidget(table_view)
        table_dialog.setLayout(table_dialog_layout)
        table_dialog.setWindowTitle("Distances of Structure Alignment")
        table_dialog.show()

    # Image
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
        cmd.scene(key="auto", action="update")

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


if __name__ == '__main__':
    app = QApplication(sys.argv)
    with open('styles/styles_alt.css', 'r', encoding="utf-8") as file:
        style = file.read()
        # Set the stylesheet of the application
        app.setStyleSheet(style)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())
