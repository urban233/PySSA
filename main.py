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

import sys
import os
import shutil
import webbrowser
import time
import logging
import PyQt5.QtCore

import utils.project_utils
import utils.settings_utils
from pathlib import Path
from PyQt5 import QtSvg
from PyQt5.QtWidgets import QHBoxLayout
from pymol import Qt
from pymol import cmd
from utils import structure_analysis_utils
from dialogs import dialog_about
from dialogs import DialogSettingsPdbPreparation
from utils import project_constants
from utils import tools
from utils import gui_utils
from uiForms.auto.auto_main_window import Ui_MainWindow
from pymolproteintools import core

# setup logger
logging.basicConfig(level=logging.DEBUG)
global_var_project_dict = {0: utils.project_utils.Project()}


class MainWindow(Qt.QtWidgets.QMainWindow):
    """This class contains all information about the MainWindow in the
    application

    """
    target_dir = f"{os.path.expanduser('~')}/Documents/data"
    renderer = ""
    SETTINGS = project_constants.SETTINGS

    def __init__(self, *args, **kwargs):
        """Constructor

        Args:
            args
            kwargs
        """
        super().__init__(*args, **kwargs)

        # build ui object
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        # sets up the status bar
        self.status_bar = Qt.QtWidgets.QStatusBar()
        self.setStatusBar(self.status_bar)

        # sets up settings.xml
        if not os.path.exists(self.SETTINGS):
            settings = utils.settings_utils.SettingsXml(self.SETTINGS)
            settings.create_settings_xml_file()
        settings = utils.settings_utils.SettingsXml(self.SETTINGS)
        self.tmp_settings = settings.load_xml_in_memory()

        # safeguard workspace path
        sg_1 = tools.safeguard_filepath_xml(self.tmp_settings, 'workspacePath', 'value')
        # safeguard pdb path
        sg_2 = tools.safeguard_filepath_xml(self.tmp_settings, 'pdbPath', 'value')
        # safeguard zip path
        sg_3 = tools.safeguard_filepath_xml(self.tmp_settings, 'zipPath', 'value')
        # safeguard cycles value
        sg_4 = tools.safeguard_numerical_value_xml(self.tmp_settings, 'cyclesValue', 'value', 'int')
        # safeguard cutoff value
        sg_5 = tools.safeguard_numerical_value_xml(self.tmp_settings, 'cutoffValue', 'value', 'float')

        if sg_1 is False or sg_2 is False or sg_3 is False or sg_4 is False or sg_5 is False:
            self.status_bar.showMessage("The settings.xml is corrupted! Please fix this issue first!")
            gui_utils.error_dialog_settings("The settings.xml is corrupted! Please fix this issue first!",
                                            "Check the log for more info.")
        else:
            self.workspace_path = utils.settings_utils.SettingsXml.get_path(self.tmp_settings,
                                                                            "workspacePath",
                                                                            "value")
            self.workspace = Qt.QtWidgets.QLabel(f"Current Workspace: {self.workspace_path}")

        # sets up defaults
        # Prediction + Analysis
        self.ui.txt_prediction_chain_ref.setEnabled(False)
        self.ui.txt_prediction_chain_model.setEnabled(False)
        self.ui.btn_prediction_start.setEnabled(False)
        self.ui.btn_prediction_cancel.setEnabled(False)
        self.ui.progress_bar_prediction.setProperty("value", 0)

        # Single Analysis
        self.ui.txt_analysis_chain_ref.setEnabled(False)
        self.ui.txt_analysis_chain_model.setEnabled(False)
        self.ui.btn_analysis_start.setEnabled(False)
        self.ui.progress_bar_analysis.setProperty("value", 0)
        # Batch
        self.ui.txt_batch_chain_ref.setEnabled(False)
        self.ui.txt_batch_chain_model.setEnabled(False)
        self.ui.btn_batch_start.setEnabled(False)
        self.ui.progress_bar_batch.setProperty("value", 0)

        # fills combo boxes
        # combo box Representation
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

        # connect elements with function
        # menu connections
        self.ui.action_file_open.triggered.connect(self.open)
        self.ui.action_file_save_as.triggered.connect(self.save_as)
        self.ui.action_file_quit.triggered.connect(self.quit_app)
        self.ui.action_file_restore_settings.triggered.connect(self.restore_settings)
        self.ui.action_wizard_prepare_model_pdbs.triggered.connect(
            self.prepare_model_pdbs)
        self.ui.actSettingsPdbPreparation.triggered.connect(
            self.open_settings_pdb_preparation)
        self.ui.action_settings_edit_all.triggered.connect(self.open_settings_global)
        self.ui.action_display_current_workspace.triggered.connect(self.display_workspace_path)
        self.ui.action_display_project_path.triggered.connect(self.display_project_path)
        self.ui.action_help_docs.triggered.connect(self.open_documentation)
        self.ui.action_help_docs_pdf.triggered.connect(self.open_documentation_pdf)
        self.ui.action_help_about.triggered.connect(self.open_about)

        # Prediction + Analysis
        # buttons connections
        self.ui.btn_prediction_load_reference.clicked.connect(self.load_reference_for_prediction)
        self.ui.btn_prediction_start.clicked.connect(self.predict)

        # checkbox
        self.ui.cb_prediction_chain_info.stateChanged.connect(
            self.enable_chain_information_input_for_prediction)

        # text fields connections
        self.ui.txt_prediction_project_name.textChanged.connect(
            self.check_prediction_if_txt_prediction_project_name_is_filled)
        self.ui.txt_prediction_load_reference.textChanged.connect(
            self.check_prediction_if_txt_prediction_load_reference_is_filled)

        # Home/Single Analysis
        # buttons connections
        self.ui.btn_analysis_load_reference.clicked.connect(self.load_reference_for_analysis)
        self.ui.btn_analysis_load_model.clicked.connect(self.load_model_for_analysis)
        self.ui.btn_analysis_start.clicked.connect(self.start_process)

        # checkbox connections
        self.ui.cb_analysis_chain_info.stateChanged.connect(
            self.enable_chain_information_input_for_analysis)

        # text fields connections
        self.ui.txt_analysis_project_name.textChanged.connect(
            self.check_analysis_if_txt_analysis_project_name_is_filled)
        self.ui.txt_analysis_chain_ref.textChanged.connect(
            self.check_analysis_if_txt_analysis_chain_ref_is_filled)
        self.ui.txt_analysis_chain_model.textChanged.connect(
            self.check_analysis_if_txt_analysis_chain_model_is_filled)
        self.ui.txt_prediction_chain_ref.textChanged.connect(
            self.check_prediction_if_txt_prediction_chain_ref_is_filled)
        self.ui.txt_prediction_chain_model.textChanged.connect(
            self.check_prediction_if_txt_prediction_chain_model_is_filled)

        # Home/Batch
        # buttons connections
        self.ui.btn_batch_load_reference.clicked.connect(self.load_reference_for_batch)
        self.ui.btn_batch_load_model.clicked.connect(self.load_model_for_batch)
        self.ui.btn_batch_start.clicked.connect(self.start_process_batch)

        # checkbox connections
        self.ui.cb_batch_chain_info.stateChanged.connect(
            self.enable_chain_information_input_for_batch)

        # text fields connections
        self.ui.txt_batch_job_name.textChanged.connect(
            self.check_batch_if_txt_batch_job_name_is_filled)
        self.ui.txt_batch_chain_ref.textChanged.connect(
            self.check_batch_if_txt_batch_chain_ref_is_filled)
        self.ui.txt_batch_chain_model.textChanged.connect(
            self.check_batch_if_txt_batch_chain_model_is_filled)

        # Results
        # buttons connections
        self.ui.btn_view_struct_alignment.clicked.connect(self.display_structure_alignment)
        self.ui.btn_view_distance_plot.clicked.connect(self.display_distance_plot)
        self.ui.btn_view_distance_histogram.clicked.connect(self.display_distance_histogram)
        self.ui.btn_view_distance_table.clicked.connect(self.display_distance_table)
        self.ui.btn_view_interesting_region.clicked.connect(self.display_interesting_region)

        # Image
        # buttons connections
        self.ui.btn_save_image.clicked.connect(self.save_image)
        self.ui.btn_preview_image.clicked.connect(self.preview_image)

        # combo box connections
        self.ui.box_representation.activated.connect(self.show_representation)
        self.ui.box_bg_color.activated.connect(self.choose_bg_color)
        self.ui.box_renderer.activated.connect(self.choose_renderer)
        self.ui.box_ray_trace_mode.activated.connect(self.choose_ray_trace_mode)
        self.ui.box_ray_texture.activated.connect(self.choose_ray_texture)

        # checkbox connections
        self.ui.cb_transparent_bg.stateChanged.connect(self.decide_transparent_bg)

        # creates tooltips
        # Home/Single Analysis
        # for buttons
        self.ui.btn_analysis_load_reference.setToolTip("Open reference pdb file")
        self.ui.btn_analysis_load_model.setToolTip("Open model pdb file")
        self.ui.btn_analysis_start.setToolTip("Start analysis process")

        # for text fields
        self.ui.txt_analysis_load_reference.setToolTip("Reference file path")
        self.ui.txt_analysis_load_model.setToolTip("Model file path")
        self.ui.txt_analysis_chain_ref.setToolTip("Enter chain(s) of reference")
        self.ui.txt_analysis_chain_model.setToolTip("Enter chain(s) of model")

        # for checkbox
        self.ui.cb_analysis_chain_info.setToolTip("Enable input of chains")

        # for statusbar
        self.status_bar.setToolTip("Status information: Current process")

        # Home/Batch tab
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

        # Image tab
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

        # Shortcuts menubar

        # setting additional parameters
        self.setWindowTitle("PySSA v0.1.0")

    def __check_start_possibility(self):
        """This function is used to determine if the Start button can be
        enabled for the single analysis.

        """
        i = 0
        j = 0
        if len(str(self.ui.txt_analysis_project_name.text())) > 0:
            i += 1
        if len(str(self.ui.txt_analysis_load_reference.text())) > 0:
            i += 1
        if len(str(self.ui.txt_analysis_load_model.toPlainText())) > 0:
            i += 1
        if self.ui.cb_analysis_chain_info.isChecked() and len(self.ui.txt_analysis_chain_ref.text()) > 0 and len(
                self.ui.txt_analysis_chain_model.text()) > 0:
            i += 1
        if self.ui.cb_analysis_chain_info.isChecked():
            j = 4
        else:
            j = 3
        if i == j:
            self.ui.btn_analysis_start.setEnabled(True)
            with open('styles/styles_start_button_ready.css', 'r') as file:
                button_style = file.read()

                # Set the stylesheet of the application
                self.ui.btn_analysis_start.setStyleSheet(button_style)
        else:
            self.ui.btn_analysis_start.setEnabled(False)
            with open('styles/styles_start_button_not_ready.css', 'r') as file:
                button_style = file.read()

                # Set the stylesheet of the application
                self.ui.btn_analysis_start.setStyleSheet(button_style)

    def __check_start_possibility_batch(self):
        """This function is used to determine if the Start button can be
        enabled for the single analysis.

        """
        i = 0
        j = 0
        if len(str(self.ui.txt_batch_job_name.text())) > 0:
            i += 1
        if len(str(self.ui.txt_batch_load_reference.text())) > 0:
            i += 1
        if len(str(self.ui.txt_batch_load_model.toPlainText())) > 0:
            i += 1
        if self.ui.cb_batch_chain_info.isChecked() and len(self.ui.txt_batch_chain_ref.text()) > 0 and len(
                self.ui.txt_batch_chain_model.text()) > 0:
            i += 1
        if self.ui.cb_batch_chain_info.isChecked():
            j = 4
        else:
            j = 3
        if i == j:
            self.ui.btn_batch_start.setEnabled(True)
            with open('styles/styles_start_button_ready.css', 'r') as file:
                button_style = file.read()

                # Set the stylesheet of the application
                self.ui.btn_batch_start.setStyleSheet(button_style)
        else:
            self.ui.btn_batch_start.setEnabled(False)
            with open('styles/styles_start_button_not_ready.css', 'r') as file:
                button_style = file.read()

                # Set the stylesheet of the application
                self.ui.btn_batch_start.setStyleSheet(button_style)

    def __check_start_possibility_prediction(self):
        """This function is used to determine if the Start button can be
        enabled for the single analysis.

        """
        i = 0
        j = 0
        if len(str(self.ui.txt_prediction_project_name.text())) > 0:
            i += 1
        if len(str(self.ui.txt_prediction_load_reference.text())) > 0:
            i += 1
        if self.ui.cb_prediction_chain_info.isChecked() and len(self.ui.txt_prediction_chain_ref.text()) > 0 and len(
                self.ui.txt_prediction_chain_model.text()) > 0:
            i += 1
        if self.ui.cb_prediction_chain_info.isChecked():
            j = 3
        else:
            j = 2
        if i == j:
            self.ui.btn_prediction_start.setEnabled(True)
            with open('styles/styles_start_button_ready.css', 'r') as file:
                button_style = file.read()

                # Set the stylesheet of the application
                self.ui.btn_prediction_start.setStyleSheet(button_style)
        else:
            self.ui.btn_prediction_start.setEnabled(False)
            with open('styles/styles_start_button_not_ready.css', 'r') as file:
                button_style = file.read()

                # Set the stylesheet of the application
                self.ui.btn_prediction_start.setStyleSheet(button_style)

    # @SLOT
    # Menu
    def open(self):
        """This function opens a project.xml file and fills the input boxes with the right values

        """
        # TODO: * it's broken and needs to be fixed
        attribute = "value"
        # open file dialog
        try:
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self,
                                                                 "Open project file",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "Plugin Project File (*.xml)")
            if file_name == ("", ""):
                # TODO: add status bar and logger message
                print("No file has been selected.")
                return
            global global_var_project_dict
            global_var_project_dict[0].load_project(file_name[0])
            # add filename to project list (results tab)
            self.ui.project_list.addItem(global_var_project_dict[0].get_project_name())
            self.ui.project_list.setCurrentRow(0)
            # fill combo box of interesting regions
            results_path = global_var_project_dict[0].get_results_path()
            dir_content = os.listdir(f"{results_path}/images/interesting_regions")
            for tmp_file in dir_content:
                self.ui.cb_interesting_regions.addItem(tmp_file)

            #
            # project_file = tools.ProjectXml(file_name[0])
            # tmp_project_file = project_file.load_xml_in_memory()
            #
            # if project_file.get_path(tmp_project_file, "predictionDone", attribute) == "False":
            #     self.ui.txt_analysis_project_name.setText(
            #         project_file.get_path(tmp_project_file,
            #                               "projectName",
            #                               attribute))
            #     self.ui.txt_analysis_load_reference.setText(
            #         project_file.get_path(tmp_project_file,
            #                               "reference",
            #                               attribute))
            #     if len(project_file.get_path(tmp_project_file, "referenceChains",
            #                                  attribute)) > 0:
            #         self.ui.cb_analysis_chain_info.setChecked(True)
            #         self.ui.txt_analysis_chain_ref.setText(project_file.get_path(
            #             tmp_project_file, "referenceChains", attribute))
            #         self.ui.txt_analysis_chain_model.setText(project_file.get_path(
            #             tmp_project_file, "modelChains", attribute))
            #     results_path = project_file.get_path(tmp_project_file, "results", attribute)
            #     session_file_name = ""
            #     for file in os.listdir(f"{results_path}/sessions"):
            #         session_file_name = file
            #     cmd.load(f"{results_path}/sessions/{session_file_name}")
            # else:
            #     self.ui.txt_prediction_project_name.setText(project_file.get_path(tmp_project_file,
            #                                                                       "projectName",
            #                                                                       attribute))
            #     self.ui.txt_prediction_load_reference.setText(project_file.get_path(tmp_project_file,
            #                                                                         "reference",
            #                                                                         attribute))
            #     if len(project_file.get_path(tmp_project_file, "referenceChains",
            #                                  attribute)) > 0:
            #         self.ui.cb_prediction_chain_info.setChecked(True)
            #         self.ui.txt_prediction_chain_ref.setText(project_file.get_path(
            #             tmp_project_file, "referenceChains", attribute))
            #         self.ui.txt_prediction_chain_model.setText(project_file.get_path(
            #             tmp_project_file, "modelChains", attribute))
            #     results_path = project_file.get_path(tmp_project_file, "results", attribute)
            #     session_file_name = ""
            #     for file in os.listdir(f"{results_path}/sessions"):
            #         session_file_name = file
            #     cmd.load(f"{results_path}/sessions/{session_file_name}")
        except FileNotFoundError:
            print("File could not be opened.")
        except ValueError:
            print("No file has been selected.")

    def save_as(self):
        """This function saves all information in the current tab

        """
        # TODO: fix saving process
        gui_utils.save_project_xml(self, self.status_bar)

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

    def prepare_model_pdbs(self):
        """This function extracts and moves the .pdb files directly from the
        download directory to another directory

        """
        try:
            # path to download directory
            tmp_path = f"{os.path.expanduser('~')}/Downloads"
            os.chdir(tmp_path)
            # path to store all pdb files (should be defined in the Settings tab!)
            target_dir = f"{os.path.expanduser('~')}/Documents/data"

            # tmp_list contains a list with all zips starting with prediction
            tmp_list = tools.filter_prediction_zips(tmp_path)

            # create temporary dir to handle unzipping
            tmp_dir = f"{tmp_path}/tmp"
            if not os.path.exists(tmp_dir):
                os.mkdir(tmp_dir)

            for file in tmp_list:
                tools.extract_and_move_model_pdb(tmp_path, tmp_dir, file, target_dir)

            shutil.rmtree(tmp_dir)
            self.status_bar.showMessage("Preparing the .pdb files was successful.")
        except Exception:
            self.status_bar.showMessage("Preparing the .pdb files failed!")

    @staticmethod
    def open_settings_global():
        """This function open the dialog for the global settings.

        """
        tools.open_global_settings()

    @staticmethod
    def open_settings_pdb_preparation():
        """This function opens the dialog for the pdb Preparation settings.

        """
        dialog = DialogSettingsPdbPreparation.DialogSettingsPdbPreparation()
        dialog.exec_()

    def display_workspace_path(self):
        """This function displays the current workspace path in the statusbar of the main window.

        """
        self.status_bar.showMessage(f"Current Workspace: {self.workspace_path}")

    def display_project_path(self):
        """This function displays the current project path in the statusbar of the main window.

        """
        try:
            global global_var_project_dict
            if global_var_project_dict[0].get_project_path() == "":
                raise ValueError
            project_path = global_var_project_dict[0].get_project_path()
            self.status_bar.showMessage(f"Current Project Path: {project_path}")
        except ValueError:
            tools.quick_log_and_display("info", "No project loaded yet.", self.status_bar,
                                        "No project loaded yet.")

    @staticmethod
    def open_documentation():
        """This function opens the official plugin documentation as HTML page.

        """
        webbrowser.open_new(f"file://{os.getcwd()}/docs/pymol_plugin/build/html/index.html")

    @staticmethod
    def open_documentation_pdf():
        """This function opens the official plugin documentation as PDF.

        """
        webbrowser.open_new(
            f"file://{os.getcwd()}/docs/pymol_plugin/build/latex/pyssa-python-pluginforsequencetostructureanalysis.pdf")

    @staticmethod
    def open_about():
        """This function opens the about dialog.

        """
        dialog = dialog_about.DialogAbout()
        dialog.exec_()

    # Prediction + Analysis
    def load_reference_for_prediction(self):
        """This function loads a reference for the analysis part

        """
        try:
            # open file dialog
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
            if file_name == ("", ""):
                raise ValueError("No file has been selected.")
            # display path in text box
            self.ui.txt_prediction_load_reference.setText(str(file_name[0]))
            self.status_bar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_prediction()
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the reference failed!")
        except ValueError:
            print("No file has been selected.")

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

    def check_prediction_if_txt_prediction_chain_ref_is_filled(self):
        """This function checks if any chains are in the text field for the
        reference.

        """
        self.__check_start_possibility_prediction()

    def check_prediction_if_txt_prediction_chain_model_is_filled(self):
        """This function checks if any chains are in the text field for the
        model.

        """
        self.__check_start_possibility_prediction()

    def check_prediction_if_txt_prediction_project_name_is_filled(self):
        """This function checks if the project name text field is filled.

        """
        self.__check_start_possibility_prediction()

    def check_prediction_if_txt_prediction_load_reference_is_filled(self):
        """This function checks if a reference pdb file is selected.

        """
        self.__check_start_possibility_prediction()

    def predict(self):
        """This function opens a webbrowser with a colab notebook, to run the
        prediction. In addition, it runs the entire analysis after the
        prediction.

        """
        self.status_bar.removeWidget(self.workspace)
        self.status_bar.showMessage("Checking user input ...")

        # check if a prediction is already finished
        if os.path.isfile(project_constants.FULL_FILENAME_PREDICTION_ZIP):
            self.status_bar.showMessage(
                f"Warning! | Current Workspace: {self.workspace_path}")
            check = gui_utils.warning_message_prediction_exists(
                f"The prediction is here: {project_constants.FULL_FILENAME_PREDICTION_ZIP} ",
                project_constants.FULL_FILENAME_PREDICTION_ZIP)
            if not check:
                return

        # project = utils.project_utils.Project()
        # project_name = self.ui.txt_prediction_project_name.text()
        # project_name_with_underscores = project_name.replace(" ", "_")
        # project_folder_path = Path(f"{self.workspace_path}/{project_name_with_underscores}")
        # self.__create_project_folder(project_folder_path)
        #
        # # setting values for the project xml file
        # project.set_project_name(project_name_with_underscores)
        # # TODO: * implement input check with QValidator and Regex
        # if len(self.ui.txt_prediction_load_reference.text()) == 4:
        #     project.set_pdb_id(self.ui.txt_prediction_load_reference.text())
        # project.set_pdb_file(self.ui.txt_prediction_load_reference.text())
        # if self.ui.txt_prediction_chain_ref.text() != "":
        #     project.set_ref_chains(self.ui.txt_prediction_chain_ref.text())
        # if self.ui.txt_prediction_chain_model.text() != "":
        #     project.set_model_chains(self.ui.txt_prediction_chain_model.text())
        # project.set_results_path(f"{project_folder_path}/results")
        #
        # # creates a pdb folder
        # self.__create_directory(Path(f"{project_folder_path}"), "pdb")
        # project_pdb_path = f"{project_folder_path}/pdb"
        # # creates a tmp folder to unzip the prediction
        # self.__create_directory(Path(f"{project_folder_path}"), "tmp")
        # project_tmp_path = f"{project_folder_path}/tmp"
        # # creates results folder
        # self.__create_directory(Path(f"{project_folder_path}"), "results")
        # project_results_path = f"{project_folder_path}/results"

        # gets reference filename and filepath
        if len(self.ui.txt_prediction_load_reference.text()) == 4:
            tmp_protein = core.protein(self.ui.txt_prediction_load_reference.text(),
                                       export_data_dir=project_pdb_path)
            tmp_protein.clean_pdb_file()
            REFERENCE_OBJ_NAME = self.ui.txt_prediction_load_reference.text()
            REFERENCE_DIR = project_pdb_path
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

        # waiting for the colab notebook to finish
        # TODO:
        #   * implement cancel button for "Abort Prediction"
        self.ui.btn_prediction_cancel.setEnabled(True)
        archive = "prediction.zip"
        source_path = Path(f"{os.path.expanduser('~')}/Downloads")
        FILE_NAME = source_path / archive
        # flag = False
        # while flag == False:
        #     print("AlphaFold is still running ...")
        #     time.sleep(5)
        #     # time.sleep(120)
        while os.path.isfile(FILE_NAME) is False:
            print("AlphaFold is still running ...")
            # time.sleep(5)
            time.sleep(20)
            # time.sleep(120)

        # ----------------------------------------------------------------- #
        # start of the analysis algorithm
        self.status_bar.showMessage("Protein structure analysis started ...")

        # extracts and moves the prediction.pdb to the workspace/pdb folder
        tools.extract_and_move_model_pdb(
            str(source_path), project_tmp_path, archive, project_pdb_path)

        # removes tmp dir
        # os.remove(project_tmp_path)

        # gets model filename and filepath
        PREDICTION_NAME = tools.get_prediction_file_name(project_pdb_path)
        full_model_file_path = f"{project_pdb_path}/{PREDICTION_NAME[0]}"
        model_file_info = Qt.QtCore.QFileInfo(full_model_file_path)
        MODEL_OBJ_NAME = model_file_info.baseName()
        MODEL_DIR = model_file_info.canonicalPath()

        project = utils.project_utils.Project()
        try:
            project_pdb_path, project_results_path = project.create_project(self.workspace_path, full_model_file_path,
                                                                            MODEL_OBJ_NAME,
                                                                            self.ui.txt_prediction_project_name,
                                                                            self.ui.txt_prediction_load_reference,
                                                                            self.ui.txt_prediction_chain_ref,
                                                                            self.ui.txt_prediction_chain_model,
                                                                            "")
        except IsADirectoryError:
            return
        except TypeError:
            return
        # create the protein object for the reference
        reference_protein: list[core.protein] = [core.protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]

        # create model protein object
        model_proteins: list[core.protein] = [core.protein(MODEL_OBJ_NAME, MODEL_DIR)]
        # sets the filepath of the model in the project xml file
        project.set_pdb_models([full_model_file_path])
        project.create_xml_file(f"{str(project_folder_path)}/{project_name_with_underscores}.xml")
        export_dir = project_results_path
        structure_analysis = structure_analysis_utils.StructureAnalysis(reference_protein, model_proteins,
                                                                        project.get_ref_chains().split(","),
                                                                        project.get_model_chains().split(","),
                                                                        export_dir)
        structure_analysis.create_selection_for_proteins(structure_analysis.ref_chains,
                                                         structure_analysis.reference_protein)
        structure_analysis.create_selection_for_proteins(structure_analysis.model_chains,
                                                         structure_analysis.model_proteins)
        structure_analysis.do_analysis_in_pymol(structure_analysis.create_protein_pairs(),
                                                self.status_bar, self.ui.progress_bar_prediction)

    # Single Analysis
    def load_reference_for_analysis(self):
        """This function opens a file dialog to choose a .pdb file as
        reference and displays the path in a text box

        """
        try:
            # open file dialog
            fileName = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                Qt.QtCore.QDir.homePath(),
                                                                "PDB Files (*.pdb)")
            # display path in text box
            if fileName == ("", ""):
                raise ValueError("No file has been selected.")
            # display path in text box
            self.ui.txt_analysis_load_reference.setText(str(fileName[0]))
            self.status_bar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_prediction()
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the reference failed!")
        except ValueError:
            print("No file has been selected.")
            self.__check_start_possibility()
        except Exception:
            self.status_bar.showMessage("Loading the reference failed!")

    def load_model_for_analysis(self):
        """This function opens a file dialog to choose a .pdb file as
        model and displays the path in a text box

        """
        try:
            # open file dialog
            fileName = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Model",
                                                                Qt.QtCore.QDir.homePath(),
                                                                "PDB Files (*.pdb)")
            # display path in text box
            self.ui.txt_analysis_load_model.setText(str(fileName[0]))
            self.status_bar.showMessage("Loading the model was successful.")
            self.__check_start_possibility()
        except Exception:
            self.status_bar.showMessage("Loading the model failed!")

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
        protein structure comparison.

        """
        self.status_bar.removeWidget(self.workspace)
        self.ui.btn_analysis_start.setEnabled(False)
        self.status_bar.showMessage("Protein structure analysis started ...")

        model_full_filepath = self.ui.txt_analysis_load_model.toPlainText()
        model_file_info = Qt.QtCore.QFileInfo(model_full_filepath)
        MODEL_OBJ_NAME = model_file_info.baseName()
        MODEL_DIR = model_file_info.canonicalPath()

        project = utils.project_utils.Project()
        try:
            project_pdb_path, project_results_path = project.create_project(self.workspace_path, model_full_filepath,
                                                                            MODEL_OBJ_NAME,
                                                                            self.ui.txt_analysis_project_name,
                                                                            self.ui.txt_analysis_load_reference,
                                                                            self.ui.txt_analysis_chain_ref,
                                                                            self.ui.txt_analysis_chain_model,
                                                                            "")
        except IsADirectoryError:
            return
        except TypeError:
            return
        # gets reference filename and filepath
        if len(self.ui.txt_analysis_load_reference.text()) == 4:
            tmp_protein = core.protein(self.ui.txt_analysis_load_reference.text(),
                                       export_data_dir=project_pdb_path)
            tmp_protein.clean_pdb_file()
            REFERENCE_OBJ_NAME = self.ui.txt_analysis_load_reference.text()
            REFERENCE_DIR = project_pdb_path
        else:
            ref_file_info = Qt.QtCore.QFileInfo(self.ui.txt_analysis_load_reference.text())
            REFERENCE_OBJ_NAME = ref_file_info.baseName()
            REFERENCE_DIR = ref_file_info.canonicalPath()

        reference_protein: list[core.protein] = [core.protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]
        model_proteins: list[core.protein] = [core.protein(MODEL_OBJ_NAME, MODEL_DIR)]
        export_dir = project_results_path
        structure_analysis = structure_analysis_utils.StructureAnalysis(reference_protein, model_proteins,
                                                                        project.get_ref_chains().split(","),
                                                                        project.get_model_chains().split(","),
                                                                        export_dir)
        structure_analysis.create_selection_for_proteins(structure_analysis.ref_chains,
                                                         structure_analysis.reference_protein)
        structure_analysis.create_selection_for_proteins(structure_analysis.model_chains,
                                                         structure_analysis.model_proteins)
        structure_analysis.do_analysis_in_pymol(structure_analysis.create_protein_pairs(),
                                                self.status_bar, self.ui.progress_bar_analysis)

    # Batch
    def load_reference_for_batch(self):
        """This function opens a file dialog to choose a .pdb file as
        reference and displays the path in a text box

        """
        try:
            # open file dialog
            fileName = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                Qt.QtCore.QDir.homePath(),
                                                                "PDB Files (*.pdb)")
            # display path in text box
            if fileName == ("", ""):
                raise ValueError("No file has been selected.")
            # display path in text box
            self.ui.txt_batch_load_reference.setText(str(fileName[0]))
            self.status_bar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_batch()
        except Exception:
            self.status_bar.showMessage("Loading the reference failed!")

    def load_model_for_batch(self):
        """This function loads multiple files as models.

        """
        try:
            # open file dialog
            file_names = Qt.QtWidgets.QFileDialog.getOpenFileNames(self, "Open Models", self.target_dir,
                                                                   "PDB Files (*.pdb)")
            # display path in text box
            for file in file_names[0]:
                self.ui.txt_batch_load_model.append(str(file))
            self.status_bar.showMessage("Loading the models was successful.")
            self.__check_start_possibility_batch()
        except Exception:
            self.status_bar.showMessage("Loading the models failed!")

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
        protein structure comparison.

        """
        self.status_bar.removeWidget(self.workspace)
        self.status_bar.showMessage("Checking user input ...")

        # creates job folder
        job = utils.project_utils.Project()
        job_name = self.ui.txt_batch_job_name.text()
        job_name_with_underscores = job_name.replace(" ", "_")
        job_folder_path = Path(f"{self.workspace_path}/{job_name_with_underscores}")
        self.__create_project_folder(job_folder_path)

        job.set_job_name(job_name_with_underscores)
        raw_models_input = self.ui.txt_batch_load_model.toPlainText()
        models = raw_models_input.split("\n")
        job.set_pdb_models(models)
        job.create_xml_file(f"{str(job_folder_path)}/{job_name_with_underscores}.xml")

        for model in models:
            model_file_info = Qt.QtCore.QFileInfo(model)
            MODEL_OBJ_NAME = model_file_info.baseName()
            MODEL_DIR = model_file_info.canonicalPath()

            project = utils.project_utils.Project()
            try:
                project_pdb_path, project_results_path = project.create_project(self.workspace_path, model, MODEL_OBJ_NAME,
                                                                                "", self.ui.txt_batch_load_reference,
                                                                                self.ui.txt_batch_chain_ref,
                                                                                self.ui.txt_batch_chain_model,
                                                                                job_name_with_underscores)
            except IsADirectoryError:
                return
            except TypeError:
                return
            # gets reference filename and filepath
            if len(self.ui.txt_batch_load_reference.text()) == 4:
                tmp_protein = core.protein(self.ui.txt_batch_load_reference.text(),
                                           export_data_dir=project_pdb_path)
                tmp_protein.clean_pdb_file()
                REFERENCE_OBJ_NAME = self.ui.txt_batch_load_reference.text()
                REFERENCE_DIR = project_pdb_path
            else:
                ref_file_info = Qt.QtCore.QFileInfo(self.ui.txt_batch_load_reference.text())
                REFERENCE_OBJ_NAME = ref_file_info.baseName()
                REFERENCE_DIR = ref_file_info.canonicalPath()
            reference_protein: list[core.protein] = [core.protein(REFERENCE_OBJ_NAME, REFERENCE_DIR)]
            model_proteins: list[core.protein] = [core.protein(MODEL_OBJ_NAME, MODEL_DIR)]
            export_dir = project_results_path
            structure_analysis = structure_analysis_utils.StructureAnalysis(reference_protein, model_proteins,
                                                                            project.get_ref_chains().split(","),
                                                                            project.get_model_chains().split(","),
                                                                            export_dir)
            structure_analysis.create_selection_for_proteins(structure_analysis.ref_chains,
                                                             structure_analysis.reference_protein)
            structure_analysis.create_selection_for_proteins(structure_analysis.model_chains,
                                                             structure_analysis.model_proteins)
            # TODO: * fix analysis with chains (hint: error selection (<pymolproteintools object>))
            structure_analysis.do_analysis_in_pymol(structure_analysis.create_protein_pairs(),
                                                    self.status_bar, self.ui.progress_bar_batch)

    # Results
    def display_structure_alignment(self):
        """This function opens a window which displays the image of the structure alignment.

        """
        png_dialog = Qt.QtWidgets.QDialog(self)
        label = Qt.QtWidgets.QLabel(self)
        global global_var_project_dict
        file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
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

    def display_distance_plot(self):
        """This function opens a window which displays the distance plot.

        """
        item = self.ui.project_list.selectedItems()
        if item is None:
            raise ValueError

        svg_dialog = Qt.QtWidgets.QDialog(self)
        svg_dialog_layout = QHBoxLayout()
        viewer = QtSvg.QSvgWidget()
        global global_var_project_dict
        model_file_name = global_var_project_dict[self.ui.project_list.currentRow()].get_model_filename()
        file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
        viewer.load(f"{file_path}/plots/distance_plot/distance_plot_{model_file_name}.svg")
        viewer.show()
        svg_dialog.setWindowTitle(f"Distance Plot of {model_file_name}")
        svg_dialog_layout.addWidget(viewer)
        svg_dialog.setLayout(svg_dialog_layout)
        svg_dialog.show()

    def display_distance_histogram(self):
        """This function opens a window which displays the distance histogram.

        """
        item = self.ui.project_list.selectedItems()
        if item is None:
            raise ValueError

        svg_dialog = Qt.QtWidgets.QDialog(self)
        svg_dialog_layout = QHBoxLayout()
        viewer = QtSvg.QSvgWidget()
        global global_var_project_dict
        model_file_name = global_var_project_dict[self.ui.project_list.currentRow()].get_model_filename()
        file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
        viewer.load(f"{file_path}/plots/distance_histogram/distance_histogram_{model_file_name}.svg")
        viewer.show()
        svg_dialog.setWindowTitle(f"Distance Histogram of {model_file_name}")
        svg_dialog_layout.addWidget(viewer)
        svg_dialog.setLayout(svg_dialog_layout)
        svg_dialog.show()

    def display_interesting_region(self):
        """This function displays an image of an interesting region.

        """
        png_dialog = Qt.QtWidgets.QDialog(self)
        label = Qt.QtWidgets.QLabel(self)
        global global_var_project_dict
        file_path = global_var_project_dict[self.ui.project_list.currentRow()].get_results_path()
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
        table_dialog = Qt.QtWidgets.QDialog(self)
        table_view = Qt.QtWidgets.QTableView()

        table_dialog_layout = QHBoxLayout()
        table_dialog_layout.addWidget(table_view)
        table_dialog.setLayout(table_dialog_layout)
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

    def update_scene(self):
        """This function updates the current selected PyMOL scene.

        TODO: create function
        """
        return

    def save_scene(self):
        """This function saves the current view as a new PyMOL scene.

        TODO: create function
        """
        return

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
            full_file_name = save_dialog.getSaveFileName(caption="Save Image",
                                                         filter="Image (*.png)")
            self.status_bar.showMessage("Creating ray-traced image ...")
            cmd.ray(2400, 2400, renderer=int(self.renderer))
            cmd.png(full_file_name[0], dpi=300)
            self.status_bar.showMessage("Finished image creation.")
        else:
            save_dialog = Qt.QtWidgets.QFileDialog()
            full_file_name = save_dialog.getSaveFileName(caption="Save Image",
                                                         filter="Image (*.png)")
            self.status_bar.showMessage("Creating draw image ...")
            cmd.draw(2400, 2400)
            cmd.png(full_file_name[0], dpi=300)
            self.status_bar.showMessage("Finished image creation.")


# comment out if run within pymol
if __name__ == '__main__':
    # Start point of the application
    app = Qt.QtWidgets.QApplication(sys.argv)
    # Open the qss styles file and read in the css-alike styling code
    with open('styles/styles.css', 'r', encoding="utf-8") as file:
        style = file.read()

        # Set the stylesheet of the application
        app.setStyleSheet(style)
    mainWindow = MainWindow()
    mainWindow.show()
    app.exec_()
