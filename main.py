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

from PyQt5 import QtSvg
from PyQt5.QtWidgets import QHBoxLayout
from pymol import Qt

import sys
import os
import shutil
import webbrowser
import time
from dialogs import dialog_image
from dialogs import dialog_settings_global
from dialogs import dialog_about
from dialogs import dialog_finished
from dialogs import DialogSettingsPdbPreparation
from dialogs import DialogWarningPredictionProject
from dialogs import DialogWarningPredictionZip
from utils import constants
from utils import tools
from utils import gui_utils
from uiForms.auto.auto_main_window import Ui_MainWindow
from pymolproteintools import core
from pymolproteintools import graphics
import matplotlib.pyplot as plt
from pymol import cmd


class MainWindow(Qt.QtWidgets.QMainWindow):
    """This class contains all information about the MainWindow in the
    application

    """
    target_dir = f"{os.path.expanduser('~')}/Documents/data"
    renderer = ""
    SETTINGS_DIR = constants.SETTINGS_DIR

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

        # sets up settings.xml
        if not os.path.exists(self.SETTINGS_DIR):
            settings = tools.SettingsXml(self.SETTINGS_DIR)
            settings.create_settings_xml_file()
        settings = tools.SettingsXml(self.SETTINGS_DIR)
        self.tmpSettings = settings.load_xml_in_memory()

        # sets up the status bar
        self.statusBar = Qt.QtWidgets.QStatusBar()
        self.setStatusBar(self.statusBar)
        self.workspacePath = tools.SettingsXml.get_path(self.tmpSettings,
                                                        "workspacePath",
                                                        "value")
        self.workspace = Qt.QtWidgets.QLabel(f"Current Workspace: {self.workspacePath}")
        # self.statusBar.addPermanentWidget(self.workspace)

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
        COMBO_BOX_DEFAULT = ""
        # combo box Representation
        self.ui.box_representation.addItem(COMBO_BOX_DEFAULT)
        self.ui.box_representation.addItem("cartoon")
        self.ui.box_representation.addItem("ribbon")
        # combo box BgColor
        self.ui.box_bg_color.addItem(COMBO_BOX_DEFAULT)
        self.ui.box_bg_color.addItem("black")
        self.ui.box_bg_color.addItem("white")
        # combo box Renderer
        self.ui.box_renderer.addItem(COMBO_BOX_DEFAULT)
        self.ui.box_renderer.addItem("default renderer")
        self.ui.box_renderer.addItem("PyMOL internal renderer")
        # combo box RayTraceMode
        self.ui.box_ray_trace_mode.addItem(COMBO_BOX_DEFAULT)
        self.ui.box_ray_trace_mode.addItem("normal color")
        self.ui.box_ray_trace_mode.addItem("normal color + black outline")
        self.ui.box_ray_trace_mode.addItem("black outline only")
        self.ui.box_ray_trace_mode.addItem("quantized color + black outline")

        # combo box Ray Texture
        self.ui.box_ray_texture.addItem(COMBO_BOX_DEFAULT)
        self.ui.box_ray_texture.addItem("None")
        self.ui.box_ray_texture.addItem("Matte 1")
        self.ui.box_ray_texture.addItem("Matte 2")
        self.ui.box_ray_texture.addItem("Swirl 1")
        self.ui.box_ray_texture.addItem("Swirl 2")
        self.ui.box_ray_texture.addItem("Fiber")

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
        self.statusBar.setToolTip("Status information: Current process")

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
        self.statusBar.setToolTip("Status information: Current process")

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

    # private functions
    def __create_directory(self, parentPath, dirName):
        """This function creates a directory with a given path and directory name

        Args:
            parentPath:
                parent path where the new directory should be created
            dirName:
                name of the new directory
        """
        newDir = f"{parentPath}/{dirName}"
        if not os.path.exists(newDir):
            os.mkdir(newDir)

    def __create_project_folder(self, DialogWarningPredictionProject):
        """This function creates a project folder.

        Args:
            DialogWarningPredictionProject:
                entire dialog which warns the user that the project folder
                has been already created
        """
        projectName = self.ui.txt_prediction_project_name.text()
        projectPath = f"{self.workspacePath}/{projectName}"
        # check if the project folder already exists
        if os.path.exists(projectPath):
            self.statusBar.showMessage(
                f"Warning! | Current Workspace: {self.workspacePath}")
            dialog = DialogWarningPredictionProject
            dialog.exec_()
            # breaks out of function
            self.statusBar.clearMessage()
            return None
        else:
            os.mkdir(projectPath)
            return projectName, projectPath

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
        if self.ui.cb_analysis_chain_info.isChecked() and len(self.ui.txt_analysis_chain_ref.text()) > 0 and len(self.ui.txt_analysis_chain_model.text()) > 0:
            i += 1
        if self.ui.cb_analysis_chain_info.isChecked():
            j = 4
        else:
            j = 3
        if i == j:
            self.ui.btn_analysis_start.setEnabled(True)
        else:
            self.ui.btn_analysis_start.setEnabled(False)

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
        if self.ui.cb_batch_chain_info.isChecked() and len(self.ui.txt_batch_chain_ref.text()) > 0 and len(self.ui.txt_batch_chain_model.text()) > 0:
            i += 1
        if self.ui.cb_batch_chain_info.isChecked():
            j = 4
        else:
            j = 3
        if i == j:
            self.ui.btn_batch_start.setEnabled(True)
        else:
            self.ui.btn_batch_start.setEnabled(False)

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
        if self.ui.cb_prediction_chain_info.isChecked() and len(self.ui.txt_prediction_chain_ref.text()) > 0 and len(self.ui.txt_prediction_chain_model.text()) > 0:
            i += 1
        if self.ui.cb_prediction_chain_info.isChecked():
            j = 3
        else:
            j = 2
        if i == j:
            self.ui.btn_prediction_start.setEnabled(True)
        else:
            self.ui.btn_prediction_start.setEnabled(False)

    # @SLOT
    # Menu
    def open(self):
        ATTRIBUTE = "value"
        # open file dialog
        try:
            fileName = Qt.QtWidgets.QFileDialog.getOpenFileName(self,
                                                                "Open project file",
                                                                Qt.QtCore.QDir.homePath(),
                                                                "Plugin Project File (project.xml)")
            if fileName == ("",""):
                raise ValueError("No file has been selected.")

            projectFile = tools.ProjectXml(fileName[0])
            tmpProjectFile = projectFile.load_xml_in_memory()

            if projectFile.get_path(tmpProjectFile, "predictionDone", ATTRIBUTE) == "False":
                self.ui.txt_analysis_project_name.setText(
                    projectFile.get_path(tmpProjectFile,
                                         "projectName",
                                         ATTRIBUTE))
                self.ui.txt_analysis_load_reference.setText(
                    projectFile.get_path(tmpProjectFile,
                                         "reference",
                                         ATTRIBUTE))
                if len(projectFile.get_path(tmpProjectFile, "referenceChains",
                                            ATTRIBUTE)) > 0:
                    self.ui.cb_analysis_chain_info.setChecked(True)
                    self.ui.txt_analysis_chain_ref.setText(projectFile.get_path(
                        tmpProjectFile, "referenceChains", ATTRIBUTE))
                    self.ui.txt_analysis_chain_model.setText(projectFile.get_path(
                        tmpProjectFile, "modelChains", ATTRIBUTE))
                resultsPath = projectFile.get_path(tmpProjectFile, "results",
                                                   ATTRIBUTE)
                sessionFileName = ""
                for file in os.listdir(f"{resultsPath}/sessions"):
                    sessionFileName = file
                cmd.load(f"{resultsPath}/sessions/{sessionFileName}")
            else:
                self.ui.txt_prediction_project_name.setText(projectFile.get_path(tmpProjectFile,
                                                                    "projectName",
                                                                                 ATTRIBUTE))
                self.ui.txt_prediction_load_reference.setText(projectFile.get_path(tmpProjectFile,
                                                                  "reference",
                                                                                   ATTRIBUTE))
                if len(projectFile.get_path(tmpProjectFile, "referenceChains",
                                            ATTRIBUTE)) > 0:
                    self.ui.cb_prediction_chain_info.setChecked(True)
                    self.ui.txt_prediction_chain_ref.setText(projectFile.get_path(
                        tmpProjectFile, "referenceChains", ATTRIBUTE))
                    self.ui.txt_prediction_chain_model.setText(projectFile.get_path(
                        tmpProjectFile, "modelChains", ATTRIBUTE))
                resultsPath = projectFile.get_path(tmpProjectFile, "results",
                                                   ATTRIBUTE)
                sessionFileName = ""
                for file in os.listdir(f"{resultsPath}/sessions"):
                    sessionFileName = file
                cmd.load(f"{resultsPath}/sessions/{sessionFileName}")
        except FileNotFoundError:
            print("File could not be opened.")
        except ValueError:
            print("No file has been selected.")

    def save_as(self):
        ATTRIBUTE = "value"
        try:
            fileName = Qt.QtWidgets.QFileDialog.getOpenFileName(self,
                                                                "Open project file",
                                                                Qt.QtCore.QDir.homePath(),
                                                                "Plugin Project File (project.xml)")
            if fileName == ("", ""):
                raise ValueError

            projectFile = tools.ProjectXml(fileName[0])
            tmpProjectFile = projectFile.load_xml_in_memory()
            resultsPath = projectFile.get_path(tmpProjectFile, "results", ATTRIBUTE)
            sessionFileName = ""
            for file in os.listdir(f"{resultsPath}/sessions"):
                sessionFileName = file
            cmd.save(f"{resultsPath}/sessions/{sessionFileName}")
        except ValueError:
            print("No file has been selected.")
        except FileNotFoundError:
            print("File not found!")

    def restore_settings(self):
        """This function deletes the old settings.xml and creates a new one,
        with the default values.

        TODO:
            * change fixed path to variable path
        """
        default_settings = tools.SettingsXml("/home/matt/Documents/settings.xml")
        default_settings.create_settings_xml_file()

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
            tmpPath = f"{os.path.expanduser('~')}/Downloads"
            os.chdir(tmpPath)
            # path to store all pdb files (should be defined in the Settings tab!)
            target_dir = f"{os.path.expanduser('~')}/Documents/data"

            # tmpList contains a list with all zips starting with prediction
            tmpList = tools.filter_prediction_zips(tmpPath)

            # create temporary dir to handle unzipping
            tmpDir = f"{tmpPath}/tmp"
            if not os.path.exists(tmpDir):
                os.mkdir(tmpDir)

            for file in tmpList:
                tools.extract_and_move_model_pdb(tmpPath, tmpDir, file, target_dir)

            shutil.rmtree(tmpDir)
            self.statusBar.showMessage("Preparing the .pdb files was successful.")
        except Exception:
            self.statusBar.showMessage("Preparing the .pdb files failed!")

    def open_settings_global(self):
        """This function open the dialog for the global settings.

        """
        dialog = dialog_settings_global.DialogSettingsGlobal()
        dialog.exec_()

    def open_settings_pdb_preparation(self):
        """This function opens the dialog for the pdb Preparation settings.

        """
        dialog = DialogSettingsPdbPreparation.DialogSettingsPdbPreparation()
        dialog.exec_()

    def display_workspace_path(self):
        self.statusBar.showMessage(f"Current Workspace: {self.workspacePath}")

    def display_project_path(self):
        """

        TODO:
            * implement displaying project path function
        """
        self.statusBar.showMessage(f"Current Project Path: ... ")

    def open_documentation(self):
        """This function opens the official plugin documentation.

        """
        webbrowser.open_new("docs/pymol_plugin/build/html/index.html")

    def open_documentation_pdf(self):
        webbrowser.open_new("docs/pymol_plugin/build/latex/pymol-plugin.pdf")

    def open_about(self):
        dialog = dialog_about.DialogAbout()
        dialog.exec_()

    # Prediction + Analysis
    def load_reference_for_prediction(self):
        try:
            # open file dialog
            fileName = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Reference",
                                                                Qt.QtCore.QDir.homePath(),
                                                                "PDB Files (*.pdb)")
            if fileName == ("", ""):
                raise ValueError("No file has been selected.")
            # display path in text box
            self.ui.txt_prediction_load_reference.setText(str(fileName[0]))
            self.statusBar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_prediction()
        except FileNotFoundError:
            self.statusBar.showMessage("Loading the reference failed!")
        except ValueError:
            print("No file has been selected.")

    def enable_chain_information_input_for_prediction(self):
        """This function enables the text boxes to enter the chains for the
        reference and the model

        """
        try:
            self.ui.txt_prediction_chain_ref.setEnabled(self.ui.cb_prediction_chain_info.checkState())
            self.ui.txt_prediction_chain_model.setEnabled(self.ui.cb_prediction_chain_info.checkState())
            self.statusBar.showMessage("Enter the chain information.")
            self.__check_start_possibility_prediction()
        except Exception:
            self.statusBar.showMessage("Unexpected Error.")

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
        self.__check_start_possibility_prediction()

    def check_prediction_if_txt_prediction_load_reference_is_filled(self):
        self.__check_start_possibility_prediction()

    def predict(self):
        """This function opens a webbrowser with a colab notebook, to run the
        prediction. In addition it runs the entire analysis after the
        prediction.

        """
        self.statusBar.removeWidget(self.workspace)
        self.statusBar.showMessage("Checking user input ...")

        # check if a prediction is already finished
        FILE_NAME = constants.FULL_FILENAME_PREDICTION_ZIP
        if os.path.isfile(FILE_NAME) == True:
            self.statusBar.showMessage(
                f"Warning! | Current Workspace: {self.workspacePath}")
            dialog = DialogWarningPredictionZip.DialogWarningPredictionZip()
            dialog.exec_()

        projectName, projectPath = gui_utils.create_project_folder(
            self.ui.txt_prediction_project_name, self.workspacePath, self.statusBar,
            DialogWarningPredictionProject.DialogWarningPredictionProject())

        # TODO:
        #   * the code block below can be subbed by the function
        #   createProjectFolder

        # creates a project folder
        # projectName = self.ui.txt_prediction_project_name.text()
        # projectPath = f"{self.workspacePath}/{projectName}"
        # # check if the project folder already exists
        # if os.path.exists(projectPath):
        #     self.statusBar.showMessage(
        #         f"Warning! | Current Workspace: {self.workspacePath}")
        #     dialog = DialogWarningPredictionProject.DialogWarningPredictionProject()
        #     dialog.exec_()
        #     # breaks out of function
        #     self.statusBar.clearMessage()
        #     self.statusBar.addPermanentWidget(self.workspace)
        #     return None
        # else:
        #     os.mkdir(projectPath)

        # TODO:
        #   * the two code blocks below can be subbed by the function
        #   createDirectory
        # creates a pdb folder
        projectPdbDir = f"{projectPath}/pdb"
        if not os.path.exists(projectPdbDir):
            os.mkdir(projectPdbDir)
        # creates a tmp folder to unzip the prediction
        tmpProjectDir = f"{projectPath}/tmp"
        if not os.path.exists(tmpProjectDir):
            os.mkdir(tmpProjectDir)

        # check if input is valid
        # gets reference filename and filepath
        if len(self.ui.txt_prediction_load_reference.text()) == 4:
            tmpREFERENCE_OBJ_NAME = self.ui.txt_prediction_load_reference.text()
            tmpProtein = core.protein(tmpREFERENCE_OBJ_NAME,
                                      export_data_dir=projectPdbDir)
            tmpProtein.clean_pdb_file()
            REFERENCE_OBJ_NAME = tmpREFERENCE_OBJ_NAME
            REFERENCE_DIR = projectPdbDir
        else:
            refFileInfo = Qt.QtCore.QFileInfo(
                self.ui.txt_prediction_load_reference.text())
            REFERENCE_OBJ_NAME = refFileInfo.baseName()
            REFERENCE_DIR = refFileInfo.canonicalPath()

        # creates results folder
        resultsPath = f"{projectPath}/Results"
        os.mkdir(resultsPath)

        # TODO:
        #   * the code block below can be subbed by the function
        #   setValuesInProjectXml
        # see function setValuesInProjectXml!!
        # fullProjectFileName = f"{projectPath}/project.xml"
        # projectFile = tools.ProjectXml(fullProjectFileName)
        # projectFile.create_settings_xml_file()
        # try:
        #     tmpProjectFile = projectFile.load_xml_in_memory()
        #     projectFile.set_value(tmpProjectFile, "projectName", "value",
        #                           projectName)
        #     projectFile.set_value(tmpProjectFile, "predictionDone", "value",
        #                           "True")
        #     projectFile.set_value(tmpProjectFile, "reference", "value",
        #                           f"{REFERENCE_DIR}/{REFERENCE_OBJ_NAME}")
        #     projectFile.set_value(tmpProjectFile, "referenceChains", "value",
        #                           self.ui.txt_prediction_chain_ref.text())
        #     projectFile.set_value(tmpProjectFile, "modelChains", "value",
        #                           self.ui.txt_prediction_chain_model.text())
        #     projectFile.set_value(tmpProjectFile, "results", "value",
        #                           resultsPath)
        # except FileNotFoundError:
        #     print("Project file could not be loaded!")

        # starting the default web browser to display the colab notebook
        self.statusBar.showMessage("Opening Google colab notebook ...")
        if self.ui.action_settings_model_w_off_colab_notebook.isChecked():
            webbrowser.open_new(constants.OFFICIAL_NOTEBOOK_URL)
        else:
            webbrowser.open_new(constants.NOTEBOOK_URL)

        # waiting for the colab notebook to finish
        # TODO:
        #   * implement cancel button for "Abort Prediction"
        self.ui.btn_prediction_cancel.setEnabled(True)
        flag = False
        while flag == False:
            print("AlphaFold is still running ...")
            time.sleep(5)
            # time.sleep(120)
        # while os.path.isfile(FILE_NAME) == False:
        #     print("AlphaFold is still running ...")
        #     time.sleep(5)
        #     # time.sleep(120)

        # start of the analysis algorithm
        self.statusBar.showMessage("Protein structure analysis started ...")

        archive = "prediction.zip"
        sourcePath = f"{os.path.expanduser('~')}/Downloads"
        # extracts and moves the prediction.pdb to the workspace/pdb folder
        tools.extract_and_move_model_pdb(
            sourcePath, tmpProjectDir, archive, projectPdbDir)

        # removes tmp dir
        os.remove(tmpProjectDir)

        # gets model filename and filepath
        PREDICTION_NAME = tools.get_prediction_file_name(projectPdbDir)
        fullModelFilePath = f"{projectPdbDir}/{PREDICTION_NAME[0]}"
        modelFileInfo = Qt.QtCore.QFileInfo(fullModelFilePath)
        MODEL_OBJ_NAME = modelFileInfo.baseName()
        MODEL_DIR = modelFileInfo.canonicalPath()

        projectFile.set_value(tmpProjectFile, "model", "value",
                              fullModelFilePath)
        projectFile.save_xml_file(tmpProjectFile)

        try:
            # gets chain information for the reference
            refChains = self.ui.txt_prediction_chain_ref.text().split(",")

            # gets chain information for the model
            modelChains = self.ui.txt_prediction_chain_model.text().split(",")

            # define a parent directory where all results will be saved
            # EXPORT_SUBDIR: str = f"{os.path.expanduser('~')}/anaconda3/envs/pymol_plugin/lib/python3.9/site-packages/pmg_tk/startup/pymol_plugin/data/results"
            EXPORT_SUBDIR = resultsPath

            # define parameters for the align command
            CYCLES: int = 1
            CUTOFF: float = 1.0
            # define a filename for the rmsd and aligned amino acids CSV file
            FILE_NAME: str = "rmsd_and_aligned_AA_form_alignment"

            # constant variables for plotting
            FIGURE_SIZE = (11.0, 6.0)
            # defines the atoms of which the distance gets computed
            DISTANCE_LABEL: str = "$\\alpha$-C"

            # define variable filenames for the different results
            ALIGNMENT_FILE_NAME: str = f"alignment_file_{MODEL_OBJ_NAME}"
            DISTANCE_FILE_NAME: str = f"distance_file_{MODEL_OBJ_NAME}"
            SESSION_FILE_NAME: str = f"session_file_{MODEL_OBJ_NAME}"

            # comment the below code in, if you wish to have a fixed y-axis for ALL plots
            # LIMIT_Y: tuple[float, float] = (0, 7)

            # create the protein object for the reference
            reference_protein: core.protein = core.protein(REFERENCE_OBJ_NAME,
                                                           REFERENCE_DIR)
            # sets selection for any chain combination of reference
            refSelection = ""
            seperator = ", "
            tmpList = []

            for chain in refChains:
                tmpSelection = f"/{reference_protein.molecule_object}//{chain}//CA"
                tmpList.append(tmpSelection)
                refSelection = seperator.join(tmpList)

            reference_protein.set_selection(refSelection)
            self.ui.progress_bar_prediction.setProperty("value", 5)
            # sets selection for any chain combination of model
            # create model protein object
            model_protein: core.protein = core.protein(MODEL_OBJ_NAME,
                                                       MODEL_DIR)
            modelSelection = ""
            tmpList = []
            for chain in modelChains:
                tmpSelection = f"/{model_protein.molecule_object}//{chain}//CA"
                tmpList.append(tmpSelection)
                modelSelection = seperator.join(tmpList)

            model_protein.set_selection(modelSelection)
            self.ui.progress_bar_prediction.setProperty("value", 10)
            # create protein pair object
            protein_pair = core.proteinpair(reference_protein, model_protein,
                                            EXPORT_SUBDIR)

            # --- Begin of PyMOL session --- #
            # reinitialize pymol session
            self.statusBar.showMessage("Reinitialize pymol session ...")

            cmd.reinitialize()
            self.statusBar.showMessage(
                "Finished reinitializing pymol session. | Load proteins ...")

            # load both proteins in the pymol session
            protein_pair.load_protein_pair()
            self.ui.progress_bar_prediction.setProperty("value", 15)
            self.statusBar.showMessage(
                "Finished loading proteins. | Color proteins ...")
            print(f"Finished loading proteins.")

            # color protein with default colors; ref: green, model: blue
            protein_pair.color_protein_pair()
            self.ui.progress_bar_prediction.setProperty("value", 20)
            self.statusBar.showMessage(
                "Finished coloring proteins. | Align proteins ...")
            print(f"Finished coloring proteins.")

            # do the structure alignment
            align_results = protein_pair.align_protein_pair(CYCLES, CUTOFF,
                                                            ALIGNMENT_FILE_NAME)
            self.ui.progress_bar_prediction.setProperty("value", 25)
            self.statusBar.showMessage(
                "Finished aligning proteins. | Calculate distances ...")
            print(f"Finished aligning proteins.")

            # do the distance calculation
            distance_results = protein_pair.calculate_distance_between_ca_atoms(
                ALIGNMENT_FILE_NAME)
            protein_pair.export_distance_between_ca_atoms(distance_results)
            self.ui.progress_bar_prediction.setProperty("value", 40)
            self.statusBar.showMessage(
                "Finished calculating distances. | Create distance plot ...")
            print(f"Finished distance calculations.")

            # create an instance of the graphics class
            graphics_instance: graphics.graphics = graphics.graphics(
                protein_pair,
                distance_results,
                FIGURE_SIZE)
            # create distance plot
            fig = graphics_instance.create_distance_plot(DISTANCE_LABEL, CUTOFF)
            self.ui.progress_bar_prediction.setProperty("value", 45)
            self.statusBar.showMessage(
                "Finished creating distance plot. | Create distance histogram ...")
            print(f"Finished creation of distance plot.")

            # save distance plot
            if not os.path.exists(f"{protein_pair.results_dir}/plots"):
                os.mkdir(f"{protein_pair.results_dir}/plots")
            if not os.path.exists(
                    f"{protein_pair.results_dir}/plots/distance_plot"):
                os.mkdir(f"{protein_pair.results_dir}/plots/distance_plot")
            plt.savefig(f"{protein_pair.results_dir}/plots/distance_plot/"
                        f"distance_plot_{MODEL_OBJ_NAME}.svg")
            plt.close(fig)

            # create distance histogram
            graphics_instance.create_distance_histogram()
            self.ui.progress_bar_prediction.setProperty("value", 50)
            self.statusBar.showMessage(
                "Finished creating distance histogram. | "
                "Take image of structure alignment ...")
            print(f"Finished creation of distance histogram.")

            # save distance histogram
            if not os.path.exists(f"{protein_pair.results_dir}/plots"):
                os.mkdir(f"{protein_pair.results_dir}/plots")
            if not os.path.exists(
                    f"{protein_pair.results_dir}/plots/distance_histogram"):
                os.mkdir(f"{protein_pair.results_dir}/plots/distance_histogram")
            plt.savefig(f"{protein_pair.results_dir}/plots/distance_histogram"
                        f"/distance_histogram_{MODEL_OBJ_NAME}.svg")

            # take image of whole structure alignment
            graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME,
                                                         "cartoon",
                                                         "test")
            self.ui.progress_bar_prediction.setProperty("value", 70)
            self.statusBar.showMessage(
                f"Finished taking image of structure alignment. | Take images "
                f"of interesting regions (within {CUTOFF} angstrom)")
            print(f"Finished creation of image which shows the whole structure "
                  f"alignment.")

            # take image of interesting regions
            graphics_instance.take_image_of_interesting_regions(3.0,
                                                                "interesting_region",
                                                                opaque_background=1)
            self.ui.progress_bar_prediction.setProperty("value", 90)
            self.statusBar.showMessage(
                f"Finished taking images of interesting regions. | "
                f"Save PyMOL session")
            print(
                f"Finished creation of images which show interesting regions with"
                f"a distance greater than {CUTOFF} angstrom.")

            # color residues by distance
            # graphics_instance.color_by_distance(ALIGNMENT_FILE_NAME)
            # print(f"Finished coloring of prediction with color_by_distance.")
            # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
            #                                              "coloredByRMSD")
            # graphics_instance.create_gif()
            # print(f"Finished creation of gif which shows the whole predicted "
            #       f"structure colored by distance.")

            # save pymol session
            protein_pair.save_session_of_protein_pair(SESSION_FILE_NAME)
            self.statusBar.showMessage(
                f"Finished saving PyMOL session. | "
                f"Protein structure analysis has successfully finished.")
            print(f"Finished saving of pymol session file.")
            self.statusBar.showMessage("")
            self.ui.progress_bar_prediction.setProperty("value", 100)
            dialog = dialog_finished.DialogFinished()
            dialog.exec_()
            self.ui.btn_prediction_start.setEnabled(True)
        except Exception:
            self.statusBar.showMessage("Protein structure analysis has failed!")

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
            self.statusBar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_prediction()
        except FileNotFoundError:
            self.statusBar.showMessage("Loading the reference failed!")
        except ValueError:
            print("No file has been selected.")
            self.__check_start_possibility()
        except Exception:
            self.statusBar.showMessage("Loading the reference failed!")

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
            self.statusBar.showMessage("Loading the model was successful.")
            self.__check_start_possibility()
        except Exception:
            self.statusBar.showMessage("Loading the model failed!")

    def enable_chain_information_input_for_analysis(self):
        """This function enables the text boxes to enter the chains for the
        reference and the model

        """
        try:
            self.ui.txt_analysis_chain_ref.setEnabled(self.ui.cb_analysis_chain_info.checkState())
            self.ui.txt_analysis_chain_model.setEnabled(self.ui.cb_analysis_chain_info.checkState())
            self.statusBar.showMessage("Enter the chain information.")
            self.__check_start_possibility()
        except Exception:
            self.statusBar.showMessage("Unexpected Error.")

    def check_analysis_if_txt_analysis_project_name_is_filled(self):
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
        self.statusBar.removeWidget(self.workspace)
        self.ui.btn_analysis_start.setEnabled(False)
        self.statusBar.showMessage("Protein structure analysis started ...")

        # creates a project folder
        projectName = self.ui.txt_analysis_project_name.text()
        projectPath = f"{self.workspacePath}/{projectName}"
        # check if the project folder already exists
        if os.path.exists(projectPath):
            self.statusBar.showMessage(
                f"Warning! | Current Workspace: {self.workspacePath}")
            dialog = DialogWarningPredictionProject.DialogWarningPredictionProject()
            dialog.exec_()
            # breaks out of function
            self.statusBar.clearMessage()
            self.statusBar.addPermanentWidget(self.workspace)
            return None
        else:
            os.mkdir(projectPath)

        # creates a pdb folder
        projectPdbDir = f"{projectPath}/pdb"
        if not os.path.exists(projectPdbDir):
            os.mkdir(projectPdbDir)

        # creates results folder
        resultsPath = f"{projectPath}/Results"
        os.mkdir(resultsPath)

        try:
            # gets reference filename and filepath
            # refFileInfo = Qt.QtCore.QFileInfo(self.ui.txt_analysis_load_reference.text())
            # REFERENCE_OBJ_NAME = refFileInfo.baseName()
            # REFERENCE_DIR = refFileInfo.canonicalPath()

            # check if input is valid
            # gets reference filename and filepath
            if len(self.ui.txt_analysis_load_reference.text()) == 4:
                tmpREFERENCE_OBJ_NAME = self.ui.txt_analysis_load_reference.text()
                tmpProtein = core.protein(tmpREFERENCE_OBJ_NAME,
                                          export_data_dir=projectPdbDir)
                tmpProtein.clean_pdb_file()
                REFERENCE_OBJ_NAME = tmpREFERENCE_OBJ_NAME
                REFERENCE_DIR = projectPdbDir
            else:
                refFileInfo = Qt.QtCore.QFileInfo(
                    self.ui.txt_analysis_load_reference.text())
                REFERENCE_OBJ_NAME = refFileInfo.baseName()
                REFERENCE_DIR = refFileInfo.canonicalPath()

            # gets model filename and filepath
            modelFileInfo = Qt.QtCore.QFileInfo(self.ui.txt_analysis_load_model.toPlainText())
            MODEL_OBJ_NAME = modelFileInfo.baseName()
            MODEL_DIR = modelFileInfo.canonicalPath()

            fullProjectFileName = f"{projectPath}/project.xml"
            projectFile = tools.ProjectXml(fullProjectFileName)
            projectFile.create_settings_xml_file()

            tmpProjectFile = projectFile.load_xml_in_memory()
            projectFile.set_value(tmpProjectFile, "projectName", "value",
                                  projectName)
            projectFile.set_value(tmpProjectFile, "predictionDone", "value",
                                  "False")
            projectFile.set_value(tmpProjectFile, "reference", "value",
                                  f"{REFERENCE_DIR}/{REFERENCE_OBJ_NAME}")
            projectFile.set_value(tmpProjectFile, "referenceChains",
                                  "value",
                                  self.ui.txt_analysis_chain_ref.text())
            projectFile.set_value(tmpProjectFile, "modelChains", "value",
                                  self.ui.txt_analysis_chain_model.text())
            projectFile.set_value(tmpProjectFile, "results", "value",
                                  resultsPath)

            # gets chain information for the reference
            refChains = self.ui.txt_analysis_chain_ref.text().split(",")

            # gets chain information for the model
            modelChains = self.ui.txt_analysis_chain_model.text().split(",")

            # define a parent directory where all results will be saved
            # EXPORT_SUBDIR: str = f"{os.path.expanduser('~')}/anaconda3/envs/pymol_plugin/lib/python3.9/site-packages/pmg_tk/startup/pymol_plugin/data/results"
            EXPORT_SUBDIR: str = resultsPath

            # define parameters for the align command
            CYCLES: int = 1
            CUTOFF: float = 1.0
            # define a filename for the rmsd and aligned amino acids CSV file
            FILE_NAME: str = "rmsd_and_aligned_AA_form_alignment"

            # constant variables for plotting
            FIGURE_SIZE = (11.0, 6.0)
            # defines the atoms of which the distance gets computed
            DISTANCE_LABEL: str = "$\\alpha$-C"

            # define variable filenames for the different results
            ALIGNMENT_FILE_NAME: str = f"alignment_file_{MODEL_OBJ_NAME}"
            DISTANCE_FILE_NAME: str = f"distance_file_{MODEL_OBJ_NAME}"
            SESSION_FILE_NAME: str = f"session_file_{MODEL_OBJ_NAME}"

            # comment the below code in, if you wish to have a fixed y-axis for ALL plots
            # LIMIT_Y: tuple[float, float] = (0, 7)

            # create the protein object for the reference
            reference_protein: core.protein = core.protein(REFERENCE_OBJ_NAME,
                                                           REFERENCE_DIR)
            # sets selection for any chain combination of reference
            refSelection = ""
            seperator = ", "
            tmpList = []

            for chain in refChains:
                tmpSelection = f"/{reference_protein.molecule_object}//{chain}//CA"
                tmpList.append(tmpSelection)
                refSelection = seperator.join(tmpList)

            reference_protein.set_selection(refSelection)
            self.ui.progress_bar_analysis.setProperty("value", 5)
            # sets selection for any chain combination of model
            # create model protein object
            model_protein: core.protein = core.protein(MODEL_OBJ_NAME,
                                                       MODEL_DIR)
            modelSelection = ""
            tmpList = []
            for chain in modelChains:
                tmpSelection = f"/{model_protein.molecule_object}//{chain}//CA"
                tmpList.append(tmpSelection)
                modelSelection = seperator.join(tmpList)

            model_protein.set_selection(modelSelection)
            self.ui.progress_bar_analysis.setProperty("value", 10)
            # create protein pair object
            protein_pair = core.proteinpair(reference_protein, model_protein,
                                            EXPORT_SUBDIR)

            # --- Begin of PyMOL session --- #
            # reinitialize pymol session
            self.statusBar.showMessage("Reinitialize pymol session ...")

            cmd.reinitialize()
            self.statusBar.showMessage(
                "Finished reinitializing pymol session. | Load proteins ...")

            # load both proteins in the pymol session
            protein_pair.load_protein_pair()
            self.ui.progress_bar_analysis.setProperty("value", 15)
            self.statusBar.showMessage(
                "Finished loading proteins. | Color proteins ...")
            print(f"Finished loading proteins.")

            # color protein with default colors; ref: green, model: blue
            protein_pair.color_protein_pair()
            self.ui.progress_bar_analysis.setProperty("value", 20)
            self.statusBar.showMessage(
                "Finished coloring proteins. | Align proteins ...")
            print(f"Finished coloring proteins.")

            # do the structure alignment
            align_results = protein_pair.align_protein_pair(CYCLES, CUTOFF,
                                                            ALIGNMENT_FILE_NAME)
            self.ui.progress_bar_analysis.setProperty("value", 25)
            self.statusBar.showMessage(
                "Finished aligning proteins. | Calculate distances ...")
            print(f"Finished aligning proteins.")

            # do the distance calculation
            distance_results = protein_pair.calculate_distance_between_ca_atoms(
                ALIGNMENT_FILE_NAME)
            protein_pair.export_distance_between_ca_atoms(distance_results)
            self.ui.progress_bar_analysis.setProperty("value", 40)
            self.statusBar.showMessage(
                "Finished calculating distances. | Create distance plot ...")
            print(f"Finished distance calculations.")

            # create an instance of the graphics class
            graphics_instance: graphics.graphics = graphics.graphics(protein_pair,
                                                                     distance_results,
                                                                     FIGURE_SIZE)
            # create distance plot
            fig = graphics_instance.create_distance_plot(DISTANCE_LABEL, CUTOFF)
            self.ui.progress_bar_analysis.setProperty("value", 45)
            self.statusBar.showMessage(
                "Finished creating distance plot. | Create distance histogram ...")
            print(f"Finished creation of distance plot.")

            # save distance plot
            if not os.path.exists(f"{protein_pair.results_dir}/plots"):
                os.mkdir(f"{protein_pair.results_dir}/plots")
            if not os.path.exists(f"{protein_pair.results_dir}/plots/distance_plot"):
                os.mkdir(f"{protein_pair.results_dir}/plots/distance_plot")
            plt.savefig(f"{protein_pair.results_dir}/plots/distance_plot/"
                        f"distance_plot_{MODEL_OBJ_NAME}.svg")
            plt.close(fig)

            # create distance histogram
            graphics_instance.create_distance_histogram()
            self.ui.progress_bar_analysis.setProperty("value", 50)
            self.statusBar.showMessage(
                "Finished creating distance histogram. | "
                "Take image of structure alignment ...")
            print(f"Finished creation of distance histogram.")

            # save distance histogram
            if not os.path.exists(f"{protein_pair.results_dir}/plots"):
                os.mkdir(f"{protein_pair.results_dir}/plots")
            if not os.path.exists(
                    f"{protein_pair.results_dir}/plots/distance_histogram"):
                os.mkdir(f"{protein_pair.results_dir}/plots/distance_histogram")
            plt.savefig(f"{protein_pair.results_dir}/plots/distance_histogram"
                        f"/distance_histogram_{MODEL_OBJ_NAME}.svg")

            # take image of whole structure alignment
            graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME,
                                                         "cartoon",
                                                         "test")
            self.ui.progress_bar_analysis.setProperty("value", 70)
            self.statusBar.showMessage(
                f"Finished taking image of structure alignment. | Take images "
                f"of interesting regions (within {CUTOFF} angstrom)")
            print(f"Finished creation of image which shows the whole structure "
                  f"alignment.")

            # take image of interesting regions
            graphics_instance.take_image_of_interesting_regions(3.0,
                                                                "interesting_region",
                                                                opaque_background=1)
            self.ui.progress_bar_analysis.setProperty("value", 90)
            self.statusBar.showMessage(
                f"Finished taking images of interesting regions. | "
                f"Save PyMOL session")
            print(f"Finished creation of images which show interesting regions with"
                  f"a distance greater than {CUTOFF} angstrom.")

            # color residues by distance
            # graphics_instance.color_by_distance(ALIGNMENT_FILE_NAME)
            # print(f"Finished coloring of prediction with color_by_distance.")
            # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
            #                                              "coloredByRMSD")
            # graphics_instance.create_gif()
            # print(f"Finished creation of gif which shows the whole predicted "
            #       f"structure colored by distance.")

            # save pymol session
            protein_pair.save_session_of_protein_pair(SESSION_FILE_NAME)
            self.statusBar.showMessage(
                f"Finished saving PyMOL session. | "
                f"Protein structure analysis has successfully finished.")
            print(f"Finished saving of pymol session file.")
            self.statusBar.showMessage("")
            self.ui.progress_bar_analysis.setProperty("value", 100)
            dialog = dialog_finished.DialogFinished()
            dialog.exec_()
            self.ui.btn_analysis_start.setEnabled(True)
        except Exception:
            self.statusBar.showMessage("Protein structure analysis has failed!")

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
            self.statusBar.showMessage("Loading the reference was successful.")
            self.__check_start_possibility_batch()
        except Exception:
            self.statusBar.showMessage("Loading the reference failed!")

    def load_model_for_batch(self):
        """This function loads multiple files as models.

        """
        try:
            # open file dialog
            fileNames = Qt.QtWidgets.QFileDialog.getOpenFileNames(self, "Open Models",
                                                                  self.target_dir,
                                                                  "PDB Files (*.pdb)")
            # display path in text box
            for file in fileNames[0]:
                self.ui.txt_batch_load_model.append(str(file))
            self.statusBar.showMessage("Loading the models was successful.")
            self.__check_start_possibility_batch()
        except Exception:
            self.statusBar.showMessage("Loading the models failed!")

    def enable_chain_information_input_for_batch(self):
        """This function enables the text boxes to enter the chains for the
        reference and the model

        """
        try:
            self.ui.txt_batch_chain_ref.setEnabled(
                self.ui.cb_batch_chain_info.checkState())
            self.ui.txt_batch_chain_model.setEnabled(
                self.ui.cb_batch_chain_info.checkState())
            self.statusBar.showMessage("Enter the chain information.")
            self.__check_start_possibility_batch()
        except Exception:
            self.statusBar.showMessage("Unexpected Error.")

    def check_batch_if_txt_batch_job_name_is_filled(self):
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
        self.statusBar.removeWidget(self.workspace)

        # creates a project folder
        projectName = self.ui.txt_prediction_project_name.text()
        projectPath = f"{self.workspacePath}/{projectName}"
        # check if the project folder already exists
        if os.path.exists(projectPath):
            self.statusBar.showMessage(
                f"Warning! | Current Workspace: {self.workspacePath}")
            dialog = DialogWarningPredictionProject.DialogWarningPredictionProject()
            dialog.exec_()
            # breaks out of function
            self.statusBar.clearMessage()
            self.statusBar.addPermanentWidget(self.workspace)
            return None
        else:
            os.mkdir(projectPath)

        # creates results folder
        resultsPath = f"{projectPath}/Results"
        os.mkdir(resultsPath)

        # gets names of all selected models
        modelFileName = self.ui.txt_batch_load_model.toPlainText()
        # creates a list with the full filename
        modelsFullName = modelFileName.split(sep="\n")
        models = []
        modelsDir = []
        # creates a list with the filename only, without the file extension
        for model in modelsFullName:
            tmpModel, tmpModelDir = tools.get_file_path_and_name(model)
            models.append(tmpModel)
            modelsDir.append(tmpModelDir)

        self.ui.btn_analysis_start.setEnabled(False)
        self.statusBar.showMessage("Protein structure analysis started ...")
        try:
            # gets reference filename and filepath
            REFERENCE_OBJ_NAME, REFERENCE_DIR = tools.get_file_path_and_name(
                self.ui.txtRefPath_batch.toPlainText())

            # gets chain information for the reference & model
            refChains = self.ui.txt_batch_chain_ref.text().split(",")
            modelChains = self.ui.txt_batch_chain_model.text().split(",")

            # define a parent directory where all results will be saved
            EXPORT_SUBDIR: str = resultsPath

            # define parameters for the align command
            CYCLES: int = 1
            CUTOFF: float = 1.0

            # constant variables for plotting
            FIGURE_SIZE = (11.0, 6.0)
            # defines the atoms of which the distance gets computed
            DISTANCE_LABEL: str = "$\\alpha$-C"

            # comment the below code in, if you wish to have a fixed y-axis for ALL plots
            # LIMIT_Y: tuple[float, float] = (0, 7)

            # create the protein object for the reference
            reference_protein: core.protein = core.protein(REFERENCE_OBJ_NAME,
                                                           REFERENCE_DIR)
            # sets selection for any chain combination of reference
            refSelection = reference_protein.create_selection_from_chains(refChains)

            reference_protein.set_selection(refSelection)
            self.ui.progress_bar_analysis.setProperty("value", 5)

            i = 0
            for MODEL_OBJ_NAME in models:
                # define variable filenames for the different results
                ALIGNMENT_FILE_NAME: str = f"alignment_file_{MODEL_OBJ_NAME}"
                DISTANCE_FILE_NAME: str = f"distance_file_{MODEL_OBJ_NAME}"
                SESSION_FILE_NAME: str = f"session_file_{MODEL_OBJ_NAME}"
                MODEL_DIR = modelsDir[i]

                # sets selection for any chain combination of model
                # create model protein object
                model_protein: core.protein = core.protein(MODEL_OBJ_NAME,
                                                           MODEL_DIR)
                modelSelection = ""
                tmpList = []
                for chain in modelChains:
                    tmpSelection = f"/{model_protein.molecule_object}//{chain}//CA"
                    tmpList.append(tmpSelection)
                    modelSelection = seperator.join(tmpList)

                model_protein.set_selection(modelSelection)
                self.ui.progress_bar_analysis.setProperty("value", 10)
                # create protein pair object
                protein_pair = core.proteinpair(reference_protein, model_protein,
                                                EXPORT_SUBDIR)

                # --- Begin of PyMOL session --- #
                # reinitialize pymol session
                self.statusBar.showMessage("Reinitialize pymol session ...")

                cmd.reinitialize()
                self.statusBar.showMessage(
                    "Finished reinitializing pymol session. | Load proteins ...")

                # load both proteins in the pymol session
                protein_pair.load_protein_pair()
                self.ui.progress_bar_analysis.setProperty("value", 15)
                self.statusBar.showMessage(
                    "Finished loading proteins. | Color proteins ...")
                print(f"Finished loading proteins.")

                # color protein with default colors; ref: green, model: blue
                protein_pair.color_protein_pair()
                self.ui.progress_bar_analysis.setProperty("value", 20)
                self.statusBar.showMessage(
                    "Finished coloring proteins. | Align proteins ...")
                print(f"Finished coloring proteins.")

                # do the structure alignment
                align_results = protein_pair.align_protein_pair(CYCLES, CUTOFF,
                                                                ALIGNMENT_FILE_NAME)
                self.ui.progress_bar_analysis.setProperty("value", 25)
                self.statusBar.showMessage(
                    "Finished aligning proteins. | Calculate distances ...")
                print(f"Finished aligning proteins.")

                # do the distance calculation
                distance_results = protein_pair.calculate_distance_between_ca_atoms(
                    ALIGNMENT_FILE_NAME)
                protein_pair.export_distance_between_ca_atoms(distance_results)
                self.ui.progress_bar_analysis.setProperty("value", 40)
                self.statusBar.showMessage(
                    "Finished calculating distances. | Create distance plot ...")
                print(f"Finished distance calculations.")

                # create an instance of the graphics class
                graphics_instance: graphics.graphics = graphics.graphics(protein_pair,
                                                                         distance_results,
                                                                         FIGURE_SIZE)
                # create distance plot
                fig = graphics_instance.create_distance_plot(DISTANCE_LABEL, CUTOFF)
                self.ui.progress_bar_analysis.setProperty("value", 45)
                self.statusBar.showMessage(
                    "Finished creating distance plot. | Create distance histogram ...")
                print(f"Finished creation of distance plot.")

                # save distance plot
                if not os.path.exists(f"{protein_pair.results_dir}/plots"):
                    os.mkdir(f"{protein_pair.results_dir}/plots")
                if not os.path.exists(f"{protein_pair.results_dir}/plots/distance_plot"):
                    os.mkdir(f"{protein_pair.results_dir}/plots/distance_plot")
                plt.savefig(f"{protein_pair.results_dir}/plots/distance_plot/"
                            f"distance_plot_{MODEL_OBJ_NAME}.svg")
                plt.close(fig)

                # create distance histogram
                graphics_instance.create_distance_histogram()
                self.ui.progress_bar_analysis.setProperty("value", 50)
                self.statusBar.showMessage(
                    "Finished creating distance histogram. | "
                    "Take image of structure alignment ...")
                print(f"Finished creation of distance histogram.")

                # save distance histogram
                if not os.path.exists(f"{protein_pair.results_dir}/plots"):
                    os.mkdir(f"{protein_pair.results_dir}/plots")
                if not os.path.exists(
                        f"{protein_pair.results_dir}/plots/distance_histogram"):
                    os.mkdir(f"{protein_pair.results_dir}/plots/distance_histogram")
                plt.savefig(f"{protein_pair.results_dir}/plots/distance_histogram"
                            f"/distance_histogram_{MODEL_OBJ_NAME}.svg")

                # take image of whole structure alignment
                graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME,
                                                             "cartoon",
                                                             "test")
                self.ui.progress_bar_analysis.setProperty("value", 70)
                self.statusBar.showMessage(
                    f"Finished taking image of structure alignment. | Take images "
                    f"of interesting regions (within {CUTOFF} angstrom)")
                print(f"Finished creation of image which shows the whole structure "
                      f"alignment.")

                # take image of interesting regions
                graphics_instance.take_image_of_interesting_regions(3.0,
                                                                    "interesting_region",
                                                                    opaque_background=1)
                self.ui.progress_bar_analysis.setProperty("value", 90)
                self.statusBar.showMessage(
                    f"Finished taking images of interesting regions. | "
                    f"Save PyMOL session")
                print(f"Finished creation of images which show interesting regions with"
                      f"a distance greater than {CUTOFF} angstrom.")

                # color residues by distance
                # graphics_instance.color_by_distance(ALIGNMENT_FILE_NAME)
                # print(f"Finished coloring of prediction with color_by_distance.")
                # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
                #                                              "coloredByRMSD")
                # graphics_instance.create_gif()
                # print(f"Finished creation of gif which shows the whole predicted "
                #       f"structure colored by distance.")

                # save pymol session
                protein_pair.save_session_of_protein_pair(SESSION_FILE_NAME)
                self.statusBar.showMessage(
                    f"Finished saving PyMOL session. | "
                    f"Protein structure analysis has successfully finished.")
                print(f"Finished saving of pymol session file.")
                self.statusBar.showMessage("")
                self.ui.progress_bar_analysis.setProperty("value", 100)
                dialog = dialog_finished.DialogFinished()
                dialog.exec_()
                self.ui.btn_analysis_start.setEnabled(True)

                i += 1
        except Exception:
            self.statusBar.showMessage("Protein structure analysis has failed!")

    # Results
    def display_structure_alignment(self):
        dialog = dialog_image.DialogImage()
        dialog.exec_()

    def display_distance_plot(self):
        app2 = Qt.QtWidgets.QDialog(self)
        appLayout = QHBoxLayout()
        viewer = QtSvg.QSvgWidget()
        viewer.load(
            "/home/matt/Documents/test_pymol/Results/plots/distance_plot/distance_plot_selected_prediction_1.svg")
        viewer.setWindowTitle("Distance Plot")
        viewer.show()
        appLayout.addWidget(viewer)
        app2.setLayout(appLayout)
        app2.exec()

    # Image
    def show_representation(self):
        """This function sets the representation.

        """
        if self.ui.box_representation.currentIndex() == 0:
            print("Please select a representation.")
            self.statusBar.showMessage("Please select a representation.")
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
            self.statusBar.showMessage("Please select a background color.")
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
            self.statusBar.showMessage("Please select a renderer.")
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
            self.statusBar.showMessage("Please select a Ray-Trace-Mode.")
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
            self.statusBar.showMessage("Please select a Ray Texture.")
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

    def preview_image(self):
        """This function previews the image

        """
        if self.ui.cb_ray_tracing.isChecked():
            self.statusBar.showMessage("Preview ray-traced image ...")
            cmd.ray(2400, 2400, renderer=int(self.renderer))
            self.statusBar.showMessage("Finished preview of ray-traced image.")
        else:
            self.statusBar.showMessage("Preview draw image ...")
            cmd.draw(2400, 2400)
            self.statusBar.showMessage("Finished preview of drawn image.")

    def save_image(self):
        """This function saves the image as a png file.

        """
        if self.ui.cb_ray_tracing.isChecked():
            saveDialog = Qt.QtWidgets.QFileDialog()
            fullFileName = saveDialog.getSaveFileName(caption="Save Image",
                                                      filter="Image (*.png)")
            self.statusBar.showMessage("Creating ray-traced image ...")
            cmd.ray(2400, 2400, renderer=int(self.renderer))
            cmd.png(fullFileName[0], dpi=300)
            self.statusBar.showMessage("Finished image creation.")
        else:
            saveDialog = Qt.QtWidgets.QFileDialog()
            fullFileName = saveDialog.getSaveFileName(caption="Save Image",
                                                      filter="Image (*.png)")
            self.statusBar.showMessage("Creating draw image ...")
            cmd.draw(2400, 2400)
            cmd.png(fullFileName[0], dpi=300)
            self.statusBar.showMessage("Finished image creation.")


# comment out if run within pymol
if __name__ == '__main__':
    # Start point of the application
    app = Qt.QtWidgets.QApplication(sys.argv)
    mainWindow = MainWindow()
    mainWindow.show()
    app.exec_()
