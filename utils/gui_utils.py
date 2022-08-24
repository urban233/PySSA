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

import os
from PyQt5.QtWidgets import QMessageBox
from utils import tools


def create_directory(parentPath, dirName):
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


def create_project_folder(txtProjectNamePrediction, workspacePath, statusBar,
                          DialogWarningPredictionProject):
    """This function creates a project folder.

    Args:
        txtProjectNamePrediction:
            textbox which contains the project name
        workspacePath:
            path of the workspace, where the project should be saved
        statusBar:
            status bar object (e.g. self.statusbar)
        DialogWarningPredictionProject:
            entire dialog which warns the user that the project folder
            has been already created
    """
    projectName = txtProjectNamePrediction.text()
    projectPath = f"{workspacePath}/{projectName}"
    # check if the project folder already exists
    if os.path.exists(projectPath):
        statusBar.showMessage(
            f"Warning! | Current Workspace: {workspacePath}")
        dialog = DialogWarningPredictionProject
        dialog.exec_()
        # breaks out of function
        statusBar.clearMessage()
        return None
    else:
        os.mkdir(projectPath)
        return projectName, projectPath


def set_values_in_project_xml(projectPath, projectName, REFERENCE_DIR,
                              REFERENCE_OBJ_NAME, txtChainRefPrediction,
                              txtChainModelPrediction, resultsPath):
    """This function sets specific values in the project xml file

    Args:
        projectPath:
            path where the project is stored
        projectName:
            name of the project
        REFERENCE_DIR:
            path where the reference pdb file is stored
        REFERENCE_OBJ_NAME:
            name of the reference
        txtChainRefPrediction:
            textbox which contains the chain information of the reference
        txtChainModelPrediction:
            textbox which contains the chain information of the model
        resultsPath:
            path where all results are stored
    """
    fullProjectFileName = f"{projectPath}/project.xml"
    projectFile = tools.ProjectXml(fullProjectFileName)
    projectFile.create_settings_xml_file()
    try:
        tmpProjectFile = projectFile.load_xml_in_memory()
        projectFile.set_value(tmpProjectFile, "projectName", "value",
                              projectName)
        projectFile.set_value(tmpProjectFile, "predictionDone", "value",
                              "True")
        projectFile.set_value(tmpProjectFile, "reference", "value",
                              f"{REFERENCE_DIR}/{REFERENCE_OBJ_NAME}")
        projectFile.set_value(tmpProjectFile, "referenceChains", "value",
                              txtChainRefPrediction.text())
        projectFile.set_value(tmpProjectFile, "modelChains", "value",
                              txtChainModelPrediction.text())
        projectFile.set_value(tmpProjectFile, "results", "value",
                              resultsPath)
    except FileNotFoundError:
        print("Project file could not be loaded!")


def error_dialog(message, message_detail):
    """This function creates an error dialog, which can be customized.

    Args:
        message:
            text which should get displayed in the dialog
        message_detail:
            additional information
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("Error")
    msg.setInformativeText(message)
    msg.setDetailedText(message_detail)
    msg.setWindowTitle("Error")
    msg.setStandardButtons(QMessageBox.Abort)
    msg.exec_()
