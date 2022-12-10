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
"""Module for functions which reduce code duplicates in the main module"""

import os
import shutil
from pathlib import Path

from PyQt5.QtGui import QIcon
from pymol import Qt
from PyQt5.QtWidgets import QMessageBox
from pyssa.gui.utilities import tools


def fill_combo_box(combo_box, item_list):
    """This function fills a pyqt combobox

    Args:
        combo_box (QtWidgets.QComboBox):
             pyqt combo box which should be filled
        item_list (list):
            list of items which should be placed in the combo box
    """
    for item in item_list:
        combo_box.addItem(item)


def create_directory(parent_path, dir_name):
    """This function creates a directory with a given path and directory name

    Args:
        parent_path:
            parent path where the new directory should be created
        dir_name:
            name of the new directory
    """
    new_dir = f"{parent_path}/{dir_name}"
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)


def choose_directory(self, txt_box_dir):
    """This function is for choosing a filepath with the
    qt-file dialog.

    Example:
        txt_box_dir: self.ui.txt_workspace_dir

    Args:
        self:
            main_window object
        txt_box_dir:
            qt textbox object which holds the current path
    """
    current_file_path = txt_box_dir.text()
    new_file_path = Qt.QtWidgets.QFileDialog.getExistingDirectory(
        self, "Open Directory", current_file_path,
        options=Qt.QtWidgets.QFileDialog.ShowDirsOnly)
    if new_file_path == "":
        txt_box_dir.setText(current_file_path)
    else:
        txt_box_dir.setText(new_file_path)


def critical_message(message, message_detail):
    """This function creates a basic critical message box, which can be customized.

    Args:
        message:
            text which should get displayed in the dialog
        message_detail:
            additional information
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("Critical")
    msg.setInformativeText(message)
    msg.setDetailedText(message_detail)
    msg.setWindowTitle("Critical")
    msg.setStandardButtons(QMessageBox.Abort)
    msg.exec_()


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


def error_dialog_settings(message, message_detail):
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

    open_global_settings_button = msg.addButton("Open Settings", QMessageBox.ActionRole)
    restore_settings_button = msg.addButton("Restore Settings", QMessageBox.ActionRole)
    # TODO:
    #  * Should a help function be implemented?
    # help_button = msg.addButton("Help", QMessageBox.ActionRole)
    msg.exec_()

    # button logic
    if msg.clickedButton() == restore_settings_button:
        tools.restore_default_settings()
    elif msg.clickedButton() == open_global_settings_button:
        tools.open_global_settings()
    # elif msg.clickedButton() == help_button:
    #     webbrowser.open_new("docs/pymol_plugin/build/html/index.html")


def warning_dialog_restore_settings(message, message_detail):
    """This function creates a warning dialog, which can be customized.

    Args:
        message:
            text which should get displayed in the dialog
        message_detail:
            additional information
    Returns (bool):
        returns a bool value whether to restore the settings or not
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Warning")
    msg.setInformativeText(message)
    msg.setDetailedText(message_detail)
    msg.setWindowTitle("Warning")
    yes_button = msg.addButton("Yes", QMessageBox.ActionRole)
    no_button = msg.addButton("No", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == yes_button:
        return True
    elif msg.clickedButton() == no_button:
        return False


def warning_message_prediction_exists(message_detail, path):
    """This function creates a warning message, which can be customized.

    Args:
        message_detail:
            additional information
        path:
            path where the prediction is stored
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Warning")
    msg.setInformativeText("A prediction already exists! Please delete it or "
                           "move it to another location.")
    msg.setDetailedText(message_detail)
    msg.setWindowTitle("Warning")
    cancel_button = msg.addButton("Cancel", QMessageBox.ActionRole)
    delete_button = msg.addButton("Delete", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == cancel_button:
        return False
    if msg.clickedButton() == delete_button:
        os.remove(Path(path))
        return True


def warning_message_project_exists(project_name, message_detail, path):
    """This function creates a warning message, which can be customized.

    Args:
        project_name:
            name of the active project
        message_detail:
            additional information
        path:
            path where the prediction is stored
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Warning")
    msg.setInformativeText(f"A project with the name {project_name} already exists! "
                           f"Please delete it or move it to another location.")
    msg.setDetailedText(message_detail)
    msg.setWindowTitle("Warning")
    cancel_button = msg.addButton("Cancel", QMessageBox.ActionRole)
    delete_button = msg.addButton("Delete", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == cancel_button:
        return False
    if msg.clickedButton() == delete_button:
        shutil.rmtree(Path(path))
        return True


def warning_message_project_gets_deleted() -> bool:
    """This function creates a warning message, which can be customized.

    Args:
        project_name:
            name of the active project
        message_detail:
            additional information
        path:
            path where the prediction is stored
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Warning")
    msg.setInformativeText(f"Are you sure you want to delete this project?")
    msg.setWindowTitle("Warning")
    cancel_button = msg.addButton("Cancel", QMessageBox.ActionRole)
    ok_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == cancel_button:
        return False
    if msg.clickedButton() == ok_button:
        return True


def warning_message_protein_gets_deleted() -> bool:
    """This function creates a warning message, which can be customized.

    Args:
        project_name:
            name of the active project
        message_detail:
            additional information
        path:
            path where the prediction is stored
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Warning")
    msg.setInformativeText(f"Are you sure you want to delete this protein?")
    msg.setWindowTitle("Warning")
    cancel_button = msg.addButton("Cancel", QMessageBox.ActionRole)
    ok_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == cancel_button:
        return False
    if msg.clickedButton() == ok_button:
        return True


def error_project_data_is_invalid(path) -> bool:
    """This function creates an error message, if the project data is invalid.

    Args:
        path:
            path where the data is missing
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("Error")
    msg.setInformativeText(f"Your project data is corrupted.")
    msg.setDetailedText(f"Data is missing in the following path: {path}")
    msg.setWindowTitle("Error")
    # cancel_button = msg.addButton("Cancel", QMessageBox.ActionRole)
    ok_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    # if msg.clickedButton() == cancel_button:
    #     return False
    if msg.clickedButton() == ok_button:
        return True


def warning_switch_pymol_session(message_detail) -> bool:
    """This function creates a warning message box, to inform the
    user to save the current active pymol session.

    Args:
        message_detail:
            detailed message string
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Warning")
    msg.setInformativeText(f"Do you want to save the current PyMOL session?")
    msg.setDetailedText(message_detail)
    msg.setWindowTitle("Warning")
    no_button = msg.addButton("No", QMessageBox.ActionRole)
    yes_button = msg.addButton("Yes", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == no_button:
        return False
    if msg.clickedButton() == yes_button:
        return True


def warning_prediction_is_running(dialog_obj):
    """This function creates a warning message box, to inform the
    user to save the current active pymol session.

    Args:
        message_detail:
            detailed message string
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Warning")
    msg.setInformativeText("The AlphaFold prediction is running ...")
    msg.setWindowTitle("Prediction is running.")
    abort_button = msg.addButton("Abort", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == abort_button:
        dialog_obj.close()
        msg.close()
        return False, msg
    return True, msg


def warning_prediction_is_finished(dialog_obj):
    """This function creates a warning message box, to inform the
    user to save the current active pymol session.

    Args:
        message_detail:
            detailed message string
    """
    # closes the previous message box
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
    msg.setText("Finish")
    msg.setInformativeText("The AlphaFold prediction is finished.")
    msg.setWindowTitle("Prediction is finished.")
    finish_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == finish_button:
        dialog_obj.close()
        msg.close()


def error_prediction_progress_lost() -> bool:
    """This function creates an error message, if the project data is invalid.

    Args:
        path:
            path where the data is missing
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText("Error")
    msg.setInformativeText(f"Your prediction data was not saved!")
    # msg.setDetailedText(f"Due to an error the prediction process was aborted.")
    msg.setWindowTitle("Error")
    ok_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == ok_button:
        return True


def hide_gui_elements(gui_elements: list):
    """This function hides gui elements

    Args:
        gui_elements:
            a list of pyqt gui elements
    """
    for gui_element in gui_elements:
        gui_element.hide()


def show_gui_elements(gui_elements: list):
    """This function shows gui elements

    Args:
        gui_elements:
            a list of pyqt gui elements
    """
    for gui_element in gui_elements:
        gui_element.show()


def disable_text_box(text_box, text_box_label):
    """This function disables a text box and grays out the corresponding label

    Args:
        text_box:
            pyqt line edit
        text_box_label:
            pyqt label which describes the text box
    """
    text_box_label.setStyleSheet("color: #E1E1E1")
    text_box.setStyleSheet("background-color: white")
    text_box.setEnabled(False)


def enable_text_box(text_box, text_box_label):
    """This function enables a text box and colors the corresponding label black

    Args:
        text_box:
            pyqt line edit
        text_box_label:
            pyqt label which describes the text box
    """
    text_box_label.setStyleSheet("color: black")
    text_box.setStyleSheet("background-color: white")
    text_box.setEnabled(True)
