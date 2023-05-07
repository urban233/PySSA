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
import typing
import os
import pathlib
import shutil
from pathlib import Path
from pymol import Qt
from PyQt5.QtWidgets import QMessageBox
from PyQt5 import QtCore
from pyssa.util import tools

if typing.TYPE_CHECKING:
    from pyssa.internal.data_structures import project


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
    current_file_path = pathlib.Path(txt_box_dir.text())
    new_file_path = pathlib.Path(Qt.QtWidgets.QFileDialog.getExistingDirectory(
        self, "Open Directory", str(current_file_path),
        options=Qt.QtWidgets.QFileDialog.ShowDirsOnly))
    if new_file_path == "":
        txt_box_dir.setText(str(current_file_path))
    else:
        txt_box_dir.setText(str(new_file_path))


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


def error_dialog_settings(message, message_detail, settings_obj):
    """This function creates an error dialog, which can be customized.

    Args:
        message:
            text which should get displayed in the dialog
        message_detail:
            additional information
        settings_obj:
            settings object of the MainWindow class
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Critical)
    msg.setText(message)
    msg.setInformativeText(message_detail)
    #msg.setDetailedText()
    msg.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
    msg.setWindowTitle("Error")

    #open_global_settings_button = msg.addButton("Open Settings", QMessageBox.ActionRole)
    restore_settings_button = msg.addButton("Restore Settings", QMessageBox.ActionRole)
    msg.exec_()

    # button logic
    if msg.clickedButton() == restore_settings_button:
        tools.restore_default_settings(settings_obj)
    #elif msg.clickedButton() == open_global_settings_button:
        #tools.open_global_settings()
    # elif msg.clickedButton() == help_button:
    #     webbrowser.open_new("docs/pyssa/build/html/index.html")


def warning_dialog_restore_settings(message):
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
    msg.setText(message)
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
    msg.setText(f"Are you sure you want to delete this project?")
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
    msg.setText(f"Are you sure you want to delete this protein?")
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
        if gui_element is not None:
            gui_element.hide()


def show_gui_elements(gui_elements: list):
    """This function shows gui elements

    Args:
        gui_elements:
            a list of pyqt gui elements
    """
    for gui_element in gui_elements:
        if gui_element is not None:
            gui_element.show()


def manage_gui_visibility(gui_elements_to_show: list, gui_elements_to_hide: list):
    """This function is a combination of "show_gui_elements" and "hide_gui_elements" to quickly
    manage the visibility of gui elements

    Args:
        gui_elements_to_show:
            list which contains the gui elements which should be displayed
        gui_elements_to_hide:
            list which contains the gui elements which should be hidden
    """
    show_gui_elements(gui_elements_to_show)
    hide_gui_elements(gui_elements_to_hide)


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


def get_prediction_name_and_seq_from_table(table) -> list[tuple[str, str]]:
    """This function gets the names and sequences of the table which stores the predictions to run

    Args:
        table:

    Returns:
        list of tuples with name and sequence
    """
    # list which consists of tuples of the protein name and protein sequence
    predictions: list[tuple[str, str]] = []
    for i in range(table.rowCount()):
        tmp_name = table.verticalHeaderItem(i).text()
        tmp_seq = table.item(i, 1).text()
        predictions.append((tmp_name, tmp_seq))
    return predictions


def write_fasta_file_from_table():
    pass


def fill_list_view_with_protein_names(app_project: 'project.Project', list_view_project_proteins):
    for tmp_protein in app_project.proteins:
        list_view_project_proteins.addItem(tmp_protein.get_molecule_object())
