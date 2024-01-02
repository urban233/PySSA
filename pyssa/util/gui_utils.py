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
"""Module for functions which reduce code duplicates in the main module."""
import typing
import os
import pathlib

from pymol import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QMessageBox
from PyQt5 import QtCore
from PyQt5 import QtWidgets

from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.util import tools, constants
from pyssa.util.void import rvoid

if typing.TYPE_CHECKING:
    from pyssa.internal.data_structures import project
    from pyssa.internal.data_structures import settings


def fill_combo_box(combo_box: QtWidgets.QComboBox, item_list: list) -> None:
    """This function fills a pyqt combobox.

    Args:
        combo_box (QtWidgets.QComboBox):
             pyqt combo box which should be filled
        item_list (list):
            list of items which should be placed in the combo box
    """
    for item in item_list:
        combo_box.addItem(item)


def create_directory(parent_path: pathlib.Path, dir_name: str) -> None:
    """This function creates a directory with a given path and directory name.

    Args:
        parent_path:
            parent path where the new directory should be created
        dir_name:
            name of the new directory
    """
    new_dir = f"{str(parent_path)}/{dir_name}"
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)


def choose_directory(self, txt_box_dir: QtWidgets.QLineEdit) -> None:  # noqa: ANN001
    """This function is for choosing a filepath with the QFileDialog.

    Example:
        txt_box_dir: self.ui.txt_workspace_dir

    Args:
        self:
            main_window object
        txt_box_dir:
            qt textbox object which holds the current path
    """
    current_file_path = pathlib.Path(txt_box_dir.text())
    new_file_path = pathlib.Path(
        Qt.QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Open Directory",
            str(current_file_path),
            options=Qt.QtWidgets.QFileDialog.ShowDirsOnly,
        ),
    )
    rvoid(os.access(new_file_path, os.W_OK))
    if new_file_path != pathlib.Path(".") and os.access(new_file_path, os.W_OK):
        txt_box_dir.setText(str(new_file_path))
    elif new_file_path != pathlib.Path(".") and not os.access(new_file_path, os.W_OK):
        basic_boxes.ok(
            "Permission error",
            "You do not have write permissions for this directroy. Please choose another one.",
            QMessageBox.Warning,
        )
        txt_box_dir.setText(str(current_file_path))
    else:
        txt_box_dir.setText(str(current_file_path))


def no_internet_dialog() -> None:
    """Displays a message box which informs the user that there is no internet connection."""
    msg = QMessageBox()
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setIcon(QMessageBox.Critical)
    msg.setText("You do not have a working internet connection but it is necessary for this operation!")
    msg.setWindowTitle("Internet connection")
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec_()


def no_internet_dialog_with_custom_msg(message: str) -> None:
    """Displays a custom message box which informs the user that there is no internet connection.

    Args:
        message: a custom message for the dialog.
    """
    msg = QMessageBox()
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setIcon(QMessageBox.Critical)
    msg.setText(message)
    msg.setWindowTitle("Internet connection")
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec_()


def error_dialog_settings(message: str, message_detail: str, settings_obj: "settings.Settings") -> None:
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
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setIcon(QMessageBox.Critical)
    msg.setText(message)
    msg.setInformativeText(message_detail)
    msg.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
    msg.setWindowTitle("Error")

    restore_settings_button = msg.addButton("Restore Settings", QMessageBox.ActionRole)
    msg.exec_()

    # button logic
    if msg.clickedButton() == restore_settings_button:
        tools.restore_default_settings(settings_obj)


def warning_dialog_restore_settings(message: str) -> bool:
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
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setIcon(QMessageBox.Warning)
    msg.setText(message)
    msg.setWindowTitle("Warning")
    yes_button = msg.addButton("Yes", QMessageBox.ActionRole)
    no_button = msg.addButton("No", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == yes_button:
        return True
    elif msg.clickedButton() == no_button:  # noqa: RET505
        return False
    return False


def warning_message_project_gets_deleted() -> bool:
    """This function creates a warning message, which can be customized."""
    msg = QMessageBox()
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Are you sure you want to delete this project?")
    msg.setWindowTitle("Warning")
    cancel_button = msg.addButton("Cancel", QMessageBox.ActionRole)
    ok_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == cancel_button:
        return False
    if msg.clickedButton() == ok_button:
        return True
    return False


def warning_message_protein_gets_deleted() -> bool:
    """This function creates a warning message, which can be customized."""
    msg = QMessageBox()
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Are you sure you want to delete this protein?")
    msg.setWindowTitle("Warning")
    cancel_button = msg.addButton("Cancel", QMessageBox.ActionRole)
    ok_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == cancel_button:
        return False
    if msg.clickedButton() == ok_button:
        return True
    return False


def warning_prediction_is_finished(dialog_obj) -> None:  # noqa: ANN001
    """This function creates a warning message box, to inform the user to save the current active pymol session.

    Args:
        dialog_obj: the web interface object.
    """
    # closes the previous message box
    msg = QMessageBox()
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
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
    """Creates an error message if the project data is invalid."""
    msg = QMessageBox()
    msg.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setIcon(QMessageBox.Critical)
    msg.setText("Error")
    msg.setInformativeText("Your prediction data was not saved!")
    # msg.setDetailedText(f"Due to an error the prediction process was aborted.")
    msg.setWindowTitle("Error")
    ok_button = msg.addButton("OK", QMessageBox.ActionRole)
    msg.exec_()
    # button logic
    if msg.clickedButton() == ok_button:
        return True
    return False


def hide_gui_elements(gui_elements: list) -> None:
    """This function hides gui elements.

    Args:
        gui_elements:
            a list of pyqt gui elements
    """
    for gui_element in gui_elements:
        if gui_element is not None:
            gui_element.hide()


def show_gui_elements(gui_elements: list) -> None:
    """This function shows gui elements.

    Args:
        gui_elements:
            a list of pyqt gui elements
    """
    for gui_element in gui_elements:
        if gui_element is not None:
            gui_element.show()


def manage_gui_visibility(gui_elements_to_show: list, gui_elements_to_hide: list) -> None:
    """Manages a combination of "show_gui_elements" and "hide_gui_elements" manage the visibility of gui elements.

    Args:
        gui_elements_to_show:
            list which contains the gui elements which should be displayed
        gui_elements_to_hide:
            list which contains the gui elements which should be hidden
    """
    show_gui_elements(gui_elements_to_show)
    hide_gui_elements(gui_elements_to_hide)


def disable_text_box(text_box: QtWidgets.QLineEdit, text_box_label: QtWidgets.QLabel) -> None:
    """This function disables a text box and grays out the corresponding label.

    Args:
        text_box:
            pyqt line edit
        text_box_label:
            pyqt label which describes the text box
    """
    text_box_label.setStyleSheet("color: #E1E1E1")
    text_box.setStyleSheet("background-color: white")
    text_box.setEnabled(False)


def enable_text_box(text_box: QtWidgets.QLineEdit, text_box_label: QtWidgets.QLabel) -> None:
    """This function enables a text box and colors the corresponding label black.

    Args:
        text_box:
            pyqt line edit
        text_box_label:
            pyqt label which describes the text box
    """
    text_box_label.setStyleSheet("color: black")
    text_box.setStyleSheet("background-color: white")
    text_box.setEnabled(True)


def fill_list_view_with_protein_names(
    app_project: "project.Project",
    list_view_project_proteins: QtWidgets.QListWidget,
) -> None:
    """Fills a list with protein names.

    Args:
        app_project: the app project.
        list_view_project_proteins: a QListWidget for the protein names.
    """
    for tmp_protein in app_project.proteins:
        list_view_project_proteins.addItem(tmp_protein.get_molecule_object())


def fill_list_view_with_protein_pair_names(
    app_project: "project.Project",
    list_view_project_proteins: QtWidgets.QListWidget,
) -> None:
    """Fills a list with protein pair names.

    Args:
        app_project: the app project.
        list_view_project_proteins: a QListWidget for the protein pair names.
    """
    for tmp_protein_pair in app_project.protein_pairs:
        list_view_project_proteins.addItem(tmp_protein_pair.name)


def setup_standard_block_box(block_box: QMessageBox, window_title: str, msg_text: str) -> QMessageBox:
    """Sets up the default values for a given block box.

    Args:
        block_box: the block box to configure.
        window_title: the window title for the block box.
        msg_text: the message for the block box.
    """
    block_box.setStandardButtons(QMessageBox.NoButton)
    block_box.setIcon(QMessageBox.Information)
    block_box.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(block_box)
    block_box.setWindowTitle(window_title)
    block_box.setText(msg_text)
    return block_box
