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
"""Module for functions which can be used across the entire project."""
import os
import pathlib
import shutil
import logging
import fnmatch
import zmq
import pygetwindow
import subprocess
from typing import TYPE_CHECKING
from urllib.error import URLError
from urllib.request import urlopen

from PyQt5.QtWidgets import QMessageBox
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from pyssa.internal.data_structures.data_classes import current_session
from pyssa.internal.portal import pymol_io
from pyssa.internal.data_structures import settings
from pyssa.internal.data_structures.data_classes import basic_protein_info
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.io_pyssa import filesystem_io
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.util import constants

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project


def check_internet_connectivity() -> bool:
    """Checks the connection to the internet."""
    timeout: float = 3
    try:
        urlopen("https://www.google.com", timeout=timeout)
    except URLError:
        return False
    else:
        return True


def create_directory(parent_path: str, dir_name: str) -> None:
    """This function creates a directory with a given path and directory name.

    Args:
        parent_path: parent path where the new directory should be created
        dir_name: name of the new directory
    """
    new_dir = pathlib.Path(f"{parent_path}/{dir_name}")
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)


def restore_default_settings(settings_obj: settings.Settings) -> None:
    """This function creates a settings.xml file which is filled with the pre-defined values."""
    settings_obj.restore_settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)


def quick_log_and_display(
    log_type: str,
    log_message: str,
    status_bar: QtWidgets.QStatusBar,
    status_bar_message: str,
) -> None:
    """This function is used to quickly log, print and display a message.

    Args:
        log_type (str):
            defines the log level
        log_message (str):
            the message which gets logged and printed to the console
        status_bar (Qt.QtWidgets.QStatusBar):
            the actual statusbar object
        status_bar_message (str):
            the message which gets displayed through the statusbar
    """
    if log_type == "info":
        logging.info(log_message)
    elif log_type == "warning":
        logging.warning(log_message)
    elif log_type == "error":
        logging.error(log_message)
    elif log_type == "critical":
        logging.critical(log_message)
    elif log_type == "debug":
        logging.debug(log_message)
    else:
        raise ValueError
    print(log_message)
    status_bar.showMessage(status_bar_message)


def scan_workspace_for_valid_projects(
    workspace_path: pathlib.Path,
    list_new_projects: QtWidgets.QListWidget,
) -> list[str]:
    """Scans the workspace for valid projects and adds them to a QListWidget.

    Args:
        workspace_path: the path to the current workspace.
        list_new_projects: a pyqt list widget.
    """
    directory_content = os.listdir(str(workspace_path))
    directory_content.sort()
    xml_files = [file for file in directory_content if file.endswith(".xml")]
    for tmp_project in xml_files:
        list_new_projects.addItem(tmp_project)
    return xml_files


def scan_workspace_for_non_duplicate_proteins(
    workspace_path: pathlib.Path,
) -> tuple[list[basic_protein_info], list[str]]:
    """This function scans the workspace directory for protein structures and eliminates all duplicates.

    Args:
        workspace_path: the path of the current workspace.

    Returns:
        dict which contains all proteins without duplicates
    """
    """Var: workspace_proteins is a list which contains all proteins from all projects in the workspace"""
    protein_names = []
    protein_infos: list[basic_protein_info] = []

    for tmp_project_file in os.listdir(workspace_path):
        """Var: project_proteins is a list which contains all proteins from a single project"""
        if not os.path.isdir(pathlib.Path(f"{workspace_path}/{tmp_project_file}")):
            xml_deserializer = filesystem_io.XmlDeserializer(pathlib.Path(f"{workspace_path}/{tmp_project_file}"))
            project_name = str(tmp_project_file).replace(".xml", "")
            for tmp_protein in xml_deserializer.xml_root.iter(element_names.PROTEIN):
                molecule_object = tmp_protein.attrib[attribute_names.PROTEIN_MOLECULE_OBJECT]
                if molecule_object not in protein_names:
                    protein_names.append(molecule_object)
                    protein_infos.append(
                        basic_protein_info.BasicProteinInfo(
                            molecule_object,
                            tmp_protein.attrib[attribute_names.ID],
                            project_name,
                        ),
                    )
    return protein_infos, protein_names


def scan_project_for_valid_proteins(
    project_path: pathlib.Path,
    list_view_project_proteins=None,  # noqa: ANN001
) -> list[str]:
    """Scans the project for valid proteins.

    Args:
        project_path: the path to the project.
        list_view_project_proteins: a pyqt list widget
    """
    directory = "pdb"
    project_proteins: list[str] = os.listdir(pathlib.Path(f"{project_path}/{directory}"))
    pattern = "*.pdb"
    # iterates over possible project directories
    if list_view_project_proteins is not None:
        for protein in project_proteins:
            if fnmatch.fnmatch(protein, pattern):
                list_view_project_proteins.addItem(protein)
    return project_proteins


def switch_page(
    stacked_widget: QtWidgets.QStackedWidget,
    lbl_page_title: QtWidgets.QLabel,
    index: int,
    text: str,
) -> None:
    """This function switches a given stackedWidget page.

    Args:
        stacked_widget:
            QStackedWidget
        lbl_page_title:
            Label which holds the page title
        index:
            the stacked widget page index
        text:
            the page title which should get displayed

    """
    stacked_widget.setCurrentIndex(index)
    lbl_page_title.setText(text)


def validate_protein_name(
    txt_for_protein_name: QtWidgets.QTextEdit,
    lbl_for_status_protein_name: QtWidgets.QLabel,
    btn_next: QtWidgets.QPushButton,
) -> None:
    """This function validates the input of the protein name in real-time.

    Args:
        txt_for_protein_name: a line edit widget which ii used for entering the protein name
        lbl_for_status_protein_name: a label which gives feedback if the protein name is legal
        btn_next: a push button which is used for the next step
    """
    # set color for lineEdit
    lbl_for_status_protein_name.setStyleSheet("""color: #FC5457; font-size: 10px;""")
    txt_for_protein_name.setStyleSheet("color: #f44336;" "background-color: white;")
    if len(txt_for_protein_name.text()) == 0:
        lbl_for_status_protein_name.setText("")
        btn_next.setEnabled(False)
        styles.color_button_not_ready(btn_next)
    elif len(txt_for_protein_name.text()) > 20:
        lbl_for_status_protein_name.setText("Project name is too long (max. 20 characters).")
        btn_next.setEnabled(False)
        styles.color_button_not_ready(btn_next)
    else:
        regex = QtCore.QRegularExpression()
        regex.setPattern("(([a-z])|([A-Z])|([0-9])|(-)|(_)){0,20}")
        validator = QtGui.QRegularExpressionValidator(regex)
        for i in range(len(txt_for_protein_name.text())):
            result = validator.validate(txt_for_protein_name.text(), i)
            if result[0] > 0:
                txt_for_protein_name.setStyleSheet("color: #000000;" "background-color: white;")
                lbl_for_status_protein_name.setText("")
                btn_next.setEnabled(True)
                styles.color_button_ready(btn_next)
            else:
                txt_for_protein_name.setStyleSheet("color: #f44336;" "background-color: white;")
                lbl_for_status_protein_name.setText("Invalid character.")
                btn_next.setEnabled(False)
                styles.color_button_not_ready(btn_next)
                return


def validate_protein_sequence(
    txt_protein_sequence: QtWidgets.QTextEdit,
    lbl_status_protein_sequence: QtWidgets.QLabel,
    btn_next: QtWidgets.QPushButton,
) -> None:
    """This function validates the input of the protein sequence in real-time.

    Args:
        txt_protein_sequence:a line edit widget which is used to enter the protein sequence
        lbl_status_protein_sequence: a label which gives feedback if the protein sequence is legal or not
        btn_next: a button which is used to get to the next step
    """
    # set color for lineEdit
    lbl_status_protein_sequence.setStyleSheet("""color: #FC5457; font-size: 10px;""")
    txt_protein_sequence.setStyleSheet("color: #f44336;" "background-color: white;")
    if len(txt_protein_sequence.toPlainText()) == 0:
        lbl_status_protein_sequence.setText("")
        btn_next.setEnabled(False)
        styles.color_button_not_ready(btn_next)
    else:
        regex = QtCore.QRegularExpression()
        regex.setPattern("(([A])|([C-I])|([K-N])|([P-T])|([V-W])|([Y]))+")
        validator = QtGui.QRegularExpressionValidator(regex)
        for i in range(len(txt_protein_sequence.toPlainText())):
            result = validator.validate(txt_protein_sequence.toPlainText(), i)
            if result[0] > 0:
                txt_protein_sequence.setStyleSheet("color: #000000;" "background-color: white;")
                lbl_status_protein_sequence.setText("")
                btn_next.setEnabled(True)
                styles.color_button_ready(btn_next)
            else:
                txt_protein_sequence.setStyleSheet("color: #f44336;" "background-color: white;")
                lbl_status_protein_sequence.setText("Invalid character.")
                btn_next.setEnabled(False)
                styles.color_button_not_ready(btn_next)
                return


def clean_scratch_folder() -> None:
    """Deletes the scratch folder and creates a new one."""
    shutil.rmtree(constants.SCRATCH_DIR)
    os.mkdir(constants.SCRATCH_DIR)


def ask_to_save_pymol_session(
    app_project: "project.Project",
    the_current_session: current_session.CurrentSession,
    app_settings: "settings.Settings",
) -> None:
    """This function asks to save the current pymol session.

    Args:
        app_project: the app project.
        the_current_session: the current pymol session, as data class object
        app_settings: the app settings.
    """
    if app_settings.ask_save_pymol_session == 1:
        session_msg = basic_boxes.yes_or_no(
            "PyMOL Session",
            "Do you want to save your current pymol session?",
            QMessageBox.Information,
        )
        if session_msg is True:
            the_current_session.session = pymol_io.convert_pymol_session_to_base64_string(the_current_session.name)
            app_project.save_pymol_session(the_current_session)
