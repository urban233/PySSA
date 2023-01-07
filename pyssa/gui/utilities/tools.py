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
"""Module for functions which can be used across the entire project"""

import os
import pathlib
import shutil
import logging
import fnmatch
from Bio import PDB
from pathlib import Path

from pymol import Qt
from pymol import cmd
from PyQt5 import QtGui

from gui.data_structures import settings
from pyssa.gui.ui.dialogs import dialog_settings_global
from pyssa.gui.utilities import styles
from pyssa.gui.utilities import constants
from pyssa.gui.data_structures.data_classes import protein_info


def create_directory(parent_path, dir_name) -> None:
    """This function creates a directory with a given path and directory name

    Args:
        parent_path:
            parent path where the new directory should be created
        dir_name:
            name of the new directory
    """
    new_dir = pathlib.Path(f"{parent_path}/{dir_name}")
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)


def extract_and_move_model_pdb(source_dir, tmp_dir, archive, target_dir) -> None:
    """This function extracts the prediction.zip archive from the AlphaFold
    colab notebook and moves the selected_prediction.pdb to a different
    directory

    Args:
        source_dir (str):
            directory where the prediction.zip archives are stored
        tmp_dir (str):
            directory where the prediction.zip archives should be extracted
        archive (str):
            name of the specific prediction.zip archive WITH the file extension
        target_dir (str):
            directory where all .pdb files should be stored
    """
    # const var for the .pdb file name
    PREDICTION_NAME = "selected_prediction"
    FOLDER_NAME = "prediction"
    # extract the zip archive to the temporary directory
    shutil.unpack_archive(f"{source_dir}/{archive}", tmp_dir)

    # creates a list which contains all files from the target directory
    tmp_list = os.listdir(target_dir)
    # gets the length of the list from above
    index = len(tmp_list)
    # splits the zip archive to get name and extension
    archive_name, extension = os.path.splitext(archive)

    # creates new and unique file name for the .pdb file
    new_filename = f"{tmp_dir}/{FOLDER_NAME}/{PREDICTION_NAME}_{index}.pdb"
    # renaming of the old file
    os.rename(f"{tmp_dir}/{FOLDER_NAME}/{PREDICTION_NAME}.pdb", new_filename)
    # move the .pdb file to the new target directory
    shutil.move(new_filename, target_dir)
    shutil.rmtree(f"{tmp_dir}/{FOLDER_NAME}")


def filter_prediction_zips(path) -> list[str]:
    """This function filters a list of file names to a new list which contains
    only files starting with "prediction"

    Args:
        path (str):
            directory path which should be filtered
    Returns:
        prediction_list:
            prediction_list which is a list of all filenames which starts with
            "prediction"
    """
    prediction_list = []
    # iterate through all file names in the given path
    for file in os.listdir(path):
        # check if the file name starts with "prediction"
        if file.startswith("prediction"):
            prediction_list.append(file)
    return prediction_list


def get_prediction_file_name(path) -> list[str]:
    """This function returns a list with all prediction names.

    Args:
        path:
            path where the predictions are
    Returns:
        prediction_list:
            a list of all filenames which starts with "prediction"
    """
    prediction_list = []
    # iterate through all file names in the given path
    for file in os.listdir(path):
        # check if the file name starts with "prediction"
        if file.startswith("selected"):
            prediction_list.append(file)
    return prediction_list


def get_file_path_and_name(text) -> tuple[str, str]:
    """This function gets the file name and path from a QFileDialog

    Args:
        text (str):
            string of the text box

    Returns:
        file_name:
            name of the file with file extension
        path:
            absolute path where the file is stored
    """
    file_info = Qt.QtCore.QFileInfo(text)
    file_name = file_info.baseName()
    path = file_info.canonicalPath()
    return file_name, path


def safeguard_filepath_xml(settings_file, xml_tag, xml_attribute):
    """This function is a safeguard for a filepath from a xml file.
    IMPORTANT: the parameters needs to be in singe quotes!

    Args:
        settings_file:
            actual settings.xml file
        xml_tag (str):
            the tag used in the xml file
        xml_attribute (str):
            the attribute used in the xml file
    Returns (bool):
        Returns a boolean value which tells if the safeguard test failed
        or succeeded

    TODO: needs to be fixed for the new settings class
    """
    if not os.path.exists(SettingsXml.get_path(settings_file, xml_tag, xml_attribute)):
        logging.critical(f"The {xml_tag} does NOT exist!")
        logging.debug(f"{xml_tag}: {SettingsXml.get_path(settings_file, xml_tag, xml_attribute)}")
        return False
    elif os.path.exists(SettingsXml.get_path(settings_file, xml_tag, xml_attribute)):
        logging.info(
            f"This is the {xml_tag}: {SettingsXml.get_path(settings_file, xml_tag, xml_attribute)}")
        return True


def safeguard_numerical_value_xml(settings_file, xml_tag, xml_attribute, data_type):
    """This function is a safeguard for a numerical value from a xml file.
    IMPORTANT: the parameters needs to be in singe quotes!

    Args:
        settings_file:
            actual settings.xml file
        xml_tag (str):
            the tag used in the xml file
        xml_attribute (str):
            the attribute used in the xml file
        data_type (str):
            defines the used data type for the element like int or float
    Returns (bool):
        Returns a boolean value which tells if the safeguard test failed
        or succeeded

    TODO: needs to be fixed for the new settings class
    """
    if data_type == "int":
        if int(SettingsXml.get_path(settings_file, xml_tag, xml_attribute)) < 0:
            logging.critical(f"The {xml_tag} is negative!")
            logging.debug(f"{xml_tag}: {SettingsXml.get_path(settings_file, xml_tag, xml_attribute)}")
            return False
        elif int(SettingsXml.get_path(settings_file, xml_tag, xml_attribute)) >= 0:
            logging.info(f"This is the {xml_tag}: {SettingsXml.get_path(settings_file, xml_tag, xml_attribute)}")
            return True
    elif data_type == "float":
        if float(SettingsXml.get_path(settings_file, xml_tag, xml_attribute)) < 0:
            logging.critical(f"The {xml_tag} is negative!")
            logging.debug(f"{xml_tag}: {SettingsXml.get_path(settings_file, xml_tag, xml_attribute)}")
            return False
        elif float(SettingsXml.get_path(settings_file, xml_tag, xml_attribute)) >= 0:
            logging.info(f"This is the {xml_tag}: {SettingsXml.get_path(settings_file, xml_tag, xml_attribute)}")
            return True


def restore_default_settings(settings_obj: settings.Settings) -> None:
    """This function creates a settings.xml file which is filled with the
    pre-defined values.

    """
    settings_obj.restore_settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)


def quick_log_and_display(log_type: str, log_message: str,
                          status_bar: Qt.QtWidgets.QStatusBar, status_bar_message: str) -> None:
    """This function is used to quickly log, print and display a message

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


def check_results_for_integrity(workspace_path, project_name) -> (bool, Path):
    """This function checks if the results data is complete or not

        Returns:
            False if data is incomplete and therefore invalid
                Path of the data which is missing
            True if data is complete and therefore valid
    """
    results_dir_struct = [
        Path(f"{workspace_path}/{project_name}"),
        Path(f"{workspace_path}/{project_name}/pdb"),
        Path(f"{workspace_path}/{project_name}/results"),
        Path(f"{workspace_path}/{project_name}/results/alignment_files"),
        Path(f"{workspace_path}/{project_name}/results/distance_csv"),
        Path(f"{workspace_path}/{project_name}/results/images"),
        Path(f"{workspace_path}/{project_name}/results/images/interesting_regions"),
        Path(f"{workspace_path}/{project_name}/results/plots"),
        Path(f"{workspace_path}/{project_name}/results/plots/distance_histogram"),
        Path(f"{workspace_path}/{project_name}/results/plots/distance_plot"),
        Path(f"{workspace_path}/{project_name}/results/sessions/"),
    ]

    is_pdb_path_empty: bool = False
    for path in results_dir_struct:
        content = os.listdir(path)
        if not content:
            if path == results_dir_struct[2]:
                is_pdb_path_empty = True
        else:
            if is_pdb_path_empty:
                # data is incomplete and therefore invalid
                return False, path

    # data is complete and therefore valid
    return True, "Data is valid"


def scan_workspace_for_valid_projects(workspace_path, list_new_projects):
    workspace_projects: list[str] = os.listdir(workspace_path)
    valid_directories = []
    # iterates over possible project directories
    for directory in workspace_projects:
        try:
            directory_content = os.listdir(f"{workspace_path}/{directory}")
            # iterates over the content in a single project directory
            for content in directory_content:
                if content == "project.json":
                    valid_directories.append(directory)
        except NotADirectoryError:
            print(f"This: {directory} is not a directory.")

    valid_directories.sort()
    for project in valid_directories:
        list_new_projects.addItem(project)
    return valid_directories


def scan_workspace_for_non_duplicate_proteins(valid_projects: list, current_project_name: str,
                                              workspace_path: str, list_widget: Qt.QtWidgets.QListWidget) -> tuple[dict, list]:
    """This function scans the workspace directory for protein structures and eliminates all duplicates

    Args:
        valid_projects (list):
            a list of all projects within the workspace
        current_project_name (str):
            name of the currently loaded project
        workspace_path (str):
            path of the current workspace
        list_widget (Qt.QtWidgets.QListWidget)
            list widget which is needed to temporarily store the results from the function "scan_project_for_valid_proteins"

    Returns:
        dict which contains all proteins without duplicates
    """
    """Var: workspace_proteins is a list which contains all proteins from all projects in the workspace"""
    workspace_proteins = []
    protein_names = []
    protein_tuples_notation = []
    for valid_project in valid_projects:
        if valid_project != current_project_name:
            """Var: project_proteins is a list which contains all proteins from a single project"""
            project_proteins = scan_project_for_valid_proteins(f"{workspace_path}/{valid_project}", list_widget)
            list_widget.clear()
            for protein in project_proteins:
                tmp_protein = protein_info.ProteinInfo(protein, f"{workspace_path}/{valid_project}/pdb/{protein}")
                workspace_proteins.append(tmp_protein)
                if tmp_protein.name not in protein_names:
                    protein_names.append(tmp_protein.name)
    # this for-loop is necessary for the creation of the protein dictionary
    for protein in workspace_proteins:
        protein_tuples_notation.append(protein.get_tuple_notation())
    protein_dict = dict(protein_tuples_notation)
    return protein_dict, protein_names


def scan_project_for_valid_proteins(project_path, list_view_project_proteins=None):
    directory = "pdb"
    project_proteins: list[str] = os.listdir(f"{project_path}/{directory}")
    valid_proteins = []
    pattern = "*.pdb"
    # iterates over possible project directories
    if list_view_project_proteins is not None:
        for protein in project_proteins:
            if fnmatch.fnmatch(protein, pattern):
                list_view_project_proteins.addItem(protein)
    return project_proteins


def switch_page(stackedWidget: Qt.QtWidgets.QStackedWidget, lbl_page_title: Qt.QtWidgets.QLabel, index: int, text: str) -> None:
    """This function switches a given stackedWidget page.

    Args:
        stackedWidget:
            QStackedWidget
        lbl_page_title:
            Label which holds the page title
        index:
            the stacked widget page index
        text:
            the page title which should get displayed

    """
    stackedWidget.setCurrentIndex(index)
    lbl_page_title.setText(text)


def get_sequence_from_pdb_file(file_path):
    # You can use a dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    parser_pdb = PDB.PDBParser()
    structure = parser_pdb.get_structure("struct", file_path)
    sequence_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                try:
                    sequence_list.append(d3to1[residue.resname])
                except KeyError:
                    print(f"Residue {residue} is not a valid residue!")
    sequence = ''.join(sequence_list)
    return sequence


def remove_pdb_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)
    else:
        print(f"There is no protein in this project under: {file_path}")


def add_chains_from_pdb_file_to_list(project_path, protein_filename, list_widget=None) -> list:
    cmd.load(f"{project_path}/pdb/{protein_filename}", object="tmp_protein")
    tmp_chains: list = cmd.get_chains("tmp_protein")
    if list_widget is not None:
        list_widget.clear()
        for chain in tmp_chains:
            list_widget.addItem(chain)
    cmd.reinitialize()
    return tmp_chains


def validate_project_name(list_of_projects, txt_for_project_name, lbl_for_status_project_name,
                          btn_for_next_step, cb_for_add_reference=None):
    """This function validates the input of the project name in real-time

    Args:
        list_of_projects:
            list widget which holds all projects from the workspace
        txt_for_project_name:
            line edit widget which is used to enter the project name
        lbl_for_status_project_name:
            label which is used to give feedback if the input is legal or not
        btn_for_next_step:
            push button which is used to execute either the next step or to create a project
        cb_for_add_reference (optional):
            checkbox widget which is used to add a reference

    """
    if list_of_projects.currentItem() is not None:
        list_of_projects.currentItem().setSelected(False)
    # set color for lineEdit
    txt_for_project_name.setStyleSheet("background-color: #FC5457")
    if len(txt_for_project_name.text()) == 0:
        lbl_for_status_project_name.setText("")
        if cb_for_add_reference is not None:
            cb_for_add_reference.setCheckable(False)
            cb_for_add_reference.setStyleSheet("color: #E1E1E1;")
        btn_for_next_step.setEnabled(False)
        styles.color_button_not_ready(btn_for_next_step)
        return
    elif len(txt_for_project_name.text()) > 20:
        lbl_for_status_project_name.setText("Project name is too long (max. 20 characters).")
        btn_for_next_step.setEnabled(False)
        styles.color_button_not_ready(btn_for_next_step)
        return
    else:
        regex = Qt.QtCore.QRegularExpression()
        regex.setPattern("(([a-z])|([A-Z])|([0-9])|(-)|(_)){0,20}")
        validator = QtGui.QRegularExpressionValidator(regex)
        for i in range(len(txt_for_project_name.text())):
            result = validator.validate(txt_for_project_name.text(), i)
            if result[0] > 0:
                txt_for_project_name.setStyleSheet("background-color: #33C065")
                lbl_for_status_project_name.setText("")
                if cb_for_add_reference is not None:
                    cb_for_add_reference.setCheckable(True)
                    cb_for_add_reference.setStyleSheet("color: black;")
                btn_for_next_step.setEnabled(True)
                styles.color_button_ready(btn_for_next_step)
            else:
                txt_for_project_name.setStyleSheet("background-color: #FC5457")
                lbl_for_status_project_name.setText("Invalid character.")
                if cb_for_add_reference is not None:
                    cb_for_add_reference.setCheckable(False)
                    cb_for_add_reference.setStyleSheet("color: #E1E1E1;")
                btn_for_next_step.setEnabled(False)
                styles.color_button_not_ready(btn_for_next_step)
                return
        item = list_of_projects.findItems(txt_for_project_name.text(),
                                          Qt.QtCore.Qt.MatchContains |
                                          Qt.QtCore.Qt.MatchExactly)
        if len(item) != 0:
            list_of_projects.setCurrentItem(item[0])
            txt_for_project_name.setStyleSheet("background-color: #FC5457")
            lbl_for_status_project_name.setText("Project name already exists.")
            if cb_for_add_reference is not None:
                cb_for_add_reference.setCheckable(False)
                cb_for_add_reference.setStyleSheet("color: #E1E1E1;")
            btn_for_next_step.setEnabled(False)
            styles.color_button_not_ready(btn_for_next_step)


def validate_search_input(list_for_projects, txt_for_search, lbl_for_status_search, txt_for_selected_project=None):
    """This function validates the input of the project name in real-time

    Args:
        list_for_projects:
            list widget where all projects from the workspace are stored
        txt_for_search:
            line edit widget which ii used for entering the search term
        lbl_for_status_search:
            label which gives feedback
        txt_for_selected_project:
            line edit widget which is used to display the selected project

    """
    if list_for_projects.currentItem() is not None:
        list_for_projects.currentItem().setSelected(False)
    # set color for lineEdit
    txt_for_search.setStyleSheet("background-color: white")
    if len(txt_for_search.text()) == 0:
        lbl_for_status_search.setText("")
        if txt_for_selected_project is not None:
            txt_for_selected_project.setText("")
        return
    else:
        item = list_for_projects.findItems(txt_for_search.text(),
                                            Qt.QtCore.Qt.MatchContains |
                                            Qt.QtCore.Qt.MatchExactly)
        if len(item) != 0:
            list_for_projects.setCurrentItem(item[0])
            lbl_for_status_search.setText("")
            if txt_for_selected_project is not None:
                txt_for_selected_project.setText(list_for_projects.currentItem().text())
        else:
            txt_for_search.setStyleSheet("background-color: #FC5457")
            lbl_for_status_search.setText("Project name does not exists.")
            if txt_for_selected_project is not None:
                txt_for_selected_project.setText("")


def validate_protein_name(txt_for_protein_name, lbl_for_status_protein_name,
                          btn_next):
    """This function validates the input of the protein name in real-time

    Args:
        txt_for_protein_name:
            line edit widget which ii used for entering the protein name
        lbl_for_status_protein_name:
            label which gives feedback if the protein name is legal
        btn_next:
            push button which is used for the next step

    """
    # set color for lineEdit
    txt_for_protein_name.setStyleSheet("background-color: #FC5457")
    if len(txt_for_protein_name.text()) == 0:
        lbl_for_status_protein_name.setText("")
        btn_next.setEnabled(False)
        styles.color_button_not_ready(btn_next)
        return
    elif len(txt_for_protein_name.text()) > 20:
        lbl_for_status_protein_name.setText("Project name is too long (max. 20 characters).")
        btn_next.setEnabled(False)
        styles.color_button_not_ready(btn_next)
        return
    else:
        regex = Qt.QtCore.QRegularExpression()
        regex.setPattern("(([a-z])|([A-Z])|([0-9])|(-)|(_)){0,20}")
        validator = QtGui.QRegularExpressionValidator(regex)
        for i in range(len(txt_for_protein_name.text())):
            result = validator.validate(txt_for_protein_name.text(), i)
            if result[0] > 0:
                txt_for_protein_name.setStyleSheet("background-color: #33C065")
                lbl_for_status_protein_name.setText("")
                btn_next.setEnabled(True)
                styles.color_button_ready(btn_next)
            else:
                txt_for_protein_name.setStyleSheet("background-color: #FC5457")
                lbl_for_status_protein_name.setText("Invalid character.")
                btn_next.setEnabled(False)
                styles.color_button_not_ready(btn_next)
                return


def validate_protein_sequence(txt_protein_sequence, lbl_status_protein_sequence,
                              btn_next):
    """This function validates the input of the protein sequence in real-time

    Args:
        txt_protein_sequence:
            line edit widget which is used to enter the protein sequence
        lbl_status_protein_sequence:
            label which gives feedback if the protein sequence is legal or not
        btn_next:
            button which is used to get to the next step

    """
    # set color for lineEdit
    txt_protein_sequence.setStyleSheet("background-color: #FC5457")
    if len(txt_protein_sequence.text()) == 0:
        lbl_status_protein_sequence.setText("")
        btn_next.setEnabled(False)
        styles.color_button_not_ready(btn_next)
        return
    else:
        regex = Qt.QtCore.QRegularExpression()
        regex.setPattern("(([A])|([C-I])|([K-N])|([P-T])|([V-W])|([Y]))+")
        validator = QtGui.QRegularExpressionValidator(regex)
        for i in range(len(txt_protein_sequence.text())):
            result = validator.validate(txt_protein_sequence.text(), i)
            if result[0] > 0:
                txt_protein_sequence.setStyleSheet("background-color: #33C065")
                lbl_status_protein_sequence.setText("")
                btn_next.setEnabled(True)
                styles.color_button_ready(btn_next)
            else:
                txt_protein_sequence.setStyleSheet("background-color: #FC5457")
                lbl_status_protein_sequence.setText("Invalid character.")
                btn_next.setEnabled(False)
                styles.color_button_not_ready(btn_next)
                return


def clean_scratch_folder():
    shutil.rmtree(constants.SCRATCH_DIR)
    os.mkdir(constants.SCRATCH_DIR)

# def create_histogram(results_hashtable):
#     y: np.ndarray = results_hashtable.get("distance")
#
#     # max distance value
#     max_distance = np.amax(y)
#
#     # calculate figure size for y direction
#     y_size: float = len(np.arange(0, max_distance, 0.25)) / 1.2
#     FIGURE_SIZE: (float, float) = (11.0, y_size)
#     # create an empty figure with no Axes
#     fig = plt.figure()
#     # create a figure with a single Axes
#     fig, ax = plt.subplots(figsize=FIGURE_SIZE)
#     # creates a basic histogram
#     counts, bins, patches = ax.hist(y, bins=np.arange(0, max_distance, 0.25), orientation="horizontal")
#     # sets the label for the x-axis
#     ax.set_xlabel("Frequency of $\\alpha$-C atoms distance")
#     # sets the label for the y-axis
#     ax.set_ylabel("Distance [$\mathring{A}$ngstrom]")
#     # set coordinates for y-axis label
#     ax.yaxis.set_label_coords(-0.12, 0.5)
#     # hide y-ticks through empty list
#     ax.set_yticks([])
#     # create label bin position
#     bins_centers = 0.15 * np.diff(bins) + bins[:-1]
#
#     i: int = 0
#     for count, y in zip(counts, bins_centers):
#         # define bin label
#         bin_name: str = f"{round(bins[i], 2)} - {round(bins[i + 1], 2)}"
#         # set bin label through annotation
#         ax.annotate(bin_name, xy=(0, y), xytext=(-70, y), textcoords="offset points")
#         i += 1
#     # sets grid
#     ax.grid(True, axis="both")
