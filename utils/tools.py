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
import shutil
import logging
from pathlib import Path

from pymol import Qt
from dialogs import dialog_settings_global
from utils import global_utils


def create_directory(parent_path, dir_name) -> None:
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


def restore_default_settings() -> None:
    """This function creates a settings.xml file which is filled with the
    pre-defined values.

    """
    global_utils.global_var_settings_obj.restore_settings()
    global_utils.global_var_settings_obj.save_settings_to_xml()


def open_global_settings() -> None:
    """This function opens the global settings dialog

    """
    dialog = dialog_settings_global.DialogSettingsGlobal()
    if not dialog.ERROR:
        dialog.exec_()
    else:
        logging.error("Settings dialog cannot be opened due to an error.")


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


def scan_workspace_for_vaild_projects(workspace_path, list_new_projects):
    workspace_projects: list[str] = os.listdir(workspace_path)
    valid_directories = []
    # iterates over possible project directories
    for directory in workspace_projects:
        try:
            directory_content = os.listdir(f"{workspace_path}/{directory}")
            # iterates over the content in a single project directory
            for content in directory_content:
                if content == "project.xml":
                    valid_directories.append(directory)
        except NotADirectoryError:
            print(f"This: {directory} is not a directory.")

    valid_directories.sort()
    for project in valid_directories:
        list_new_projects.addItem(project)

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
