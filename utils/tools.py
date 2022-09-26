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
import shutil
import logging

import numpy as np
from matplotlib import pyplot as plt

import utils.settings_utils
from pymol import Qt
from dialogs import dialog_settings_global
from xml.dom import minidom
from utils import project_constants
from utils.settings_utils import SettingsXml


def create_histogram(results_hashtable):
    y: np.ndarray = results_hashtable.get("distance")

    # max distance value
    max_distance = np.amax(y)

    # calculate figure size for y direction
    y_size: float = len(np.arange(0, max_distance, 0.25)) / 1.2
    FIGURE_SIZE: (float, float) = (11.0, y_size)
    # create an empty figure with no Axes
    fig = plt.figure()
    # create a figure with a single Axes
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)
    # creates a basic histogram
    counts, bins, patches = ax.hist(y, bins=np.arange(0, max_distance, 0.25), orientation="horizontal")
    # sets the label for the x-axis
    ax.set_xlabel("Frequency of $\\alpha$-C atoms distance")
    # sets the label for the y-axis
    ax.set_ylabel("Distance [$\mathring{A}$ngstrom]")
    # set coordinates for y-axis label
    ax.yaxis.set_label_coords(-0.12, 0.5)
    # hide y-ticks through empty list
    ax.set_yticks([])
    # create label bin position
    bins_centers = 0.15 * np.diff(bins) + bins[:-1]

    i: int = 0
    for count, y in zip(counts, bins_centers):
        # define bin label
        bin_name: str = f"{round(bins[i], 2)} - {round(bins[i + 1], 2)}"
        # set bin label through annotation
        ax.annotate(bin_name, xy=(0, y), xytext=(-70, y), textcoords="offset points")
        i += 1
    # sets grid
    ax.grid(True, axis="both")


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


def extract_and_move_model_pdb(source_dir, tmp_dir, archive, target_dir):
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
    tmpList = os.listdir(target_dir)
    # gets the length of the list from above
    index = len(tmpList)
    # splits the zip archive to get name and extension
    archive_name, extension = os.path.splitext(archive)

    # creates new and unique file name for the .pdb file
    new_filename = f"{tmp_dir}/{FOLDER_NAME}/{PREDICTION_NAME}_{index}.pdb"
    # renaming of the old file
    os.rename(f"{tmp_dir}/{FOLDER_NAME}/{PREDICTION_NAME}.pdb", new_filename)
    # move the .pdb file to the new target directory
    shutil.move(new_filename, target_dir)
    shutil.rmtree(f"{tmp_dir}/{FOLDER_NAME}")


def filter_prediction_zips(path):
    """This function filters a list of file names to a new list which contains
    only files starting with "prediction"

    Args:
        path (str):
            directory path which should be filtered
    Returns:
        list:
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


def get_prediction_file_name(path):
    """This function returns a list with all prediction names.

    Args:
        path:
            path where the predictions are
    Returns:

    """
    prediction_list = []
    # iterate through all file names in the given path
    for file in os.listdir(path):
        # check if the file name starts with "prediction"
        if file.startswith("selected"):
            prediction_list.append(file)
    return prediction_list


def get_file_path_and_name(text):
    """This function gets the file name and path from a QFileDialog

    Args:
        text (str):
            string of the text box
    """
    file_info = Qt.QtCore.QFileInfo(text)
    OBJ_NAME = file_info.baseName()
    DIR = file_info.canonicalPath()
    return OBJ_NAME, DIR


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


def restore_default_settings():
    """This function creates a settings.xml file which is filled with the
    pre-defined values.

    """
    default_settings = utils.settings_utils.SettingsXml(project_constants.SETTINGS)
    default_settings.create_settings_xml_file()


def open_global_settings():
    """This function opens the global settings dialog

    """
    dialog = dialog_settings_global.DialogSettingsGlobal()
    if not dialog.ERROR:
        dialog.exec_()
    else:
        logging.error("Settings dialog cannot be opened due to an error.")


def quick_log_and_display(log_type: str, log_message: str,
                          status_bar: Qt.QtWidgets.QStatusBar, status_bar_message: str):
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


# TODO: class needs to be tested
class XmlFile:
    pass

    @staticmethod
    def create_xml_file(tag_list, element_list, save_path_of_xml, attribute="value"):
        """This function create the settings xml with the format:

        Args:
            tag_list (list):
                list which contains the xml file tags
            element_list (list):
                list which contains the elements
            save_path_of_xml (str):
                the complete filepath where the xml file should be created
            attribute (str):
                defines the attribute to access the elements


        <?xml version="1.0" ?>
        <root>
            <workspacePath DEFAULT_ATTRIBUTE=DEFAULT_WORKSPACE_PATH/>
            <pdbpath DEFAULT_ATTRIBUTE=DEFAULT_PDB_PATH>
            <zippath DEFAULT_ATTRIBUTE=DEFAULT_ZIP_PATH/>
            <cyclesValue DEFAULT_ATTRIBUTE=DEFAULT_CYCLES_VALUE/>
            <cutoffValue DEFAULT_ATTRIBUTE=DEFAULT_CUTOFF_VALUE/>
        </root>
        """
        if len(tag_list) == 0:
            raise ValueError("The tag_list is empty!")
        if len(element_list) == 0:
            raise ValueError("The element_list is empty!")
        if len(tag_list) != len(element_list):
            raise ValueError("The tag_list and element_list have not the same"
                             "number of elements!")

        root = minidom.Document()
        root_node = root.createElement("root")
        root.appendChild(root_node)

        index = 0
        for tag in tag_list:
            xml_tag = root.createElement(tag)
            xml_tag.setAttribute(attribute, element_list[index])
            root_node.appendChild(xml_tag)
            index += 1

        if not os.path.exists(save_path_of_xml):
            os.mkdir(save_path_of_xml)
        # save xml file to filesystem
        with open(save_path_of_xml, "w") as f:
            f.write(root.toprettyxml())

    @staticmethod
    def load_xml_in_memory(save_path_of_xml):
        """This function loads a xml file into the memory.

        Args:
            save_path_of_xml (str):
                the complete filepath where the xml file should be created
        Note:
            This function should be used once to load the xml file into the
            memory.
        """
        path_as_string = str(save_path_of_xml)
        return minidom.parse(path_as_string)

    @staticmethod
    def get_path(xml_file, tag, attribute):
        """This functions returns the value of the path node.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            tag (str):
                e.g. pdbpath or zippath node
            attribute (str):
                e.g. name
        """
        path_name = xml_file.getElementsByTagName(tag)
        path = path_name[0].getAttribute(attribute)
        return path

    @staticmethod
    def set_value(xml_file, tag, attribute, value):
        """This function changes a specific value in the xml file

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            tag (str):
                 e.g. pdbpath or zippath node
            attribute (str):
                 e.g. name
            value (str):
                 new value which should be set to the attribute
        """
        path_name = xml_file.getElementsByTagName(tag)
        path_name[0].setAttribute(attribute, value)

    @staticmethod
    def save_xml_file(xml_file, save_path_of_xml):
        """This function saves the opened xml file.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            save_path_of_xml (str):
                the complete filepath where the xml file should be created
        """
        with open(save_path_of_xml, "w") as f:
            f.write(xml_file.toxml())
