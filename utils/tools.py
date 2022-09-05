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
from typing import Any

from pymol import Qt

import os
import shutil
import logging
from utils import tools
from dialogs import dialog_settings_global
from xml.dom import minidom
from datetime import date
from utils import constants
from pathlib import Path


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
    default_settings = tools.SettingsXml(constants.SETTINGS)
    default_settings.create_settings_xml_file()


def open_global_settings():
    """This function opens the global settings dialog

    """
    dialog = dialog_settings_global.DialogSettingsGlobal()
    if not dialog.ERROR:
        dialog.exec_()
    else:
        logging.error("Settings dialog cannot be opened due to an error.")


class SettingsXml:
    root = minidom.Document()

    def __init__(self, pathToXml):
        self.pathToXml = pathToXml

    def create_settings_xml_file(self):
        """This function create the settings xml with the format:

        <?xml version="1.0" ?>
        <root>
            <workspacePath DEFAULT_ATTRIBUTE=DEFAULT_WORKSPACE_PATH/>
            <pdbpath DEFAULT_ATTRIBUTE=DEFAULT_PDB_PATH>
            <zippath DEFAULT_ATTRIBUTE=DEFAULT_ZIP_PATH/>
            <cyclesValue DEFAULT_ATTRIBUTE=DEFAULT_CYCLES_VALUE/>
            <cutoffValue DEFAULT_ATTRIBUTE=DEFAULT_CUTOFF_VALUE/>
        </root>
        """
        DEFAULT_WORKSPACE_PATH = str(Path(f"{os.path.expanduser('~')}/Documents"))
        DEFAULT_PDB_PATH = str(Path(f"{os.path.expanduser('~')}/Documents"))
        DEFAULT_ZIP_PATH = str(Path(f"{os.path.expanduser('~')}/Downloads"))
        DEFAULT_CYCLES_VALUE = "1"
        DEFAULT_CUTOFF_VALUE = "1.0"

        DEFAULT_ATTRIBUTE = "value"

        root_node = self.root.createElement("root")
        self.root.appendChild(root_node)

        workspace_path_node = self.root.createElement("workspacePath")
        # init node/attribute with default values
        workspace_path_node.setAttribute(DEFAULT_ATTRIBUTE,
                                         DEFAULT_WORKSPACE_PATH)
        root_node.appendChild(workspace_path_node)

        pdb_path_node = self.root.createElement("pdbPath")
        # init node/attribute with default values
        pdb_path_node.setAttribute(DEFAULT_ATTRIBUTE, DEFAULT_PDB_PATH)
        root_node.appendChild(pdb_path_node)

        zip_path_node = self.root.createElement("zipPath")
        # init node/attribute with default values
        zip_path_node.setAttribute(DEFAULT_ATTRIBUTE, DEFAULT_ZIP_PATH)
        root_node.appendChild(zip_path_node)

        cycles_value_node = self.root.createElement("cyclesValue")
        # init node/attribute with default values
        cycles_value_node.setAttribute(DEFAULT_ATTRIBUTE, DEFAULT_CYCLES_VALUE)
        root_node.appendChild(cycles_value_node)

        cutoff_value_node = self.root.createElement("cutoffValue")
        # init node/attribute with default values
        cutoff_value_node.setAttribute(DEFAULT_ATTRIBUTE, DEFAULT_CUTOFF_VALUE)
        root_node.appendChild(cutoff_value_node)

        # if os.path.exists(f"{self.pathToXml}"):
        #    os.remove(self.pathToXml)
        if not os.path.exists(constants.SETTINGS_DIR):
            os.mkdir(constants.SETTINGS_DIR)
        # save xml file to filesystem
        with open(self.pathToXml, "w") as f:
            f.write(self.root.toprettyxml())

    def load_xml_in_memory(self):
        """This function loads a xml file into the memory.

        Note:
            This function should be used once to load the xml file into the
            memory.
        """
        path_as_string = str(self.pathToXml)
        return minidom.parse(path_as_string)

    @staticmethod
    def get_path(xml_file, Tag, attribute):
        """This functions returns the value of the path node.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            Tag (str):
                e.g. pdbpath or zippath node
            attribute (str):
                e.g. name
        """
        path_name = xml_file.getElementsByTagName(Tag)
        path = path_name[0].getAttribute(attribute)
        return path

    @staticmethod
    def set_value(xml_file, Tag, attribute, value):
        """This function changes a specific value in the xml file

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            Tag (str):
                 e.g. pdbpath or zippath node
            attribute (str):
                 e.g. name
            value (str):
                 new value which should be set to the attribute
        """
        pathName = xml_file.getElementsByTagName(Tag)
        pathName[0].setAttribute(attribute, value)

    def save_xml_file(self, xml_file):
        """This function saves the opened xml file. The path of the class will
        be used as save path.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
        """
        with open(self.pathToXml, "w") as f:
            f.write(xml_file.toxml())


class ProjectXml:
    root = minidom.Document()

    def __init__(self, pathToXml):
        self.pathToXml = pathToXml

    def create_project_xml_file(self):
        """This function create the settings xml with the format:

        <?xml version="1.0" ?>
        <root>
            <date DEFAULT_ATTRIBUTE=DATE/>
            <projectName DEFAULT_ATTRIBUTE=PROJECT_NAME/>
            <predictionDone DEFAULT_ATTRIBUTE=PREDICTION_DONE/>
            <reference DEFAULT_ATTRIBUTE=REFERENCE/>
            <model DEFAULT_ATTRIBUTE=REFERENCE/>
            <referenceChains DEFAULT_ATTRIBUTE=REFERENCE/>
            <modelChains DEFAULT_ATTRIBUTE=REFERENCE/>
            <results DEFAULT_ATTRIBUTE=REFERENCE/>
        </root>
        """
        DATE = str(date.today())
        PROJECT_NAME = "projectname"
        PREDICTION_DONE = "False"
        REFERENCE = "reference"
        MODEL = "model"
        REFERENCE_CHAINS = ""
        MODEL_CHAINS = ""
        RESULTS = "/path/to/results/directory"
        DEFAULT_ATTRIBUTE = "value"

        root_node = self.root.createElement("root")
        self.root.appendChild(root_node)

        date_node = self.root.createElement("date")
        # init node/attribute with default values
        date_node.setAttribute(DEFAULT_ATTRIBUTE,
                               DATE)
        root_node.appendChild(date_node)

        project_name_node = self.root.createElement("projectName")
        # init node/attribute with default values
        project_name_node.setAttribute(DEFAULT_ATTRIBUTE, PROJECT_NAME)
        root_node.appendChild(project_name_node)

        prediction_done_node = self.root.createElement("predictionDone")
        # init node/attribute with default values
        prediction_done_node.setAttribute(DEFAULT_ATTRIBUTE, PREDICTION_DONE)
        root_node.appendChild(prediction_done_node)

        reference_node = self.root.createElement("reference")
        # init node/attribute with default values
        reference_node.setAttribute(DEFAULT_ATTRIBUTE, REFERENCE)
        root_node.appendChild(reference_node)

        model_node = self.root.createElement("model")
        # init node/attribute with default values
        model_node.setAttribute(DEFAULT_ATTRIBUTE, MODEL)
        root_node.appendChild(model_node)

        reference_chains_node = self.root.createElement("referenceChains")
        # init node/attribute with default values
        reference_chains_node.setAttribute(DEFAULT_ATTRIBUTE, REFERENCE_CHAINS)
        root_node.appendChild(reference_chains_node)

        model_chains_node = self.root.createElement("modelChains")
        # init node/attribute with default values
        model_chains_node.setAttribute(DEFAULT_ATTRIBUTE, MODEL_CHAINS)
        root_node.appendChild(model_chains_node)

        results_node = self.root.createElement("results")
        # init node/attribute with default values
        results_node.setAttribute(DEFAULT_ATTRIBUTE, RESULTS)
        root_node.appendChild(results_node)

        if os.path.exists(f"{self.pathToXml}"):
            os.remove(self.pathToXml)
        # save xml file to filesystem
        with open(self.pathToXml, "w") as f:
            f.write(self.root.toprettyxml())

    def load_xml_in_memory(self):
        """This function loads a xml file into the memory.

        Note:
            This function should be used once to load the xml file into the
            memory.
        """
        return minidom.parse(self.pathToXml)

    @staticmethod
    def get_path(xml_file, Tag, attribute):
        """This functions returns the value of the path node.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            Tag (str):
                e.g. pdbpath or zippath node
            attribute (str):
                e.g. name
        """
        path_name = xml_file.getElementsByTagName(Tag)
        path = path_name[0].getAttribute(attribute)
        return path

    @staticmethod
    def set_value(xml_file, Tag, attribute, value):
        """This function changes a specific value in the xml file

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
            Tag (str):
                 e.g. pdbpath or zippath node
            attribute (str):
                 e.g. name
            value (str):
                 new value which should be set to the attribute
        """
        pathName = xml_file.getElementsByTagName(Tag)
        pathName[0].setAttribute(attribute, value)

    def save_xml_file(self, xml_file):
        """This function saves the opened xml file. The path of the class will
        be used as save path.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
        """
        with open(self.pathToXml, "w") as f:
            f.write(xml_file.toxml())


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


class Project:
    """This class is for the projects used in the plugin
    Args:
        project_name:
            name of the project
        pdb_file:
            path of the loaded pdb file
        pdb_id:
            name of the PDB ID
        ref_chains:
            chains which are used for the reference protein
        model_chains:
            chains which are used for the model protein
        """
    project_name = ""
    pdb_file = ""  # either file or id!
    pdb_id = ""
    ref_chains = ""
    model_chains = ""


    def __setattr__(self, name: str, value: Any) -> None:
        super().__setattr__(name, value)

    def __getattribute__(self, name: str) -> Any:
        return super().__getattribute__(name)
