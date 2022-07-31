from pymol import Qt

import os
import shutil
from xml.dom import minidom
from datetime import date


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
    os.rename(f"{tmp_dir}/{FOLDER_NAME}/{PREDICTION_NAME}.pdb",new_filename)
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
        DEFAULT_WORKSPACE_PATH = "/home/$USER/Documents"
        DEFAULT_PDB_PATH = "/home/$USER/Documents"
        DEFAULT_ZIP_PATH = "/home/$USER/Downloads"
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


class ProjectXml:
    root = minidom.Document()

    def __init__(self, pathToXml):
        self.pathToXml = pathToXml

    def create_settings_xml_file(self):
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