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

"""Module for handling everything related to the settigs.xml"""

import os
from pathlib import Path
from xml.dom import minidom
from utils import project_constants


class SettingsXml:
    """This class is for the handling of the settings.xml file.

    """
    root = minidom.Document()

    def __init__(self, path_to_xml):
        self.path_to_xml = path_to_xml

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
        if not os.path.exists(project_constants.SETTINGS_DIR):
            os.mkdir(project_constants.SETTINGS_DIR)
        # save xml file to filesystem
        with open(self.path_to_xml, "w", encoding="utf-8") as file:
            file.write(self.root.toprettyxml())

    def load_xml_in_memory(self):
        """This function loads a xml file into the memory.

        Note:
            This function should be used once to load the xml file into the
            memory.
        """
        path_as_string = str(self.path_to_xml)
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

    def save_xml_file(self, xml_file):
        """This function saves the opened xml file. The path of the class will
        be used as save path.

        Args:
            xml_file:
                the xml file which comes from the function load_xml_in_memory
        """
        with open(self.path_to_xml, "w", encoding="utf-8") as file:
            file.write(xml_file.toxml())
