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
"""Module for handling everything related to the settings.xml"""

import os
from pathlib import Path
from xml.etree import ElementTree
from xml.dom import minidom


class Settings:
    """This class is for the handling of the settings.xml file.
    Var:
        workspace_path:
            path of the workspace
        model_path:
            path where models are stored
        prediction_path:
            path where the models from predictions are initially stored
        cycles:
            number of structure alignment cycles
        cutoff:
            cutoff of the structure alignment
    """
    def __init__(self, dir_settings: str, filename: str) -> None:
        """Constructor

        Args:
            dir_settings:
                directory where the settings.xml is stored
            filename:
                name of the settings.xml
        """
        self._workspace_path = str(Path(f"{os.path.expanduser('~')}/Documents"))
        self._prediction_path = str(Path(f"{os.path.expanduser('~')}/Downloads"))
        self._cycles: int = 0
        self._cutoff: float = 1.0
        self._app_launch = 1
        self._dir_settings: str = dir_settings
        self._filename: str = filename

        # uses default values
        # self.restore_settings()
        if not os.path.exists(f"{self._dir_settings}/{self._filename}"):
            self._app_launch = 0
            self.construct_xml_element_pairs()
            self.save_settings_to_xml()

    def construct_xml_element_pairs(self):
        settings = [
            [
                "workspace_path", self._workspace_path
            ],
            [
                "prediction_path", self._prediction_path
            ],
            [
                "cycles", self._cycles
            ],
            [
                "cutoff", self._cutoff
            ],
            [
                "first_app_launch", self._app_launch
            ],
        ]
        return settings

    def get_workspace_path(self) -> str:
        """This function gets the value of the workspace_path variable

        Returns (str):
            workspace_path
        """
        return self._workspace_path

    def set_workspace_path(self, value) -> None:
        """This function gets the value of the workspace_path variable

        """
        self._workspace_path = value

    def get_prediction_path(self) -> str:
        """This function gets the value of the prediction_path variable

        Returns (str):
            prediction_path
        """
        return self._prediction_path

    def set_prediction_path(self, value) -> None:
        """This function gets the value of the prediction_path variable

        """
        self._prediction_path = value

    def get_cycles(self) -> int:
        """This function gets the value of the cycles variable

        Returns (int):
            cycles
        """
        return int(self._cycles)

    def set_cycles(self, value) -> None:
        """This function gets the value of the cycles variable

        """
        self._cycles = value

    def get_cutoff(self) -> float:
        """This function gets the value of the cutoff variable

        Returns (float):
            cutoff
        """
        return float(self._cutoff)

    def set_cutoff(self, value) -> None:
        """This function gets the value of the cutoff variable

        """
        self._cutoff = value

    def get_app_launch(self) -> int:
        """This function gets the value of the app_launch variable

       Returns (int):
           app_launch
       """
        return int(self._app_launch)

    def set_app_launch(self, value) -> None:
        """This function gets the value of the app_launch variable

        """
        self._app_launch = value

    def save_settings_to_xml(self):
        root = ElementTree.Element("settings")
        settings = self.construct_xml_element_pairs()

        for setting in settings:
            tmp_xml_element = ElementTree.Element(setting[0])
            tmp_xml_element.text = str(setting[1])
            root.append(tmp_xml_element)

        # creates a pretty xml file
        xml_as_string = minidom.parseString(ElementTree.tostring(root)).toprettyxml(indent="   ")
        with open(f"{self._dir_settings}/{self._filename}", "w") as f:
            f.write(xml_as_string)

    def load_settings_from_xml(self):
        tree = ElementTree.parse(f"{self._dir_settings}/{self._filename}")
        root = tree.getroot()
        child_content = []
        for child in root:
            child_content.append(child.text)
        self.set_workspace_path(str(child_content[0]))
        self.set_prediction_path(str(child_content[1]))
        self.set_cycles(int(child_content[2]))
        self.set_cutoff(float(child_content[3]))

    def delete_settings_xml(self) -> None:
        os.remove(f"{self._dir_settings}/{self._filename}")

        #settings = self.construct_xml_element_pairs()
        # settings[0][1] = str(child_content[0])
        # settings[1][1] = str(child_content[1])
        # settings[2][1] = int(child_content[2])
        # settings[3][1] = float(child_content[3])

    # def save_settings_to_xml(self) -> None:
    #     """This function saves the xml file to the default location which is defined in the
    #     settings object
    #
    #     """
    #     DEFAULT_ATTRIBUTE = "value"
    #
    #     root = minidom.Document()
    #     root_node = root.createElement("root")
    #     root.appendChild(root_node)
    #
    #     workspace_path_node = root.createElement("workspacePath")
    #     # init node/attribute with default values
    #     workspace_path_node.setAttribute(DEFAULT_ATTRIBUTE,
    #                                      self._workspace_path)
    #     root_node.appendChild(workspace_path_node)
    #
    #     pdb_path_node = root.createElement("pdbPath")
    #     # init node/attribute with default values
    #     pdb_path_node.setAttribute(DEFAULT_ATTRIBUTE, self._model_path)
    #     root_node.appendChild(pdb_path_node)
    #
    #     zip_path_node = root.createElement("zipPath")
    #     # init node/attribute with default values
    #     zip_path_node.setAttribute(DEFAULT_ATTRIBUTE, self._prediction_path)
    #     root_node.appendChild(zip_path_node)
    #
    #     cycles_value_node = root.createElement("cyclesValue")
    #     # init node/attribute with default values
    #     cycles_value_node.setAttribute(DEFAULT_ATTRIBUTE, str(self._cycles))
    #     root_node.appendChild(cycles_value_node)
    #
    #     cutoff_value_node = root.createElement("cutoffValue")
    #     # init node/attribute with default values
    #     cutoff_value_node.setAttribute(DEFAULT_ATTRIBUTE, str(self._cutoff))
    #     root_node.appendChild(cutoff_value_node)
    #
    #     if not os.path.exists(self._dir_settings):
    #         os.mkdir(self._dir_settings)
    #     # save xml file to filesystem
    #     with open(f"{self._dir_settings}/{self._filename}", "w", encoding="utf-8") as file:
    #         file.write(root.toprettyxml())
    #         file.close()
    #
    # def load_settings_from_xml(self) -> None:
    #     """This function loads the information from the xml file into the xml object
    #
    #     """
    #     xml_file = minidom.parse(f"{self._dir_settings}/{self._filename}")
    #     self.set_workspace_path(xml_file.getElementsByTagName('workspacePath')[0].getAttribute('value'))
    #     self.set_model_path(xml_file.getElementsByTagName('pdbPath')[0].getAttribute('value'))
    #     self.set_prediction_path(xml_file.getElementsByTagName('zipPath')[0].getAttribute('value'))
    #     self.set_cycles(xml_file.getElementsByTagName('cyclesValue')[0].getAttribute('value'))
    #     self.set_cutoff(xml_file.getElementsByTagName('cutoffValue')[0].getAttribute('value'))

    # def restore_settings(self):
    #     """This function resets the settings to the default values
    #
    #     """
    #     self._workspace_path = str(Path(f"{os.path.expanduser('~')}/Documents"))
    #     self._model_path = str(Path(f"{os.path.expanduser('~')}/Documents"))
    #     self._prediction_path = str(Path(f"{os.path.expanduser('~')}/Downloads"))
    #     self._cycles: int = 0
    #     self._cutoff: float = 1.0
