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
import pyssa.gui.data_structures.project

import os
from pathlib import Path
from xml.etree import ElementTree
from xml.dom import minidom


class Settings:
    """This class is for the handling of the settings.xml file.
        Var:
            workspace_path:
                path of the workspace
            prediction_path:
                path where the models from predictions are initially stored
            cycles:
                number of structure alignment cycles
            cutoff:
                cutoff of the structure alignment
        """
    #dir_settings: str, filename: str


    def __init__(self, ) -> None:
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

        self.settings = [
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
        ]

        #self._dir_settings: str = dir_settings
        #self._filename: str = filename
        # uses default values
        # if not os.path.exists(f"{self._dir_settings}/{self._filename}"):
        #     self.save_settings_to_xml()


    def save_settings_to_xml(self):
        root = ElementTree.Element("settings")

        for setting in self.settings:
            tmp_xml_element = ElementTree.Element(setting[0])
            tmp_xml_element.text = str(setting[1])
            root.append(tmp_xml_element)

        # creates a pretty xml file
        xml_as_string = minidom.parseString(ElementTree.tostring(root)).toprettyxml(indent="   ")
        with open("test_settings.xml", "w") as f:
            f.write(xml_as_string)

    def load_settings_from_xml(self, xml_file):
        tree = ElementTree.parse(xml_file)
        root = tree.getroot()
        child_content = []
        for child in root:
            child_content.append(child.text)
        self.settings[0][1] = str(child_content[0])
        self.settings[1][1] = str(child_content[1])
        self.settings[2][1] = int(child_content[2])
        self.settings[3][1] = float(child_content[3])


if __name__ == '__main__':
    settings = Settings()
    settings.save_settings_to_xml()
    settings.load_settings_from_xml("test_settings.xml")



if __name__ == '__main__':
    test_project = pyssa.gui.data_structures.project.Project()
    test_project.__setattr__("project_name", "bmp2")
    print(test_project.__getattribute__("project_name"))
