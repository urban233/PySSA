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
"""Module for the project class"""

import os
import platform
from datetime import datetime
from pathlib import Path
from xml.etree import ElementTree
from xml.dom import minidom
from pyssa.gui.utilities import gui_utils


class Project:
    """This class is for the projects used in the plugin
    Var:
        project_name:
            name of the project
        ref_chains:
            chains which are used for the reference Protein
        model_chains:
            chains which are used for the model Protein
        results_path:
            path where all results are stored
    """
    _date = datetime.now()
    _operating_system = platform.system()
    _workspace = ""
    _project_name: str = ""
    _ref_chains: str = "no_chains_selected"
    _model_chains: str = "no_chains_selected"
    _folder_paths: list[Path] = []
    _session_file_name: str = "session_file_model_s.pse"

    def __init__(self, project_name, workspace_path) -> None:
        """Constructor

        Args:
            project_name:
                name of the project
            workspace_path:
                path of the workspace
        """
        project_name_with_underscores = project_name.replace(" ", "_")
        self._project_name = project_name_with_underscores
        self._workspace = workspace_path

        self._folder_paths = [
            Path(f"{workspace_path}/{project_name_with_underscores}"),
            Path(f"{workspace_path}/{project_name_with_underscores}/pdb"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/alignment_files"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/distance_csv"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/images"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/images/interesting_regions"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/plots"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/plots/distance_histogram"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/plots/distance_plot"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/sessions/"),
        ]

    def create_project_tree(self) -> bool:
        """This function creates the directory structure for a new project.

        """
        # check if the project folder already exists
        if os.path.exists(self.get_project_path()):
            detailed_message = f"The project is located under: {self.get_project_path()}."
            flag = gui_utils.warning_message_project_exists(self._project_name, detailed_message,
                                                            self.get_project_path())
            if flag is False:
                return False

        for folder in self._folder_paths:
            os.mkdir(folder)

    def get_session_file(self) -> str:
        """This function gets the value of the session_file variable

        Returns (str):
            session_file as complete path with file name
        """
        return f"{self._folder_paths[10]}/{self._session_file_name}"

    def get_project_name(self) -> str:
        """This function gets the value of the project_name variable

        Returns (str):
            project_name
        """
        return self._project_name

    def set_project_name(self, value: str) -> None:
        """This function gets the value of the project_name variable

        """
        self._project_name = value

    def get_ref_chains(self) -> str:
        """This function gets the value of the ref_chains variable

        Returns (str):
            project_name
        """
        return self._ref_chains

    def set_ref_chains(self, value: str) -> None:
        """This function gets the value of the ref_chains variable

        """
        self._ref_chains = value

    def get_model_chains(self) -> str:
        """This function gets the value of the model_chains variable

        Returns (str):
            project_name
        """
        return self._model_chains

    def set_model_chains(self, value: str) -> None:
        """This function gets the value of the model_chains variable

        """
        self._model_chains = value

    def get_project_path(self) -> str:
        """This function returns the project path of the project

        Returns:
            the project path
        """
        return str(self._folder_paths[0])

    def get_pdb_path(self) -> str:
        """This function returns the pdb path of the project

        Returns:
            the pdb path
        """
        return str(self._folder_paths[1])

    def construct_xml_element_pairs(self):
        project = [
            [
                "date", self._date
            ],
            [
                "operating_system", self._operating_system
            ],
            [
                "workspace", self._workspace
            ],
            [
                "project_name", self._project_name
            ],
            [
                "reference_chains", self._ref_chains
            ],
            [
                "model_chains", self._model_chains
            ],
        ]
        return project

    def save_project_to_xml(self):
        root = ElementTree.Element("project")
        settings = self.construct_xml_element_pairs()

        for setting in settings:
            tmp_xml_element = ElementTree.Element(setting[0])
            tmp_xml_element.text = str(setting[1])
            root.append(tmp_xml_element)

        # creates a pretty xml file
        xml_as_string = minidom.parseString(ElementTree.tostring(root)).toprettyxml(indent="   ")
        with open(f"{self._folder_paths[0]}/project.xml", "w") as f:
            f.write(xml_as_string)

    # def create_xml_file(self, attribute="value") -> None:
    #     """This function create the settings xml with the format:
    #
    #     Args:
    #         attribute (optional str):
    #             defines the attribute to access the elements
    #
    #     <?xml version="1.0" ?>
    #     <root>
    #         <project_name DEFAULT_ATTRIBUTE=project_name_value/>
    #         <pdb_file DEFAULT_ATTRIBUTE=pdb_file_path>
    #         <pdb_id DEFAULT_ATTRIBUTE=pdb_id_value/>
    #         <pdb_models DEFAULT_ATTRIBUTE=pdb_models_value/>
    #         <ref_chains DEFAULT_ATTRIBUTE=ref_chains_value/>
    #         <model_chains DEFAULT_ATTRIBUTE=model_chains_value/>
    #         <results_path DEFAULT_ATTRIBUTE=results_path_value/>
    #     </root>
    #     """
    #     root = minidom.Document()
    #     root_node = root.createElement("root")
    #     root.appendChild(root_node)
    #
    #     job_name_path_node = root.createElement("job_name")
    #     job_name_path_node.setAttribute(attribute, self._job_name)
    #     root_node.appendChild(job_name_path_node)
    #
    #     project_name_path_node = root.createElement("project_name")
    #     project_name_path_node.setAttribute(attribute, self._project_name)
    #     root_node.appendChild(project_name_path_node)
    #
    #     pdb_file_path_node = root.createElement("pdb_file")
    #     pdb_file_path_node.setAttribute(attribute, self._pdb_file)
    #     root_node.appendChild(pdb_file_path_node)
    #
    #     pdb_id_path_node = root.createElement("pdb_id")
    #     pdb_id_path_node.setAttribute(attribute, self._pdb_id)
    #     root_node.appendChild(pdb_id_path_node)
    #
    #     pdb_models_path_node = root.createElement("pdb_model")
    #     pdb_models_path_node.setAttribute(attribute, self.get_pdb_models())
    #     root_node.appendChild(pdb_models_path_node)
    #
    #     ref_chains_path_node = root.createElement("ref_chains")
    #     ref_chains_path_node.setAttribute(attribute, self._ref_chains)
    #     root_node.appendChild(ref_chains_path_node)
    #
    #     model_chains_path_node = root.createElement("model_chains")
    #     model_chains_path_node.setAttribute(attribute, self._model_chains)
    #     root_node.appendChild(model_chains_path_node)
    #
    #     results_path_path_node = root.createElement("results_path")
    #     results_path_path_node.setAttribute(attribute, self.get_results_path())
    #     root_node.appendChild(results_path_path_node)
    #
    #     with open(f"{self.get_project_path()}/{self.get_project_name()}.xml", "w", encoding="utf-8") as file:
    #         file.write(root.toprettyxml())
