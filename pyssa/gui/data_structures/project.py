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
import json
import os
import pathlib
import platform
from datetime import datetime
from pathlib import Path
from xml.etree import ElementTree
from xml.dom import minidom
from pyssa.gui.utilities import gui_utils
from pyssa.pymol_protein_tools import protein_pair
from pyssa.pymol_protein_tools import protein


class Project:
    """This class is for the projects used in the plugin

    Attributes:
        _project_name:
            The name of the project.
        _workspace:
            The path of the projects' workspace.
        project_path:
            The path where the project is located.
        _date:
            The creation date of the project.
        _operating_system:
            The operating system, where the project was created on.
        _folder_paths:
            A list of all major paths of the project.
        _session_file_name:
            The name of the pymol session.
        proteins:
            A list of protein.Protein objects which are used within the project.
        protein_pairs:
            A list of protein_pair.ProteinPair objects which are used within the project.

    """
    def __init__(self, project_name: str, workspace_path: pathlib.Path) -> None:
        """Constructor

        Args:
            project_name:
                name of the project
            workspace_path:
                path of the workspace
        """
        self._project_name: str = project_name
        self._workspace: pathlib.Path = workspace_path

        self.project_path: pathlib.Path = Path(f"{workspace_path}/{project_name}")
        self._date = str(datetime.now())
        self._operating_system = platform.system()
        # self._workspace_path: pathlib.Path = ""
        self._folder_paths: list[Path] = []
        self._session_file_name: str = "session_file_model_s.pse"
        self.proteins: list[protein.Protein] = []
        self.protein_pairs: list[protein_pair.ProteinPair] = []
        self.create_folder_paths()

    def add_new_protein(self, protein_name) -> None:
        """This function adds a new protein to the project.

        Args:
            protein_name:
                name of the protein
        """
        new_protein = protein.Protein(protein_name)
        self.proteins.append(new_protein)

    def add_existing_protein(self, value_protein: protein.Protein) -> None:
        """This function adds an existing protein object to the project.

        Args:
            protein:
                name of the existing protein object
        """
        self.proteins.append(value_protein)

    def add_protein_pair(self, value_protein_pair: protein_pair.ProteinPair) -> None:
        """This function adds an existing protein_pair object to the project.

        Args:
            protein_pair:
                name of the existing protein_pairs object
        """
        self.protein_pairs.append(value_protein_pair)

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

    def create_folder_paths(self):
        self._folder_paths = [
            Path(f"{self._workspace}/{self._project_name}"),
            Path(f"{self._workspace}/{self._project_name}/pdb"),
            Path(f"{self._workspace}/{self._project_name}/results"),
            Path(f"{self._workspace}/{self._project_name}/.objects/"),
            Path(f"{self._workspace}/{self._project_name}/.objects/proteins"),
            Path(f"{self._workspace}/{self._project_name}/.objects/protein_pairs"),
            # Path(f"{workspace_path}/{self._project_name}/results/alignment_files"),
            # Path(f"{workspace_path}/{self._project_name}/results/distance_csv"),
            # Path(f"{workspace_path}/{self._project_name}/results/images"),
            # Path(f"{workspace_path}/{self._project_name}/results/images/interesting_regions"),
            # Path(f"{workspace_path}/{self._project_name}/results/plots"),
            # Path(f"{workspace_path}/{self._project_name}/results/plots/distance_histogram"),
            # Path(f"{workspace_path}/{self._project_name}/results/plots/distance_plot"),
            # Path(f"{workspace_path}/{self._project_name}/results/sessions/"),
        ]

    def get_session_file(self) -> str:
        """This function gets the value of the session_file variable

        Returns (str):
            session_file as complete path with file name
        """
        return f"{self._folder_paths[10]}/{self._session_file_name}"

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

    def get_results_path(self):
        """This function returns the results path of the project

        Returns:
            the results path
        """
        return str(self._folder_paths[2])

    def get_objects_path(self):
        """This function returns the objects path of the project

        Returns:
            the object path
        """
        return str(self._folder_paths[3])

    def get_objects_proteins_path(self):
        """This function returns the objects path of the project

                Returns:
                    the object path
                """
        return str(self._folder_paths[4])

    def get_objects_protein_pairs_path(self):
        """This function returns the objects path of the project

                Returns:
                    the object path
                """
        return str(self._folder_paths[5])

    def get_number_of_proteins(self):
        return len(self.proteins)

    def get_project_name(self):
        return self._project_name

    def set_folder_paths(self, value):
        self._folder_paths = value

    def serialize_project(self, filepath, filename) -> None:
        """This function serialize the protein object

        """
        if not os.path.exists(filepath):
            print(f"The filepath: {filepath} does not exists!")
            return
        # self._folder_paths.__str__()
        project_dict = self.__dict__
        protein_names = []
        for tmp_protein in self.proteins:
            protein_names.append(tmp_protein.molecule_object)
        protein_pair_names = []
        for tmp_protein_pair in self.protein_pairs:
            protein_pair_names.append(tmp_protein_pair.name)
        update = {
            '_workspace': str(self._workspace),
            'project_path': str(self.project_path),
            '_folder_paths': [str(self._folder_paths[0]), str(self._folder_paths[1]), str(self._folder_paths[2]), str(self._folder_paths[3]), str(self._folder_paths[4]), str(self._folder_paths[5])],
            'proteins': protein_names,
            'protein_pairs': protein_pair_names,
        }
        project_dict.update(update)
        print(project_dict)
        project_file = open(f"{filepath}\\{filename}.json", "w", encoding="utf-8")
        json.dump(project_dict, project_file, indent=4)

    @staticmethod
    def deserialize_project(filepath):
        """This function constructs the protein object from
        the json file

        Returns:
            a complete project object deserialized from a json file
        """
        try:
            project_obj_file = open(pathlib.Path(f"{filepath}/project.json"), "r", encoding="utf-8")
        except FileNotFoundError:
            print(f"There is no valid json file under: {filepath}")
            return
        project_dict = json.load(project_obj_file)
        tmp_project: Project = Project(project_dict.get("_project_name"), project_dict.get("_workspace"))
        tmp_project.set_folder_paths(project_dict.get("_folder_paths"))
        tmp_project.project_path = project_dict.get("project_path")
        # adding proteins to projects object
        if len(os.listdir(tmp_project.get_objects_proteins_path())) > 0:
            for tmp_protein_file in os.listdir(tmp_project.get_objects_proteins_path()):
                tmp_protein = protein.Protein.deserialize_protein(pathlib.Path(f"{tmp_project.get_objects_proteins_path()}/{tmp_protein_file}"))
                tmp_project.add_existing_protein(tmp_protein)
        # adding protein pairs to projects object
        if len(os.listdir(tmp_project.get_objects_protein_pairs_path())) > 0:
            for tmp_protein_pair_file in os.listdir(tmp_project.get_objects_protein_pairs_path()):
                tmp_protein_pair = protein_pair.ProteinPair.deserialize_protein_pair(pathlib.Path(f"{tmp_project.get_objects_protein_pairs_path()}/{tmp_protein_pair_file}"))
                tmp_project.add_protein_pair(tmp_protein_pair)
        return tmp_project

    def search_protein(self, protein_name):
        for tmp_protein in self.proteins:
            if tmp_protein.molecule_object == protein_name:
                return tmp_protein
        print(f"No matching protein with the name {protein_name} found.")

    def search_protein_pair(self, protein_pair_name):
        """This function searches all protein_pairs within the project and returns true if the project contains the pair

        Args:
            protein_pair_name:
                name of the pair to search
        Returns:
            True: is in the project
            False: is NOT in the project
        """
        for tmp_protein_pair in self.protein_pairs:
            if tmp_protein_pair.name == protein_pair_name:
                return True
        print(f"No matching protein with the name {protein_pair_name} found.")
        return False

    def get_specific_protein_pair(self, protein_pair_name):
        """This function gets a specific protein_pair by name from the project

        Args:
            protein_pair_name:
                name of the pair to get
        Returns:
            the actual protein pair
        """
        for tmp_protein_pair in self.protein_pairs:
            if tmp_protein_pair[2].name == protein_pair_name:
                return tmp_protein_pair[2]
        print(f"No matching protein with the name {protein_pair_name} found.")

    def get_specific_protein_pair_tuple(self, protein_pair_name):
        """This function gets a specific protein_pair by name from the project

        Args:
            protein_pair_name:
                name of the pair to get
        Returns:
            the actual protein pair
        """
        for tmp_protein_pair in self.protein_pairs:
            if tmp_protein_pair[2].name == protein_pair_name:
                return tmp_protein_pair
        print(f"No matching protein with the name {protein_pair_name} found.")